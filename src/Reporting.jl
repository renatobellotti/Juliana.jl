using DataFrames
using LaTeXStrings
using Printf
using PrettyTables
using Statistics


function homogeneity(dose, mask)
    return Juliana.dvh_d(dose, mask, 95f0) - Juliana.dvh_d(dose, mask, 5f0)
end


function constraint_fulfillment_df(reports)
    # Make sure the constraint fulfillments are the same for all reports.
    @assert all([length(r.constraint_fulfillments) == length(reports[1].constraint_fulfillments) for r in reports])
    if length(reports) > 1
        for (a, b) in zip([r.constraint_fulfillments for r in reports]...)
            @assert a[1] == b[1]
        end
    end

    # Build a table that summarises the constraint fulfillments.
    table = Vector{Dict{String, Union{String, Real}}}(undef, length(reports[1].constraint_fulfillments))
    constraints = [entry[1] for entry in reports[1].constraint_fulfillments]
    for (i, constraint) in enumerate(constraints)
        quantity = ""
        if Juliana.is_maximum_constraint(constraint)
            quantity = "Dmax"
        elseif constraint.kind == Juliana.constraint_mean
            quantity = "Dmean"
        elseif constraint.kind == Juliana.constraint_dvh_d
            quantity = "D$(@sprintf("%.2f", constraint.volume*100))%"
        else
            @error "Unknown constraint type $(constraint.kind); $(constraint)"
        end

        table[i] = Dict{String, Union{String, Real}}(
            "structure_name" => constraint.structure_name,
            "priority" => "$(constraint.priority)",
            "quantity" => quantity,
            "threshold" => constraint.dose,
        )

        for report in reports
            achieved = report.constraint_fulfillments[i][2]
            table[i]["$(report.label)"] = @sprintf("%.1f", achieved)
            table[i]["fulfilled_$(report.label)"] = achieved <= constraint.dose ? "yes" : "no"
        end
    end
    
    # Build a DataFrame and sort the columns.
    ordered_columns = [
        "structure_name",
        "priority",
        "quantity",
        "threshold",
    ]

    for report in reports
        push!(ordered_columns, "$(report.label)")
    end

#     for report in reports
#         push!(ordered_columns, "fulfilled_$(report.label)")
#     end

    return DataFrame(table)[!, ordered_columns]
end


function target_coverage_df(reports)
    @assert all([length(r.target_coverages) == length(reports[1].target_coverages) for r in reports])

    table = Vector{Dict}(undef, length(reports))

    for (i, report) in enumerate(reports)
        table[i] = Dict(
            "label" => report.label,
            "Dmax (overall)" => report.Dmax,
        )

        for coverage in report.target_coverages
            table[i]["Dmin_$(coverage.target_name) [Gy]"] = coverage.Dmin
            table[i]["V95_$(coverage.target_name) [%]"] = coverage.V95
            table[i]["Dmax_$(coverage.target_name) [Gy]"] = coverage.Dmax
            table[i]["HI_$(coverage.target_name)"] = coverage.HI_D5_D95
        end
    end

    ordered_columns = [
        "label",
        "Dmax (overall)",
    ]

    for coverage in reports[1].target_coverages
        push!(ordered_columns, "Dmin_$(coverage.target_name) [Gy]")
        push!(ordered_columns, "V95_$(coverage.target_name) [%]")
        push!(ordered_columns, "Dmax_$(coverage.target_name) [Gy]")
        push!(ordered_columns, "HI_$(coverage.target_name)")
    end

    df_target = DataFrame(table)[!, ordered_columns]
    
    df_new = DataFrame(collect(Matrix(df_target)[:, 2:end]'), df_target[!, "label"])
    df_new[!, "label"] = names(df_target)[2:end]
    
    ordered_columns = ["label"]
    for report in reports
        push!(ordered_columns, report.label)
    end
    df_new = df_new[!, ordered_columns]
    
    for j in 2:size(df_new, 2)
        for i in 1:size(df_new, 1)
            df_new[i, j] = @sprintf("%.1f", df_new[i, j])
        end
    end
    
    df_new.label = convert.(LaTeXString, df_new.label)
    
    return df_new
end


function reports_to_pdf(reports, work_dir; paths_to_plots=[])
    df_constraints = Juliana.constraint_fulfillment_df(reports)
    df_target = Juliana.target_coverage_df(reports)

    latex_report_path = "$(work_dir)/output_file.tex"

    tmp_target = copy(df_target)
    tmp_target.label = replace.(tmp_target.label, "_" => "\\_")

    tmp_constraints = copy(df_constraints)
    tmp_constraints.structure_name = replace.(tmp_constraints.structure_name, "_" => "\\_")
    rename!(tmp_constraints, replace.(names(tmp_constraints), "_" => "\\_"))
    tmp_constraints.priority = replace.(tmp_constraints.priority, "_" => "\\_")


    open(latex_report_path, "w") do file
        write(file, raw"\documentclass[a4paper]{article}")
        write(file, "\n")
        write(file, raw"\usepackage[margin=2cm]{geometry}")
        write(file, "\n")
        write(file, raw"\usepackage{adjustbox}")
        write(file, "\n")
        write(file, raw"\usepackage{booktabs}")
        write(file, "\n")
        write(file, raw"\usepackage[dvipsnames]{xcolor}")
        write(file, "\n")
        write(file, "\n")
        write(file, raw"\begin{document}")
        write(file, "\n")

        # Target metrics table.
        pretty_table(file, tmp_target, backend=Val(:latex), tf=tf_latex_booktabs)
        write(file, "\n")

        # OAR metrics table.
        write(file, "\n")
        highlighters = Vector{Union{Tuple, LatexHighlighter}}(undef, 0)
        # c == column
        for c in 1:length(reports)
            push!(
                highlighters,
                # Highlight red if the column does not fulfill the threshold
                LatexHighlighter((data, i, j) -> (j == 4+c) && (parse(Float32, data[i, j]) > data[i, 4]), ["color{red}","textbf"]),
            )
            push!(
                highlighters,
                # Highlight red if the column does not fulfill the threshold
                LatexHighlighter((data, i, j) -> (j == 4+c) && (parse(Float32, data[i, j]) <= data[i, 4]), ["color{ForestGreen}","textbf"]),
            )
        end
        constraint_table = pretty_table(
            String,
            tmp_constraints,
            backend=Val(:latex), tf=tf_latex_booktabs,
            highlighters=Tuple(highlighters),
        )
        constraint_table = replace(
            constraint_table,
            raw"\begin{tabular}" => raw"\begin{adjustbox}{width=\textwidth}\begin{tabular}",
            raw"\end{tabular}" => raw"\end{tabular}\end{adjustbox}",
        )
        write(file, constraint_table)
        write(file, "\n")

        for (i, path) in enumerate(paths_to_plots)
            write(file, raw"\includegraphics[width=0.45\textwidth]{" * path * raw"}")
            if mod(i, 2) == 0
                write(file, raw"\\\\")
            end
            write(file, "\n")
        end

        write(file, raw"\end{document}")
    end
    
    build_dir = "$(work_dir)/build"
    mkpath(build_dir)
    read(`pdflatex -output-directory=$(build_dir) $latex_report_path`, String)
    
    return "$(build_dir)/output_file.pdf"
end


###################################
# Reports for robustness analysis.
###################################
function constraint_fulfillment_df(report::RobustReport)
    reports = report.reports
    # Make sure the constraint fulfillments are the same for all reports.
    @assert all([length(r.constraint_fulfillments) == length(reports[1].constraint_fulfillments) for r in reports])
    if length(reports) > 1
        for (a, b) in zip([r.constraint_fulfillments for r in reports]...)
            @assert a[1] == b[1]
        end
    end

    metrics = (
        ("minimum", minimum),
        ("median", median),
        ("mean", mean),
        ("maximum", maximum),
    )

    # Build a table that summarises the constraint fulfillments.
    table = Vector{Dict{String, Union{String, Real}}}(undef, length(reports[1].constraint_fulfillments))
    constraints = [entry[1] for entry in reports[1].constraint_fulfillments]
    for (i, constraint) in enumerate(constraints)
        quantity = ""
        if Juliana.is_maximum_constraint(constraint)
            quantity = "Dmax"
        elseif constraint.kind == Juliana.constraint_mean
            quantity = "Dmean"
        elseif constraint.kind == Juliana.constraint_dvh_d
            quantity = "D$(@sprintf("%.2f", constraint.volume*100))%"
        else
            @error "Unknown constraint type $(constraint.kind); $(constraint)"
        end

        table[i] = Dict{String, Union{String, Real}}(
            "structure_name" => constraint.structure_name,
            "priority" => "$(constraint.priority)",
            "quantity" => quantity,
            "threshold" => constraint.dose,
        )

        values = Vector{Float32}(undef, length(reports))
        for (j, report) in enumerate(reports)
            achieved = report.constraint_fulfillments[i][2]
            values[j] = achieved
        end

        for (label, f) in metrics
            metric = f(values)
            table[i][label] = @sprintf("%.1f", metric)
            table[i]["fulfilled_$(label)"] = metric <= constraint.dose ? "yes" : "no"
        end
    end

    # Build a DataFrame and sort the columns.
    ordered_columns = [
        "structure_name",
        "priority",
        "quantity",
        "threshold",
    ]

    for (label, _) in metrics
        push!(ordered_columns, label)
    end

    return DataFrame(table)[!, ordered_columns]
end


function target_coverage_df(report::RobustReport)
    reports = report.reports
    @assert all([length(r.target_coverages) == length(reports[1].target_coverages) for r in reports])

    n_targets = length(reports[1].target_coverages)
    n_scenarios = length(report.reports)

    dfs = Dict()
    for i in 1:n_targets
        target_name = reports[1].target_coverages[i].target_name

        values_Dmin = Vector{Float32}(undef, n_scenarios)
        values_Dmax = Vector{Float32}(undef, n_scenarios)
        values_V95 = Vector{Float32}(undef, n_scenarios)
        values_HI = Vector{Float32}(undef, n_scenarios)

        for j in 1:n_scenarios
            values_Dmin[j] = reports[j].target_coverages[i].Dmin
            values_Dmax[j] = reports[j].target_coverages[i].Dmax
            values_V95[j] = reports[j].target_coverages[i].V95
            values_HI[j] = reports[j].target_coverages[i].HI_D5_D95
        end
        
        table = []
        
        # Target Dmin.
        push!(table, Dict(
            "target name" => target_name,
            "metric" => "Dmin",
            "minimum" => minimum(values_Dmin),
            "median" => median(values_Dmin),
            "mean" => mean(values_Dmin),
            "maximum" => maximum(values_Dmin),
        ))
        # Target Dmax.
        push!(table, Dict(
            "target name" => target_name,
            "metric" => "Dmax",
            "minimum" => minimum(values_Dmax),
            "median" => median(values_Dmax),
            "mean" => mean(values_Dmax),
            "maximum" => maximum(values_Dmax),
        ))
        # Target V95.
        push!(table, Dict(
            "target name" => target_name,
            "metric" => "V95",
            "minimum" => minimum(values_V95),
            "median" => median(values_V95),
            "mean" => mean(values_V95),
            "maximum" => maximum(values_V95),
        ))
        # Target HI.
        push!(table, Dict(
            "target name" => target_name,
            "metric" => "HI",
            "minimum" => minimum(values_HI),
            "median" => median(values_HI),
            "mean" => mean(values_HI),
            "maximum" => maximum(values_HI),
        ))

        df = DataFrame(table)[!, ["target name", "metric", "minimum", "median", "mean", "maximum"]]
        df[:, ["minimum", "median", "mean", "maximum"]] .= round.(df[:, ["minimum", "median", "mean", "maximum"]], digits=1)
        dfs[target_name] = df
    end

    return dfs
end


function reports_to_pdf(report::Juliana.RobustReport, work_dir; paths_to_plots=[])
    df_constraints = constraint_fulfillment_df(report)
    df_targets = target_coverage_df(report)

    latex_report_path = "$(work_dir)/output_file.tex"

#     tmp_target = copy(df_target)
#     tmp_target.label = replace.(tmp_target.label, "_" => "\\_")

    tmp_constraints = copy(df_constraints)
    tmp_constraints.structure_name = replace.(tmp_constraints.structure_name, "_" => "\\_")
    rename!(tmp_constraints, replace.(names(tmp_constraints), "_" => "\\_"))
    tmp_constraints.priority = replace.(tmp_constraints.priority, "_" => "\\_")


    open(latex_report_path, "w") do file
        write(file, raw"\documentclass[a4paper]{article}")
        write(file, "\n")
        write(file, raw"\usepackage[margin=2cm]{geometry}")
        write(file, "\n")
        write(file, raw"\usepackage{adjustbox}")
        write(file, "\n")
        write(file, raw"\usepackage{booktabs}")
        write(file, "\n")
        write(file, raw"\usepackage[dvipsnames]{xcolor}")
        write(file, "\n")
        write(file, "\n")
        write(file, raw"\begin{document}")
        write(file, "\n")

        # Target metrics table.
        for (_, df) in df_targets
            df = copy(df)
            df[:, "target name"] .= replace.(df[:, "target name"], "_" => "\\_")
            pretty_table(file, df, backend=Val(:latex), tf=tf_latex_booktabs)
            write(file, "\n")
        end

        # OAR metrics table.
        write(file, "\n")
        highlighters = Vector{Union{Tuple, LatexHighlighter}}(undef, 0)
        # c == column
        # -4: Each column that is not in ["structure_name", "priority", "quantity", "threshold"]
        # is a result.
        n_results = size(df_constraints, 2)-4
        for c in 1:n_results
            push!(
                highlighters,
                # Highlight red if the column does not fulfill the threshold
                LatexHighlighter((data, i, j) -> (j == 4+c) && (parse(Float32, data[i, j]) > data[i, 4]), ["color{red}","textbf"]),
            )
            push!(
                highlighters,
                # Highlight red if the column does not fulfill the threshold
                LatexHighlighter((data, i, j) -> (j == 4+c) && (parse(Float32, data[i, j]) <= data[i, 4]), ["color{ForestGreen}","textbf"]),
            )
        end
        constraint_table = pretty_table(
            String,
            tmp_constraints,
            backend=Val(:latex), tf=tf_latex_booktabs,
            highlighters=Tuple(highlighters),
        )
        constraint_table = replace(
            constraint_table,
            raw"\begin{tabular}" => raw"\begin{adjustbox}{width=\textwidth}\begin{tabular}",
            raw"\end{tabular}" => raw"\end{tabular}\end{adjustbox}",
        )
        write(file, constraint_table)
        write(file, "\n")

        for (i, path) in enumerate(paths_to_plots)
            write(file, raw"\includegraphics[width=0.45\textwidth]{" * path * raw"}")
            if mod(i, 2) == 0
                write(file, raw"\\\\")
            end
            write(file, "\n")
        end

        write(file, raw"\end{document}")
    end
    
    build_dir = "$(work_dir)/build"
    mkpath(build_dir)
    read(`pdflatex -output-directory=$(build_dir) $latex_report_path`, String)
    
    return "$(build_dir)/output_file.pdf"
end