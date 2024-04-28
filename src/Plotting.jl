using Colors, ColorSchemes
using Printf
using PyPlot


function plot_dvh(
    evaluation_doses::AbstractArray{T, 1},
    dose_distributions,
    # key is the label, the value is the structure and its color
    structures;
    linewidth=1,
    legendfontsize=14,
    ticklabelfontsize=16,
    axislabelfontsize=18) where {T}
    fig, ax = subplots(figsize=(8 * 1.3, 4.5))

    n_doses = length(dose_distributions)
    @assert n_doses in [1, 2, 3]

    handles = [
        PyPlot.matplotlib.lines.Line2D([0], [0], color=:black, label=dose_distributions[1][1], linewidth=linewidth),
    ]
    if n_doses >= 2
        push!(
            handles,
            PyPlot.matplotlib.lines.Line2D([0], [0], color=:black, linestyle=:dashed, label=dose_distributions[2][1], linewidth=linewidth),
        )
    end
    if n_doses == 3
        push!(
            handles,
            PyPlot.matplotlib.lines.Line2D([0], [0], color=:black, linestyle=:dotted, label=dose_distributions[3][1], linewidth=linewidth),
        )
    end

    for (label, (structure, color)) in structures
        volumes = Juliana.dvh_v(dose_distributions[1][2], structure.mask, evaluation_doses)
        stop_index = findfirst(volumes .== 0)
        stop_index = isnothing(stop_index) ? length(volumes) : stop_index
        nonzero = 1:stop_index
        x = evaluation_doses[nonzero]
        y = volumes[nonzero]
        ax.plot(x, y, color=color, linewidth=linewidth)

        if n_doses >= 2
            volumes = Juliana.dvh_v(dose_distributions[2][2], structure.mask, evaluation_doses)
            stop_index = findfirst(volumes .== 0)
            stop_index = isnothing(stop_index) ? length(volumes) : stop_index
            nonzero = 1:stop_index
            x = evaluation_doses[nonzero]
            y = volumes[nonzero]
            ax.plot(x, y, linestyle=:dashed, color=color, linewidth=linewidth)
        end

        if n_doses == 3
            volumes = Juliana.dvh_v(dose_distributions[3][2], structure.mask, evaluation_doses)
            stop_index = findfirst(volumes .== 0)
            stop_index = isnothing(stop_index) ? length(volumes) : stop_index
            nonzero = 1:stop_index
            x = evaluation_doses[nonzero]
            y = volumes[nonzero]
            ax.plot(x, y, linestyle=:dotted, color=color, linewidth=linewidth)
        end

        push!(handles, PyPlot.matplotlib.lines.Line2D([0], [0], color=color, label=label, linewidth=linewidth))
    end

    ax.set_xlabel("Dose [Gy]", fontsize=axislabelfontsize, fontweight="bold")
    ax.set_ylabel("Volume [%]", fontsize=axislabelfontsize, fontweight="bold")
    ax.tick_params(labelsize=ticklabelfontsize)
    ax.set_xlim([0, evaluation_doses[end]])
    ax.set_ylim([0, 100])
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)

    ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.3, 1.0), fontsize=legendfontsize);
    fig.tight_layout()
            
    return fig, ax
end


function plot_distributions(ct,
                            dose_distributions,
                            structures,
                            x_start, x_end,
                            y_start, y_end,
                            iz,
                            target_dose,
                            colorscheme;
                            labelsize=20)
    @assert length(dose_distributions) in (0, 1, 2)
    
    fig, axes = nothing, nothing
    if length(dose_distributions) == 2
        fig, axes = subplots(
            1, 3,
            figsize=(8 * 1.3, 4.5),
            gridspec_kw=Dict("width_ratios" => [0.04, 1, 1]),
            layout="constrained",
        )
    else
        fig, axes = subplots(
            1, 2,
            figsize=(8, 4.5),
            gridspec_kw=Dict("width_ratios" => [0.04, 1]),
            layout="constrained",
        )
    end

    img = nothing

    if length(dose_distributions) == 0
        dose_distributions = [("", zeros(size(ct.data)...))]
    end

    for (ax, (label, dose_distribution)) in zip(axes[2:end], dose_distributions)
        ax.imshow(
            ct.data[x_start:x_end, y_start:y_end, iz]',
            vmin=-1000,
            vmax=3095,
            cmap="gray",
        );
        img = ax.imshow(
            dose_distribution[x_start:x_end, y_start:y_end, iz]' / target_dose * 100,
            vmin=0,
            vmax=110,
            cmap=ColorMap(colorscheme),
        )

        # Add contours.
        for (label, (structure, color)) in structures
            points = structure.points
            points = points[points[:, 3] .≈ (iz-1)*ct.grid.spacing[3], :]
            for p in Juliana.split_by_contour_index(points)
                p = p[:, 1:3]
                p = (p' .- ct.grid.origin) ./ ct.grid.spacing
                p = (p[1:2, :] .- [x_start, y_start]) .+ 1
                ax.add_patch(PyPlot.matplotlib.patches.Polygon(
                    p',
                    fill=false,
                    color=color,
                    linewidth=3,
                ))
            end
        end

        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines["bottom"].set_visible(false)
        ax.spines["top"].set_visible(false)
        ax.spines["left"].set_visible(false)
        ax.spines["right"].set_visible(false)

        ax.text(3, 3, label, color=:white, fontsize=labelsize, fontweight=:bold, va="top")
    end

    # Add colorbar.
    axes[1].clear()
    cbar = fig.colorbar(img, ticks=[0, 20, 40, 60, 80, 95, 100, 105], cax=axes[1])
    cbar.ax.tick_params(labelsize=16)
    axes[1].set_ylabel("Dose [%]", fontsize=18, fontweight=:bold)

    #fig.subplots_adjust(wspace=0.01)
    #fig.tight_layout()

    return fig, axes
end


function build_colorscheme(max_to_plot=110)
    HOTSPOT_INDEX = convert(Int64, floor(105 * 100 / max_to_plot) + 1)

    my_colorscheme = convert.(RGBA, colorschemes[:jet1].colors)
    new = copy(my_colorscheme)
    for (i, c) in enumerate(my_colorscheme)
        if i <= 15
            new[i] = RGBA(c.r, c.g, c.b, 0)
        else
            new[i] = RGBA(c.r, c.g, c.b, 0.6)
        end

        if i > HOTSPOT_INDEX
            new[i] = RGBA(100 / 256, 0, 100 / 256, 0.6)
        end
    end
    return new
end


function plot_depth_dose_curve(depth_dose_curves::Juliana.LookUpTable, energy::Real, ax)
    i = argmin(abs.(depth_dose_curves.energies .- energy))
    E = depth_dose_curves.energies[i]

    x = [depth_dose_curves.x0[i] + j * depth_dose_curves.dx[i] for j in 1:size(depth_dose_curves.table, 2)]
    y = depth_dose_curves.table[i, :]
    
    label = @sprintf("E = %6.2f [MeV]", energy / 1000)
    ax.plot(x, y, label=label)
    ax.legend(prop=Dict(
        "family" => "monospace",
    ))
end

function plot_depth_dose_curve(depth_dose_curves::Juliana.LookUpTable, energy::Real)
    fig, ax = subplots()

    ax.set_xlabel("Water-equivalent depth [cm]", fontsize=16)
    ax.set_ylabel("Value [a. u.]", fontsize=16)
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    plot_depth_dose_curve(depth_dose_curves, energy, ax)
    
    return fig, ax
end
