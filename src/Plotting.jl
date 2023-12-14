using Colors, ColorSchemes
using PyPlot


function plot_dvh(
    evaluation_doses::AbstractArray{T, 1},
    dose_distributions,
    # key is the label, the value is the structure and its color
    structures) where {T}
    fig, ax = subplots(figsize=(8 * 1.3, 4.5))


    handles = [
        PyPlot.matplotlib.lines.Line2D([0], [0], color=:black, label=dose_distributions[1][1]),
        PyPlot.matplotlib.lines.Line2D([0], [0], color=:black, linestyle=:dashed, label=dose_distributions[2][1]),
        PyPlot.matplotlib.lines.Line2D([], [], visible=false, label="")
    ]

    for (label, (structure, color)) in structures
        volumes = Juliana.dvh_v(dose_distributions[1][2], structure.mask, evaluation_doses)
        nonzero = 1:findfirst(volumes .== 0)
        x = evaluation_doses[nonzero]
        y = volumes[nonzero]
        ax.plot(x, y, color=color)

        volumes = Juliana.dvh_v(dose_distributions[2][2], structure.mask, evaluation_doses)
        nonzero = 1:findfirst(volumes .== 0)
        x = evaluation_doses[nonzero]
        y = volumes[nonzero]
        ax.plot(x, y, linestyle=:dashed, color=color)

        push!(handles, PyPlot.matplotlib.lines.Line2D([0], [0], color=color, label=label))
    end

    ax.set_xlabel("Dose [Gy]", fontsize=18)
    ax.set_ylabel("Volume [%]", fontsize=18)
    ax.tick_params(labelsize=16)
    ax.set_xlim([0, evaluation_doses[end]])
    ax.set_ylim([0, 100])
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)

    ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.3, 1.0), fontsize=14);
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
                            colorscheme)
    fig, axes = subplots(
        1, 3,
        figsize=(8 * 1.3, 4.5),
        gridspec_kw=Dict("width_ratios" => [1, 1, 0.04]),
    )

    img = nothing

    if length(dose_distributions) == 0
        dose_distributions = [("", zeros(size(ct.data)...))]
    end

    for (ax, (label, dose_distribution)) in zip(axes, dose_distributions)
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
            points = structure.points[:, 1:3]
            points = points[points[:, 3] .≈ (iz-1)*ct.grid.spacing[3], :]

            points = (points' .- ct.grid.origin) ./ ct.grid.spacing
            points = (points[1:2, :] .- [x_start, y_start]) .+ 1
            ax.add_patch(PyPlot.matplotlib.patches.Polygon(
                points',
                fill=false,
                color=color,
                linewidth=3,
            ))
        end

        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines["bottom"].set_visible(false)
        ax.spines["top"].set_visible(false)
        ax.spines["left"].set_visible(false)
        ax.spines["right"].set_visible(false)

        ax.text(3, 3, label, color=:white, fontsize=20, fontweight=:bold, va="top")
    end

    # Add colorbar.
    axes[3].clear()
    cbar = fig.colorbar(img, ticks=[0, 20, 40, 60, 80, 95, 100, 105], cax=axes[3])
    cbar.ax.tick_params(labelsize=16)
    axes[3].set_ylabel("Dose [%]", fontsize=18, fontweight=:bold)

    fig.tight_layout()

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
