module QoPlots

export plot_dm, plot_dm_square_size_legend, plot_bloch

include("phase_cmap.jl")

using Colors
using Printf
using PyCall
using PyPlot
using QuantumOpticsBase

rgb(c) = convert(RGB, c) |> a -> (a.r, a.g, a.b)

function brighten(color, factor)
    luv = convert(Luv, color)
    rgb(Luv(luv.l * clamp(factor, 0.0, 1.0), luv.u, luv.v))
end

function draw_square(x, y, area, color, ax; angle_tick=false, angle=0.0)
    hs = area^0.5 / 2
    xcorners = [x - hs, x + hs, x + hs, x - hs]
    ycorners = [y - hs, y - hs, y + hs, y + hs]

    square = ax.fill(xcorners, ycorners, color=rgb(color))
    for a in square
        a.set_clip_on(false)
    end

    if angle_tick && hs > 0.04
        len = hs * 0.8
        ax.add_line(plt.Line2D([x, x + len * cos(angle)],
            [y, y + len * sin(angle)], color=brighten(color, 0.6), linewidth=0.8))
    end
end

function colour(arg)
    fraction = arg / 2π + 1.5
    idx = round(Int, (fraction % 1) * (length(cmap_kovesi_c6) - 1)) + 1
    cmap_kovesi_c6[idx]
end
#colour(arg) = RGB(cmocean.cm.phase((arg / 2π + 0.63) % 1)[1:3]...)

function draw_phase_ring(x, y, radius, ax)
    mc = pyimport("matplotlib.collections")
    mp = pyimport("matplotlib.patches")

    patches = []
    steps = 256
    degree_each = 360 / steps
    for i in 1:steps
        push!(patches, mp.Wedge((x, y), radius,
            (i - 1) * degree_each, i * degree_each, width=radius / 2))
    end
    p = mc.PatchCollection(patches)
    p.set_color([rgb(colour(i / steps * 2π)) for i in 0:(steps - 1)])
    ax.add_collection(p)

    padded = 1.1 * radius
    color = rgb(convert(Luv, colour(0)) |> a-> Luv(a.l, 0, 0))
    ax.text(x + padded, y, L"$+$", ha="left", va="center", size=6, color=color)
    ax.text(x - padded, y, L"$-$", ha="right", va="center", size=6, color=color)
    ax.text(x, y + padded, L"$\mathrm{i}$", ha="center", va="bottom", size=6, color=color)
    ax.text(x, y - padded, L"$-\mathrm{i}$", ha="center", va="top", size=6, color=color)
end

function computational_basis_labels(dims)
    [join("$v" for v in t) for t in Base.product((0:(d - 1) for d in dims)...)][:]
end


function plot_dm_square_size_legend(x, y, scale, ax; direction=:horizontal, num_dims=4, textsize=8)
    legend_color = convert(Luv, colour(0)) |> a-> Luv(a.l, 0, 0)

    if num_dims == 2
        areas = [0.5, 0.01]
    else
        # Decide whether to show the 0.5 or 0.001 square depending on the scale, i.e.
        # the data values in the corresponding plot.
        # (The cutoff is just chosen based on what empirically looked good).
        show_0p5 = scale < 1.72
        areas = [0.5, 0.25, 0.1, 0.01, 0.001][show_0p5 ? (1:end - 1) : (2:end)]
    end

    next_square_pos = [x, y]
    if direction == :horizontal
        inc_dir = [1, 0]
    elseif direction == :vertical
        inc_dir = [0, -1]
    else
        error("Unexpected legend direction: '$(direction)'")
    end

    for (i, area) in enumerate(areas)
        draw_square(next_square_pos..., scale * area, legend_color, ax)
        size_ = sqrt(scale * area)

        ax.text((next_square_pos - [0, 0.06 + size_ / 2])..., "\$ $area\$",
            horizontalalignment="center",
            verticalalignment="top",
            color=rgb(legend_color),
            size=textsize)

        next_square_pos += (0.8size_ + 0.55) * inc_dir
    end

    draw_phase_ring((next_square_pos + 0.3 * inc_dir)..., 0.25, ax)

    ax.axis("equal")
end

function plot_dm(ρ::AbstractArray; xlabels=nothing, ylabels=nothing, show_legend=true,
                 show_numbers=false, hide_zeros=false, scale=nothing, textsize=8, ax=nothing)
    num_rows, num_cols = size(ρ)

    if ax === nothing
        figsize = [2, 2] * num_rows / 2
        if show_legend
            figsize += [0.0, 0.5]
        end
        _, ax = subplots(1, 1, figsize=figsize, constrained_layout=true)
    end

    if xlabels === nothing
        xlabels = []
    end
    if ylabels === nothing
        ylabels = []
    end

    ax.axis("equal")

    if length(xlabels) == 0 && length(ylabels) == 0
        ax.axis("off")
    end
    ax.set_frame_on(false)

    #ax.fill([0, num_cols, num_cols, 0], [0, 0, num_rows, num_rows], color="w")

    if scale === nothing
        max_magn = maximum(abs.(ρ))
        scale = 1 / (2 * max_magn)
    end

    path_effects = pyimport("matplotlib.patheffects")

    for x in 1:num_cols
        for y in 1:num_rows
            z = ρ[y, x]

            if !show_numbers || x >= y
                # Draw square.
                if abs(z) > 0.0
                    draw_square(x, num_rows + 1 - y, abs(z) * scale,
                        colour(angle(z)), ax, angle=angle(z), angle_tick=(x != y),)
                end
            end

            if show_numbers && x <= y && (!hide_zeros || abs(z) >= 1e-3)
                # Draw value text.
                if x == y
                    # On the diagonal, don't display phase term.
                    abstext = @sprintf "\$ %.3f \$" abs(z)
                    t = ax.text(x, num_rows + 1 - y - 0.023, abstext,
                        horizontalalignment="center",
                        verticalalignment="center",
                        color=brighten(colour(0.0), 0.65),
                        size=textsize)
                    t.set_path_effects([
                        path_effects.Stroke(linewidth=textsize/4, foreground="white", alpha=0.3),
                        path_effects.Normal()])
                else
                    abstext = @sprintf "\$ %.3f \\cdot \$" abs(z)
                    ax.text(x, num_rows + 1 - y + 0.02, abstext,
                        horizontalalignment="center",
                        verticalalignment="bottom",
                        color=brighten(colour(0.0), 0.65),
                        size=textsize)

                    phasetext = @sprintf "\$ \\mathrm{e}^{\\mathrm{i} %.3f \\pi} \$" (mod(angle(z), 2π) / pi)
                    ax.text(x, num_rows + 1 - y - 0.04, phasetext,
                        horizontalalignment="center",
                        verticalalignment="top",
                        color=brighten(colour(0.0), 0.65),
                        size=textsize)
                end
            end
        end
    end

    ticks = collect(1.0:1.0:size(ρ)[1])

    xa = ax.xaxis
    xa.set_ticks(ticks)
    if length(xlabels) > 0
        xa.set_ticklabels(xlabels)
        ax.xaxis.tick_top()
    end

    ya = ax.yaxis
    ya.set_ticks(ticks)
    if length(ylabels) > 0
        ya.set_ticklabels(ylabels[end:-1:1])
    end

    if show_legend
        plot_dm_square_size_legend(0.25, -0.05, scale, ax, num_dims=num_cols, textsize=textsize)
        ax.spines["left"].set_position(("data", 0.4))
    end

    ax
end

function plot_dm(ρ::DenseOpType; xlabels=nothing, ylabels=nothing, kwargs...)
    if xlabels === nothing && ylabels === nothing
        ls = computational_basis_labels(ρ.basis_l.shape)
        xlabels = ["\$\\langle$l|\$" for l in ls]
        ylabels = ["\$|$l\\rangle\$" for l in ls]
    end
    plot_dm(ρ.data; xlabels=xlabels, ylabels=ylabels, kwargs...)
end

function plot_bloch(ρs)
    qutip = pyimport("qutip")
    b = qutip.Bloch()
    pos(ρ) = [real(expect(f(SpinBasis(1//2)), ρ)) for f in [sigmax, sigmay, sigmaz]]
    b.add_vectors([pos(ρ) for ρ in ρs])
    b
end

end # module
