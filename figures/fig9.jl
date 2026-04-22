include("plotting_env_and_load.jl")
## Figure nine
# Batchelor lengths and resolution
files = (l_R_ρ_105_dT2_diagnostics, nl_R_ρ_105_dT2_diagnostics,
         l_R_ρ_105_dT1_diagnostics, nl_R_ρ_105_dT1_diagnostics,
         l_R_ρ_105_dT05_diagnostics, nl_R_ρ_105_dT05_diagnostics)
fig = Figure(size = (800, 800))
ax = [Axis(fig[i, 1],
          title = "Batchelor lengths and grid resolution",
          xlabel = "time (mins)",
          ylabel = "Length (log10(mm))") for i ∈ 1:3]
ts_length = 3:240
linestyles = (:solid, :dash)

for (i, file) ∈ enumerate(files)

    ls = i % 2 == 0 ? linestyles[2] : linestyles[1]
    jldopen(file) do f
        j = i ≤ 2 ? 1 : i > 4 ? 3 : 2
        Ba = 1e3 * f["Ba"][ts_length]
        lines!(ax[j], log10.(Ba), label = all_labels[i], linestyle = ls, color = Makie.wong_colors()[i])

        j > 1 ? ax[j].title = "" : nothing
        j < 3 ? hidexdecorations!(ax[j], grid = false, ticks = false) : nothing
    end

        jldopen(file) do f
            Ba = 1e3 * f["Ba"][ts_length]
            spcaings = vcat(diff(f["dims"]["x_caa"]),
                            diff(f["dims"]["y_aca"]),
                            diff(f["dims"]["z_aac"]))
            Δ = minimum(spcaings) * 1e3
            println(Δ)
            j = i ≤ 2 ? 1 : i > 4 ? 3 : 2
            i % 2 == 0 ? lines!(ax[j], log10.(fill(Δ, length(Ba))), color = :red, linestyle = :dot,
                                label = "Δ") : nothing
        end

end
files = (l_fbc_diagnostics, nl_fbc_diagnostics)
fbc_labels = ("II′", "V′")
fbc_colors = (Makie.wong_colors()[7], RGBAf(0.5, 0.2, 0.8, 0.5))
for (i, file) ∈ enumerate(files)

    jldopen(file) do f
        Ba = 1e3 * f["Ba"][ts_length]
        ls = i % 2 == 0 ? linestyles[2] : linestyles[1]
        lines!(ax[2], log10.(Ba), label = fbc_labels[i], linestyle = ls, color = fbc_colors[i])
    end

end
for i ∈ 1:3
    axislegend(ax[i], position = :rt, nbanks = 3)
end
linkyaxes!(ax[2], ax[1])
linkyaxes!(ax[3], ax[1])
linkxaxes!(ax[2], ax[1])
linkxaxes!(ax[3], ax[1])
xlims!(ax[1], 2, 240)
# panel labels
panel_labels = ["(a)", "(b)", "(c)"]
for (i, a) ∈ enumerate(ax)
    text!(a, 0, 1,
          text = panel_labels[i],
          font = :bold,
          align = (:left, :top),
          offset = (35, 2),
          space = :relative,
          fontsize = 24
          )
end
##
save("fig9.png", fig)
