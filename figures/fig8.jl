include("plotting_env_and_load.jl")
## Figure eight
# Flux boundary conditions
fig = Figure(size = (800, 800))
files = (l_fbc_diagnostics, nl_fbc_diagnostics)
labels = ("II′","V′")
titles = ("(a) Experiment II′",
            "(b) Experiment V′")
ax = [Axis(fig[i, 1]) for i ∈ 1:3]
ts_range = 1:420
for (i, file) ∈ enumerate(files)

    f = jldopen(file)
    ΔS = f["ΔS"][ts_range]
    ΔT = f["ΔT"][ts_range]
    Rᵨ = f["R_ρ"][ts_range]
     t = f["dims"]["time"][ts_range] ./ 60
    z✶ = f["dims"]["z✶"][f["S_interface_idx"][ts_range]]
    close(f)

    lines!(ax[i], t, ΔS ./ ΔS[1], color = :blue, label = L"ΔS / ΔS_{0}")
	lines!(ax[i], t, ΔT ./ ΔT[1], color = :red, label = L"ΔΘ / ΔΘ_{0}")
	ylims!(ax[i], 0, 1.1)
    ax[i].title = titles[i]
	ax2 = Axis(fig[i, 1],
                yaxisposition = :right,
			    yticklabelcolor = :dodgerblue,
			    rightspinecolor = :dodgerblue,
				ylabelcolor = :dodgerblue,
			    ytickcolor = :dodgerblue,
				ylabel = L"R_{\rho}")
	lines!(ax2, t, Rᵨ, color = :dodgerblue, label = L"R_{\rho}")
    linkxaxes!(ax2, ax[i])
    hidexdecorations!(ax2)
	axislegend(ax[i], merge = true, position = :rb)
    lines!(ax[3], t[1:end], z✶, label = labels[i])
    ax[3].xlabel = "time (min)"
    ax[3].ylabel = L"$z*$ (m)"
            ax[3].title = "(c) Salinity interface height"
end

files = (l_R_ρ_105_dT1_diagnostics, nl_R_ρ_105_dT1_diagnostics)
labels_nf = ("II", "V")
for (i, file) ∈ enumerate(files)
    jldopen(file) do f
        t = f["dims"]["time"]
        z✶ = f["dims"]["z✶"][f["S_interface_idx"]]
        Rᵨ = f["R_ρ"]
        lines!(ax[3], t[2:end] ./ 60, z✶, label = labels_nf[i],
                linestyle = :dash, color = Makie.wong_colors()[i])
    end
end

hidexdecorations!(ax[1], grid = false, ticks = false)
hidexdecorations!(ax[2], grid = false, ticks = false)
fig
# linear fit
linfit_range = 250:420
Rᵨ_mean_l = jldopen(l_fbc_diagnostics) do f
    mean(f["R_ρ"][linfit_range])
end
f = jldopen(l_fbc_diagnostics)
z✶ = f["dims"]["z✶"][f["S_interface_idx"][linfit_range]]
t = f["dims"]["time"][linfit_range]
Rᵨ_mean_nl = mean(f["R_ρ"][linfit_range])
A = [ones(length(t)) t]
a, b = A \ z✶
lines!(ax[3], t./60, a .+ t .* b, label = "Linear fit", linestyle = :dot, linewidth = 3)
Legend(fig[4, :], ax[3], orientation = :horizontal, nbanks = 2)
fig
##
save("fig8.png", fig)
