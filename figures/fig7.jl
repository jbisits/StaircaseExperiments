include("plotting_env_and_load.jl")
## Figure seven
fig = Figure(size = (800, 1000))
# density ratio and salinity midpoint
ax_interface = Axis(fig[1, 1],
                    title = "(a) Salinity interface height",
                    xlabel = "time (min)",
                    ylabel = L"$z*$ (m)")
ax_R_rho = Axis(fig[2, 1],
                title = "(b) Density ratio",
                xlabel = "time ( min)",
                ylabel = L"R_{\rho}")

for (i, file) ∈ enumerate(files)
    jldopen(file) do f
        t = f["dims"]["time"][1:ts_length]
        z✶ = f["dims"]["z✶"][f["S_interface_idx"]][1:ts_length]
        Rᵨ = f["R_ρ"][1:ts_length]
        ls = i % 2 == 0 ? :dash : :solid
        lines!(ax_interface, t ./ 60, z✶, label = all_labels[i], linestyle = ls)
        lines!(ax_R_rho, t ./ 60, Rᵨ, label = all_labels[i], linestyle = ls)
    end
end
hidexdecorations!(ax_R_rho, grid = false, ticks = false)
hidexdecorations!(ax_interface, grid = false, ticks = false)
linkxaxes!(ax_R_rho, ax_interface)

# energetics
files = (l_R_ρ_105_dT2_diagnostics, nl_R_ρ_105_dT2_diagnostics,
         l_R_ρ_105_dT1_diagnostics, nl_R_ρ_105_dT1_diagnostics,
         l_R_ρ_105_dT05_diagnostics, nl_R_ρ_105_dT05_diagnostics)
ax = [Axis(fig[i, 1], xlabel = "time (min)", ylabel = L"\mathcal{APE}") for i ∈ 3:4]
ax[1].title = "(c) Available potential energy"
ax[2].title = "(d) Time derivateive of APE, experiments I and IV"
ts_length = 240
j = 1
linestyles = (:solid, :dash)
for (i, f) ∈ enumerate(files)

    jldopen(f) do ds

        t = ds["dims"]["time"][1:ts_length]
        Ep = ds["Ep"][1:ts_length]
        Eb = ds["Eb"][1:ts_length]
        ape = Ep .- Eb
        ape0 = ape[1]
        ape .-= ape0
        ape ./= ape0 # (ape - ape₀) ./ ape₀
        Δt = diff(t)
        dₜape = diff(ape) ./ Δt
        ls = i % 2 == 0 ? linestyles[2] : linestyles[1]
        lines!(ax[1], 1:ts_length, ape, linestyle = ls, label = all_labels[i])
        if i < 3
            lines!(ax[2], 2:ts_length, dₜape, linestyle = ls)
            ax[2].ylabel = L"d$\mathcal{APE}$/dt (s)"
        end

    end
end
hidexdecorations!(ax[1], grid = false, ticks = false)
linkxaxes!(ax[1], ax[2])

Legend(fig[5, 1], ax_R_rho, "Experiment", orientation = :horizontal, nbanks = 2)
##
save("fig7.png", fig)

