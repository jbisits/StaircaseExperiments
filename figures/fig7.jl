include("plotting_env_and_load.jl")
## Figure seven
# energetics
files = (l_R_ρ_105_dT2_diagnostics, nl_R_ρ_105_dT2_diagnostics,
         l_R_ρ_105_dT1_diagnostics, nl_R_ρ_105_dT1_diagnostics,
         l_R_ρ_105_dT05_diagnostics, nl_R_ρ_105_dT05_diagnostics)
fig = Figure(size = (800, 600))
ax = [Axis(fig[i, 1], xlabel = "time (min)", ylabel = "Non-dimensional energy") for i ∈ 1:2]
ax[1].title = "(a) Available potential energy"
ax[2].title = "(b) Time derivateive of APE, I and IV"
ts_length = 240
j = 1
linestyles = (:solid, :dash)
for (i, f) ∈ enumerate(files)

    # j = i ∈ (1, 2) ? 1 : 2
    jldopen(f) do ds

        t = ds["dims"]["time"][1:ts_length]
        Δt = diff(t)
        Ep0 = ds["Ep"][1]
        Epfinal = ds["Eb"][ts_length]
        Ep = ds["Ep"][1:ts_length] .- Ep0
        Eb = ds["Eb"][1:ts_length] .- Ep0
        dₜEp = diff(Ep) ./ Δt
        dₜEb = diff(Eb) ./ Δt
        ape = Ep .- Eb
        dₜape = diff(ape) ./ Δt
        ek = diff(ds["∫Eₖ"][1:ts_length]) ./ Δt
        # lines!(ax, 2:120, Ep, label = "Ep, "*labels[i])
        # lines!(ax, 2:120, Eb, label = "Eb, "*labels[i])
        ls = i % 2 == 0 ? linestyles[2] : linestyles[1]
        lines!(ax[1], 1:ts_length, ape, linestyle = ls, label = all_labels[i])
        if i < 3
            lines!(ax[2], 2:ts_length, dₜape, linestyle = ls)
            ax[2].ylabel = "Non-dimensional rate"
        end

    end
end
hidexdecorations!(ax[1], grid = false, ticks = false)
hidexdecorations!(ax[2], grid = false, ticks = false)
linkxaxes!(ax[1], ax[2])
axislegend(ax[1], orientation = :horizontal, nbanks = 2, position = :rb)

# density ratio and salinity midpoint
files = (l_R_ρ_105_dT2_diagnostics, nl_R_ρ_105_dT2_diagnostics,
         l_R_ρ_105_dT1_diagnostics, nl_R_ρ_105_dT1_diagnostics,
         l_R_ρ_105_dT05_diagnostics, nl_R_ρ_105_dT05_diagnostics)
ax_R_rho = Axis(fig[3, 1],
                title = "(c) Density ratio",
                xlabel = "time ( min)",
                ylabel = L"R_{\rho}")
ax_interface = Axis(fig[4, 1],
                    title = "(d) Salinity interface height",
                    xlabel = "time (min)",
                    ylabel = L"$z*$ (m)")
ts_range = 1:240
for (i, file) ∈ enumerate(files)
    jldopen(file) do f
        t = f["dims"]["time"][ts_range]
        z✶ = f["dims"]["z✶"][f["S_interface_idx"]][ts_range]
        Rᵨ = f["R_ρ"][ts_range]
        ls = i % 2 == 0 ? :dash : :solid
        lines!(ax_interface, t ./ 60, z✶, label = all_labels[i], linestyle = ls)
        lines!(ax_R_rho, t ./ 60, Rᵨ, label = all_labels[i], linestyle = ls)
    end
end
hidexdecorations!(ax_R_rho, grid = false, ticks = false)
linkxaxes!(ax_R_rho, ax_interface)
axislegend(ax_R_rho, position = :lt, orientation = :horizontal, nbanks = 2)
##
save("fig7_timeseries.png", fig)

