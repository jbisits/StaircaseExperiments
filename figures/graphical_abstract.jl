include("plotting_env_and_load.jl")
## Graphical abstract
files = (l_R_ρ_105_dT2_diagnostics, nl_R_ρ_105_dT2_diagnostics)
         # l_R_ρ_105_dT1_diagnostics, nl_R_ρ_105_dT1_diagnostics,
         # l_R_ρ_105_dT05_diagnostics, nl_R_ρ_105_dT05_diagnostics)
fig = Figure(size = (800, 400))
axprofile = Axis(fig[1, 1], xlabel = L"$σ_{0}′$ (kgm$^{-3}$)", ylabel = "z (m)")
hidedecorations!(axprofile)
axσ = [Axis(fig[1, j], xlabel = "time (min)", ylabel = "z (m)") for j ∈ 2:3]
σ_colorrange_dT2 = jldopen(files[1]) do ds
        extrema(ds["σ_ha"][:, 3] .- mean(ds["σ_ha"][:, 1]))
end
colorrange = σ_colorrange_dT2
linestyles = (:solid, :dash)
ts_length = 240
for (i, f) ∈ enumerate(files)

    jldopen(f) do ds
        _z = ds["dims/z_aac"]
        σ_ha′ = ds["σ_ha"][:, 1:ts_length ]' .- mean(ds["σ_ha"][:, 1])'
        hm = heatmap!(axσ[i], 1:ts_length, _z, σ_ha′; colorrange, colormap = :diff)
        hidedecorations!(axσ[i])

        # profile
        t = ds["dims/time"]
        profile_time = 4
        ls = i % 2 == 0 ? linestyles[2] : linestyles[1]
        if i ∈ (1, 2)
            lines!(axprofile, σ_ha′[profile_time, :], _z, color = Makie.wong_colors()[i],
                    linestyle = ls, label = all_labels[i])
        end
    end
end
fig
##
save("graphical_abstract.png", fig)
