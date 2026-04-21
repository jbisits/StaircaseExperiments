include("plotting_env_and_load.jl")
## Figure six
# Density hovmollers for Rρ = 1.05
files = (l_R_ρ_105_dT2_diagnostics, nl_R_ρ_105_dT2_diagnostics,
         l_R_ρ_105_dT1_diagnostics, nl_R_ρ_105_dT1_diagnostics,
         l_R_ρ_105_dT05_diagnostics, nl_R_ρ_105_dT05_diagnostics)
fig = Figure(size = (800, 1000))
axprofile = [Axis(fig[j, 1], xlabel = L"$σ_{0}′$ (kgm$^{-3}$)", ylabel = "z (m)") for j ∈ 1:3]
axσ = [Axis(fig[j, i], xlabel = "time (min)", ylabel = "z (m)") for i ∈ 2:3, j ∈ 1:3]
σ_colorrange_dT2 = jldopen(files[1]) do ds
        extrema(ds["σ_ha"][:, 3] .- mean(ds["σ_ha"][:, 1]))
end
σ_colorrange_dT1 = jldopen(files[3]) do ds
        extrema(ds["σ_ha"][:, 3] .- mean(ds["σ_ha"][:, 1]))
end
σ_colorrange_dT05 = jldopen(files[5]) do ds
        extrema(ds["σ_ha"][:, 3] .- mean(ds["σ_ha"][:, 1]))
end
colorranges = (σ_colorrange_dT2, σ_colorrange_dT1, σ_colorrange_dT05)
linestyles = (:solid, :dash)
ts_length = 240
for (i, f) ∈ enumerate(files)

    jldopen(f) do ds
        crange = if i ∈ 1:2
                    colorranges[1]
                elseif i ∈ 3:4
                    colorranges[2]
                elseif i ∈ 5:6
                    colorranges[3]
                end
        _z = ds["dims/z_aac"]
        σ_ha′ = ds["σ_ha"][:, 1:ts_length ]' .- mean(ds["σ_ha"][:, 1])'
        # if i == 5
        #     σ_ha′ .+= 0.001
        # end
        hm = heatmap!(axσ[i], 1:ts_length, _z, σ_ha′; colorrange = crange, colormap = :diff)
        axσ[i].title = "Experiment " * all_labels[i]
        hideydecorations!(axσ[i], ticks = false)
        if i % 2 == 0
            Colorbar(fig[Int(i / 2), 4], hm, label = L"$σ_{0}′$ (kgm$^{-3}$)")
        end
        if i < 5
            hidexdecorations!(axσ[i], ticks = false)
        end

        # profile
        t = ds["dims/time"]
        profile_time = 4
        ls = i % 2 == 0 ? linestyles[2] : linestyles[1]
        if i ∈ (1, 2)
            lines!(axprofile[1], σ_ha′[profile_time, :], _z, color = Makie.wong_colors()[i],
                    linestyle = ls, label = all_labels[i])
            if i % 2 == 0
                axprofile[1].title = "t = $(round(Int, t[profile_time]/60))min"
            end
            hidexdecorations!(axprofile[1], ticks = false, grid = false)
        elseif i ∈ (3, 4)
            lines!(axprofile[2], σ_ha′[profile_time, :], _z, color = Makie.wong_colors()[i],
                    linestyle = ls, label = all_labels[i])
            if i % 2 == 0
                axprofile[2].title = "t = $(round(Int, t[profile_time]/60))min"
            end
            hidexdecorations!(axprofile[2], ticks = false, grid = false)
        elseif i ∈ (5, 6)
            lines!(axprofile[3], σ_ha′[profile_time, :], _z, color = Makie.wong_colors()[i],
                    linestyle = ls, label = all_labels[i])
            if i % 2 == 0
                # xlims!(axprofile[3], (-0.02, 0.02))
                axprofile[3].title = "t = $(round(Int, t[profile_time]/60))min"
                # axprofile[3].xticklabelrotation = π / 4
            end
        end
    end
end
linkxaxes!(axprofile[2], axprofile[1])
linkxaxes!(axprofile[3], axprofile[1])
for i ∈ 1:3
    axislegend(axprofile[i])
end
fig
##
save("fig6_density_hovs.png", fig)
