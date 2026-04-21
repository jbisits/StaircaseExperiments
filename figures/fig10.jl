include("plotting_env_and_load.jl")
## Figure ten
# initial evlotuion
nl_files = (nl_R_ρ_105_dT2_diagnostics, nl_R_ρ_105_dT1_diagnostics, nl_R_ρ_105_dT05_diagnostics)
l_files = (l_R_ρ_105_dT2_diagnostics, l_R_ρ_105_dT1_diagnostics, l_R_ρ_105_dT05_diagnostics)
merged_files = (l_files, nl_files)
eos_type = (leos, nleos)
ordered_labels = (("I", "II", "III"), ("IV", "V", "VI"))
fig = Figure(size = (900, 900))
for (j, files) ∈ enumerate(merged_files)
    ax = [Axis(fig[j, i],
                xlabel = L"$S$ (gkg$^{-1}$)",
                ylabel = L"$Θ$ (°C)",
                title = "Experiment " * ordered_labels[j][i]) for i ∈ eachindex(files)]
    file = jldopen(files[1])
    S = file["S_ha"][:, 4]
    Θ = file["T_ha"][:, 4]
    N = 100
    S_range = range(extrema(S)..., length=N)
    T_range = range(extrema(Θ)..., length=N)
    eos_vec = fill(eos_type[j], (N, N))
    ρ_grid = total_density.(T_range' .* ones(length(T_range)),
                            S_range .* ones(length(S_range))',
                            fill(0, (N, N)), eos_vec)
    close(file)
    for (i, f) ∈ enumerate(files)
        file = jldopen(f)
        S = file["S_ha"][:, 4]
        Θ = file["T_ha"][:, 4]
        σ = file["σ_ha"][:, 4]
        z = file["dims/z_aac"]
        t = file["dims/time"]
        κₛ = file["attrib/κₛ (m²s⁻¹)"]
        κₜ = file["attrib/κₜ (m²s⁻¹)"]
        id = file["attrib/interface_depth"]
        close(file)

        ΔS = S[end] - S[1]
        ΔΘ = Θ[end] - Θ[1]
        _t = 5000
        S_model = erf_tracer_solution.(z, S[1], ΔS, κₛ, _t, id)
        Θ_model = erf_tracer_solution.(z, Θ[1], ΔΘ, κₜ, _t, id)
        σ_model = total_density.(Θ_model, S_model, fill(0, length(S_model)), fill(eos_type[j], length(S_model)))
        ρmodel_shallow = σ_model[end]
        ρmodel_deep = σ_model[1]
        ρmodel_max = maximum(σ_model)
        ρmodel_min = minimum(σ_model)
        R_Δρ_model = round(abs(ρmodel_min - ρmodel_shallow) / abs(ρmodel_max - ρmodel_deep), digits = 3)

        ρ_shallow = σ[end]
        ρ_deep = σ[1]
        ρ_max = maximum(σ)
        ρ_min = minimum(σ)
        R_Δρ_sim = round(abs(ρ_min - ρ_shallow) / abs(ρ_max - ρ_deep), digits = 3)

        contour!(ax[i], S_range, T_range, ρ_grid, levels = [ρmodel_shallow, ρmodel_deep, ρmodel_max, ρmodel_min],
                color = Makie.wong_colors()[3], label = "Isopycnals", linestyle = :dot)
        lines!(ax[i], S_model, Θ_model, label = "1D model", linewidth = 3)
        lines!(ax[i], S, Θ; label = "Simulation output", color = Makie.wong_colors()[2],
                linestyle = :dash, linewidth = 3)
        text!(ax[i], 34.7, -1.2, text = L"1D model $\delta_{\mathrm{mol}} =$ %$(R_Δρ_model)", align = (:right, :top))
        text!(ax[i], 34.7, -1.4, text = L"Simulation $\delta_{\mathrm{mol}} =$ %$(R_Δρ_sim)", align = (:right, :top))
        ax[i].subtitle = L"$t~=$ %$(round(t[4]/60, digits = 1)) min"
        if i > 1
                hideydecorations!(ax[i], grid = false, ticks = false)
                linkyaxes!(ax[i], ax[1])
                linkxaxes!(ax[i], ax[1])
        end
        if j == 1
            hidexdecorations!(ax[i], grid = false, ticks = false)
        end
    end
end
Legend(fig[3, :], ax[1], orientation = :horizontal)
fig
##
save("fig10_ST_simluation.png", fig)
