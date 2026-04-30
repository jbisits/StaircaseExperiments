include("plotting_env_and_load.jl")
## Figure four
# Density asymmetry heatmap
Sᵤ_range = range(33.54, Sₗ, length = 400)
Θᵤ_range = range(-1.6, Θₗ, length = 400)
Rᵨ_leos = Array{Float64}(undef, length(Θᵤ_range), length(Sᵤ_range))
Rᵨ_nleos = similar(Rᵨ_leos)
σ₀_nonlinear_max = similar(Rᵨ_leos)
σ₀_nonlinear_min = similar(Rᵨ_leos)
σ₀_linear_max = similar(Rᵨ_leos)
σ₀_linear_min = similar(Rᵨ_leos)
σ₀ᵘ_leos = similar(Rᵨ_leos)
σ₀ᵘ_nleos = similar(Rᵨ_leos)
σ₀ˡ_nleos = gsw_rho(Sₗ, Θₗ, 0)
σ₀ˡ_leos = total_density(Θₗ, Sₗ, 0, leos)
τ = 0.1
Δρ′_nlinear = similar(Rᵨ_leos)
δ(C) = 1 + 2*C + 2*sqrt(C^2 + C)
δ_leos = similar(Rᵨ_leos)
δ_nleos = similar(Rᵨ_leos)
for (j, _Sᵤ) ∈ enumerate(Sᵤ_range)

    _κₛ = τ * κₜ
    salinity = [_Sᵤ, Sₗ]
    _ΔS = _Sᵤ - Sₗ
    for (i, _Θᵤ) ∈ enumerate(Θᵤ_range)

        temperature = [_Θᵤ, Θₗ]
        _ΔΘ = _Θᵤ - Θₗ

        Rᵨ_leos[i, j] = compute_R_ρ(salinity, temperature, interface_depth, leos)
        Rᵨ_nleos[i, j] = compute_R_ρ(salinity, temperature, interface_depth, nleos)

        S = erf_tracer_solution.(z, Sₗ, _ΔS, _κₛ, t, interface_depth)
        T = erf_tracer_solution.(z, Θₗ, _ΔΘ, κₜ, t, interface_depth)
        σ₀_nonlinear = gsw_rho.(S, T, 0)
        σ₀_nonlinear_max[i, j] = maximum(σ₀_nonlinear)
        σ₀_nonlinear_min[i, j] = minimum(σ₀_nonlinear)

        σ₀_linear = total_density.(T, S, 0, leos_vec)
        σ₀_linear_max[i, j] = maximum(σ₀_linear)
        σ₀_linear_min[i, j] = minimum(σ₀_linear)

        σ₀ᵘ_nleos[i, j] = gsw_rho(_Sᵤ, _Θᵤ, 0)
        σ₀ᵘ_leos[i, j] = total_density(_Θᵤ, _Sᵤ, 0, leos)

        slope = _ΔΘ / _ΔS
        S_mix = range(salinity..., length = 1000)
        Θ_mix = Θₗ .+ slope .* (S_mix .- Sₗ)
        ρ_linear = total_density.(Θ_mix, S_mix, 0, fill(leos, length(S_mix)))
        ρ_l_max = maximum(ρ_linear)
        ρ_nlinear = gsw_rho.(S_mix, Θ_mix, 0)
        ρ_nl_max = maximum(ρ_nlinear)
        Δρ_linear = abs(ρ_linear[1] - ρ_linear[end])
        Δρ_nlinear = abs(ρ_nlinear[1] - ρ_nlinear[end])
        Δρ′_linear = ρ_l_max - ρ_linear[end]
        Δρ′_nlinear[i, j] = ρ_nl_max - ρ_nlinear[end]
        δ_leos[i, j] = δ(Δρ′_linear / Δρ_linear)
        δ_nleos[i, j] = δ(Δρ′_nlinear[i, j] / (ρ_nlinear[end] - ρ_nlinear[1]))
    end

end
Δσ_lower_nonlinear = abs.(σ₀_nonlinear_max .- σ₀ˡ_nleos)
Δσ_upper_nonlinear = abs.(σ₀_nonlinear_min .- σ₀ᵘ_nleos)
Δσ_nonlinear = Δσ_upper_nonlinear ./ Δσ_lower_nonlinear

Δσ_lower_linear = abs.(σ₀_linear_max .- σ₀ˡ_leos)
Δσ_upper_linear = abs.(σ₀_linear_min .- σ₀ᵘ_leos)
Δσ_linear = Δσ_upper_linear ./ Δσ_lower_linear

replace!(x -> x > 10 || x < 1 ? NaN : x, Rᵨ_nleos)
for (i, c) ∈ enumerate(eachcol(reverse(Rᵨ_nleos, dims = 1)))
    for j ∈ eachindex(c)
        Δσ_nonlinear[i, j] = isnan(Rᵨ_nleos[i, j]) ? NaN : Δσ_nonlinear[i, j]
        δ_nleos[i, j] = isnan(Rᵨ_nleos[i, j]) ? NaN : δ_nleos[i, j]
    end
end
replace!(x -> x > 10 || x < 1 ? NaN : x, Rᵨ_leos)
for (i, c) ∈ enumerate(eachcol(reverse(Rᵨ_leos, dims = 1)))
    for j ∈ eachindex(c)
        Δσ_linear[i, j] = isnan(Rᵨ_leos[i, j]) ? NaN : Δσ_linear[i, j]
        δ_leos[i, j] = isnan(Rᵨ_leos[i, j]) ? NaN : δ_leos[i, j]
    end
end

ΔΘ = Θᵤ_range .- Θₗ
ΔS = Sᵤ_range .- Sₗ

arctic_obs = [[-0.04], [-0.014]]
ΔΘ_expts = [-2, -1, -0.5]
ΔS_expts_linear = [-0.12, -0.06, -0.03]
ΔS_expts_nlinear = [-0.12, -0.069, -0.037]

# 0.95 asymmetry
find_dmol = findall(Δσ_nonlinear .> 0.95)
ΔΘ_asym = [ΔΘ[find_dmol[i][2]] for i ∈ eachindex(find_dmol)]

##
markersize = 15
linear_colour = :black
nlinear_colour = linear_colour
fig = Figure(size = (800, 1000))
ax_lRᵨ = Axis(fig[1, 1], xlabel = "ΔΘ (°C)", ylabel = "ΔS (gkg⁻¹)", title = "Linear eos")
hm = heatmap!(ax_lRᵨ, ΔΘ, ΔS, Rᵨ_leos, colormap = :amp)
scatter!(ax_lRᵨ, arctic_obs...; color = :blue, markersize)
scatter!(ax_lRᵨ, ΔΘ_expts, ΔS_expts_linear;
         color = linear_colour, markersize, marker = linear_expt_markers)
ax_nlRᵨ = Axis(fig[1, 2], xlabel = "ΔΘ (°C)", ylabel = "ΔS (gkg⁻¹)", title = "Nonlinear eos")
hm_nlRᵨ = heatmap!(ax_nlRᵨ, ΔΘ, ΔS, Rᵨ_nleos, colormap = :amp)
scatter!(ax_nlRᵨ, arctic_obs...; markersize, color = :blue)
scatter!(ax_nlRᵨ, ΔΘ_expts, ΔS_expts_nlinear;
         color = nlinear_colour, markersize, marker = nlinear_expt_markers)
hideydecorations!(ax_nlRᵨ, grid = false, ticks = false)
hidexdecorations!(ax_nlRᵨ, grid = false, ticks = false)
hidexdecorations!(ax_lRᵨ, grid = false, ticks = false)
Colorbar(fig[1, 3], hm_nlRᵨ, label = L"R_{\rho}")

colorrange = (minimum(Δσ_nonlinear[.!isnan.(Δσ_nonlinear)]), 1)
ax_lR_Δρ = Axis(fig[2, 1], xlabel = "ΔΘ (°C)", ylabel = "ΔS (gkg⁻¹)")
hm = heatmap!(ax_lR_Δρ, ΔΘ, ΔS, Δσ_linear; colorrange, colormap = :batlow)
scatter!(ax_lR_Δρ, arctic_obs...; color = :blue, markersize)
scatter!(ax_lR_Δρ, ΔΘ_expts, ΔS_expts_linear;
         color = linear_colour, markersize, marker = linear_expt_markers)
ax_nlR_Δρ = Axis(fig[2, 2], xlabel = "ΔΘ (°C)", ylabel = "ΔS (gkg⁻¹)")
hideydecorations!(ax_nlR_Δρ, grid = false, ticks = false)
hm_nlR_Δρ = heatmap!(ax_nlR_Δρ, ΔΘ, ΔS, Δσ_nonlinear; colorrange, colormap = :batlow)
scatter!(ax_nlR_Δρ, arctic_obs...; markersize, color = :blue)
scatter!(ax_nlR_Δρ, ΔΘ_expts, ΔS_expts_nlinear;
         color = nlinear_colour, markersize, marker = nlinear_expt_markers)
hidexdecorations!(ax_lR_Δρ, grid = false, ticks = false)
hidexdecorations!(ax_nlR_Δρ, grid = false, ticks = false)
Colorbar(fig[2, 3], hm_nlR_Δρ, label = L"\delta_{\mathrm{mol}}")

colorrange = (0, 1)
custom_colors = [:grey75, :tomato]
cat_cmap = cgrad(custom_colors, 2, categorical=true)
δ_leos_cat = ifelse.(δ_leos .== 1, 0, δ_leos)
δ_nleos_cat = ifelse.(δ_nleos .== 1, 0, δ_nleos)
ax_δ_linear = Axis(fig[3, 1], xlabel = "ΔΘ (°C)", ylabel = "ΔS (gkg⁻¹)")
# hm_δ_linear = heatmap!(ax_δ_linear, ΔΘ, ΔS, δ_leos; colorrange, colormap = :turbid, highclip = :red, lowclip = :orange)
hm_δ_linear = heatmap!(ax_δ_linear, ΔΘ, ΔS, δ_leos_cat; colorrange, colormap = cat_cmap)
scatter!(ax_δ_linear, arctic_obs...; markersize, color = :blue)
scatter!(ax_δ_linear, ΔΘ_expts, ΔS_expts_linear; color = linear_colour,
         label = linear_expt_labels, markersize, marker = linear_expt_markers)
ax_δ_nlinear = Axis(fig[3, 2], xlabel = "ΔΘ (°C)", ylabel = "ΔS (gkg⁻¹)")
# hm_δ_nlinear = heatmap!(ax_δ_nlinear, ΔΘ, ΔS, δ_nleos; colorrange, colormap = :turbid, highclip = :red, lowclip = :orange)
hm_δ_nlinear = heatmap!(ax_δ_nlinear, ΔΘ, ΔS, δ_nleos_cat; colorrange, colormap = cat_cmap)
scatter!(ax_δ_nlinear, arctic_obs...; markersize, color = :blue, label = "Arctic interfaces")
scatter!(ax_δ_nlinear, ΔΘ_expts, ΔS_expts_nlinear; color = nlinear_colour,
         label = nlinear_expt_labels, markersize, marker = nlinear_expt_markers)
hideydecorations!(ax_δ_nlinear, grid = false, ticks = false)
Colorbar(fig[3, 3], hm_δ_nlinear, ticks=([0.25, 0.75], [L"\delta_{\mathrm{turb}} = 1",L"\delta_{\mathrm{turb}} > 1"]), ticklabelrotation = π/2,
         ticksvisible = false)

# panel labels
panel_labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]
axs = [ax_lRᵨ, ax_lR_Δρ, ax_δ_linear, ax_nlRᵨ, ax_nlR_Δρ, ax_δ_nlinear]
for (i, a) ∈ enumerate(axs)
    text!(a, 1, 0,
          text = panel_labels[i],
          font = :bold,
          align = (:right, :bottom),
          offset = (-4, 2),
          space = :relative,
          fontsize = 24
          )
end
# Legend
legend_markers = [MarkerElement(color = :black, marker = m; markersize)
                  for m ∈ vcat(linear_expt_markers, nlinear_expt_markers)]
legend_expts = vcat(linear_expt_labels, nlinear_expt_labels)
Legend(fig[4, 1], legend_markers, legend_expts, "DNS Experiments",
       orientation = :horizontal, nbanks = 2)
lmarker_2 = [MarkerElement(color = :blue, marker = :circle; markersize)]
lobs = ["Timmermans et al. (2008)"]
Legend(fig[4, 2], lmarker_2, lobs, "Observations",
       orientation = :horizontal)
##
save("fig4.png", fig)
