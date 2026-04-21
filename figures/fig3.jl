include("plotting_env_and_load.jl")
## Figure three
# τ impact on asymmetry
Sᵤ_range = (34.51, 34.54, 34.58)
ΔS_range = Sᵤ_range .- Sₗ
κₛ_range = range(1e-9, 8e-8, length = 100)
τ_range = κₛ_range ./ κₜ
σ₀_nonlinear_max = Array{Float64}(undef, length(κₛ_range), length(Sᵤ_range))
σ₀_nonlinear_min = similar(σ₀_nonlinear_max)
σ₀_linear_max = similar(σ₀_nonlinear_max)
σ₀_linear_min = similar(σ₀_nonlinear_max)
Rᵨ_leos = Array{Float64}(undef, length(κₛ_range))
Rᵨ_nleos = similar(Rᵨ_leos)
σ₀ˡ_nleos = similar(σ₀_nonlinear_max)
σ₀ᵘ_nleos = similar(σ₀_nonlinear_max)
σ₀ˡ_leos = similar(σ₀_nonlinear_max)
σ₀ᵘ_leos = similar(σ₀_nonlinear_max)
t = 10000 # seconds, this value is not stricly important as maximum density does not increase
for (j, _Sᵤ) ∈ enumerate(Sᵤ_range)
    Rᵨ_leos[j] = compute_R_ρ([_Sᵤ, Sₗ], [Θᵤ, Θₗ], interface_depth, leos)
    Rᵨ_nleos[j] = compute_R_ρ([_Sᵤ, Sₗ], [Θᵤ, Θₗ], interface_depth, nleos)
    for i ∈ eachindex(κₛ_range)
        S = erf_tracer_solution.(z, Sₗ, ΔS_range[j], κₛ_range[i], t, interface_depth)
        T = erf_tracer_solution.(z, Θₗ, ΔΘ, κₜ, t, interface_depth)
        σ₀_nonlinear = gsw_rho.(S, T, 0)
        σ₀_nonlinear_max[i, j] = maximum(σ₀_nonlinear)
        σ₀_nonlinear_min[i, j] = minimum(σ₀_nonlinear)

        σ₀_linear = total_density.(T, S, 0, leos_vec)
        σ₀_linear_max[i, j] = maximum(σ₀_linear)
        σ₀_linear_min[i, j] = minimum(σ₀_linear)
        σ₀ˡ_nleos[i, j] = gsw_rho(Sₗ, Θₗ, 0)
        σ₀ᵘ_nleos[i, j] = gsw_rho(_Sᵤ, Θᵤ, 0)
        σ₀ˡ_leos[i, j] = total_density(Θₗ, Sₗ, 0, leos)
        σ₀ᵘ_leos[i, j] = total_density(Θᵤ, _Sᵤ, 0, leos)
    end
end
Δσ_lower_nonlinear = abs.(σ₀_nonlinear_max .- σ₀ˡ_nleos)
Δσ_upper_nonlinear = abs.(σ₀_nonlinear_min .- σ₀ᵘ_nleos)
Δσ_nonlinear = Δσ_upper_nonlinear ./ Δσ_lower_nonlinear

Δσ_lower_linear = abs.(σ₀_linear_max .- σ₀ˡ_leos)
Δσ_upper_linear = abs.(σ₀_linear_min .- σ₀ᵘ_leos)
Δσ_linear = Δσ_upper_linear ./ Δσ_lower_linear

find_τ_01 = findfirst(τ_range .> 0.1)
τ_range[find_τ_01]

fig = Figure(size = (800, 400))
linestyle = [:solid, :dash, :dot, :dashdot]
ax1 = Axis(fig[1, 1], title = "Asymmetry parameter variation", titlefont = :bold, xlabel = L"τ", ylabel = L"\delta_{\mathrm{mol}}")
for i ∈ eachindex(Sᵤ_range)
    lines!(ax1, τ_range, Δσ_linear[:, i], linestyle = linestyle[i], color = Makie.wong_colors()[1],
           label = L"$ρ_{\mathrm{linear}}\text{, }R_{\rho} =$ %$(round(Rᵨ_leos[i], digits = 2))")
end
for i ∈ eachindex(Sᵤ_range)
    lines!(ax1, τ_range, Δσ_nonlinear[:, i], linestyle = linestyle[i], color = Makie.wong_colors()[2],
           label = L"$ρ_{\mathrm{nonlinear}}\text{, }R_{\rho}$ = %$(round(Rᵨ_nleos[i], digits = 2))")
end
# axislegend(ax1, position = :rc)
Legend(fig[2, 1], ax1, orientation = :horizontal, nbanks = 3)
fig
##
save("fig3_tau_asymmetry.png", fig)
