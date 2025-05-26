using CairoMakie, StaircaseShenanigans
using StaircaseShenanigans: compute_R_ρ
using JLD2, GibbsSeaWater, StatsBase, ColorSchemes
using SeawaterPolynomials: TEOS10EquationOfState, total_density, haline_contraction, thermal_expansion
using SpecialFunctions: erf

cd("figures")
## EOS's and other constants for the paper
ρ₀ = gsw_rho(34.7, 0.5, 0.5)
leos = CustomLinearEquationOfState(-0.5, 34.64, reference_density = ρ₀)
nleos = TEOS10EquationOfState(reference_density = ρ₀)
erf_tracer_solution(z, Cₗ::Number, ΔC::Number, κ::Number, t, interface_depth) =
    Cₗ + 0.5 * ΔC * (1 + erf((z - interface_depth) / sqrt(4 * κ * t)))

Sᵤ, Θᵤ = 34.56, -1.5 # cabbeling expt from project two
Sₗ, Θₗ = 34.7, 0.5
ΔS = Sᵤ - Sₗ
ΔΘ = Θᵤ - Θₗ
κₛ, κₜ = 1e-9, 1e-7
Nz = 1400
z = range(-0.5, 0, length = Nz) # range for density profile
leos_vec = fill(leos, Nz)
nleos_vec = fill(nleos, Nz)
interface_depth = -0.25
t = 5000

## Load output
# when done..

## Figure theme
markersize = 10
publication_theme = Theme(font="CMU Serif", fontsize = 20,
                          Axis=(titlesize = 22,
                                xlabelsize = 20, ylabelsize = 20,
                                xgridstyle = :dash, ygridstyle = :dash,
                                xtickalign = 0, ytickalign = 0,
                                yticksize = 6.5, xticksize = 6.5),
                          Legend=(framecolor = (:black, 0.5),
                                  backgroundcolor = (:white, 0.5),
                                  labelsize = 20),
                          Colorbar=(ticksize=12,
                                    tickalign=1,
                                    spinewidth=0.5))
new_theme = merge(theme_latexfonts(), publication_theme)
set_theme!(new_theme)

## Figure
# Asymmetric effects of non-linear eos

S₀ᵘ, S₀ˡ = fill(Sᵤ, Int(Nz / 2)), fill(Sₗ, Int(Nz / 2))
Θ₀ᵘ, Θ₀ˡ = fill(Θᵤ, Int(Nz / 2)), fill(Θₗ, Int(Nz / 2))
linear_σ₀ᵘ = total_density.(Θ₀ᵘ, S₀ᵘ, 0, leos_vec[1:Int(Nz / 2)])
linear_σ₀ˡ = total_density.(Θ₀ˡ, S₀ˡ, 0, leos_vec[Int(Nz / 2)+1:end])
nlinear_σ₀ᵘ = gsw_rho.(S₀ᵘ, Θ₀ᵘ, 0)
nlinear_σ₀ˡ = gsw_rho.(S₀ˡ, Θ₀ˡ, 0)
S = erf_tracer_solution.(z, Sₗ, ΔS, κₛ, t, interface_depth)
T = erf_tracer_solution.(z, Θₗ, ΔΘ, κₜ, t, interface_depth)
τ = κₛ / κₜ
linear_σ₀ = total_density.(T, S, 0, leos_vec)
linear_σ₀′ = linear_σ₀ .- mean(vcat(linear_σ₀ᵘ, linear_σ₀ˡ))
max_linear_σ₀′, max_linear_idx = findmax(linear_σ₀′)
min_linear_σ₀′, min_linear_idx = findmin(linear_σ₀′)
# linear_σ₀′_norm = 2 * (linear_σ₀′ .- min_linear_σ₀′) ./ (max_linear_σ₀′ - min_linear_σ₀′) .- 1
linear_σ₀ᵘ′ = linear_σ₀ᵘ .- mean(vcat(linear_σ₀ᵘ, linear_σ₀ˡ))
linear_σ₀ˡ′ = linear_σ₀ˡ .- mean(vcat(linear_σ₀ᵘ, linear_σ₀ˡ))
# nlinear_σ₀ = total_density.(T, S, 0, nleos_vec)
nlinear_σ₀ = gsw_rho.(S, T, 0)
nlinear_σ₀′ = nlinear_σ₀ .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))
max_nlinear_σ₀′, max_nlinear_idx = findmax(nlinear_σ₀′)
min_nlinear_σ₀′, min_nlinear_idx = findmin(nlinear_σ₀′)
nlinear_σ₀ᵘ′ = nlinear_σ₀ᵘ .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))
nlinear_σ₀ˡ′ = nlinear_σ₀ˡ .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))
# nlinear_σ₀′_norm = 2 * (nlinear_σ₀′ .- min_nlinear_σ₀′) ./ (max_nlinear_σ₀′ .- min_nlinear_σ₀′) .- 1

# Salinity and temperature panel
fig = Figure(size = (800, 600))
ax = [Axis(fig[1, i]) for i ∈ 1:2]
lines!(ax[1], S₀ˡ, z[1:Int(Nz/2)]; color = (:blue, 0.5), label = "Initial salinity")
lines!(ax[1], S₀ᵘ, z[Int(Nz/2)+1:end]; color = (:blue, 0.5))
lines!(ax[1], S, z; color = :blue, linestyle = :dash, label = L"Salinity, $t = t^{*}$")
xlims!(ax[1], 34.54, 34.72)
ax[1].xlabel = "Salinity (gkg⁻¹)"
ax[1].xlabelcolor = :blue
ax[1].xticklabelcolor = :blue
ax[1].ylabel = "z (m)"
axT = Axis(fig[1, 1];
            xaxisposition = :top,
            xticklabelcolor = :red,
            xlabel = "Θ (°C)",
            xlabelcolor = :red,
            title = "(a) Temperature and salinity profiles")
lines!(axT, Θ₀ˡ, z[1:Int(Nz/2)]; color = (:red, 0.5), label = "Initial temperature")
lines!(axT, Θ₀ᵘ, z[Int(Nz/2)+1:end]; color = (:red, 0.5))
lines!(axT, T, z; color = :red, linestyle = :dash, label = L"Temperature, $t = t^{*}$")
axislegend(axT)
axislegend(ax[1], position = :lb)

# density panel
lines!(ax[2], linear_σ₀ˡ′, z[1:Int(Nz/2)]; color = (:black, 0.5), label = "Initial density")
lines!(ax[2], linear_σ₀ᵘ′, z[Int(Nz/2)+1:end]; color = (:black, 0.5))
lines!(ax[2], linear_σ₀′, z; linestyle = :dash, label = L"$\rho_{\mathrm{linear}}$, $t = t^{*}$")
# lines!(ax[2], nlinear_σ₀ˡ′, z[1:Int(Nz/2)]; color = (:green, 0.5), label = L"$\rho_{\mathrm{nonlinear}}, t = 83$")
# lines!(ax[2], nlinear_σ₀ᵘ′, z[Int(Nz/2)+1:end]; color = (:green, 0.5))
lines!(ax[2], nlinear_σ₀′, z; linestyle = :dash, label = L"$\rho_{\mathrm{nonlinear}}$, $t = t^{*}$")
ax[2].title = "(b) Density profiles"
ax[2].xlabel = L"$σ_{0}'$ (kgm⁻³)"
# scatter!(ax[2], max_linear_σ₀′, z[max_linear_idx])
# scatter!(ax[2], min_linear_σ₀′, z[min_linear_idx])
# scatter!(ax[2], max_nlinear_σ₀′, z[max_linear_idx])
# scatter!(ax[2], min_nlinear_σ₀′, z[min_linear_idx])
hideydecorations!(ax[2], grid = false)
axislegend(ax[2])
linkyaxes!(ax[1], ax[2])
# colsize!(fig.layout, 1, Relative(3/5))
fig
##
save("S_T_sigma_profiles.png", fig)

## Figure
# Profiles in salinity-temperature space
σ_grad = get(ColorSchemes.dense, range(0.2, 1, length = 4))
fig = Figure(size = (700, 600))
eos_title = ["(a) Linear eos", "(b) Nonlinear eos"]
ax = [Axis(fig[i, 1], title = eos_title[i], xlabel = "S (gkg⁻¹)", ylabel = "Θ (°C)") for i ∈ 1:2]
N = 2000
S_range, Θ_range = range(minimum(S), maximum(S), length = N), range(minimum(T), maximum(T), length = N)
S_grid, Θ_grid = ones(N) .* S_range', ones(N)' .* Θ_range
# Linear
linear_σ₀_grid = total_density.(Θ_grid, S_grid, 0, fill(leos, size(S_grid)))  .- mean(vcat(linear_σ₀ᵘ, linear_σ₀ˡ))
linear_σ₀_min = minimum(linear_σ₀) .- mean(vcat(linear_σ₀ᵘ, linear_σ₀ˡ))
linear_σ₀_lower = total_density(Θₗ, Sₗ, 0, leos) .- mean(vcat(linear_σ₀ᵘ, linear_σ₀ˡ))
linear_σ₀_upper = total_density(Θᵤ, Sᵤ, 0, leos) .- mean(vcat(linear_σ₀ᵘ, linear_σ₀ˡ))
linear_σ₀_max = maximum(linear_σ₀) .- mean(vcat(linear_σ₀ᵘ, linear_σ₀ˡ))

cn = contour!(ax[1], S_range, Θ_range, linear_σ₀_grid';
              levels = [linear_σ₀_min, linear_σ₀_upper, linear_σ₀_lower, linear_σ₀_max],
            #   colormap = :dense,
              color = σ_grad,
              labels = true,
              linewidth = 0.5,
              label = "Isopycnals",
              labelsize = 16)

lines!(ax[1], S, T, linear_σ₀; label = "S-T profile", color = :tomato, linestyle = :dash)
S_minmax = [S[min_linear_idx], S[max_linear_idx]]
T_minmax = [T[min_linear_idx], T[max_linear_idx]]

scatter!(ax[1], S[1], T[1]; markersize, label = "Initial lower density", color = σ_grad[3])
scatter!(ax[1], S[end], T[end]; markersize, label = "Initial upper density", color = σ_grad[2])
scatter!(ax[1], S_minmax[1], T_minmax[1]; markersize, label = "Minimum density", color = σ_grad[1])
scatter!(ax[1], S_minmax[2], T_minmax[2]; markersize, label = "Maximum density", color = σ_grad[end])

hidexdecorations!(ax[1], grid = false, ticks = false)
# axislegend(ax[1], position = :lt)
# Nonlinear
# nlinear_σ₀_grid = total_density.(Θ_grid, S_grid, 0, fill(nleos, size(S_grid)))
# nlinear_σ₀_lower = total_density(Θₗ, Sₗ, 0, nleos)
# nlinear_σ₀_upper = total_density(Θᵤ, Sᵤ, 0, nleos)
nlinear_σ₀_grid = gsw_rho.(S_grid, Θ_grid, 0) .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))
nlinear_σ₀_lower = gsw_rho(Sₗ, Θₗ, 0) .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))
nlinear_σ₀_upper = gsw_rho(Sᵤ, Θᵤ, 0) .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))
nlinear_σ₀_min = minimum(nlinear_σ₀) .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))
nlinear_σ₀_max = maximum(nlinear_σ₀) .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))

lines!(ax[2], S, T, nlinear_σ₀; label = "S-T profile", color = :tomato, linestyle = :dash)
S_minmax = [S[min_nlinear_idx], S[max_nlinear_idx]]
T_minmax = [T[min_nlinear_idx], T[max_nlinear_idx]]

cn = contour!(ax[2], S_range, Θ_range, nlinear_σ₀_grid';
              levels = [nlinear_σ₀_min, nlinear_σ₀_upper, nlinear_σ₀_lower, nlinear_σ₀_max],
            #   colormap = :dense,
              labels = true,
              color = σ_grad,
              linewidth = 0.5,
              label = "Isopycnals",
              labelsize = 16)
scatter!(ax[2], S_minmax[1], T_minmax[1]; markersize, label = "Minimum density", color = σ_grad[1])
scatter!(ax[2], S[end], T[end]; markersize, label = "Initial upper density", color = σ_grad[2])
scatter!(ax[2], S[1], T[1]; markersize, label = "Initial lower density", color = σ_grad[3])
scatter!(ax[2], S_minmax[2], T_minmax[2]; markersize, label = "Maximum density", color = σ_grad[end])

Legend(fig[3, :], ax[2], orientation = :horizontal, nbanks = 2)
fig
##
save("S_T_sigma_ST_space_2panel.png", fig)

## Figure
# Density asymmetry
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

# Rᵨ_leos = compute_R_ρ([Sᵤ, Sₗ], temperature, interface_depth, leos)
# Rᵨ_nleos = compute_R_ρ([Sᵤ, Sₗ], temperature, interface_depth, nleos)
fig = Figure(size = (900, 600))
linestyle = [:solid, :dash, :dot, :dashdot]
ax1 = Axis(fig[1, 1], title = L"(a) Asymmetry due to $\tau$", titlefont = :bold, xlabel = L"τ", ylabel = L"R_{\Delta\rho}")
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

# Figure
# density ratio vs Δρ in each layer for different τ
τ = (0.01, 0.05, 0.1)
Sᵤ_range = range(33.52, 34.59, length = 100)
temperature = [Θᵤ, Θₗ]
Rᵨ_leos = Array{Float64}(undef, length(Sᵤ_range), length(τ))
Rᵨ_nleos = similar(Rᵨ_leos)
σ₀_nonlinear_max = similar(Rᵨ_leos)
σ₀_nonlinear_min = similar(Rᵨ_leos)
σ₀_linear_max = similar(Rᵨ_leos)
σ₀_linear_min = similar(Rᵨ_leos)
σ₀ᵘ_leos = similar(Rᵨ_leos)
σ₀ᵘ_nleos = similar(Rᵨ_leos)
σ₀ˡ_nleos = gsw_rho(Sₗ, Θₗ, 0)
σ₀ˡ_leos = total_density(Θₗ, Sₗ, 0, leos)
for j ∈ eachindex(τ)

    _κₛ = τ[j] * κₜ
    for (i, _Sᵤ) ∈ enumerate(Sᵤ_range)

        salinity = [_Sᵤ, Sₗ]
        _ΔS = _Sᵤ - Sₗ

        Rᵨ_leos[i, j] = compute_R_ρ(salinity, temperature, interface_depth, leos)
        Rᵨ_nleos[i, j] = compute_R_ρ(salinity, temperature, interface_depth, nleos)

        S = erf_tracer_solution.(z, Sₗ, _ΔS, _κₛ, t, interface_depth)
        T = erf_tracer_solution.(z, Θₗ, ΔΘ, κₜ, t, interface_depth)
        σ₀_nonlinear = gsw_rho.(S, T, 0)
        σ₀_nonlinear_max[i, j] = maximum(σ₀_nonlinear)
        σ₀_nonlinear_min[i, j] = minimum(σ₀_nonlinear)

        σ₀_linear = total_density.(T, S, 0, leos_vec)
        σ₀_linear_max[i, j] = maximum(σ₀_linear)
        σ₀_linear_min[i, j] = minimum(σ₀_linear)

        σ₀ᵘ_nleos[i, j] = gsw_rho(_Sᵤ, Θᵤ, 0)
        σ₀ᵘ_leos[i, j] = total_density(Θᵤ, _Sᵤ, 0, leos)
    end

end
Δσ_lower_nonlinear = abs.(σ₀_nonlinear_max .- σ₀ˡ_nleos)
Δσ_upper_nonlinear = abs.(σ₀_nonlinear_min .- σ₀ᵘ_nleos)
Δσ_nonlinear = Δσ_upper_nonlinear ./ Δσ_lower_nonlinear

Δσ_lower_linear = abs.(σ₀_linear_max .- σ₀ˡ_leos)
Δσ_upper_linear = abs.(σ₀_linear_min .- σ₀ᵘ_leos)
Δσ_linear = Δσ_upper_linear ./ Δσ_lower_linear

Rᵨ_cab = compute_R_ρ([34.551, Sₗ], temperature, interface_depth, nleos)
# fig = Figure(size = (500, 500))
ax2 = Axis(fig[1, 2], title = L"(b) Asymmetry due to $R_{\rho}$", xlabel = L"R_{\rho}", ylabel =  L"R_{\Delta\rho}")
linestyle = [:solid, :dash, :dot, :dashdot]
for i ∈ eachindex(τ)
    lines!(ax2, Rᵨ_leos[:, i], Δσ_linear[:, i]; color = Makie.wong_colors()[1], linestyle = linestyle[i], label = L"$ρ_{\mathrm{linear}}\text{, }\tau =$ %$(round((τ[i]), digits = 2))")
end
for i ∈ eachindex(τ)
    lines!(ax2, Rᵨ_nleos[:, i], Δσ_nonlinear[:, i]; color = Makie.wong_colors()[2], linestyle = linestyle[i], label = L"$ρ_{\mathrm{nonlinear}}\text{, }\tau =$ %$(round((τ[i]), digits = 2))")
end
# vlines!(ax2, Rᵨ_cab, label = "Rᵨ_cab", linestyle = :dash, color = :red)
# vlines!(ax2, 1.22, linestyle = :dash)
# vlines!(ax2, 1.23, linestyle = :dash)
linkyaxes!(ax1, ax2)
hideydecorations!(ax2, grid = false, ticks = false)
# axislegend(ax2, position = :rb, orientation = :horizontal, nbanks = 3)
Legend(fig[2, 2], ax2, orientation = :horizontal, nbanks = 3)
fig
##
save("density_asymmetry.png", fig)
