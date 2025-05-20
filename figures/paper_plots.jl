using CairoMakie, StaircaseShenanigans
using JLD2, GibbsSeaWater, StatsBase #, ColorSchemes if needed
using SeawaterPolynomials: TEOS10EquationOfState, total_density, haline_contraction, thermal_expansion
using SpecialFunctions: erf

cd("figures")
## EOS's and other constants for the paper
ρ₀ = gsw_rho(34.7, 0.5, 0.5)
leos = CustomLinearEquationOfState(-0.5, 34.64, reference_density = ρ₀)
nleos = TEOS10EquationOfState(reference_density = ρ₀)
erf_tracer_solution(z, Cₗ::Number, ΔC::Number, κ::Number, t, interface_depth) =
    Cₗ + 0.5 * ΔC * (1 + erf((z - interface_depth) / sqrt(4 * κ * t)))

Sᵤ, Θᵤ = 34.58, -1.5 # cabbeling expt from project two
Sₗ, Θₗ = 34.7, 0.5
ΔS = Sᵤ - Sₗ
ΔΘ = Θᵤ - Θₗ
κₛ, κₜ = 1e-9, 1e-7
z = range(-0.5, 0, length = 1400) # range for density profile
leos_vec = fill(leos, length(z))
nleos_vec = fill(nleos, length(z))
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

S = erf_tracer_solution.(z, Sₗ, ΔS, κₛ, t, interface_depth)
T = erf_tracer_solution.(z, Θₗ, ΔΘ, κₜ, t, interface_depth)
τ = κₛ / κₜ
linear_σ₀ = total_density.(T, S, 0, leos_vec)
linear_σ₀′ = linear_σ₀ .- mean(linear_σ₀)
max_linear_σ₀′, max_linear_idx = findmax(linear_σ₀′)
min_linear_σ₀′, min_linear_idx = findmin(linear_σ₀′)
linear_σ₀′_norm = 2 * (linear_σ₀′ .- min_linear_σ₀′) ./ (max_linear_σ₀′ - min_linear_σ₀′) .- 1

# nlinear_σ₀ = total_density.(T, S, 0, nleos_vec)
nlinear_σ₀ = gsw_rho.(S, T, 0)
nlinear_σ₀′ = nlinear_σ₀ .- mean(nlinear_σ₀)
max_nlinear_σ₀′, max_nlinear_idx = findmax(nlinear_σ₀′)
min_nlinear_σ₀′, min_nlinear_idx = findmin(nlinear_σ₀′)
nlinear_σ₀′_norm = 2 * (nlinear_σ₀′ .- min_nlinear_σ₀′) ./ (max_nlinear_σ₀′ .- min_nlinear_σ₀′) .- 1

# S_minmax = [S[min_idx], S[max_idx]]
# T_minmax = [T[min_idx], T[max_idx]]

# Salinity and temperature panel
fig = Figure(size = (800, 600))
ax = [Axis(fig[1, i]) for i ∈ 1:2]
lines!(ax[1], S, z; color = (:blue, 0.5), label = "Salinity")
xlims!(ax[1], 34.57, 34.72)
ax[1].xlabel = "Salinity (gkg⁻¹)"
ax[1].xlabelcolor = :blue
ax[1].xticklabelcolor = :blue
ax[1].ylabel = "z (m)"
axT = Axis(fig[1, 1];
            xaxisposition = :top,
            xticklabelcolor = :red,
            xlabel = "Θ (°C)",
            xlabelcolor = :red,
            title = "Temperature and salinity profiles")
lines!(axT, T, z; color = (:red, 0.5), label = "Temperature")
axislegend(axT)
axislegend(ax[1], position = :lb)

# density panel
lines!(ax[2], linear_σ₀′_norm, z; color = (:black, 0.5), label = "Linear eos")
lines!(ax[2], nlinear_σ₀′_norm, z; color = (:green, 0.5), label = "Non linear eos")
ax[2].title = "Density profiles"
ax[2].xlabel = L"$\hat{σ}_{0}$ (kgm⁻³)"
# scatter!(ax[2], max_linear_σ₀′, z[max_linear_idx])
# scatter!(ax[2], min_linear_σ₀′, z[min_linear_idx])
# scatter!(ax[2], max_nlinear_σ₀′, z[max_linear_idx])
# scatter!(ax[2], min_nlinear_σ₀′, z[min_linear_idx])
hideydecorations!(ax[2], grid = false)
axislegend(ax[2])
linkyaxes!(ax[1], ax[2])
# colsize!(fig.layout, 1, Relative(3/5))
fig
save("S_T_sigma_profiles.png", fig)
## Figure
# Profiles in salinity-temperature space
fig = Figure(size = (700, 800))
eos_title = ["Linear eos", "Nonlinear eos"]
ax = [Axis(fig[i, 1], title = eos_title[i], xlabel = "S (gkg⁻¹)", ylabel = "Θ (°C)") for i ∈ 1:2]
N = 2000
S_range, Θ_range = range(minimum(S), maximum(S), length = N), range(minimum(T), maximum(T), length = N)
S_grid, Θ_grid = ones(N) .* S_range', ones(N)' .* Θ_range
# Linear
linear_σ₀_grid = total_density.(Θ_grid, S_grid, 0, fill(leos, size(S_grid)))
linear_σ₀_min = minimum(linear_σ₀)
linear_σ₀_lower = total_density(Θₗ, Sₗ, 0, leos)
linear_σ₀_upper = total_density(Θᵤ, Sᵤ, 0, leos)
linear_σ₀_max = maximum(linear_σ₀)

cn = contour!(ax[1], S_range, Θ_range, linear_σ₀_grid';
              levels = [linear_σ₀_min, linear_σ₀_upper, linear_σ₀_lower, linear_σ₀_max],
              colormap = :dense, label = "Isopycnals")

lines!(ax[1], S, T, linear_σ₀; label = "S-T profile", color = :grey, linestyle = :dash)
S_minmax = [S[min_linear_idx], S[max_linear_idx]]
T_minmax = [T[min_linear_idx], T[max_linear_idx]]

scatter!(ax[1], S_minmax[1], T_minmax[1]; markersize, label = "Minimum density", color = :orange)
scatter!(ax[1], S_minmax[2], T_minmax[2]; markersize, label = "Maximum density", color = :magenta)

hidexdecorations!(ax[1], grid = false, ticks = false)
# axislegend(ax[1], position = :lt)
# Nonlinear
# nlinear_σ₀_grid = total_density.(Θ_grid, S_grid, 0, fill(nleos, size(S_grid)))
# nlinear_σ₀_lower = total_density(Θₗ, Sₗ, 0, nleos)
# nlinear_σ₀_upper = total_density(Θᵤ, Sᵤ, 0, nleos)
nlinear_σ₀_grid = gsw_rho.(S_grid, Θ_grid, 0)
nlinear_σ₀_lower = gsw_rho(Sₗ, Θₗ, 0)
nlinear_σ₀_upper = gsw_rho(Sᵤ, Θᵤ, 0)
nlinear_σ₀_min = minimum(nlinear_σ₀)
nlinear_σ₀_max = maximum(nlinear_σ₀)

cn = contour!(ax[2], S_range, Θ_range, nlinear_σ₀_grid';
              levels = [nlinear_σ₀_min, nlinear_σ₀_upper, nlinear_σ₀_lower, nlinear_σ₀_max],
              colormap = :dense, label = "Isopycnals")
lines!(ax[2], S, T, nlinear_σ₀; label = "S-T profile", color = :grey, linestyle = :dash)
S_minmax = [S[min_nlinear_idx], S[max_nlinear_idx]]
T_minmax = [T[min_nlinear_idx], T[max_nlinear_idx]]

scatter!(ax[2], S_minmax[1], T_minmax[1]; markersize, label = "Minimum density", color = :orange)
scatter!(ax[2], S_minmax[2], T_minmax[2]; markersize, label = "Maximum density", color = :magenta)

Legend(fig[3, :], ax[2], orientation = :horizontal, nbanks = 2)
fig
save("S_T_sigma_ST_space.png", fig)
