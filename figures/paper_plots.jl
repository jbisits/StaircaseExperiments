using StaircaseShenanigans
using StaircaseShenanigans: compute_R_ρ
using JLD2, GibbsSeaWater, StatsBase, ColorSchemes
using SeawaterPolynomials: TEOS10EquationOfState, total_density, haline_contraction, thermal_expansion
using Printf
using SpecialFunctions: erf
using CairoMakie # GLMakie better for surface plots

cd("figures")
## EOS's and other constants for the paper
ρ₀ = gsw_rho(34.7, 0.5, 0)
leos = CustomLinearEquationOfState(-0.5, 34.64, reference_density = ρ₀)
leos_func(S, Θ) = CustomLinearEquationOfState(Θ, S, reference_density = ρ₀)
nleos = TEOS10EquationOfState(reference_density = ρ₀)
erf_tracer_solution(z, Cₗ::Number, ΔC::Number, κ::Number, t, interface_depth) =
    Cₗ + 0.5 * ΔC * (1 + erf((z - interface_depth) / sqrt(4 * κ * t)))

Sᵤ, Θᵤ = 34.58, -1.5 # cabbeling expt from project two
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
# rundown experiments
rundown_path = joinpath(@__DIR__, "../single_interface_rundown/rundown_1.05/step")
nl_R_ρ_105_dT2_diagnostics  = joinpath(rundown_path, "dns_res_dT2_nonlineareos/step_diagnostics.jld2")
l_R_ρ_105_dT2_diagnostics   = joinpath(rundown_path, "dns_res_dT2_lineareos/step_diagnostics.jld2")
nl_R_ρ_105_dT1_diagnostics  = joinpath(rundown_path, "dns_res_dT1_nonlineareos/step_diagnostics.jld2")
l_R_ρ_105_dT1_diagnostics   = joinpath(rundown_path, "dns_res_dT1_lineareos/step_diagnostics.jld2")
nl_R_ρ_105_dT05_diagnostics = joinpath(rundown_path, "dns_res_dT05_nonlineareos/step_diagnostics.jld2")
l_R_ρ_105_dT05_diagnostics  = joinpath(rundown_path, "dns_res_dT05_lineareos/step_diagnostics.jld2")
# flux bc experiments
fluxbc_path = joinpath(@__DIR__, "../single_interface_fluxbcs/R_rho_1.05/deltatheta_1/")
l_fbc_diagnostics  = joinpath(fluxbc_path, "alt_lineareos/step_diagnostics.jld2")
nl_fbc_diagnostics = joinpath(fluxbc_path, "longer_nonlineareos/step_diagnostics.jld2")
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
# labels
linear_expt_labels = ["I", "II", "III"]
nlinear_expt_labels = ["IV", "V", "VI"]
all_labels = collect(Iterators.flatten(zip(linear_expt_labels, nlinear_expt_labels)))
linear_expt_markers = [:utriangle, :dtriangle, :rtriangle]
nlinear_expt_markers = [:cross, :xcross, :star5]

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
linear_σ₀ᵘ′ = linear_σ₀ᵘ .- mean(vcat(linear_σ₀ᵘ, linear_σ₀ˡ))
linear_σ₀ˡ′ = linear_σ₀ˡ .- mean(vcat(linear_σ₀ᵘ, linear_σ₀ˡ))
nlinear_σ₀ = gsw_rho.(S, T, 0)
nlinear_σ₀′ = nlinear_σ₀ .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))
max_nlinear_σ₀′, max_nlinear_idx = findmax(nlinear_σ₀′)
min_nlinear_σ₀′, min_nlinear_idx = findmin(nlinear_σ₀′)
nlinear_σ₀ᵘ′ = nlinear_σ₀ᵘ .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))
nlinear_σ₀ˡ′ = nlinear_σ₀ˡ .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))

# Salinity and temperature panel
fig = Figure(size = (800, 600))
ax = [Axis(fig[1, i]) for i ∈ 1:2]
lines!(ax[1], S₀ˡ, z[1:Int(Nz/2)]; color = (:blue, 0.5), label = "Initial salinity")
lines!(ax[1], S₀ᵘ, z[Int(Nz/2)+1:end]; color = (:blue, 0.5))
lines!(ax[1], S, z; color = :blue, linestyle = :dash, label = L"Salinity, $t = t^{*}$")
xlims!(ax[1], 34.56, 34.72)
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
lines!(ax[2], nlinear_σ₀′, z; linestyle = :dash, label = L"$\rho_{\mathrm{nonlinear}}$, $t = t^{*}$")
ax[2].title = "(b) Density profiles"
ax[2].xlabel = L"$σ_{0}'$ (kgm⁻³)"
hideydecorations!(ax[2], grid = false)
axislegend(ax[2])
linkyaxes!(ax[1], ax[2])
fig
##
save("fig1_S_T_sigma_profiles.png", fig)

## Figure
# Profiles in salinity-temperature space
# σ_grad = get(ColorSchemes.dense, range(0.2, 1, length = 4))
σ_grad = get(ColorSchemes.diff, [0.2, 0.3, 0.6, 1])
σ_grad = cgrad(σ_grad, 4, categorical = true)
fig = Figure(size = (800, 600))
eos_title = ["(a) Linear equation of state", "(b) Nonlinear equation of state"]
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

function my_contour_label_formatter(level::Real)::String
    lev_short = round(level; digits = 3)
    x = @sprintf("%.1e", lev_short)
    string(x)
end

σ_linear_levels = [linear_σ₀_min, linear_σ₀_upper, linear_σ₀_lower, linear_σ₀_max]
cn = contour!(ax[1], S_range, Θ_range, linear_σ₀_grid';
              levels = σ_linear_levels,
              colormap = σ_grad,
              colorrange = extrema(σ_linear_levels),
            #   color = σ_grad,
            #   labelformatter = my_contour_label_formatter,
            #   labels = true,
            #   labelsize = 14,
              linewidth = 0.75,
            #   linestyle = :dot,
              label = "Isopycnals")
Colorbar(fig[1, 2], colorrange = extrema(σ_linear_levels),
        colormap = cgrad(σ_grad, 4, categorical = true),
        ticks = (range(extrema(σ_linear_levels)..., length = 4), string.(round.(σ_linear_levels, digits = 3))),
        label = L"$\Delta \sigma_{0}'$ (kgm$^{-3}$)")
lines!(ax[1], S, T, linear_σ₀; label = "S-T profile", color = :tomato, linestyle = :dash)
S_minmax = [S[min_linear_idx], S[max_linear_idx]]
T_minmax = [T[min_linear_idx], T[max_linear_idx]]

scatter!(ax[1], S[1], T[1]; markersize, label = "Initial lower density", color = σ_grad[3])
scatter!(ax[1], S[end], T[end]; markersize, label = "Initial upper density", color = σ_grad[2])
scatter!(ax[1], S_minmax[1], T_minmax[1]; markersize, label = "Minimum density", color = σ_grad[1])
scatter!(ax[1], S_minmax[2], T_minmax[2]; markersize, label = "Maximum density", color = σ_grad[end])

hidexdecorations!(ax[1], grid = false, ticks = false)
nlinear_σ₀_grid = gsw_rho.(S_grid, Θ_grid, 0) .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))
nlinear_σ₀_lower = gsw_rho(Sₗ, Θₗ, 0) .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))
nlinear_σ₀_upper = gsw_rho(Sᵤ, Θᵤ, 0) .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))
nlinear_σ₀_min = minimum(nlinear_σ₀) .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))
nlinear_σ₀_max = maximum(nlinear_σ₀) .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))

lines!(ax[2], S, T, nlinear_σ₀; label = "S-T profile", color = :tomato, linestyle = :dash)
S_minmax = [S[min_nlinear_idx], S[max_nlinear_idx]]
T_minmax = [T[min_nlinear_idx], T[max_nlinear_idx]]

σ_nl_levels = [nlinear_σ₀_min, nlinear_σ₀_upper, nlinear_σ₀_lower, nlinear_σ₀_max]
cn = contour!(ax[2], S_range, Θ_range, nlinear_σ₀_grid';
              levels = σ_nl_levels,
              colormap = σ_grad,
              colorrange = (-maximum(σ_nl_levels), maximum(σ_nl_levels)),
            #   labelformatter = my_contour_label_formatter,
            #   labels = true,
            #   color = σ_grad,
                # labelsize = 14,
              linewidth = 0.75,
            #   linestyle = :dot,
              label = "Isopycnals")
Colorbar(fig[2, 2], colorrange = extrema(σ_nl_levels),
        colormap = σ_grad,
        ticks = (range(nlinear_σ₀_min, nlinear_σ₀_max, length = 4), string.(round.(σ_nl_levels, digits = 3))),
        label = L"$\Delta \sigma_{0}'$ (kgm$^{-3}$)")
scatter!(ax[2], S_minmax[1], T_minmax[1]; markersize, label = "Minimum density", color = σ_grad[1])
scatter!(ax[2], S[end], T[end]; markersize, label = "Initial upper density", color = σ_grad[2])
scatter!(ax[2], S[1], T[1]; markersize, label = "Initial lower density", color = σ_grad[3])
scatter!(ax[2], S_minmax[2], T_minmax[2]; markersize, label = "Maximum density", color = σ_grad[end])

Legend(fig[3, :], ax[2], orientation = :horizontal, nbanks = 2)
fig
##
save("fig2_S_T_sigma_ST_space_2panel.png", fig)

## Figure
# DNS flow evolution
# file = jldopen(nl_R_ρ_105)
files = (l_R_ρ_105_dT2_diagnostics, nl_R_ρ_105_dT2_diagnostics)
file = jldopen(files[2])
x = file["dims/x_caa"]
y = file["dims/y_aca"]
z = file["dims/z_aac"]
times = file["dims/time"]

x_xz = repeat(x, 1, length(z))
x_xy = repeat(x, 1, length(y))
y_xz = repeat(y, 1, length(z))
y_xy = repeat(y, 1, length(x))
z_xz = repeat(reshape(z, 1, length(z)), length(x), 1)
z_yz = repeat(reshape(z, 1, length(z)), length(x), 1)

x_xz_velocity = x[1] * ones(length(x), length(z))
x_xz_velocity_end = x[end] * ones(length(x), length(z))
y_xz_density = y[1] * ones(length(x), length(z))
y_xz_density_end = y[end] * ones(length(x), length(z))

z_xy_top = z[end] * ones(length(x), length(y))
close(file)

horizontal_ticks = LinRange(-0.025, 0.025, 5)
T_colorrange = (-1.5, 1.5)
S_colorrange = (-0.09, 0.09)
w_colorrnage =  (-0.0005, 0.0005)
fig = Figure(size = (1500, 1300))
snapshots = [[60.0 * 3, 60.0 * 6, 60 * 12, 60.0 * 24],
             [180.00000000000003, 60.0 * 6, 60 * 12, 60.0 * 24]]
for (j, path) ∈ enumerate(files)

    file = jldopen(path)

    initial_S = file["S"]["xzslice_0.0"]
    initial_T = file["T"]["yzslice_0.0"]
    slices = [(S_xz = file["S"]["xzslice_$(snapshot)"],
            T_yz = file["T"]["yzslice_$(snapshot)"],
            velocity_zmean = file["w"]["w_zmean_$(snapshot)"][:, :, 1]) for snapshot ∈ snapshots[j]]
    close(file)

    ax = [Axis3(fig[j, i],
                # title = i == 1 ? "Experiment I" : "Experiment IV",
                aspect=(1/3, 1/3, 1),
                titlefont = :regular,
                xlabel = "x (cm)",
                ylabel = "y (cm)",
                zlabel = "z (m)",
                xticks = (horizontal_ticks, string.(horizontal_ticks .* 100)),
                yticks = (horizontal_ticks, string.(horizontal_ticks .* 100)),
                xlabeloffset = 40,
                ylabeloffset = 40,
                zlabeloffset = 50,
                xlabelsize = 14,
                ylabelsize = 14,
                zlabelsize = 16,
                xticklabelsize = 12,
                yticklabelsize = 12,
                zticklabelsize = 16,
                zlabelrotation = π / 2,
                limits = ((x[1], x[end]), (y[1], y[end]), (z[1], z[end])),
                elevation = π / 6.5,
                azimuth = 1.25π,
                xspinesvisible = false,
                yspinesvisible = false,
                zspinesvisible = false,
                zgridvisible = false,
                protrusions = 20
                ) for i ∈ 1:4]

    for i ∈ eachindex(snapshots[j])
        sf_S = surface!(ax[i], x_xz, y_xz_density, z_xz; color = slices[i].S_xz .- initial_S,
                        colormap = :curl,
                        colorrange = S_colorrange)
        sf_T = surface!(ax[i], x_xz_velocity, y_xz, z_yz; color = slices[i].T_yz .- initial_T,
                        colormap = :delta, colorrange = T_colorrange,
                        backlight = 5f0, shading = FastShading)
        sf_w = surface!(ax[i], x, y, z_xy_top; color = slices[i].velocity_zmean, colormap = :balance,
                        colorrange = w_colorrnage)
        ax[i].title = j == 1 ? "I, t = $(round(snapshots[j][i] / 60)) minutes" :
                               "IV, t = $(round(snapshots[j][i] / 60)) minutes"
        if j == 2
            if i == 1
                Colorbar(fig[j+1, 1:2], sf_S, label = "S′ (gkg⁻¹)", vertical = false, flipaxis = false)
                Colorbar(fig[j+1, 3:4], sf_T, label = "Θ′ (°C)", vertical = false, flipaxis = false)
                Colorbar(fig[:, 5], sf_w, label = "w (ms⁻¹)")
            end
        end
        if i > 1
            hidezdecorations!(ax[i], ticks = false)
        end
    end
end
Label(fig[0, :], "Salinity and temperature evolution", fontsize = 22, font = :bold)
rowgap!(fig.layout, 2, Relative(0.05))
rowgap!(fig.layout, 3, Relative(0.05))
fig
##
save("fig3_S_and_T_dns_evolution.png", fig)

## Figure
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
ax1 = Axis(fig[1, 1], title = L"Asymmetry due to $\tau$", titlefont = :bold, xlabel = L"τ", ylabel = L"R_{\Delta\rho}")
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
save("fig4_tau_asymmetry.png", fig)

## Figure
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
        # δ_leos[i, j] = Δρ′_linear #/ Δρ_linear
        # δ_nleos[i, j] = Δρ′_nlinear[i, j] #/ Δρ_nlinear
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

# δ_nleos = δ_nleos ./ maximum(δ_nleos[.!isnan.(δ_nleos)]) # normalise

ΔΘ = Θᵤ_range .- Θₗ
ΔS = Sᵤ_range .- Sₗ

arctic_obs = [[-0.04], [-0.014]]
ΔΘ_expts = [-2, -1, -0.5]
ΔS_expts_linear = [-0.12, -0.06, -0.03]
ΔS_expts_nlinear = [-0.12, -0.069, -0.037]

##
markersize = 15
linear_colour = :black
nlinear_colour = linear_colour
fig = Figure(size = (800, 1000))
ax_lRᵨ = Axis(fig[1, 1], xlabel = "ΔΘ (°C)", ylabel = "ΔS (gkg⁻¹)", title = "Linear eos")
hm = heatmap!(ax_lRᵨ, ΔΘ, ΔS, Rᵨ_leos, colormap = :amp)
scatter!(ax_lRᵨ, arctic_obs...; color = :red, markersize)
scatter!(ax_lRᵨ, ΔΘ_expts, ΔS_expts_linear;
         color = linear_colour, markersize, marker = linear_expt_markers)
ax_nlRᵨ = Axis(fig[1, 2], xlabel = "ΔΘ (°C)", ylabel = "ΔS (gkg⁻¹)", title = "Nonlinear eos")
hm_nlRᵨ = heatmap!(ax_nlRᵨ, ΔΘ, ΔS, Rᵨ_nleos, colormap = :amp)
scatter!(ax_nlRᵨ, arctic_obs...; markersize, color = :red)
scatter!(ax_nlRᵨ, ΔΘ_expts, ΔS_expts_nlinear;
         color = nlinear_colour, markersize, marker = nlinear_expt_markers)
hideydecorations!(ax_nlRᵨ, grid = false, ticks = false)
hidexdecorations!(ax_nlRᵨ, grid = false, ticks = false)
hidexdecorations!(ax_lRᵨ, grid = false, ticks = false)
Colorbar(fig[1, 3], hm_nlRᵨ, label = "Rᵨ")

colorrange = (minimum(Δσ_nonlinear[.!isnan.(Δσ_nonlinear)]), 1)
ax_lR_Δρ = Axis(fig[2, 1], xlabel = "ΔΘ (°C)", ylabel = "ΔS (gkg⁻¹)")
hm = heatmap!(ax_lR_Δρ, ΔΘ, ΔS, Δσ_linear; colorrange, colormap = :batlow)
scatter!(ax_lR_Δρ, arctic_obs...; color = :red, markersize)
scatter!(ax_lR_Δρ, ΔΘ_expts, ΔS_expts_linear;
         color = linear_colour, markersize, marker = linear_expt_markers)
ax_nlR_Δρ = Axis(fig[2, 2], xlabel = "ΔΘ (°C)", ylabel = "ΔS (gkg⁻¹)")
hideydecorations!(ax_nlR_Δρ, grid = false, ticks = false)
hm_nlR_Δρ = heatmap!(ax_nlR_Δρ, ΔΘ, ΔS, Δσ_nonlinear; colorrange, colormap = :batlow)
scatter!(ax_nlR_Δρ, arctic_obs...; markersize, color = :red)
scatter!(ax_nlR_Δρ, ΔΘ_expts, ΔS_expts_nlinear;
         color = nlinear_colour, markersize, marker = nlinear_expt_markers)
hidexdecorations!(ax_lR_Δρ, grid = false, ticks = false)
hidexdecorations!(ax_nlR_Δρ, grid = false, ticks = false)
Colorbar(fig[2, 3], hm_nlR_Δρ, label = "R_Δρ")

# colorrange = (0, maximum(δ_nleos[.!isnan.(δ_nleos)]))
colorrange = (1, 10)
ax_δ_linear = Axis(fig[3, 1], xlabel = "ΔΘ (°C)", ylabel = "ΔS (gkg⁻¹)")
hm_δ_linear = heatmap!(ax_δ_linear, ΔΘ, ΔS, δ_leos; colorrange, colormap = :turbid, highclip = :red, lowclip = :orange)
scatter!(ax_δ_linear, arctic_obs...; markersize, color = :red)
scatter!(ax_δ_linear, ΔΘ_expts, ΔS_expts_linear; color = linear_colour,
         label = linear_expt_labels, markersize, marker = linear_expt_markers)
ax_δ_nlinear = Axis(fig[3, 2], xlabel = "ΔΘ (°C)", ylabel = "ΔS (gkg⁻¹)")
hm_δ_nlinear = heatmap!(ax_δ_nlinear, ΔΘ, ΔS, δ_nleos; colorrange, colormap = :turbid, highclip = :red, lowclip = :orange)
scatter!(ax_δ_nlinear, arctic_obs...; markersize, color = :red, label = "Arctic interfaces")
scatter!(ax_δ_nlinear, ΔΘ_expts, ΔS_expts_nlinear; color = nlinear_colour,
         label = nlinear_expt_labels, markersize, marker = nlinear_expt_markers)
hideydecorations!(ax_δ_nlinear, grid = false, ticks = false)
Colorbar(fig[3, 3], hm_δ_nlinear, label = "Δρ′")
legend_markers = [MarkerElement(color = :black, marker = m; markersize)
                  for m ∈ vcat(linear_expt_markers, nlinear_expt_markers)]
legend_expts = vcat(linear_expt_labels, nlinear_expt_labels)
Legend(fig[4, 1], legend_markers, legend_expts, "DNS Experiments",
       orientation = :horizontal, nbanks = 2)
lmarker_2 = [MarkerElement(color = :red, marker = :circle; markersize)]
lobs = ["Timmermans et al. (2008)"]
Legend(fig[4, 2], lmarker_2, lobs, "Observations",
       orientation = :horizontal)
fig
##
save("fig5_density_asymmetry.png", fig)
##
# RΔρ for observations
find_Θ = findall(-0.5 .<= ΔΘ .<= -0.1) # Southern Ocean
SO_vals = Δσ_nonlinear[find_Θ, :]
replace!(SO_vals, NaN => 0)
nz = findall(SO_vals .!= 0)
extrema(SO_vals[nz])

find_Θ = findall(-0.81 .<= ΔΘ .<= -0.79)
SO_vals = Δσ_nonlinear[find_Θ, :]
replace!(SO_vals, NaN => 0)
nz = findall(SO_vals .!= 0)
extrema(SO_vals[nz])

## Figure
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
        S = file["S_ha"][:, 2+i]
        Θ = file["T_ha"][:, 2+i]
        σ = file["σ_ha"][:, 2+i]
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
        text!(ax[i], 34.7, -1.2, text = L"1D model $R_{\Delta\rho} =$ %$(R_Δρ_model)", align = (:right, :top))
        text!(ax[i], 34.7, -1.4, text = L"Simulation $R_{\Delta\rho} =$ %$(R_Δρ_sim)", align = (:right, :top))
        ax[i].subtitle = L"$t~=$ %$(round(t[2+i]/60, digits = 1)) min"
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
save("fig6_ST_simluation.png", fig)

## Figure
# Density hovmollers for Rρ = 1.05
files = (l_R_ρ_105_dT2_diagnostics, nl_R_ρ_105_dT2_diagnostics,
         l_R_ρ_105_dT1_diagnostics, nl_R_ρ_105_dT1_diagnostics,
         l_R_ρ_105_dT05_diagnostics, nl_R_ρ_105_dT05_diagnostics)
fig = Figure(size = (800, 1000))
axσ = [Axis(fig[j, i], xlabel = "time (min)", ylabel = "z (m)") for i ∈ 1:2, j ∈ 1:3]
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
# axwT = [Axis(fig[2, i], xlabel = "time (min)", ylabel = "z (m)") for i ∈ eachindex(files)]
# wT_colorrange = jldopen(files[1]) do ds
#         extrema(ds["ha_wT"])
# end
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
        σ_ha′ = ds["σ_ha"][:, 1:ts_length ]' .- mean(ds["σ_ha"][:, 1])'
        hm = heatmap!(axσ[i], 1:ts_length, z, σ_ha′; colorrange = crange, colormap = :diff)
        axσ[i].title = "Experiment " * all_labels[i]
        if i % 2 == 0
            Colorbar(fig[Int(i / 2), 3], hm, label = L"$σ_{0}′$ (kgm$^{-3}$)")
            hideydecorations!(axσ[i], ticks = false)
        end
        if i < 5
            hidexdecorations!(axσ[i], ticks = false)
        end
    end

end
fig
##
save("fig7_density_hovs.png", fig)
##
## Figure
# salinity midpoint
files = (l_R_ρ_105_dT2_diagnostics, nl_R_ρ_105_dT2_diagnostics,
         l_R_ρ_105_dT1_diagnostics, nl_R_ρ_105_dT1_diagnostics,
         l_R_ρ_105_dT05_diagnostics, nl_R_ρ_105_dT05_diagnostics)
fig_interface = Figure(size = (800, 600))
ax_R_rho = Axis(fig_interface[1, 1],
                title = "(a) Density ratio",
                xlabel = "time ( min)",
                ylabel = L"R_{\rho}")
ax_interface = Axis(fig_interface[2, 1],
                    title = "(b) Salinity interface height",
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
fig_interface
##
save("salinity_interfaces.png", fig_interface)

## Figure
# energetics
files = (l_R_ρ_105_dT2_diagnostics, nl_R_ρ_105_dT2_diagnostics,
         l_R_ρ_105_dT1_diagnostics, nl_R_ρ_105_dT1_diagnostics,
         l_R_ρ_105_dT05_diagnostics, nl_R_ρ_105_dT05_diagnostics)
fig = Figure(size = (800, 600))
ax = [Axis(fig[i, 1], xlabel = "time (min)", ylabel = "Non-dimensional energy") for i ∈ 1:2]
ax[1].title = "Available potential energy"
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
        if i == 1
            hidexdecorations!(ax[i], grid = false, ticks = false)
        end

    end
end
linkxaxes!(ax[1], ax[2])
axislegend(ax[1], orientation = :horizontal, nbanks = 2, position = :rb)
fig
##
save("fig9_ape.png", fig)

## Figure
# Flux boundary conditions
fig = Figure(size = (800, 800))
files = (l_fbc_diagnostics, nl_fbc_diagnostics)
labels = ("II with flux boundary conditions","V with flux boundary conditions")
titles = ("Experiment II with flux boundary conditions",
          "Experiment V with flux boundary conditions")
ax = [Axis(fig[i, 1]) for i ∈ 1:3]
ts_range = 1:360
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
				ylabel = "Rᵨ")
	lines!(ax2, t, Rᵨ, color = :dodgerblue, label = L"R_{\rho}")
    linkxaxes!(ax2, ax[i])
    hidexdecorations!(ax2)
	axislegend(ax[i], merge = true, position = :rb)
    lines!(ax[3], t[1:end], z✶, label = labels[i])
    ax[3].xlabel = "time (min)"
    ax[3].ylabel = L"$z*$ (m)"
    ax[3].title = "Salinity interface"
end

files = (l_R_ρ_105_dT1_diagnostics, nl_R_ρ_105_dT1_diagnostics)
ts_range = 1:240
labels_nf = ("II", "V")
for (i, file) ∈ enumerate(files)
    jldopen(file) do f
        t = f["dims"]["time"][ts_range]
        z✶ = f["dims"]["z✶"][f["S_interface_idx"]][ts_range]
        Rᵨ = f["R_ρ"][ts_range]
        lines!(ax[3], t ./ 60, z✶, label = labels_nf[i],
                linestyle = :dash, color = Makie.wong_colors()[i])
    end
end

Legend(fig[4, :], ax[3], orientation = :horizontal, nbanks = 2)
hidexdecorations!(ax[1], grid = false, ticks = false)
hidexdecorations!(ax[2], grid = false, ticks = false)
fig
##
f = jldopen(nl_fbc_diagnostics)
linfit_range = 100:360
z✶ = f["dims"]["z✶"][f["S_interface_idx"][linfit_range]]
t = f["dims"]["time"][linfit_range]
A = [ones(length(t)) t]
a, b = A \ z✶
lines!(ax[3], t./60, a .+ t .* b)
fig
uₑ = b
ΔS = abs(mean(f["ΔS"][linfit_range]))
ΔT = abs(mean(f["ΔT"][linfit_range]))
Jₛ = f["FluxBCs/Jˢ"] * f["FluxBCs/top_salinity_scale"]
Jₜ = f["FluxBCs/Jᵀ"]
Ẽ = (uₑ * ΔS) / (2 * 0.25 * Jₛ + uₑ * ΔS)
Ẽ = (uₑ * ΔT) / (2 * 0.25 * Jₜ + uₑ * ΔT)
##
save("flux_bcs.png", fig)

## Calculations
# Δρ′
# I, IV
Θᵤ, Θₗ = -1.5, 0.5
Sᵤ, Sₗ = 34.58, 34.7
# II, V
Θᵤ, Θₗ = -0.5, 0.5
Sᵤ, Sₗ = 34.64, 34.7
Sᵤ, Sₗ = 34.631, 34.7
# III, VI
Θᵤ, Θₗ = 0, 0.5
Sᵤ, Sₗ = 34.67, 34.7
Sᵤ, Sₗ = 34.663, 34.7


ΔΘ_mix = Θₗ - Θᵤ
ΔS_mix = Sₗ - Sᵤ
S_mix = range(Sᵤ, Sₗ, length = 1000)
Θ_mix = Θᵤ .+ (ΔΘ_mix / ΔS_mix) * (S_mix .- Sᵤ)
ρ_mix = total_density.(Θ_mix, S_mix, fill(0, length(S_mix)), fill(nleos, length(S_mix)))
ρ_max = maximum(ρ_mix)
ρₗ = total_density(Θₗ, Sₗ, 0, nleos)
ρ_max - ρₗ

## Figure
# Batchelor lengths and resolution
files = (l_R_ρ_105_dT2_diagnostics, nl_R_ρ_105_dT2_diagnostics,
         l_R_ρ_105_dT1_diagnostics, nl_R_ρ_105_dT1_diagnostics,
         l_R_ρ_105_dT05_diagnostics, nl_R_ρ_105_dT05_diagnostics)
fig = Figure(size = (800, 800))
ax = [Axis(fig[i, 1],
          title = "Batchelor lengths and resolution",
          xlabel = "time (mins)",
          ylabel = "Length (log10(mm))") for i ∈ 1:3]
ts_length = 3:240
linestyles = (:solid, :dash)

for (i, file) ∈ enumerate(files)

    ls = i % 2 == 0 ? linestyles[2] : linestyles[1]
    jldopen(file) do f
        j = i ≤ 2 ? 1 : i > 4 ? 3 : 2
        Ba = 1e3 * f["Ba"][ts_length]
        lines!(ax[j], log10.(Ba), label = all_labels[i], linestyle = ls, color = Makie.wong_colors()[i])

        j > 1 ? ax[j].title = "" : nothing
        j < 3 ? hidexdecorations!(ax[j], grid = false, ticks = false) : nothing
    end

    jldopen(file) do f
        Ba = 1e3 * f["Ba"][ts_length]
        spcaings = vcat(diff(f["dims"]["x_caa"]),
                        diff(f["dims"]["y_aca"]),
                        diff(f["dims"]["z_aac"]))
        Δ = minimum(spcaings) * 1e3
        println(Δ)
        j = i ≤ 2 ? 1 : i > 4 ? 3 : 2
        i % 2 == 0 ? lines!(ax[j], log10.(fill(Δ, length(Ba))), color = :red, linestyle = :dot,
                            label = "Grid resolution") : nothing

        i % 2 == 0 ? axislegend(ax[j], position = :rt) : nothing
    end

end
fig
##
save("fig11_Batchelor_lengths.png", fig)
## Figure
# Energy budgets
files = (l_R_ρ_105_dT2_diagnostics, l_R_ρ_105_dT1_diagnostics, l_R_ρ_105_dT05_diagnostics,
         nl_R_ρ_105_dT2_diagnostics, nl_R_ρ_105_dT1_diagnostics, nl_R_ρ_105_dT05_diagnostics)
expt_labels_ordered = ("I", "II", "III", "IV", "V", "VI")
fig = Figure(size = (800, 1000))
ax = [Axis(fig[i, j],
          title = "Batchelor lengths and resolution",
          xlabel = "time (mins)") for i ∈ 1:3, j ∈ 1:2]
ts_length = 1:240
for (i, file) ∈ enumerate(files)

    ls = i % 2 == 0 ? linestyles[2] : linestyles[1]
    jldopen(file) do f

        Δt = diff(f["dims"]["time"][ts_length])
        dₜek = diff(f["∫Eₖ"][ts_length]) ./ Δt
        ε = 0.5 * (f["∫ε"][1:ts_length[end-1]] .+ f["∫ε"][2:ts_length[end]])
        ∫wb = 0.5 * (f["∫wb"][1:ts_length[end-1]] .+ f["∫wb"][2:ts_length[end]])
        RHS = ∫wb .- ε

        mean_abs_err = round(mean(abs.(dₜek .- RHS)), digits = 15)
        lines!(ax[i], dₜek, label = L"\rho_{0}^{-1}\frac{\mathrm{d}}{\mathrm{d}t}\int E_{k}~\mathrm{d}V")
        lines!(ax[i], RHS, label = L"\int bw~\mathrm{d}V - ε", linestyle = :dash)
        ax[i].title = expt_labels_ordered[i] * ", MAE = $(mean_abs_err)"
        i ∈ (3, 6) ? nothing : hidexdecorations!(ax[i], grid = false, ticks = false)
        i < 4 ? nothing : hideydecorations!(ax[i], grid = false, ticks = false)
    end

end
linkyaxes!(ax[4], ax[1])
linkyaxes!(ax[5], ax[2])
linkyaxes!(ax[6], ax[3])
Legend(fig[4, :], ax[1], orientation = :horizontal)
fig
##
save("fig12_energy_budgets.png", fig)

##
# δ parameter
# McDougall expt 25th Oct check
S = [0.04, 6.01]
Θ = [15.5, 30.39]
ΔS = diff(S)
Θ = gsw_ct_from_t.(S, Θ, 0)
ΔΘ = diff(Θ)
S_mix = range(S..., length = 1000)
T_mix = Θ[2] .+ (ΔΘ / ΔS) .* (S_mix .- S[2])
ρ = gsw_rho.(S, Θ, 0)
Δρ = diff(ρ)
Δρ′ = maximum(gsw_rho.(S_mix, T_mix, 0)) - gsw_rho(S[2], Θ[2], 0)

δ(C) = 1 + 2*C + 2*sqrt(C^2 + C)
δ(Δρ′/Δρ[1]) # checks out

## My experiments
#  Non-linear eos
δᵢ_nleos = Vector{Float64}(undef, 3)
# ΔΘ = -1
S = [34.58, 34.7]
ΔS = diff(S)
Θ = [-1.5, 0.5]
ΔΘ = diff(Θ)
S_mix = range(S..., length = 1000)
T_mix = Θ[2] .+ (ΔΘ / ΔS) .* (S_mix .- S[2])
ρ = gsw_rho.(S, Θ, 0)
Δρ = diff(ρ)
Δρ′ = maximum(gsw_rho.(S_mix, T_mix, 0)) - gsw_rho(S[2], Θ[2], 0)
δᵢ_nleos[1] = δ(Δρ′/Δρ[1])
# ΔΘ = -1
S = [34.631, 34.7]
ΔS = diff(S)
Θ = [-0.5, 0.5]
ΔΘ = diff(Θ)
S_mix = range(S..., length = 1000)
T_mix = Θ[2] .+ (ΔΘ / ΔS) .* (S_mix .- S[2])
ρ = gsw_rho.(S, Θ, 0)
Δρ = diff(ρ)
Δρ′ = maximum(gsw_rho.(S_mix, T_mix, 0)) - gsw_rho(S[2], Θ[2], 0)
δᵢ_nleos[2] = δ(Δρ′/Δρ[1])
# ΔΘ = -0.5
S = [34.663, 34.7]
ΔS = diff(S)
Θ = [0.0, 0.5]
ΔΘ = diff(Θ)
S_mix = range(S..., length = 1000)
T_mix = Θ[2] .+ (ΔΘ / ΔS) .* (S_mix .- S[2])
ρ = gsw_rho.(S, Θ, 0)
Δρ = diff(ρ)
Δρ′ = maximum(gsw_rho.(S_mix, T_mix, 0)) - gsw_rho(S[2], Θ[2], 0)
δᵢ_nleos[3] =  δ(Δρ′/Δρ[1])
for i ∈ eachindex(δᵢ_nleos)
    println("$(nlinear_expt_labels[i]) = $(round(δᵢ_nleos[i], digits = 2))")
end
# Linear eos
δᵢ_leos = Vector{Float64}(undef, 3)
# ΔΘ = -2
S = [34.58, 34.7]
ΔS = diff(S)
Θ = [-1.5, 0.5]
ΔΘ = diff(Θ)
S_mix = range(S..., length = 1000)
T_mix = Θ[2] .+ (ΔΘ / ΔS) .* (S_mix .- S[2])
ρ = total_density.(Θ, S, fill(0, length(2)), fill(leos, length(S)))
Δρ = diff(ρ)
Δρ′ = maximum(total_density.(T_mix, S_mix, fill(0, length(S_mix)), fill(leos, length(S_mix))) .-
                total_density(Θ[2], S[2], 0, leos))
δᵢ_leos[1] = δ(Δρ′/Δρ[1])
# ΔΘ = -1
S = [34.631, 34.7]
ΔS = diff(S)
Θ = [-0.5, 0.5]
ΔΘ = diff(Θ)
S_mix = range(S..., length = 1000)
T_mix = Θ[2] .+ (ΔΘ / ΔS) .* (S_mix .- S[2])
ρ = total_density.(Θ, S, fill(0, length(2)), fill(leos, length(S)))
Δρ = diff(ρ)
Δρ′ = maximum(total_density.(T_mix, S_mix, fill(0, length(S_mix)), fill(leos, length(S_mix))) .-
                total_density(Θ[2], S[2], 0, leos))
δᵢ_leos[2] = δ(Δρ′/Δρ[1])
# ΔΘ = -0.5
S = [34.663, 34.7]
ΔS = diff(S)
Θ = [0.0, 0.5]
ΔΘ = diff(Θ)
S_mix = range(S..., length = 1000)
T_mix = Θ[2] .+ (ΔΘ / ΔS) .* (S_mix .- S[2])
ρ = total_density.(Θ, S, fill(0, length(2)), fill(leos, length(S)))
Δρ = diff(ρ)
Δρ′ = maximum(total_density.(T_mix, S_mix, fill(0, length(S_mix)), fill(leos, length(S_mix))) .-
                total_density(Θ[2], S[2], 0, leos))
δᵢ_leos[3] = δ(Δρ′/Δρ[1])
for i ∈ eachindex(δᵢ_leos)
    println("$(linear_expt_labels[i]) = $(round(δᵢ_leos[i], digits = 2))")
end
