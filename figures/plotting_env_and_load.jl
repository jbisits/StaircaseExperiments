using StaircaseShenanigans
using StaircaseShenanigans: compute_R_ρ
using JLD2, GibbsSeaWater, StatsBase, ColorSchemes
using SeawaterPolynomials: TEOS10EquationOfState, total_density, haline_contraction, thermal_expansion
using Printf
using SpecialFunctions: erf
using CairoMakie # GLMakie better for surface plots

# Channge to plotting directory
figure_dir = @__DIR__
if pwd() != figure_dir
    cd(figure_dir)
end

## EOS's and other constants
ρ₀ = gsw_rho(34.7, 0.5, 0)
leos = CustomLinearEquationOfState(-0.5, 34.6, reference_density = ρ₀)
leos_func(S, Θ) = CustomLinearEquationOfState(Θ, S, reference_density = ρ₀)
nleos = TEOS10EquationOfState(reference_density = ρ₀)
cab_eos = RoquetEquationOfState(:Cabbeling, reference_density = ρ₀)
erf_tracer_solution(z, Cₗ::Number, ΔC::Number, κ::Number, t, interface_depth) =
    Cₗ + 0.5 * ΔC * (1 + erf((z - interface_depth) / sqrt(4 * κ * t)))
δ(C) = 1 + 2*C + 2*sqrt(C^2 + C)

Sᵤ, Θᵤ = 34.58, -1.5
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
rundown_path = joinpath(@__DIR__, "../DNS_DDC_expts/")
nl_R_ρ_105_dT2_diagnostics  = joinpath(rundown_path, "dns_res_dT2_nonlineareos/step_diagnostics.jld2")
l_R_ρ_105_dT2_diagnostics   = joinpath(rundown_path, "dns_res_dT2_lineareos/step_diagnostics.jld2")
nl_R_ρ_105_dT1_diagnostics  = joinpath(rundown_path, "dns_res_dT1_nonlineareos/step_diagnostics.jld2")
l_R_ρ_105_dT1_diagnostics   = joinpath(rundown_path, "dns_res_dT1_lineareos/step_diagnostics.jld2")
nl_R_ρ_105_dT05_diagnostics = joinpath(rundown_path, "dns_res_dT05_nonlineareos/step_diagnostics.jld2")
l_R_ρ_105_dT05_diagnostics  = joinpath(rundown_path, "dns_res_dT05_lineareos/step_diagnostics.jld2")
# flux bc experiments
fluxbc_path = joinpath(@__DIR__, "../DNS_DDC_expts/")
l_fbc_diagnostics  = joinpath(fluxbc_path, "dns_fluxbc_res_lineareos/step_diagnostics.jld2")
nl_fbc_diagnostics = joinpath(fluxbc_path, "dns_fluxbc_res_nonlineareos/step_diagnostics.jld2")

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
