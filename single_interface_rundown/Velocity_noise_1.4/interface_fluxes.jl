using JLD2, CairoMakie, GibbsSeaWater

cd(@__DIR__)
linear_interface_fluxes = load("lineareos_single_interface_360min/interface_fluxes_linear.jld2")
S_linear_mol_flux = 0.5 *(linear_interface_fluxes["S_molecular_flux"][3, 1:end-1] .+
                          linear_interface_fluxes["S_molecular_flux"][3, 2:end])
T_linear_mol_flux = 0.5 *(linear_interface_fluxes["T_molecular_flux"][3, 1:end-1] .+
                          linear_interface_fluxes["T_molecular_flux"][3, 2:end])
nonlinear_interface_fluxes = load("nonlineareos_single_interface_360min/interface_fluxes_nonlinear.jld2")
S_nonlinear_mol_flux = 0.5 *(nonlinear_interface_fluxes["S_molecular_flux"][3, 1:end-1] .+
                             nonlinear_interface_fluxes["S_molecular_flux"][3, 2:end])
T_nonlinear_mol_flux = 0.5 *(nonlinear_interface_fluxes["T_molecular_flux"][3, 1:end-1] .+
                             nonlinear_interface_fluxes["T_molecular_flux"][3, 2:end])

z✶ = range(1, 0, length = 70*70*1000) # hack because not saved
R_ρ_linear = linear_interface_fluxes["R_ρ"][2:end]
R_ρ_nonlinear = nonlinear_interface_fluxes["R_ρ"][2:end]
# The flux computation is the cumulative sum to the interface from the top to the bottom of
# the domain. This means there is an increase salinity and temperature content above the
# interface (i.e. increase in S and T in the upper layer) so we should have positive flux for
# both salinity and temperautre. This is because positive flux = flux upwards in Oceananigans.
# Note the flux depth here is regridded onto z✶ ∈ [0, 1].

## Fluxes
fig = Figure(size = (800, 600))
ax_S = Axis(fig[1, 1], ylabel = "Salinity flux", title = "Salinity flux through diffusive interface")
lines!(ax_S, linear_interface_fluxes["S_flux"][2, :], label = "Linear eos")
lines!(ax_S, nonlinear_interface_fluxes["S_flux"][2, :], label = "Nonlinear eos")
hidexdecorations!(ax_S, ticks = false)
axislegend(ax_S, position = :rb)
ax_T = Axis(fig[2, 1], xlabel = "time (min)", ylabel = "Temperature flux",
            title = "Temperature flux through interface")
lines!(ax_T, linear_interface_fluxes["T_flux"][2, :], label = "Linear eos")
lines!(ax_T, nonlinear_interface_fluxes["T_flux"][2, :], label = "Nonlinear eos")
axislegend(ax_T, position = :rb)
ax_z = Axis(fig[3, 1], xlabel = "time (min)", ylabel = "Depth (m)",
            title = "Depth of interface")
lines!(ax_z, linear_interface_fluxes["T_interface_depth"], label = "T interface linear eos")
lines!(ax_z, linear_interface_fluxes["S_interface_depth"], label = "S interface linear eos")
lines!(ax_z, nonlinear_interface_fluxes["T_interface_depth"], label = "T interface nonlinear eos")
lines!(ax_z, nonlinear_interface_fluxes["S_interface_depth"], label = "S interface nonlinear eos")
axislegend(ax_z, nbanks = 2, position = :rc)
fig

## Fₜ normalise by 4/3 flux law - horizontal line => Fₜ ∝ ΔΘ^(4/3)
Fₜ_linear = linear_interface_fluxes["T_flux"][2, :]
Fₜ_nonlinear = nonlinear_interface_fluxes["T_flux"][2, :]
λ = 0.085
ρ₀ = gsw_rho(34.7, 0.5, 0)
ν = 1e-6
κₜ = 1e-7
g = 9.81
ΔΘ = 2
Fₜᵀ = λ * κₜ * (g / (ρ₀*ν*κₜ))^(1/3) * ΔΘ^(4/3) # Turner
τ = 0.01
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], xlabel = "R_ρ", ylabel = "Fₜ / Fₜᵀ")
lines!(ax, R_ρ_linear, Fₜ_linear / Fₜᵀ, label = "linear eos")
lines!(ax, R_ρ_nonlinear, Fₜ_nonlinear / Fₜᵀ, label = "nonlinear eos")
axislegend(ax, position = :rb)
fig
## flux ratio vs R_ρ
Fₛ_linear = linear_interface_fluxes["S_flux"][2, :]
Fₛ_nonlinear = nonlinear_interface_fluxes["S_flux"][2, :]
linear_Fₛ_Fₜ = Fₛ_linear ./ Fₜ_linear
nonlinear_Fₛ_Fₜ = Fₛ_nonlinear ./ Fₜ_nonlinear
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], xlabel = "t", ylabel = "Fₛ / Fₜ")
lines!(ax, R_ρ_linear, linear_Fₛ_Fₜ, label = "linear eos")
lines!(ax, R_ρ_nonlinear, nonlinear_Fₛ_Fₜ, label = "nonlinear eos")
axislegend(ax, position = :rb)
lines!(ax, R_ρ_linear, τ^(6/5) .* R_ρ_linear)
fig
