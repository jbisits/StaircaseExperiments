using JLD2, CairoMakie

cd(@__DIR__)
linear_interface_fluxes = load("lineareos_single_interface_360min/interface_fluxes_linear.jld2")
S_linear_mol_flux = 0.5 *(linear_interface_fluxes["Smolecular_flux"][3, 1:end-1] .+
                          linear_interface_fluxes["Smolecular_flux"][3, 2:end])
T_linear_mol_flux = 0.5 *(linear_interface_fluxes["Tmolecular_flux"][3, 1:end-1] .+
                          linear_interface_fluxes["Tmolecular_flux"][3, 2:end])
nonlinear_interface_fluxes = load("nonlineareos_single_interface_360min/interface_fluxes_nonlinear.jld2")
S_nonlinear_mol_flux = 0.5 *(nonlinear_interface_fluxes["Smolecular_flux"][3, 1:end-1] .+
                             nonlinear_interface_fluxes["Smolecular_flux"][3, 2:end])
T_nonlinear_mol_flux = 0.5 *(nonlinear_interface_fluxes["Tmolecular_flux"][3, 1:end-1] .+
                             nonlinear_interface_fluxes["Tmolecular_flux"][3, 2:end])

z✶ = range(1, 0, length = 70*70*1000) # hack because not saved

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
axislegend(ax_S)
ax_T = Axis(fig[2, 1], xlabel = "time (min)", ylabel = "Temperature flux",
            title = "Temperature flux through interface")
lines!(ax_T, linear_interface_fluxes["T_flux"][2, :], label = "Linear eos")
lines!(ax_T, nonlinear_interface_fluxes["T_flux"][2, :], label = "Nonlinear eos")
axislegend(ax_T)
ax_z = Axis(fig[3, 1], xlabel = "time (min)", ylabel = "Depth (m)",
            title = "Depth of interface")
lines!(ax_z, z✶[linear_interface_fluxes["T_interface_idx"]], label = "Linear eos")
lines!(ax_z, z✶[nonlinear_interface_fluxes["T_interface_idx"]], label = "Nonlinear eos")
axislegend(ax_z)
fig
