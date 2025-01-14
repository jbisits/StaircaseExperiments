include("../S_and_T_fluxes.jl")

cd(@__DIR__)
tracers = "tracers.nc"
flux_file = "interface_fluxes_linear.jld2"

φ_interface_flux!(flux_file, tracers, :S)
φ_interface_flux!(flux_file, tracers, :T)
φ_molelcuar_flux!(flux_file, tracers, :S)
φ_molelcuar_flux!(flux_file, tracers, :T)
