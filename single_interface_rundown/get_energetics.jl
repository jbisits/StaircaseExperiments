using NCDatasets, JLD2, Statistics
using StaircaseShenanigans: potential_and_background_potential_energy!

output_path = joinpath(@__DIR__, "dns_rundown_1.05/step/lineareos_single_interface_120min")
co = joinpath(output_path, "computed_output.nc")
tracers = joinpath(output_path, "tracers.nc")
linear_energetics = joinpath(output_path, "energetics.jld2")
potential_and_background_potential_energy!(linear_energetics, co, tracers, "")

output_path = joinpath(@__DIR__, "dns_rundown_1.05/step/nonlineareos_single_interface_60min")
co = joinpath(output_path, "computed_output.nc")
tracers = joinpath(output_path, "tracers.nc")
nlinear_energetics = joinpath(output_path, "energetics.jld2")
potential_and_background_potential_energy!(nlinear_energetics, co, tracers, "")
