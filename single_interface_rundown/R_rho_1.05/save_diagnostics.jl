using StaircaseShenanigans

diagnostics_file = "diagnostics.jld2"

linear_path = joinpath(@__DIR__, "R_rho_1.05", "step", "lineareos_single_interface_480min")
tracers = joinpath(linear_path, "tracers.nc")
computed_output = joinpath(linear_path, "computed_output.nc")

save_diagnostics!(diagnostics_file, tracers, computed_output, eos = "linear")

nonlinear_path = joinpath(@__DIR__, "R_rho_1.05", "step", "nonlineareos_single_interface_480min")
tracers = joinpath(nonlinear_path, "tracers.nc")
computed_output = joinpath(nonlinear_path, "computed_output.nc")

save_diagnostics!(diagnostics_file, tracers, computed_output, eos = "nonlinear")
