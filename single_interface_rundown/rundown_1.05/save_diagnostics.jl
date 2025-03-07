using StaircaseShenanigans

initial_state = "tanh"

diagnostics_file = initial_state*"/diagnostics.jld2"

linear_path = joinpath(@__DIR__, initial_state, "lineareos_single_interface_480min")
tracers = joinpath(linear_path, "tracers.nc")
computed_output = joinpath(linear_path, "computed_output.nc")

save_diagnostics!(diagnostics_file, tracers, computed_output, group = "linear")

nonlinear_path = joinpath(@__DIR__, initial_state, "nonlineareos_single_interface_480min")
tracers = joinpath(nonlinear_path, "tracers.nc")
computed_output = joinpath(nonlinear_path, "computed_output.nc")

save_diagnostics!(diagnostics_file, tracers, computed_output, group = "nonlinear")
