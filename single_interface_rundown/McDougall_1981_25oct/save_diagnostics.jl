using StaircaseShenanigans

diagnostics_file = "diagnostics.jld2"

data_path = joinpath(@__DIR__, "nonlineareos_single_interface_480min")
tracers = joinpath(data_path, "tracers.nc")
computed_output = joinpath(data_path, "computed_output.nc")

save_diagnostics!(diagnostics_file, tracers, computed_output, group = "diags")
