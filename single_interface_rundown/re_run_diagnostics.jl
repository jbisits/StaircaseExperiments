using NCDatasets, JLD2, Statistics

# update output path to where files are located
output_path = joinpath(@__DIR__, "dns_rundown_1.05/step/nonlineareos_single_interface_60min")
co = joinpath(output_path, "computed_output.nc")
tracers = joinpath(output_path, "tracers.nc")
velocities = joinpath(output_path, "velocities.nc")

diags = joinpath(output_path, "step_diagnostics.jld2")
if isfile(diags)
    rm(diags)
end
save_diagnostics!(diags, tracers, computed_output, velocities)
