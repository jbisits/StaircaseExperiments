using StaircaseShenanigans

diagnostics_file = "diagnostics.jld2"
time_length = ("60min", "240min", "480min")
groups = ("noise", "lessnoise", "nonoise")

for (i, t) ∈ time_length
    data_path = joinpath(@__DIR__, "nonlineareos_single_interface_"*t)
    tracers = joinpath(data_path, "tracers.nc")
    computed_output = joinpath(data_path, "computed_output.nc")

    save_diagnostics!(diagnostics_file, tracers, computed_output, group = groups[i])
end
