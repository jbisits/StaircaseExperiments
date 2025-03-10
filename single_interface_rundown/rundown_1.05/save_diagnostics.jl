using StaircaseShenanigans

initial_states = ("step", "tanh")
groups = ("lineareos", "nonlineareos")

for initial_state ∈ initial_states

    for group ∈ groups

        diagnostics_file = initial_state*"/diagnostics.jld2"

        path = joinpath(@__DIR__, initial_state, group*"_single_interface_480min")
        tracers = joinpath(path, "tracers.nc")
        computed_output = joinpath(path, "computed_output.nc")

        save_diagnostics!(diagnostics_file, tracers, computed_output; group)

    end

end
