using StaircaseShenanigans

initial_states = ("step")
groups = ("lineareos", "nonlineareos")

# for initial_state ∈ initial_states
initial_state = "step"
    diagnostics_file = initial_state*"/diff_initial_h_diagnostics.jld2"
    for group ∈ groups


        path = joinpath(@__DIR__, initial_state, group*"_single_interface_60min")
        tracers = joinpath(path, "tracers.nc")
        computed_output = joinpath(path, "computed_output.nc")

        save_diagnostics!(diagnostics_file, tracers, computed_output; group)

    end

# end
