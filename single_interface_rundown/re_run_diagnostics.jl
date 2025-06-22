using StaircaseShenanigans, CairoMakie, GibbsSeaWater

# Non-linear eos
output_path = joinpath(@__DIR__, "dns_rundown_1.05_dT_-1.0/step/nonlineareos_single_interface_120min")
computed_output = joinpath(output_path, "computed_output.nc")
tracers = joinpath(output_path, "tracers.nc")
velocities = joinpath(output_path, "velocities.nc")

@info "Re-producing animations"
animate_density(computed_output, "σ",
                xslice = 17, yslice = 17, density_limit_adjustment = 0.04)
animate_tracers(tracers, xslice = 17, yslice = 17,
                S_limit_adjustment = 0.025,
                Θ_limit_adjustment = 0.5)
animate_vertical_velocity(velocities, xslice = 17, yslice = 17)

@info "Re-computing diagnostics"
diags = joinpath(output_path, "step_diagnostics.jld2")
if isfile(diags)
    rm(diags)
end
save_diagnostics!(diags, tracers, computed_output, velocities)


# Linear eos
output_path = joinpath(@__DIR__, "dns_rundown_1.05_dT_-1.0/step/lineareos_single_interface_120min")
computed_output = joinpath(output_path, "computed_output.nc")
tracers = joinpath(output_path, "tracers.nc")
velocities = joinpath(output_path, "velocities.nc")

@info "Re-producing animations"
animate_density(computed_output, "σ",
                xslice = 17, yslice = 17, density_limit_adjustment = 0.04)
animate_tracers(tracers, xslice = 17, yslice = 17,
                S_limit_adjustment = 0.025,
                Θ_limit_adjustment = 0.5)
animate_vertical_velocity(velocities, xslice = 17, yslice = 17)

@info "Re-computing diagnostics"
diags = joinpath(output_path, "step_diagnostics.jld2")
if isfile(diags)
    rm(diags)
end
save_diagnostics!(diags, tracers, computed_output, velocities)
