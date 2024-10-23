using StaircaseShenanigans

architecture = GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-9, T = 1e-7))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = -1.0)
resolution = (Nx = 100, Ny = 100, Nz = 1400)

## Setup the model
eos = CustomLinearEquationOfState(-0.5, 34.6)
model = DNSModel(architecture, domain_extent, resolution, diffusivities, eos)

## Set Initial conditions
number_of_interfaces = 1
depth_of_interfaces = [-0.5]
salinity = [34.58, 34.70]
temperature = [-1.5, 0.5]

step_ics = StepInitialConditions(model, number_of_interfaces, depth_of_interfaces, salinity, temperature)

sdns = StaircaseDNS(model, step_ics)

set_staircase_initial_conditions!(sdns)

## Build simulation
Δt = 1e-2
stop_time = 6 * 60 * 60 # seconds
save_schedule = 60  # seconds
output_path = joinpath(@__DIR__, "tracer_content_restoring/")
simulation = SDNS_simulation_setup(sdns, Δt, stop_time, save_schedule, save_computed_output!,
                                    StaircaseShenanigans.no_velocities!,
                                    S_and_T_tracer_restoring_callbacks!; output_path
                                    )

## Run
run!(simulation)
