using StaircaseShenanigans

architecture = GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-9, T = 1e-7))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = -1.0)
resolution = (Nx = 124, Ny = 124, Nz = 1400)

## Setup the model
eos = CustomLinearEquationOfState(-0.5, 34.6)
model = DNSModel(architecture, domain_extent, resolution, diffusivities, eos)

## Initial conditions
number_of_steps = 1
depth_of_steps = [-0.5]
salinity = [34.58, 34.7]
temperature = [-1.5, 0.5]

step_ics = StepInitialConditions(model, number_of_steps, depth_of_steps, salinity, temperature)

sdns = StaircaseDNS(model, step_ics)

set_staircase_initial_conditions!(sdns)

Δt = 1e-2
stop_time = 11 * 60 * 60 # seconds
save_schedule = 60  # seconds
output_path = joinpath(@__DIR__, "linear_simulation/")
simulation = SDNS_simulation_setup(sdns, Δt, stop_time, save_schedule, save_computed_output!; output_path)

run!(simulation)
