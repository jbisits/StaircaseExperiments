using StaircaseShenanigans

architecture = GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-9, T = 1e-7))
domain_extent = (Lx = 0.08, Ly = 0.08, Lz = -2.0)
resolution = (Nx = 80, Ny = 80, Nz = 2000)

## Initial conditions
number_of_steps = 5
depth_of_steps = [-0.8, -0.9, -1.0, -1.1]
salinity = [34.54, 34.60, 34.66, 34.72, 34.78]
temperature = [-1.5, -1.0, -0.5, 0.0, 0.5]

## Setup the model and initial conditions
eos = CustomLinearEquationOfState(0.0, 34.6)
model = DNSModel(architecture, domain_extent, resolution, diffusivities, eos)

step_ics = StepInitialConditions(model, number_of_steps, depth_of_steps, salinity, temperature)

sdns = StaircaseDNS(model, step_ics)

set_staircase_initial_conditions!(sdns)

## Setup the simulation
Δt = 1e-2
stop_time = 11 * 60 * 60 # seconds
save_schedule = 60  # seconds
output_path = joinpath(@__DIR__, "stability_test/")
simulation = SDNS_simulation_setup(sdns, Δt, stop_time, save_schedule, save_computed_output!; output_path)

## and go!
run!(simulation)
