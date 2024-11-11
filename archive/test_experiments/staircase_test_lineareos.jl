using StaircaseShenanigans

architecture = GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-9, T = 1e-7))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = -1.0)
resolution = (Nx = 100, Ny = 100, Nz = 1000)

## Initial conditions
number_of_steps = 4
depth_of_steps = [-0.2, -0.4, -0.6, -0.8]
salinity = [34.57, 34.60, 34.63, 34.66, 34.69]
temperature = [-1.5, -1.15, -0.8, -0.45, -0.1]

## Setup the model and initial conditions
eos = CustomLinearEquationOfState(0.0, 34.5)
model = DNSModel(architecture, domain_extent, resolution, diffusivities, eos)

step_ics = StepInitialConditions(model, number_of_steps, depth_of_steps, salinity, temperature)

sdns = StaircaseDNS(model, step_ics)

set_staircase_initial_conditions!(sdns)

## Setup the simulation
Δt = 1e-2
stop_time = 7 * 60 * 60 # seconds
save_schedule = 60  # seconds
output_path = joinpath(@__DIR__, "outputs_test/")
simulation = SDNS_simulation_setup(sdns, Δt, stop_time, save_schedule, save_computed_output!; output_path)

## and go!
run!(simulation)
