using StaircaseShenanigans

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-9, T = 1e-7))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = -1.0)
resolution = (Nx = 100, Ny = 100, Nz = 1000)

## Setup the model
model = DNSModel(architecture, domain_extent, resolution, diffusivities)

## Initial conditions
number_of_steps = 4
depth_of_steps = [-0.2, -0.4, -0.6, -0.8]
salinity = [34.58, 34.6, 34.62, 34.64, 34.66]
temperature = [-1.5, -1.0, -0.5, 0.0, 0.5]

step_ics = StepInitialConditions(number_of_steps, depth_of_steps, salinity, temperature)

sdns = StaircaseDNS(model, step_ics)

set_staircase_initial_conditions!(sdns)

Δt = 1e-2
max_Δt = 0.1
stop_time = 5 * 60 * 60 # seconds
save_schedule = 60  # seconds
checkpointer_time_interval = 30 * 60 # seconds
output_path = joinpath(@__DIR__, "outputs_test/")
simulation = SDNS_simulation_setup(sdns, Δt, max_Δt, stop_time, save_schedule; output_path)

run!(simulation)
