using StaircaseShenanigans

architecture = GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-9, T = 1e-7))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = -1.0)
resolution = (Nx = 100, Ny = 100, Nz = 1000)

## Initial conditions
number_of_steps = 3
depth_of_steps = [-0.4, -0.6]
salinity = [34.57, 34.63, 34.69]
temperature = [-1.5, -0.8, -0.1]

## Setup the model and initial conditions
eos = CustomLinearEquationOfState(0.0, 34.5)

rate = 1/(5 * 60)
mask = OuterStairMask(-0.15, -0.85)
S_target = OuterStairTargets(-0.15, 34.57, -0.85, 34.69)
S_relaxation = Relaxation(; rate, mask, target = S_target)
T_target = OuterStairTargets(-0.15, -1.5, -0.85, -0.1)
T_relaxation = Relaxation(; rate, mask, target = T_target)
relaxation = (S = S_relaxation, T = T_relaxation)

model = DNSModel(architecture, domain_extent, resolution, diffusivities, eos)
noise = 1e-6 * randn(size(model.velocities.w))
set!(model, w = noise) # to kick off, will move to constructor soon

step_ics = StepInitialConditions(model, number_of_steps, depth_of_steps, salinity, temperature)

sdns = StaircaseDNS(model, step_ics)

set_staircase_initial_conditions!(sdns)

## Setup the simulation
Δt = 1e-2
stop_time = 6 * 60 * 60 # seconds
save_schedule = 60  # seconds
output_path = joinpath(@__DIR__, "stability_test/")
simulation = SDNS_simulation_setup(sdns, Δt, stop_time, save_schedule, save_computed_output!; output_path)

## and go!
run!(simulation)
