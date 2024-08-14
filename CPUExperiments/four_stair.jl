using StaircaseShenanigans, CairoMakie

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-4, κ = (S = 1e-7, T = 1e-5))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = -1.0)
resolution = (Nx = 10, Ny = 10, Nz = 100)

## Setup the model
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
set!(model, w = noise)

## Initial conditions
number_of_steps = 4
depth_of_steps = [-0.2, -0.4, -0.6, -0.8]
## Nonlinear eos
salinity = [34.57, 34.60, 34.63, 34.66, 34.69]
temperature = [-1.5, -1.0, -0.576, -0.204, 0.133]
## Custom linear eos
salinity = [34.57, 34.60, 34.63, 34.66, 34.69]
temperature = [-1.5, -1.15, -0.8, -0.45, -0.1]
## stability with three
number_of_steps = 3
depth_of_steps = [-0.4, -0.6]
salinity = [34.57, 34.63, 34.69]
temperature = [-1.5, -0.8, -0.1]

step_ics = StepInitialConditions(model, number_of_steps, depth_of_steps, salinity, temperature)

sdns = StaircaseDNS(model, step_ics)

set_staircase_initial_conditions!(sdns)

visualise_initial_conditions(sdns, 5, 5)

visualise_initial_density(sdns, 5, 5, 0)
