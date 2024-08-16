using StaircaseShenanigans, CairoMakie

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-4, κ = (S = 1e-7, T = 1e-5))
domain_extent = (Lx = 0.08, Ly = 0.08, Lz = -2.0)
resolution = (Nx = 8, Ny = 8, Nz = 200)

## Setup the model
eos = CustomLinearEquationOfState(0.0, 34.6)

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
temperature = [-1.5, -1.0, -0.5, 0.0, 0.5]
## stability with three, linear, R_ρ = 2.06
number_of_steps = 3
depth_of_steps = [-0.95, -1.05]
salinity = [34.495, 34.63, 34.765]
temperature = [-1.5, -0.5, 0.5]
## stability with three, nonlinear, R_ρ = 2.07
number_of_steps = 3
depth_of_steps = [-0.95, -1.05]
salinity = [34.53, 34.63, 34.7665]
temperature = [-1.5, -0.5, 0.5]

step_ics = StepInitialConditions(model, number_of_steps, depth_of_steps, salinity, temperature)

sdns = StaircaseDNS(model, step_ics)

set_staircase_initial_conditions!(sdns)

visualise_initial_conditions(sdns, 5, 5)

visualise_initial_density(sdns, 5, 5, 0)
