using StaircaseShenanigans, CairoMakie

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-4, κ = (S = 1e-7, T = 1e-5))
domain_extent = (Lx = 0.08, Ly = 0.08, Lz = -1.0)
resolution = (Nx = 8, Ny = 8, Nz = 200)

## Setup the model
eos = CustomLinearEquationOfState(2.0, 34.7)

rate = 1/(5 * 60)
mask = OuterStairMask(-0.15, -0.85)
S_target = OuterStairTargets(-0.15, 34.57, -0.85, 34.69)
S_relaxation = Relaxation(; rate, mask, target = S_target)
T_target = OuterStairTargets(-0.15, -1.5, -0.85, -0.1)
T_relaxation = Relaxation(; rate, mask, target = T_target)
relaxation = (S = S_relaxation, T = T_relaxation)

model = DNSModel(architecture, domain_extent, resolution, diffusivities)

## Initial conditions

#### Four stairs
number_of_steps = 4
depth_of_steps = [-0.95, -1.05, -1.15]
## Nonlinear eos
salinity = [34.5585, 34.6, 34.65, 34.7083]
temperature = [-1.5, -1.0, -0.5, 0.0]
## Custom linear eos, R_ρ = 1.52
salinity = [34.54, 34.60, 34.66, 34.72]
temperature = [-1.5, -1.0, -0.5, 0.0]

#### Three stairs
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
## stability with three, different R_ρ
number_of_steps = 3
depth_of_steps = [-0.95, -1.05]
salinity = [34.495, 34.545, 34.765]
temperature = [-1.5, -1.0, 0.5]
## stability with three, nonlinear, R_ρ = 2.07
number_of_steps = 3
depth_of_steps = [-0.95, -1.05]
salinity = [34.53, 34.57, 34.7665]
temperature = [-1.5, -1.0, 0.5]

#### Five stairs
number_of_steps = 5
depth_of_steps = [-0.8, -0.9, -1.0, -1.1]
## Nonlinear eos
salinity = [34.5585, 34.6, 34.65, 34.7083, 34.775]
temperature = [-1.5, -1.0, -0.5, 0.0, 0.5]
## Custom linear eos, R_ρ = 1.52
salinity = [34.54, 34.60, 34.66, 34.72, 34.78]
temperature = [-1.5, -1.0, -0.5, 0.0, 0.5]

## Salt fingering
number_of_steps = 1
depth_of_steps = [-0.5]
salinity = [34.9, 34.7]
temperature = [2.8, 0.5]

step_ics = StepInitialConditions(model, number_of_steps, depth_of_steps, salinity, temperature)

sdns = StaircaseDNS(model, step_ics)

set_staircase_initial_conditions!(sdns)

visualise_initial_conditions(sdns, 5, 5)

visualise_initial_density(sdns, 5, 5, 0)
