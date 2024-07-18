using StaircaseShenanigans, CairoMakie

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-4, κ = (S = 1e-7, T = 1e-5))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = -1.0)
resolution = (Nx = 10, Ny = 10, Nz = 100)

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

z = znodes(sdns.model.grid, Center())
S = interior(sdns.model.tracers.S, 5, 5, :)
T = interior(sdns.model.tracers.T, 5, 5, :)

fig = Figure(size = (1000, 500))
axS = Axis(fig[1, 1], title = "Initial salinity", xlabel = "Absolute salinity (gkg⁻¹)", ylabel = "z (m)")
lines!(axS, S, z)

axT = Axis(fig[1, 2], title = "Initial temperature", xlabel = "Conservative temperature (°C)", ylabel = "z (m)")
lines!(axT, T, z, color = :red)
hideydecorations!(axT, grid = false)

fig
