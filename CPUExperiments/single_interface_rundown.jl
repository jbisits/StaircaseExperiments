using StaircaseShenanigans

architecture = CPU() # or GPU()
Pr = 7   # Prandtl
τ = 0.05 # diff ratio
ν = 2.5e-6 # set this get the others
diffusivities = diffusivities_from_ν(ν; τ, Pr)
diffusivities = (;ν, κ=(S=1e-7, T=1e-5))
domain_extent = (Lx = 0.01, Ly = 0.01, Lz=-1.0)
domain_topology = (x = Periodic, y = Periodic, z = Bounded)
resolution = (Nx = 5, Ny = 5, Nz = 50)
eos = CustomLinearEquationOfState(-0.5, 34.6)
model_setup = (;architecture, diffusivities, domain_extent, domain_topology, resolution, eos)

## Initial conditions
depth_of_interface = -0.5
salinity = [34.58, 34.70]
temperature = [-1.5, 0.5]
interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature,
                                    interface_smoothing = TanhInterfaceThickness(0.01, 0.02))
initial_noise = nothing

## setup model
sdns = StaircaseDNS(model_setup, interface_ics, initial_noise)

## Build simulation
stop_time = 2 * 60 * 60 # seconds
save_schedule = 60  # seconds
output_path = joinpath(@__DIR__, "output")
checkpointer_time_interval = 60 * 60 # seconds
simulation = SDNS_simulation_setup(sdns, stop_time, save_computed_output!,
                                    save_vertical_velocities!;
                                    checkpointer_time_interval,
                                    output_path, save_schedule)
## Run
run!(simulation)

## Compute density ratio
compute_R_ρ!(simulation.output_writers[:computed_output].filepath,
             simulation.output_writers[:tracers].filepath, eos)

## Produce animations
reduced_path = findlast('/', simulation.output_writers[:computed_output].filepath)
animation_path = simulation.output_writers[:computed_output].filepath[1:(reduced_path-1)]
cd(animation_path)
@info "Producing animations"
using CairoMakie
animate_density(simulation.output_writers[:computed_output].filepath, "σ", xslice = 2, yslice = 2)
animate_tracers(simulation.output_writers[:tracers].filepath, xslice = 2, yslice = 2)


# mwe for differing diffusivities
using Oceananigans, CairoMakie

## Setting S diffusivity first
grid = RectilinearGrid(CPU(),
                       topology = (Periodic, Flat, Bounded),
                       size = (5, 50),
                       x = (-0.01/2, 0.01/2),
                       z = (-1, 0))

buoyancy = SeawaterBuoyancy()
tracers = (:S, :T)

closure = ScalarDiffusivity(; ν=1e-4, κ=(S=1e-7, T=1e-5))

model = NonhydrostaticModel(; grid, buoyancy, tracers, closure)
Tᵤ, Tₗ = -1.5, 0.5
T₀(x, z) = z <= -0.5 ? Tᵤ : Tₗ
Sᵤ, Sₗ = 34.58, 34.57
S₀(x, z) = z <= -0.5 ? Sᵤ : Sₗ
set!(model, S = S₀, T = T₀)

simulation = Simulation(model, Δt=0.1, stop_time=60*60)
wizard = TimeStepWizard()
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

tracers = Dict("S" => model.tracers.S, "T" => model.tracers.T)
simulation.output_writers[:tracers] = JLD2Writer(model, tracers;
                                                filename = "tracers_κₛ_κₜ",
                                                schedule = TimeInterval(60))

run!(simulation)

S = FieldTimeSeries("tracers_κₛ_κₜ.jld2", "S")
T = FieldTimeSeries("tracers_κₛ_κₜ.jld2", "T")
x = xnodes(S.grid, Center())
z = znodes(S.grid, Center())

fig = Figure(size = (1000, 1000))
axS₀ = Axis(fig[1, 1], ylabel = "z", title = "Salinity t = 0")
heatmap!(axS₀, x, z, interior(S, :, 1, :, 1), colormap = :haline)
axT₀ = Axis(fig[1, 2], ylabel = "z", title = "Temperature t = 0")
heatmap!(axT₀, x, z, interior(T, :, 1, :, 1), colormap = :thermal)
axSend = Axis(fig[2, 1], xlabel = "x", ylabel = "z", title = "Salinity t = 61")
heatmap!(axSend, x, z, interior(S, :, 1, :, 61), colormap = :haline)
axTend = Axis(fig[2, 2], xlabel = "x", ylabel = "z", title = "Temperature t = 61")
heatmap!(axTend, x, z, interior(T, :, 1, :, 61), colormap = :thermal)
fig
save("Sdiff_set_first.png", fig)

## Setting T diffusivity first
grid = RectilinearGrid(CPU(),
                       topology = (Periodic, Flat, Bounded),
                       size = (5, 50),
                       x = (-0.01/2, 0.01/2),
                       z = (-1, 0))

buoyancy = SeawaterBuoyancy()
tracers = (:S, :T)

closure = ScalarDiffusivity(; ν=1e-4, κ=(T=1e-5, S=1e-7)) # T set first

model = NonhydrostaticModel(; grid, buoyancy, tracers, closure)
Tᵤ, Tₗ = -1.5, 0.5
T₀(x, z) = z <= -0.5 ? Tᵤ : Tₗ
Sᵤ, Sₗ = 34.58, 34.57
S₀(x, z) = z <= -0.5 ? Sᵤ : Sₗ
set!(model, S = S₀, T = T₀)

simulation = Simulation(model, Δt=0.1, stop_time=60*60)
wizard = TimeStepWizard()
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

tracers = Dict("S" => model.tracers.S, "T" => model.tracers.T)
simulation.output_writers[:tracers] = JLD2Writer(model, tracers;
                                                filename = "tracers_κₜ_κₛ",
                                                schedule = TimeInterval(60))

run!(simulation)

S = FieldTimeSeries("tracers_κₜ_κₛ.jld2", "S")
T = FieldTimeSeries("tracers_κₜ_κₛ.jld2", "T")
x = xnodes(S.grid, Center())
z = znodes(S.grid, Center())

fig = Figure(size = (1000, 1000))
axS₀ = Axis(fig[1, 1], ylabel = "z", title = "Salinity t = 0")
heatmap!(axS₀, x, z, interior(S, :, 1, :, 1), colormap = :haline)
axT₀ = Axis(fig[1, 2], ylabel = "z", title = "Temperature t = 0")
heatmap!(axT₀, x, z, interior(T, :, 1, :, 1), colormap = :thermal)
axSend = Axis(fig[2, 1], xlabel = "x", ylabel = "z", title = "Salinity t = 61")
heatmap!(axSend, x, z, interior(S, :, 1, :, 61), colormap = :haline)
axTend = Axis(fig[2, 2], xlabel = "x", ylabel = "z", title = "Temperature t = 61")
heatmap!(axTend, x, z, interior(T, :, 1, :, 61), colormap = :thermal)
fig
save("Tdiff_set_first.png", fig)
