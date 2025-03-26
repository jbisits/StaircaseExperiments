using StaircaseShenanigans

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-4, κ = (S = 1e-7, T = 1e-6))
domain_extent = (Lx = 0.01, Ly = 0.01, Lz = (-100, -99))
domain_topology = (x = Periodic, y = Periodic, z = Bounded)
resolution = (Nx = 5, Ny = 5, Nz = 50)
eos = CustomLinearEquationOfState(-0.5, 34.6)
model_setup = (;architecture, diffusivities, domain_extent, domain_topology, resolution, eos)

## Initial conditions
depth_of_interface = -99.5
salinity = [34.58, 34.70]
temperature = [-1.5, 0.5]
interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature,
                                    interface_smoothing = TanhInterfaceThickness(0.01, 0.02))
initial_noise = NoiseAtDepth([-0.27, -0.23], TracerNoise(0.0, 1e-6))

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
animate_density(simulation.output_writers[:computed_output].filepath, "σ", xslice = 25, yslice = 25)
animate_tracers(simulation.output_writers[:tracers].filepath, xslice = 25, yslice = 25)
