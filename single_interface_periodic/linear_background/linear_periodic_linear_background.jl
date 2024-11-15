using StaircaseShenanigans

architecture = GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-9, T = 1e-7))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = -1.0)
resolution = (Nx = 100, Ny = 100, Nz = 1000)
eos = CustomLinearEquationOfState(-0.5, 34.6)
model_setup = (;architecture, diffusivities, domain_extent, resolution, eos)

## Initial conditions
depth_of_interface = -0.5
salinity = [34.56, 34.70]
temperature = [-1.5, 0.5]
interface_ics = PeriodoicSingleInterfaceICs(eos, depth_of_interface, salinity, temperature, BackgroundLinear())
tracer_noise = TracerNoise(1e-8, 1e-8)

## setup model
sdns = StaircaseDNS(model_setup, interface_ics, tracer_noise)

## Build simulation
stop_time = 3 * 60 * 60 # seconds
output_path = joinpath(@__DIR__, "tracer_noise")
simulation = SDNS_simulation_setup(sdns, stop_time, save_computed_output!,
                                   save_vertical_velocities!; output_path)
## Run
run!(simulation)

## Compute density ratio
compute_R_ρ!(simulation.output_writers[:computed_output].filepath,
             simulation.output_writers[:tracers].filepath, eos)

## Produce animations
cd(output_path)
@info "Producing animations"
using CairoMakie
animate_density_anomaly(simulation.output_writers[:computed_output].filepath, "σ")
animate_tracers_anomaly(simulation.output_writers[:tracers].filepath)
animate_density(simulation.output_writers[:computed_output].filepath, "σ")
animate_tracers(simulation.output_writers[:tracers].filepath)
