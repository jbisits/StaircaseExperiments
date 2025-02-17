using StaircaseShenanigans, GibbsSeaWater

architecture = GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-9, T = 1e-7))
domain_extent = (Lx = 0.07, Ly = 0.07, Lz = -1.0)
domain_topology = (x = Periodic, y = Periodic, z = Periodic)
resolution = (Nx = 70, Ny = 70, Nz = 1000)
ρ₀ = gsw_rho(34.7, 0.5, 0)
eos = TEOS10EquationOfState(reference_density = ρ₀)
model_setup = (;architecture, diffusivities, domain_extent, domain_topology, resolution, eos)

## Initial conditions
depth_of_interface = -0.5
salinity = [34.62, 34.70]
temperature = [-0.5, 0.5]
interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature,
                                    background_state = BackgroundLinear())
# noise magnitude = 0.05ΔS, 0.05ΔΘ.
noise = (velocities = VelocityNoise(1e-2), tracers = TracerNoise(0.007, 0.1))

## setup model
sdns = StaircaseDNS(model_setup, interface_ics, noise)

## Build simulation
stop_time = 4 * 60 * 60 # seconds
output_path = joinpath(@__DIR__, "linear_background_$(round(interface_ics.R_ρ, digits = 2))")
simulation = SDNS_simulation_setup(sdns, stop_time, save_computed_output!,
                                   save_vertical_velocities!; output_path)
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
animate_density_anomaly(simulation.output_writers[:computed_output].filepath, "σ")
animate_tracers_anomaly(simulation.output_writers[:tracers].filepath)
animate_density(simulation.output_writers[:computed_output].filepath, "σ")
animate_tracers(simulation.output_writers[:tracers].filepath)
