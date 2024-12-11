using StaircaseShenanigans

restart = true

architecture = GPU()
diffusivities = (ν=1e-6, κ=(S=1e-9, T=1e-7))
domain_extent = (Lx=0.1, Ly=0.1, Lz=-1.0)
domain_topology = (x = Periodic, y = Periodic, z = Periodic)
resolution = (Nx=100, Ny=100, Nz=1000)
eos = CustomLinearEquationOfState(-0.5, 34.6)
model_setup = (;architecture, diffusivities, domain_extent, domain_topology, resolution, eos)

## Initial conditions
depth_of_interface = -0.5
salinity = [34.54, 34.70]
temperature = [-1.5, 0.5]
interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature,
                                    background_state = BackgroundTanh(50))
noise = VelocityNoise(1e-2)

## setup model
sdns = StaircaseDNS(model_setup, interface_ics, noise)

## Build simulation
stop_time = 5 * 60 * 60 # seconds
output_path = joinpath(@__DIR__, "tanh_background_velocity_noise_$(round(interface_ics.R_ρ, digits = 2))_more_stable")
checkpointer_time_interval = 60 * 60 # seconds
simulation = SDNS_simulation_setup(sdns, stop_time, save_computed_output!,
                                    save_vertical_velocities!;
                                    output_path, checkpointer_time_interval,
                                    overwrite_saved_output = restart)
## Run
# simulation.stop_time = _ * 60 * 60 # update to pickup from a checkpoint
pickup = restart ? false : joinpath(@__DIR__,
                                    "string_to_chekcpoint") # need to edit this with the actual string
run!(simulation; pickup)

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
