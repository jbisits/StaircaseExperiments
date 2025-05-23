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
interface_ics = PeriodoicSingleInterfaceICs(eos, depth_of_interface, salinity, temperature,
                                            background_state = BackgroundTanh(),
                                            interface_smoothing = Tanh)
noise = (tracers = TracerNoise(1e-5, 1e-5), velocities = VelocityNoise())

## setup model
sdns = StaircaseDNS(model_setup, interface_ics, noise)

## Build simulation
Δt = 1e-3
stop_time = 3 * 60 * 60 # seconds
save_schedule = 10  # seconds
output_path = joinpath(@__DIR__, "tracer_noise")
simulation = SDNS_simulation_setup(sdns, stop_time, save_computed_output!,
                                    StaircaseShenanigans.save_vertical_velocities!;
                                    output_path, Δt, save_schedule, max_Δt = 1e-1)
