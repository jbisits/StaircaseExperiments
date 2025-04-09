using StaircaseShenanigans, GibbsSeaWater

restart = true

architecture = GPU()
diffusivities = (ν=1e-5, κ=(S=1.4e-8, T=1.4e-6))
domain_extent = (Lx=0.07, Ly=0.07, Lz=-0.5)
domain_topology = (x = Periodic, y = Periodic, z = Bounded)
resolution = (Nx=35, Ny=35, Nz=250)
ρ₀ = gsw_rho(34.7, 0.5, 0)
eos = TEOS10EquationOfState(reference_density = ρ₀)
model_setup = (;architecture, diffusivities, domain_extent, domain_topology, resolution, eos)
# bcs from a rundown model and are an approximation/test to see if can simulate
# effect of interfaces either side.
Jᵀ = 1.86e-5
T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Jᵀ), bottom = FluxBoundaryCondition(Jᵀ))
Jˢ = 6.36e-8
S_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Jˢ), bottom = FluxBoundaryCondition(Jˢ))
boundary_conditions = (T=T_bcs, S=S_bcs)
dns_model = DNSModel(model_setup...; boundary_conditions)

## Initial conditions
depth_of_interface = -0.25
salinity = [34.56, 34.70]
temperature = [-1.5, 0.5]
interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature)

initial_noise = NoiseAtDepth([depth_of_interface-0.02, depth_of_interface+0.02], VelocityNoise(1e-5))
## setup model
sdns = StaircaseDNS(dns_model, interface_ics; initial_noise)

## Build simulation
stop_time = 2 * 60 * 60 # seconds
initial_state = interface_ics.interface_smoothing isa TanhInterfaceThickness ?  "tanh" : "step"
output_path = joinpath(@__DIR__, "fluxbcs_$(round(interface_ics.R_ρ, digits = 2))", initial_state)
save_schedule = 30
checkpointer_time_interval = 60 * 60 # seconds
Δt = 1e-4
max_Δt = 2e-3
simulation = SDNS_simulation_setup(sdns, stop_time, save_computed_output!,
                                   save_vertical_velocities!;
                                   output_path,
                                   save_schedule,
                                   checkpointer_time_interval,
                                   overwrite_saved_output = restart,
                                   max_Δt,
                                   Δt)
## Run
# simulation.stop_time = 18 * 60 * 60 # update to pickup from a checkpoint
pickup = restart ? false : readdir(simulation.output_writers[:checkpointer].dir, join = true)[1]
run!(simulation; pickup)

## Produce animations
reduced_path = findlast('/', simulation.output_writers[:computed_output].filepath)
animation_path = simulation.output_writers[:computed_output].filepath[1:(reduced_path-1)]
cd(animation_path)
@info "Producing animations"
using CairoMakie
animate_density(simulation.output_writers[:computed_output].filepath, "σ", xslice = 17, yslice = 17)
animate_tracers(simulation.output_writers[:tracers].filepath, xslice = 17, yslice = 17)
