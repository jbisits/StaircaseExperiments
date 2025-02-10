using StaircaseShenanigans, GibbsSeaWater

architecture = GPU()
diffusivities = (ν=1e-6, κ=(S=1e-9, T=1e-7))
domain_extent = (Lx=0.07, Ly=0.07, Lz=-1.0)
domain_topology = (x = Periodic, y = Periodic, z = Bounded)
resolution = (Nx=70, Ny=70, Nz=1000)
ρ₀ = gsw_rho(34.57, 0.5, 0)
eos = CustomLinearEquationOfState(-0.5, 34.6, reference_density = ρ₀)
model_setup = (;architecture, diffusivities, domain_extent, domain_topology, resolution, eos)
# bcs from a rundown model and are an approximation/test to see if can simulate
# effect of interfaces either side.
Jᵀ = 3.5e-5
T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Jᵀ), bottom = FluxBoundaryCondition(Jᵀ))
Jˢ = 2.5e-7
S_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Jˢ), bottom = FluxBoundaryCondition(Jˢ))
boundary_conditions = (T=T_bcs, S=S_bcs)
dns_model = DNSModel(model_setup...; boundary_conditions)

## Initial conditions
depth_of_interface = -0.5
salinity = [34.54, 34.70]
temperature = [-1.5, 0.5]
interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature,
                                    interface_smoothing = TanhInterfaceSteepness(100.0))
noise = VelocityNoise(1e-2)

## setup model
sdns = StaircaseDNS(dns_model, interface_ics, initial_noise = noise)

## Build simulation
stop_time = 6 * 60 * 60 # seconds
output_path = joinpath(@__DIR__, "flux_bcs_$(round(interface_ics.R_ρ, digits = 2))")
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
animate_density(simulation.output_writers[:computed_output].filepath, "σ")
animate_tracers(simulation.output_writers[:tracers].filepath)
