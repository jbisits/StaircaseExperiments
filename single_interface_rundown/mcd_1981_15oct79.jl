using StaircaseShenanigans, GibbsSeaWater

## Initial salinity and temperature from McDougall 1981 lab experiments
Sₚ = [0.03, 6.02]
T = [18.57, 31.58]
salinity = gsw_sa_from_sp.(Sₚ, 0, 149, -35)
temperature = gsw_ct_from_t.(salinity, T, 0)
# ρ₀ = gsw_rho.(salinity, temperature, 0)
# Δρ = diff(ρ₀) # = 1.1239926439408237

Sₗ = gsw_sa_from_sp(6.02, 0, 149, -35)
Θₗ = gsw_ct_from_t(Sₗ, 31.58, 0)
ρ₀ = gsw_rho(Sₗ, Θₗ, 0)
##
architecture = GPU()
diffusivities = (ν=1e-6, κ=(S=1e-9, T=1e-7))
domain_extent = (Lx=0.1, Ly=0.1, Lz=-0.5)
domain_topology = (x = Periodic, y = Periodic, z = Bounded)
resolution = (Nx=100, Ny=100, Nz=500)
eos = TEOS10EquationOfState(reference_density = ρ₀)
model_setup = (;architecture, diffusivities, domain_extent, domain_topology, resolution, eos)
dns_model = DNSModel(model_setup...)

## Initial conditions
depth_of_interface = -0.25
interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature)
# α = gsw_alpha.(salinity, temperature, 0)
# δ = ((α[1] - α[2]) / (α[1] + α[2])) * (1 / (1 - interface_ics.R_ρ)) # = 0.7662374347402843
noise = VelocityNoise(1e-2)

## setup model
sdns = StaircaseDNS(dns_model, interface_ics, initial_noise = noise)

## Build simulation
stop_time = 4 * 60 * 60 # seconds
output_path = joinpath(@__DIR__, "McDougall1981_15oct_$(round(interface_ics.R_ρ, digits = 2))")
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
