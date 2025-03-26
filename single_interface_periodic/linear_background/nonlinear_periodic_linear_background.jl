using StaircaseShenanigans, GibbsSeaWater

restart = true

architecture = GPU()
diffusivities = (ν=7e-6, κ=(S=1e-8, T=1e-6))
domain_extent = (Lx=0.1, Ly=0.1, Lz=-0.5)
domain_topology = (x = Periodic, y = Periodic, z = Periodic)
resolution = (Nx=100, Ny=100, Nz=500)
ρ₀ = gsw_rho(34.7, 0.5, 0)
eos = TEOS10EquationOfState(reference_density = ρ₀)
model_setup = (;architecture, diffusivities, domain_extent, domain_topology, resolution, eos)
model = DNSModel(model_setup...) # needed for grid

## Initial conditions
depth_of_interface = -0.25
salinity = [34.58, 34.70]
temperature = [-1.5, 0.5]

# ## stabilise with T
# Sᵤ, Sₗ = salinity
# S_range = range(Sᵤ, Sₗ, length = resolution.Nz)
# ΔS = diff(Array(salinity))[1]
# Δz = domain_extent.Lz - 0
# S_bf = BackgroundField(linear_background, parameters = (Cᵤ = Sᵤ, ΔC = ΔS, Δz = Δz, Lz = 0))

# Tᵤ, Tₗ = temperature

# z = znodes(model.grid, Center())
# p = Array(gsw_p_from_z.(z, 60))
# reverse!(p)
# ρᵤ, ρₗ = gsw_rho(Sᵤ, Tᵤ, p[1]), gsw_rho(Sₗ, Tₗ, p[end])
# ρ_range = range(ρᵤ, ρₗ, length = resolution.Nz)

# found_T = gsw_ct_from_rho.(ρ_range, S_range, p)
# found_T = [found_T[i][1] for i ∈ eachindex(found_T)]
# reverse!(found_T)
# T_background_profile = reshape(found_T, (1, 1, resolution.Nz))
# T_bf = similar(model.tracers.T)
# set!(T_bf, T_background_profile)

# background_fields = (S = S_bf, T = T_bf)
# interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature,
#                                     background_state = background_fields)

interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature,
                                    background_state = BackgroundLinear())
# noise magnitude = 0.05ΔS, 0.05ΔΘ.
# noise = (velocities = VelocityNoise(1e-4), tracers = TracerNoise(0.004, 0.05))
noise = tracers = TracerNoise(0.004, 0.05)

## setup model
sdns = StaircaseDNS(model_setup, interface_ics, noise)

## Build simulation
stop_time = 5 * 60 * 60 # seconds
output_path = joinpath(@__DIR__, "linear_background_$(round(interface_ics.R_ρ, digits = 2))")
checkpointer_time_interval = 60 * 60 # seconds
simulation = SDNS_simulation_setup(sdns, stop_time, save_computed_output!,
                                   save_vertical_velocities!; output_path,
                                   checkpointer_time_interval,
                                   overwrite_saved_output = restart)

## Run
# simulation.stop_time = 3 * 60 * 60
pickup = restart ? false : readdir(simulation.output_writers[:checkpointer].dir, join = true)[1]
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
animate_density(simulation.output_writers[:computed_output].filepath, "σ")
animate_tracers(simulation.output_writers[:tracers].filepath)
animate_density_anomaly(simulation.output_writers[:computed_output].filepath, "σ")
animate_tracers_anomaly(simulation.output_writers[:tracers].filepath)
