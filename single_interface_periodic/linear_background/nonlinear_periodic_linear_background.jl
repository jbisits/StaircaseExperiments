using StaircaseShenanigans, GibbsSeaWater
using SeawaterPolynomials: total_density

restart = true

architecture = GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-9, T = 1e-7))
domain_extent = (Lx = 0.07, Ly = 0.07, Lz = -1.0)
domain_topology = (x = Periodic, y = Periodic, z = Periodic)
resolution = (Nx = 70, Ny = 70, Nz = 1000)
ρ₀ = gsw_rho(34.7, 0.5, 0)
eos = TEOS10EquationOfState(reference_density = ρ₀)
model_setup = (;architecture, diffusivities, domain_extent, domain_topology, resolution, eos)
model = DNSModel(model_setup...) # needed for grid

## Initial conditions
depth_of_interface = -0.5
salinity = [34.58, 34.70]
temperature = [-1.5, 0.5]
interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature)

Sᵤ, Sₗ = salinity
S_range = range(Sᵤ, Sₗ, length = resolution.Nz)
ΔS = diff(Array(salinity))[1]
S_b = BackgroundField(linear_background, parameters = (Cᵤ = Sᵤ, ΔC = ΔS, Lz = abs(domain_extent.Lz)))

Tᵤ, Tₗ = temperature

ρᵤ, ρₗ = total_density(Tᵤ, Sᵤ, 0, eos), total_density(Tₗ, Sₗ, 0, eos)
ρ_range = range(ρᵤ, ρₗ, length = resolution.Nz)

Δρ = diff(Array(ρ_range))[1]

N_T = 100000
T_range = range(Tᵤ, Tₗ, length = N_T)
find = [findfirst(total_density.(T_range, S_range[n + 1], 0, fill(eos, length(N_T))) .- (ρᵤ + n * Δρ).≤ Δρ) for n in 0:resolution.Nz-1]
T_background_profile = reshape(reverse(T_range[find]), (1, 1, resolution.Nz))
T_bf = similar(model.tracers.T)
set!(T_bf, T_background_profile)

background_fields = (S = S_b, T = T_bf)

model = DNSModel(model_setup...; background_fields) # overwrite with background_fields

# noise magnitude = 0.05ΔS, 0.05ΔΘ.
noise = (velocities = VelocityNoise(1e-2), tracers = TracerNoise(0.004, 0.05))

## setup model
sdns = StaircaseDNS(model, interface_ics, noise)

## Build simulation
stop_time = 1 * 60 * 60 # seconds
output_path = joinpath(@__DIR__, "linear_background_$(round(interface_ics.R_ρ, digits = 2))")
checkpointer_time_interval = 60 * 60 # seconds
simulation = SDNS_simulation_setup(sdns, stop_time, save_computed_output!,
                                   save_vertical_velocities!; output_path,
                                   checkpointer_time_interval,
                                   overwrite_saved_output = restart)

## Run
# simulation.stop_time = 12 * 60 * 60
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
animate_density_anomaly(simulation.output_writers[:computed_output].filepath, "σ")
animate_tracers_anomaly(simulation.output_writers[:tracers].filepath)
animate_density(simulation.output_writers[:computed_output].filepath, "σ")
animate_tracers(simulation.output_writers[:tracers].filepath)
