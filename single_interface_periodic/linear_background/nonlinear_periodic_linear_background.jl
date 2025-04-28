using StaircaseShenanigans, GibbsSeaWater

restart = true

architecture = GPU()
Pr = 7   # Prandtl
τ = 0.05 # diff ratio
ν = 2.5e-6 # set this get the others
diffusivities = diffusivities_from_ν(ν; τ, Pr)
domain_extent = (Lx=0.05, Ly=0.05, Lz=-1.0)
domain_topology = (x = Periodic, y = Periodic, z = Bounded)
resolution = (Nx=50, Ny=50, Nz=500)
ρ₀ = gsw_rho(34.7, 0.5, 0)
eos = TEOS10EquationOfState(reference_density = ρ₀)
model_setup = (;architecture, diffusivities, domain_extent, domain_topology, resolution, eos)
dns_model = DNSModel(model_setup...; TD = VerticallyImplicitTimeDiscretization())

## Initial conditions
depth_of_interface = -0.5
salinity = [34.58, 34.70]
temperature = [-1.5, 0.5]

# ## get constant density gradient from linear salinity (i.e. temperature is nonlinear)
Sᵤ, Sₗ = salinity
S_range = range(Sᵤ, Sₗ, length = resolution.Nz)
ΔS = diff(Array(salinity))[1]
Δz = domain_extent.Lz - 0
S_bf = BackgroundField(linear_background, parameters = (Cᵤ = Sᵤ, ΔC = ΔS, Δz = Δz, Lz = 0))

Tᵤ, Tₗ = temperature

z = znodes(model.grid, Center())
p = Array(gsw_p_from_z.(z, 60))
reverse!(p)
ρᵤ, ρₗ = gsw_rho(Sᵤ, Tᵤ, p[1]), gsw_rho(Sₗ, Tₗ, p[end])
ρ_range = range(ρᵤ, ρₗ, length = resolution.Nz)

found_T = gsw_ct_from_rho.(ρ_range, S_range, p)
found_T = [found_T[i][1] for i ∈ eachindex(found_T)]
reverse!(found_T)
T_background_profile = reshape(found_T, (1, 1, resolution.Nz))
T_bf = similar(model.tracers.T)
set!(T_bf, T_background_profile)

background_fields = (S = S_bf, T = T_bf)
interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature,
                                    background_state = background_fields)

# interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature,
#                                     background_state = BackgroundLinear())
initial_noise = (velocities = VelocityNoise(1e-2), tracers = TracerNoise(1e-4, 1e-2))

## setup model
sdns = StaircaseDNS(model_setup, interface_ics, initial_noise)

## Build simulation
stop_time = Int(4 * 60 * 60) # seconds
initial_state = interface_ics.interface_smoothing isa TanhInterfaceThickness ?  "tanh" : "step"
output_path = joinpath(@__DIR__, "dns_rundown_$(round(interface_ics.R_ρ, digits = 2))", initial_state)
save_schedule = 60
checkpointer_time_interval = 60 * 60 # seconds
Δt = 1e-3
max_Δt = 7e-2
simulation = SDNS_simulation_setup(sdns, stop_time, save_computed_output!,
                                   save_vertical_velocities!;
                                   output_path,
                                   save_schedule,
                                   checkpointer_time_interval,
                                   overwrite_saved_output = restart,
                                   max_Δt,
                                   Δt)
## Run
restart ? nothing : simulation.stop_time = Int(10 * 60 * 60) # update to pickup from a checkpoint
pickup = restart ? false : readdir(simulation.output_writers[:checkpointer].dir, join = true)[1]
run!(simulation; pickup)

## Compute density ratio
compute_R_ρ!(simulation.output_writers[:computed_output].filepath,
             simulation.output_writers[:tracers].filepath, eos)

## Produce animations
reduced_path = findlast('/', simulation.output_writers[:computed_output].filepath)
output_path = simulation.output_writers[:computed_output].filepath[1:(reduced_path-1)]
cd(output_path)
@info "Producing animations"
using CairoMakie
animate_density(simulation.output_writers[:computed_output].filepath, "σ", xslice = 25, yslice = 25)
animate_tracers(simulation.output_writers[:tracers].filepath, xslice = 25, yslice = 25)
