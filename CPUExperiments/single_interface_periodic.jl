using StaircaseShenanigans, GibbsSeaWater

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-5, κ = (S = 1e-7, T = 1e-5))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = -1.0)
domain_topology = (x = Periodic, y = Periodic, z = Periodic)
resolution = (Nx = 5, Ny = 5, Nz = 50)
eos = TEOS10EquationOfState()#CustomLinearEquationOfState(-0.5, 34.6)
model_setup = (;architecture, diffusivities, domain_extent, domain_topology, resolution, eos)
model = DNSModel(model_setup...) # needed for grid

## Initial conditions
depth_of_interface = -0.5
salinity = [34.58, 34.70]
temperature = [-1.5, 0.5]

Sᵤ, Sₗ = salinity
S_range = range(Sᵤ, Sₗ, length = resolution.Nz)
ΔS = diff(Array(salinity))[1]
Δz = domain_extent.Lz - 0
S_b = BackgroundField(linear_background, parameters = (Cᵤ = Sᵤ, ΔC = ΔS, Δz = Δz, Lz = 0))

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

background_fields = (S = S_b, T = T_bf)
interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature,
                                    background_state = background_fields)
##
noise = (velocities = VelocityNoise(1e-2), tracers = TracerNoise(0.004, 0.05))

## setup model
sdns = StaircaseDNS(model_setup, interface_ics, noise)

## Build simulation
Δt = 1e-1
stop_time = 200 * 60 # seconds
save_schedule = 60  # seconds
output_path = joinpath(@__DIR__, "tanh_background")
simulation = SDNS_simulation_setup(sdns, stop_time, save_computed_output!,
                                    save_vertical_velocities!;
                                    output_path, max_Δt = 5, Δt, save_schedule)
## Run
run!(simulation)

## Compute density ratio
compute_R_ρ!(simulation.output_writers[:computed_output].filepath,
             simulation.output_writers[:tracers].filepath, eos)
