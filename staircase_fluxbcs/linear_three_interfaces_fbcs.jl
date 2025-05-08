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
ρ₀ = gsw_rho(34.57, 0.5, 0)
eos = CustomLinearEquationOfState(-0.5, 34.6, reference_density = ρ₀)
model_setup = (;architecture, diffusivities, domain_extent, domain_topology, resolution, eos)
Jᵀ = 2.5e-5
T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.5*Jᵀ), bottom = FluxBoundaryCondition(0.5*Jᵀ))
Jˢ = 3.2e-7
S_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Jˢ), bottom = FluxBoundaryCondition(Jˢ))
boundary_conditions = (T=T_bcs, S=S_bcs)
model = DNSModel(model_setup...; boundary_conditions, TD = VerticallyImplicitTimeDiscretization())

number_of_interfaces = 3
depth_of_interfaces = [-0.25, -0.5, -0.75]
salinity = [34.56, 34.594, 34.64, 34.7]
temperature = [-1.5, -1.0, -0.33, 0.52]
staircase_ics = StaircaseICs(model, number_of_interfaces, depth_of_interfaces, salinity, temperature)

initial_noise = (velocities = VelocityNoise(1e-2), tracers = TracerNoise(1e-4, 1e-2))

## setup model
sdns = StaircaseDNS(model, staircase_ics; initial_noise)

## Build simulation
stop_time = Int(10 * 60 * 60) # seconds
initial_state = "step" # can update if smoothing is added
output_path = joinpath(@__DIR__, "fluxbcs_$(round(staircase_ics.R_ρ[2], digits = 2))", initial_state)
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
animate_density(simulation.output_writers[:computed_output].filepath, "σ",
                xslice = 25, yslice = 25, density_limit_adjustment = 0.04)
animate_tracers(simulation.output_writers[:tracers].filepath, xslice = 25, yslice = 25,
                S_limit_adjustment = 0.02,
                Θ_limit_adjustment = 0.25)
