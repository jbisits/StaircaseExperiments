using StaircaseShenanigans, GibbsSeaWater, CairoMakie, JLD2

restart = true

architecture = GPU()
Pr = 7   # Prandtl
τ = 0.1 # diff ratio
ν = 2.5e-6 # set this get the others
diffusivities = diffusivities_from_ν(ν; τ, Pr)
domain_extent = (Lx=0.05, Ly=0.05, Lz=-0.5)
domain_topology = (x = Periodic, y = Periodic, z = Bounded)
resolution = (Nx=110, Ny=110, Nz=1100)
ρ₀ = gsw_rho(34.7, 0.5, 0)
eos = CustomLinearEquationOfState(-0.5, 34.6, reference_density = ρ₀)
model_setup = (;architecture, diffusivities, domain_extent, domain_topology, resolution, eos)
# bcs from a rundown model and are an approximation/test to see if can simulate
# effect of interfaces either side.
Jᵀ = 7.5e-6
Jˢ = 0.45 * 2.384e-7
T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Jᵀ), bottom = FluxBoundaryCondition(Jᵀ))
S_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Jˢ), bottom = FluxBoundaryCondition(Jˢ))
boundary_conditions = (T=T_bcs, S=S_bcs)
dns_model = DNSModel(model_setup...; boundary_conditions, TD = VerticallyImplicitTimeDiscretization())

## Initial conditions
depth_of_interface = -0.25
salinity = [34.64, 34.70]
Tᵤ, Tₗ = -0.5, 0.5
ΔT = Tᵤ - Tₗ
temperature = [Tᵤ, Tₗ]
interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature)

initial_noise = tracers = TracerNoise(1e-4, 1e-2)

## setup model
sdns = StaircaseDNS(dns_model, interface_ics; initial_noise)

## Build simulation
stop_time = Int(4 * 60 * 60) # seconds
initial_state = interface_ics.interface_smoothing isa TanhInterfaceThickness ?  "tanh" : "step"
output_path = joinpath(@__DIR__, "dns_fluxbcs_$(round(interface_ics.R_ρ, digits = 2))_dT_$(ΔT)", initial_state)
save_schedule = 60
checkpointer_time_interval = 60 * 60 # seconds
Δt = 1e-3
max_Δt = 1.75e-2
simulation = SDNS_simulation_setup(sdns, stop_time, save_computed_output!,
                                   save_vertical_velocities!;
                                   output_path,
                                   save_schedule,
                                   checkpointer_time_interval,
                                   overwrite_saved_output = restart,
                                   max_Δt,
                                   Δt)
## Run
# simulation.stop_time = 6 * 60 * 60 # update to pickup from a checkpoint
pickup = restart ? false : readdir(simulation.output_writers[:checkpointer].dir, join = true)[1]
run!(simulation; pickup)

## Produce animations
compute_R_ρ!(simulation.output_writers[:computed_output].filepath,
             simulation.output_writers[:tracers].filepath, (-0.2, -0.1), (-0.4, -0.3), eos)

reduced_path = findlast('/', simulation.output_writers[:computed_output].filepath)
animation_path = simulation.output_writers[:computed_output].filepath[1:(reduced_path-1)]
cd(animation_path)
@info "Producing animations"
animate_density(simulation.output_writers[:computed_output].filepath, "σ",
                xslice = 17, yslice = 17, density_limit_adjustment = 0.04)
animate_tracers(simulation.output_writers[:tracers].filepath, xslice = 17, yslice = 17,
                S_limit_adjustment = 0.025,
                Θ_limit_adjustment = 0.5)
animate_vertical_velocity(simulation.output_writers[:velocities].filepath, xslice = 17, yslice = 17)

## compute diagnostics
diags = initial_state*"_diagnostics.jld2"
if isfile(diags)
    rm(diags)
end
save_diagnostics!(diags,
                  simulation.output_writers[:tracers].filepath,
                  simulation.output_writers[:computed_output].filepath,
                  simulation.output_writers[:velocities].filepath)
jldopen(diags, "a+") do f
    f["FluxBCs/Jˢ"] = Jˢ
    f["FluxBCs/Jᵀ"] = Jᵀ
end
