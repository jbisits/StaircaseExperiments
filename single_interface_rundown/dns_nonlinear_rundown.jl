using StaircaseShenanigans, GibbsSeaWater, CairoMakie

restart = false

architecture = GPU()
Pr = 7   # Prandtl
τ = 0.1 # diff ratio
ν = 2.5e-6 # set this get the others
diffusivities = diffusivities_from_ν(ν; τ, Pr)
domain_extent = (Lx=0.05, Ly=0.05, Lz=-0.5)
domain_topology = (x = Periodic, y = Periodic, z = Bounded)
resolution = (Nx=110, Ny=110, Nz=1100) # DNS resolution
ρ₀ = gsw_rho(34.7, 0.5, 0)
eos = TEOS10EquationOfState(reference_density = ρ₀)
model_setup = (;architecture, diffusivities, domain_extent, domain_topology, resolution, eos)
dns_model = DNSModel(model_setup...; TD = VerticallyImplicitTimeDiscretization())

## Initial conditions
depth_of_interface = -0.25
salinity = [34.631, 34.70]
Tᵤ, Tₗ = -0.5, 0.5
ΔT = Tᵤ - Tₗ
temperature = [Tᵤ, Tₗ]
interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature)

# initial_noise = (velocities = VelocityNoise(1e-2), tracers = TracerNoise(1e-4, 1e-2))
initial_noise = TracerNoise(1e-4, 1e-2)
## setup model
sdns = StaircaseDNS(dns_model, interface_ics; initial_noise)

## Build simulation
stop_time = Int(2 * 60 * 60) # seconds
initial_state = interface_ics.interface_smoothing isa TanhInterfaceThickness ?  "tanh" : "step"
output_path = joinpath(@__DIR__, "dns_rundown_$(round(interface_ics.R_ρ, digits = 2))_dT_$(ΔT)", initial_state)
save_schedule = 60
checkpointer_time_interval = 60 * 60 # seconds
Δt = 1e-3
max_Δt = 1.75e-2 # DNS timesetp
simulation = SDNS_simulation_setup(sdns, stop_time, save_computed_output!,
                                   save_vertical_velocities!;
                                   output_path,
                                   save_schedule,
                                   checkpointer_time_interval,
                                   overwrite_saved_output = restart,
                                   max_Δt,
                                   Δt)
## Run
restart ? nothing : simulation.stop_time = Int(5 * 60 * 60) # update to pickup from a checkpoint
pickup = restart ? false : readdir(simulation.output_writers[:checkpointer].dir, join = true)[1]
run!(simulation; pickup)

## Compute density ratio
compute_R_ρ!(simulation.output_writers[:computed_output].filepath,
             simulation.output_writers[:tracers].filepath, (-0.2, -0.1), (-0.4, -0.3), eos)

## Produce animations
reduced_path = findlast('/', simulation.output_writers[:computed_output].filepath)
output_path = simulation.output_writers[:computed_output].filepath[1:(reduced_path-1)]
cd(output_path)
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

# # check on the output
# output_path = joinpath(@__DIR__, "rundown_1.05", "step", "smallerdt_nonlineareos")
# output = jldopen(joinpath(output_path, "step_diagnostics.jld2"))
# keys(output)
# heatmap(output["T_ha"]', colormap = :thermal)
# heatmap(output["σ_ha"]', colormap = :dense)
# heatmap(output["N²_ha"]', colormap = :balance)
