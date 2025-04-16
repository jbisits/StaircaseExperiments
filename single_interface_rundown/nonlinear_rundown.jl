using StaircaseShenanigans, GibbsSeaWater

restart = true

architecture = GPU()
Pr = 7   # Prandtl
τ = 0.1 # diff ratio
ν = 2.5e-6 # set this get the others
κₜ, κₛ = κₛ_and_κₜ_from_ν(ν; τ, Pr)
Sc = ν / κₛ
diffusivities = (;ν, κ=(S=κₛ, T=κₜ))
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
interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature)

# initial_noise = NoiseAtDepth([depth_of_interface-0.02, depth_of_interface+0.02], TracerNoise(2e-4, 0.0))
initial_noise = (velocities = VelocityNoise(1e-2), tracers = TracerNoise(1e-4, 1e-2))
## setup model
sdns = StaircaseDNS(dns_model, interface_ics; initial_noise)

## Build simulation
stop_time = Int(1 * 60 * 60) # seconds
initial_state = interface_ics.interface_smoothing isa TanhInterfaceThickness ?  "tanh" : "step"
output_path = joinpath(@__DIR__, "rundown_$(round(interface_ics.R_ρ, digits = 2))", initial_state)
save_schedule = 30
checkpointer_time_interval = 60 * 60 # seconds
Δt = 1e-3
max_Δt = 1e-2
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

## compute diagnostics
diags = "diagnostics.jld2"
if isfile(diags)
    rm(diags)
end
save_diagnostics!(diags,
                  simulation.output_writers[:tracers].filepath,
                  simulation.output_writers[:computed_output].filepath,
                  simulation.output_writers[:velocities].filepath)
