using StaircaseShenanigans, GibbsSeaWater, CairoMakie

restart = true

architecture = GPU()
Pr = 7   # Prandtl
τ = 0.05 # diff ratio
ν = 2.5e-6 # set this get the others
diffusivities = diffusivities_from_ν(ν; τ, Pr)
domain_extent = (Lx=0.05, Ly=0.05, Lz=-0.75)
domain_topology = (x = Periodic, y = Periodic, z = Bounded)
resolution = (Nx=50, Ny=50, Nz=500)
ρ₀ = gsw_rho(34.7, 0.5, 0)
eos = TEOS10EquationOfState(reference_density = ρ₀)
model_setup = (;architecture, diffusivities, domain_extent, domain_topology, resolution, eos)
model = DNSModel(model_setup...; TD = VerticallyImplicitTimeDiscretization())

number_of_interfaces = 3
depth_of_interfaces = [-0.1875, -0.375, -0.5625]
salinity = [34.58, 34.61, 34.65, 34.7]
temperature = [-1.51, -0.87, -0.19, 0.51]
staircase_ics = StaircaseICs(eos, number_of_interfaces, depth_of_interfaces, salinity, temperature)

# initial_noise = (velocities = VelocityNoise(1e-2), tracers = TracerNoise(1e-4, 1e-2))
initial_noise = TracerNoise(1e-4, 1e-2)
## setup model
sdns = StaircaseDNS(model, staircase_ics; initial_noise)

## Build simulation
stop_time = Int(4 * 60 * 60) # seconds
initial_state = "step" # can update if smoothing is added
output_path = joinpath(@__DIR__, "rundown_$(round(staircase_ics.R_ρ[2], digits = 2))", initial_state)
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
             simulation.output_writers[:tracers].filepath,
             [(-0.15, -0.05), (-0.31, -0.21), (-0.4, -0.5)],
             [(-0.31, -0.21), (-0.4, -0.5), (-0.7, -0.6)],
             eos)

## Produce animations
reduced_path = findlast('/', simulation.output_writers[:computed_output].filepath)
output_path = simulation.output_writers[:computed_output].filepath[1:(reduced_path-1)]
cd(output_path)

@info "Producing animations"
animate_density(simulation.output_writers[:computed_output].filepath, "σ",
                xslice = 25, yslice = 25, density_limit_adjustment = 0.04)
animate_tracers(simulation.output_writers[:tracers].filepath, xslice = 25, yslice = 25,
                S_limit_adjustment = 0.025,
                Θ_limit_adjustment = 0.5)
animate_vertical_velocity(simulation.output_writers[:velocities].filepath, xslice = 17, yslice = 17)

diags = initial_state*"_diagnostics.jld2"
if isfile(diags)
    rm(diags)
end
save_diagnostics!(diags,
                  simulation.output_writers[:tracers].filepath,
                  simulation.output_writers[:computed_output].filepath,
                  simulation.output_writers[:velocities].filepath)

# using JLD2, NCDatasets
# ds = NCDataset(simulation.output_writers[:computed_output].filepath)
# R_ρ = "R_rho.jld2"
# if isfile(R_ρ)
#     rm(R_ρ)
# end
# jldopen(R_ρ, "w") do f
#     f["R_ρ1"] = ds[:R_ρ1][:]
#     f["R_ρ2"] = ds[:R_ρ2][:]
#     f["R_ρ3"] = ds[:R_ρ3][:]
# end
# close(ds)
