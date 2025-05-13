using StaircaseShenanigans, GibbsSeaWater, CairoMakie

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
T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.15*Jᵀ),
                                bottom = FluxBoundaryCondition(0.25*Jᵀ))
Jˢ = 3.2e-7
S_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.15*Jˢ),
                                bottom = FluxBoundaryCondition(0.25*Jˢ))
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
stop_time = Int(16 * 60 * 60) # seconds
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
             simulation.output_writers[:tracers].filepath,
             [(-0.2, -0.1), (-0.45, -0.35), (-0.7, -0.6)],
             [(-0.45, -0.35), (-0.7, -0.6), (-0.9, -0.8)],
             eos)

## Produce animations
reduced_path = findlast('/', simulation.output_writers[:computed_output].filepath)
output_path = simulation.output_writers[:computed_output].filepath[1:(reduced_path-1)]
cd(output_path)
@info "Producing animations"
animate_density(simulation.output_writers[:computed_output].filepath, "σ",
                xslice = 25, yslice = 25, density_limit_adjustment = 0.04)
animate_tracers(simulation.output_writers[:tracers].filepath, xslice = 25, yslice = 25,
                S_limit_adjustment = 0.02,
                Θ_limit_adjustment = 0.25)

using JLD2, NCDatasets
ds = NCDataset(simulation.output_writers[:computed_output].filepath)
R_ρ = "R_rho.jld2"
if isfile(R_ρ)
    rm(R_ρ)
end
jldopen(R_ρ, "w") do f
    f["R_ρ1"] = ds[:R_ρ1][:]
    f["R_ρ2"] = ds[:R_ρ2][:]
    f["R_ρ3"] = ds[:R_ρ3][:]
end
close(ds)

# using CairoMakie, JLD2
# output_path = joinpath(@__DIR__, "nonlinear_weaker_top_fluxes")
# data = joinpath(output_path, "R_rho.jld2")
# f = jldopen(data)
# R_ρ1 = f["R_ρ1"]
# R_ρ2 = f["R_ρ2"]
# R_ρ3 = f["R_ρ3"]
# close(f)

# fig, ax = lines(R_ρ1, label = "Interface 1")
# lines!(ax, R_ρ2, label = "Interface 2")
# lines!(ax, R_ρ3, label = "Interface 3")
# ax.title = "R_ρ with fluxbcs nonlinear eos"
# ax.xlabel = "time (mins)"
# ax.ylabel = "R_ρ"
# axislegend(ax)
# fig
# save(joinpath(output_path, "R_rho.png"), fig)
