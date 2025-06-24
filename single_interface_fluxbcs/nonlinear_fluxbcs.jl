using StaircaseShenanigans, GibbsSeaWater, CairoMakie, JLD2

restart = true

architecture = GPU()
Pr = 7   # Prandtl
τ = 0.1 # diff ratio
ν = 2.5e-6 # set this get the others
diffusivities = diffusivities_from_ν(ν; τ, Pr)
domain_extent = (Lx=0.05, Ly=0.05, Lz=-0.5)
domain_topology = (x = Periodic, y = Periodic, z = Bounded)
resolution = (Nx=50, Ny=50, Nz=750)
ρ₀ = gsw_rho(34.7, 0.5, 0.0)
eos = TEOS10EquationOfState(reference_density = ρ₀)
model_setup = (;architecture, diffusivities, domain_extent, domain_topology, resolution, eos)
# bcs from a rundown model and are an approximation/test to see if can simulate
# effect of interfaces either side.
# ΔT = -1
Jᵀ = 6.4e-6
Jˢ = 1.272e-7
# # ΔT = -0.5
# Jᵀ = 2.6e-6
# Jˢ = 4.77e-8
top_salinity_scale = 0.5
T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Jᵀ), bottom = FluxBoundaryCondition(Jᵀ))
S_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(top_salinity_scale*Jˢ), bottom = FluxBoundaryCondition(Jˢ))
boundary_conditions = (T=T_bcs, S=S_bcs)
dns_model = DNSModel(model_setup...; boundary_conditions, TD = VerticallyImplicitTimeDiscretization())

## Initial conditions
depth_of_interface = -0.25
# dT = 1 + Rᵨ = 1.05
salinity = [34.631, 34.70]
Tᵤ, Tₗ = -0.5, 0.5
# # dT = 0.5 + Rᵨ = 1.05
# salinity = [34.663, 34.70]
# Tᵤ, Tₗ = 0.0, 0.5
ΔT = Tᵤ - Tₗ
temperature = [Tᵤ, Tₗ]
interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature)

# initial_noise = (velocities = VelocityNoise(1e-2), tracers = TracerNoise(1e-4, 1e-2))
initial_noise = TracerNoise(1e-4, 1e-2)
## setup model
sdns = StaircaseDNS(dns_model, interface_ics; initial_noise)

## Build simulation
stop_time = Int(6 * 60 * 60) # seconds
initial_state = interface_ics.interface_smoothing isa TanhInterfaceThickness ?  "tanh" : "step"
output_path = joinpath(@__DIR__, "fluxbcs_$(round(interface_ics.R_ρ, digits = 2))_dT_$(ΔT)", initial_state)
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

diags = initial_state*"_diagnostics.jld2"
if isfile(diags)
    rm(diags)
end
save_diagnostics!(diags,
                  simulation.output_writers[:tracers].filepath,
                  simulation.output_writers[:computed_output].filepath,
                  simulation.output_writers[:velocities].filepath)
jldopen("step_diagnostics.jld2", "a+") do f
    f["FluxesBCs/Jˢ"] = Jˢ
    f["FluxBCs/Jᵀ"] = Jᵀ
    f["FluxBCs/top_salinity_scale"] = top_salinity_scale
end

# using JLD2, NCDatasets
# ds = NCDataset(simulation.output_writers[:computed_output].filepath)
# R_ρ = "R_rho.jld2"
# if isfile(R_ρ)
#     rm(R_ρ)
# end
# jldopen(R_ρ, "w") do f
#     f["R_ρ"] = ds[:R_ρ][:]
# end
# close(ds)

# local plot of figure
# using JLD2, CairoMakie
# output_path = joinpath(@__DIR__, "R_rho_1.05", "nonlineareos", "equal_top_and_bottom_fluxes")
# data = joinpath(output_path, "R_rho.jld2")
# f = jldopen(data)
# R_ρ = f["R_ρ"]
# close(f)

# fig, ax = lines(R_ρ)
# ax.title = "R_ρ with fluxbcs nonlinear eos"
# ax.xlabel = "time (mins)"
# ax.ylabel = "R_ρ"
# fig
# save(joinpath(output_path, "R_rho.png"), fig)
