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
ρ₀ = gsw_rho(34.7, 0.5, 0)
eos = TEOS10EquationOfState(reference_density = ρ₀)
model_setup = (;architecture, diffusivities, domain_extent, domain_topology, resolution, eos)
# no-slip velocities at top and bottom
no_slip_bc = ValueBoundaryCondition(0.0)
no_slip_bcs = FieldBoundaryConditions(top = no_slip_bc, bottom = no_slip_bc)
# Dirichlet temperature
Tᵤ, Tₗ = -1.5, 0.5
T_top_bc = ValueBoundaryCondition(Tᵤ)
T_bottom_bc = ValueBoundaryCondition(Tₗ)
T_bcs = FieldBoundaryConditions(top = T_top_bc, bottom = T_bottom_bc)
# Dirichlet salinity
Sᵤ, Sₗ = 34.58, 34.70
S_top_bc = ValueBoundaryCondition(Sᵤ)
S_bottom_bc = ValueBoundaryCondition(Sₗ)
S_bcs = FieldBoundaryConditions(top = S_top_bc, bottom = S_bottom_bc)
boundary_conditions = (u=no_slip_bcs, v=no_slip_bcs, T=T_bcs, S=S_bcs)
dns_model = DNSModel(model_setup...; boundary_conditions, TD = VerticallyImplicitTimeDiscretization())

## Initial conditions
depth_of_interface = -0.5
salinity = [Sᵤ, Sₗ]
temperature = [Tᵤ, Tₗ]
dns_model = DNSModel(model_setup...; boundary_conditions, TD = VerticallyImplicitTimeDiscretization())

## Initial conditions
depth_of_interface = -0.5
salinity = [34.58, 34.70]
temperature = [-1.5, 0.5]
interface_ics = SingleInterfaceICs(eos, depth_of_interface, salinity, temperature)

initial_noise = (velocities = VelocityNoise(1e-2), tracers = TracerNoise(1e-4, 1e-2))
## setup model
sdns = StaircaseDNS(dns_model, interface_ics; initial_noise)

## Build simulation
stop_time = Int(16 * 60 * 60) # seconds
initial_state = interface_ics.interface_smoothing isa TanhInterfaceThickness ?  "tanh" : "step"
output_path = joinpath(@__DIR__, "fluxbcs_$(round(interface_ics.R_ρ, digits = 2))", initial_state)
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
             simulation.output_writers[:tracers].filepath, (-0.4, -0.2), (-0.8, -0.6), eos)

reduced_path = findlast('/', simulation.output_writers[:computed_output].filepath)
animation_path = simulation.output_writers[:computed_output].filepath[1:(reduced_path-1)]
cd(animation_path)
@info "Producing animations"
animate_density(simulation.output_writers[:computed_output].filepath, "σ",
                xslice = 17, yslice = 17, density_limit_adjustment = 0.04)
animate_tracers(simulation.output_writers[:tracers].filepath, xslice = 17, yslice = 17,
                S_limit_adjustment = 0.025,
                Θ_limit_adjustment = 0.5)

using JLD2, NCDatasets
ds = NCDataset(simulation.output_writers[:computed_output].filepath)
R_ρ = "R_rho.jld2"
if isfile(R_ρ)
    rm(R_ρ)
end
jldopen(R_ρ, "w") do f
    f["R_ρ"] = ds[:R_ρ][:]
end
close(ds)

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
