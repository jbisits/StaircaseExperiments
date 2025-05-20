using StaircaseShenanigans, GibbsSeaWater, CairoMakie

architecture = CPU()
Pr = 7   # Prandtl
τ = 0.05 # diff ratio
ν = 2.5e-6 # set this get the others
diffusivities = diffusivities_from_ν(ν; τ, Pr)
domain_extent = (Lx=0.05, Ly=0.05, Lz=-0.5)
domain_topology = (x = Periodic, y = Periodic, z = Bounded)
resolution = (Nx=5, Ny=5, Nz=1000)
ρ₀ = gsw_rho(34.7, 0.5, 0)
eos = TEOS10EquationOfState(reference_density = ρ₀)
eos = CustomLinearEquationOfState(-0.5, 34.6, reference_density = ρ₀)
model_setup = (;architecture, diffusivities, domain_extent, domain_topology, resolution, eos)
model = DNSModel(model_setup...; TD = VerticallyImplicitTimeDiscretization())

number_of_interfaces = 3
depth_of_interfaces = [-0.01, -0.25, -0.49]
#### nonlinear eos
# salinity = [34.56, 34.60, 34.64, 34.7]
# temperature = [-1.45, -0.75, -0.19, 0.51]
#### linear eos
salinity = [34.56, 34.594, 34.64, 34.7]
temperature = [-1.5, -1.0, -0.33, 0.52]
####
staircase_ics = StaircaseICs(model, number_of_interfaces, depth_of_interfaces, salinity, temperature)

initial_noise = (velocities = VelocityNoise(1e-2), tracers = TracerNoise(1e-4, 1e-2))

## setup model
sdns = StaircaseDNS(model, staircase_ics; initial_noise)

visualise_initial_conditions(sdns, 2, 2)

visualise_initial_density(sdns, 2, 2, 0)
