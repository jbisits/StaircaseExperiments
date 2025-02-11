using StaircaseShenanigans, GLMakie, JLD2

grid = RectilinearGrid(size = 100, z = (-1, 0), topology=(Flat, Flat, Periodic))

EOS = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(equation_of_state = EOS)

closure = ScalarDiffusivity(ν = 1e-5, κ = (S = 1e-7, T = 1e-5))

model = NonhydrostaticModel(grid = grid,
                            tracers = (:T, :S),
                            buoyancy = buoyancy,
                            closure = closure)

T₀ = vcat(fill(0.5, 50), fill(-1.5, 50))
S₀ = vcat(fill(34.7, 50), fill(34.56, 50))
w₀ = vcat(fill(-1, 50), fill(1, 50))

set!(model, S = S₀, T = T₀, w = w₀)

zplot = znodes(model.grid, Center())
fig = Figure(size = (1000, 600))
axS = Axis(fig[1, 1], title = "Initial salinity", xlabel = "S", ylabel = "z")
lines!(axS, interior(model.tracers.S, 1, 1, :), zplot)
axT = Axis(fig[1, 2], title = "Initial temperature", xlabel = "T", ylabel = "z")
lines!(axT, interior(model.tracers.T, 1, 1, :), zplot)
fig
# with halos
lines(parent(model.velocities.w.data)[1, 1, :], -1.03:0.01:0.02)
lines(parent(model.tracers.S.data)[1, 1, :], -1.03:0.01:0.02)
lines(parent(model.tracers.T.data)[1, 1, :], -1.03:0.01:0.02)

simulation = Simulation(model,  Δt = 0.5, stop_time = 0.5 * 10000)

simulation.output_writers[:outputs] = JLD2OutputWriter(model, (T = model.tracers.T, S = model.tracers.S,
                                                                w = model.velocities.w),
                                                        filename = "water_column.jld2",
                                                        schedule = TimeInterval(2))

run!(simulation)
##
water_column_output = joinpath(pwd(), "CPUExperiments/water_column.jld2")
S_ts = FieldTimeSeries(water_column_output, "S")
T_ts = FieldTimeSeries(water_column_output, "T")
w_ts = FieldTimeSeries(water_column_output, "w")

fig = Figure(size = (600, 500))
ax = Axis(fig[1, 1])
sl = Slider(fig[2, 1], range = 1:2501, startvalue = 1)
n = sl.value
snapshot = @lift parent(w_ts)[1, 1, :, $n]
lines!(ax, snapshot, -1.03:0.01:0.02)
fig

## With background state
S_b(z, t) = -0.1 * z
T_b(z, t) = -0.1 * z
model = NonhydrostaticModel(grid = grid,
                            tracers = (:T, :S),
                            background_fields = (T = T_b, S = S_b),
                            buoyancy = buoyancy,
                            closure = closure)

T₀ = 1e-3 * randn(size(model.tracers.T))
S₀ = 1e-3 * randn(size(model.tracers.S))
w₀ = randn(size(model.velocities.w))

set!(model, S = S₀, T = T₀, w = w₀)
S_b_field = Field(model.background_fields.tracers.S)
compute!(S_b_field)
S_total = compute!(Field(model.tracers.S + S_b_field))
lines(parent(S_total)[1, 1, :], -1.03:0.01:0.02)

simulation = Simulation(model,  Δt = 0.5, stop_time = 0.5 * 10000)

simulation.output_writers[:outputs] = JLD2OutputWriter(model, (T = model.tracers.T + model.background_fields.tracers.T,
                                                                S = model.tracers.S + model.background_fields.tracers.S,
                                                                w = model.velocities.w),
                                                        filename = "water_column_background.jld2",
                                                        schedule = TimeInterval(2))
run!(simulation)
