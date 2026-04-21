include("plotting_env_and_load.jl")
## Figure one
# Asymmetric effects of non-linear eos

S₀ᵘ, S₀ˡ = fill(Sᵤ, Int(Nz / 2)), fill(Sₗ, Int(Nz / 2))
Θ₀ᵘ, Θ₀ˡ = fill(Θᵤ, Int(Nz / 2)), fill(Θₗ, Int(Nz / 2))
linear_σ₀ᵘ = total_density.(Θ₀ᵘ, S₀ᵘ, 0, leos_vec[1:Int(Nz / 2)])
linear_σ₀ˡ = total_density.(Θ₀ˡ, S₀ˡ, 0, leos_vec[Int(Nz / 2)+1:end])
nlinear_σ₀ᵘ = gsw_rho.(S₀ᵘ, Θ₀ᵘ, 0)
nlinear_σ₀ˡ = gsw_rho.(S₀ˡ, Θ₀ˡ, 0)
S = erf_tracer_solution.(z, Sₗ, ΔS, κₛ, t, interface_depth)
T = erf_tracer_solution.(z, Θₗ, ΔΘ, κₜ, t, interface_depth)
τ = κₛ / κₜ
linear_σ₀ = total_density.(T, S, 0, leos_vec)
linear_σ₀′ = linear_σ₀ .- mean(vcat(linear_σ₀ᵘ, linear_σ₀ˡ))
max_linear_σ₀′, max_linear_idx = findmax(linear_σ₀′)
min_linear_σ₀′, min_linear_idx = findmin(linear_σ₀′)
linear_σ₀ᵘ′ = linear_σ₀ᵘ .- mean(vcat(linear_σ₀ᵘ, linear_σ₀ˡ))
linear_σ₀ˡ′ = linear_σ₀ˡ .- mean(vcat(linear_σ₀ᵘ, linear_σ₀ˡ))
nlinear_σ₀ = gsw_rho.(S, T, 0)
nlinear_σ₀′ = nlinear_σ₀ .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))
max_nlinear_σ₀′, max_nlinear_idx = findmax(nlinear_σ₀′)
min_nlinear_σ₀′, min_nlinear_idx = findmin(nlinear_σ₀′)
nlinear_σ₀ᵘ′ = nlinear_σ₀ᵘ .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))
nlinear_σ₀ˡ′ = nlinear_σ₀ˡ .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))

# Salinity and temperature panel
fig = Figure(size = (800, 600))
ax = [Axis(fig[1, i]) for i ∈ 1:2]
lines!(ax[1], S₀ˡ, z[1:Int(Nz/2)]; color = (:blue, 0.5), label = "Initial salinity")
lines!(ax[1], S₀ᵘ, z[Int(Nz/2)+1:end]; color = (:blue, 0.5))
lines!(ax[1], S, z; color = :blue, linestyle = :dash, label = L"Salinity, $t = t^{*}$")
xlims!(ax[1], 34.56, 34.72)
ax[1].xlabel = "Salinity (gkg⁻¹)"
ax[1].xlabelcolor = :blue
ax[1].xticklabelcolor = :blue
ax[1].ylabel = "z (m)"
xticks_positions = range(34.567, 34.713, length = 5)
xticks_vals = string.(round.(range(34.567, 34.713, length = 5), digits = 3))
ax[1].xticks = (xticks_positions, xticks_vals)
ax[1].xgridvisible = false
axT = Axis(fig[1, 1];
            xaxisposition = :top,
            xticklabelcolor = :red,
            xlabel = "Θ (°C)",
            xlabelcolor = :red,
            ygridvisible = false,
            title = "(a) Temperature and salinity profiles")
lines!(axT, Θ₀ˡ, z[1:Int(Nz/2)]; color = (:red, 0.5), label = "Initial temperature")
lines!(axT, Θ₀ᵘ, z[Int(Nz/2)+1:end]; color = (:red, 0.5))
lines!(axT, T, z; color = :red, linestyle = :dash, label = L"Temperature, $t = t^{*}$")
axislegend(axT)
axislegend(ax[1], position = :lb)

# density panel
lines!(ax[2], linear_σ₀ˡ′, z[1:Int(Nz/2)]; color = (:black, 0.5), label = "Initial density")
lines!(ax[2], linear_σ₀ᵘ′, z[Int(Nz/2)+1:end]; color = (:black, 0.5))
lines!(ax[2], linear_σ₀′, z; linestyle = :dash, label = L"$\rho_{\mathrm{linear}}$, $t = t^{*}$")
lines!(ax[2], nlinear_σ₀′, z; linestyle = :dash, label = L"$\rho_{\mathrm{nonlinear}}$, $t = t^{*}$")
ax[2].title = "(b) Density profiles"
ax[2].xlabel = L"$σ_{0}'$ (kgm⁻³)"
hideydecorations!(ax[2], grid = false)
axislegend(ax[2])
linkyaxes!(ax[1], ax[2])
fig
##
save("fig1_S_T_sigma_profiles.png", fig)
