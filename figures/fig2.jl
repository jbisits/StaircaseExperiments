include("plotting_env_and_load.jl")
## Figure two
# Profiles in salinity-temperature space
σ_grad = get(ColorSchemes.diff, [0.2, 0.3, 0.6, 1])
σ_grad = cgrad(σ_grad, 4, categorical = true)
fig = Figure(size = (800, 600))
eos_title = ["(a) Linear equation of state", "(b) Nonlinear equation of state"]
ax = [Axis(fig[i, 1], title = eos_title[i], xlabel = "S (gkg⁻¹)", ylabel = "Θ (°C)") for i ∈ 1:2]
N = 2000
S_range, Θ_range = range(minimum(S), maximum(S), length = N), range(minimum(T), maximum(T), length = N)
S_grid, Θ_grid = ones(N) .* S_range', ones(N)' .* Θ_range
# Linear
linear_σ₀_grid = total_density.(Θ_grid, S_grid, 0, fill(leos, size(S_grid)))  .- mean(vcat(linear_σ₀ᵘ, linear_σ₀ˡ))
linear_σ₀_min = minimum(linear_σ₀) .- mean(vcat(linear_σ₀ᵘ, linear_σ₀ˡ))
linear_σ₀_lower = total_density(Θₗ, Sₗ, 0, leos) .- mean(vcat(linear_σ₀ᵘ, linear_σ₀ˡ))
linear_σ₀_upper = total_density(Θᵤ, Sᵤ, 0, leos) .- mean(vcat(linear_σ₀ᵘ, linear_σ₀ˡ))
linear_σ₀_max = maximum(linear_σ₀) .- mean(vcat(linear_σ₀ᵘ, linear_σ₀ˡ))

function my_contour_label_formatter(level::Real)::String
    lev_short = round(level; digits = 3)
    x = @sprintf("%.1e", lev_short)
    string(x)
end

σ_linear_levels = [linear_σ₀_min, linear_σ₀_upper, linear_σ₀_lower, linear_σ₀_max]
cn = contour!(ax[1], S_range, Θ_range, linear_σ₀_grid';
              levels = σ_linear_levels,
              colormap = σ_grad,
              colorrange = extrema(σ_linear_levels),
              linewidth = 0.75,
              label = "Isopycnals")
Colorbar(fig[1, 2], colorrange = extrema(σ_linear_levels),
        colormap = cgrad(σ_grad, 4, categorical = true),
        ticks = (range(extrema(σ_linear_levels)..., length = 4), string.(round.(σ_linear_levels, digits = 3))),
        label = L"$\Delta \sigma_{0}'$ (kgm$^{-3}$)")
lines!(ax[1], S, T, linear_σ₀; label = "S-T profile", color = :tomato, linestyle = :dash)
S_minmax = [S[min_linear_idx], S[max_linear_idx]]
T_minmax = [T[min_linear_idx], T[max_linear_idx]]

scatter!(ax[1], S[1], T[1]; markersize, label = "Initial lower density", color = σ_grad[3])
scatter!(ax[1], S[end], T[end]; markersize, label = "Initial upper density", color = σ_grad[2])
scatter!(ax[1], S_minmax[1], T_minmax[1]; markersize, label = "Minimum density", color = σ_grad[1])
scatter!(ax[1], S_minmax[2], T_minmax[2]; markersize, label = "Maximum density", color = σ_grad[end])

hidexdecorations!(ax[1], grid = false, ticks = false)
nlinear_σ₀_grid = gsw_rho.(S_grid, Θ_grid, 0) .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))
nlinear_σ₀_lower = gsw_rho(Sₗ, Θₗ, 0) .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))
nlinear_σ₀_upper = gsw_rho(Sᵤ, Θᵤ, 0) .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))
nlinear_σ₀_min = minimum(nlinear_σ₀) .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))
nlinear_σ₀_max = maximum(nlinear_σ₀) .- mean(vcat(nlinear_σ₀ᵘ, nlinear_σ₀ˡ))

lines!(ax[2], S, T, nlinear_σ₀; label = "S-T profile", color = :tomato, linestyle = :dash)
S_minmax = [S[min_nlinear_idx], S[max_nlinear_idx]]
T_minmax = [T[min_nlinear_idx], T[max_nlinear_idx]]

σ_nl_levels = [nlinear_σ₀_min, nlinear_σ₀_upper, nlinear_σ₀_lower, nlinear_σ₀_max]
cn = contour!(ax[2], S_range, Θ_range, nlinear_σ₀_grid';
              levels = σ_nl_levels,
              colormap = σ_grad,
              colorrange = (-maximum(σ_nl_levels), maximum(σ_nl_levels)),
              linewidth = 0.75,
              label = "Isopycnals")
Colorbar(fig[2, 2], colorrange = extrema(σ_nl_levels),
        colormap = σ_grad,
        ticks = (range(nlinear_σ₀_min, nlinear_σ₀_max, length = 4), string.(round.(σ_nl_levels, digits = 3))),
        label = L"$\Delta \sigma_{0}'$ (kgm$^{-3}$)")
scatter!(ax[2], S_minmax[1], T_minmax[1]; markersize, label = "Minimum density", color = σ_grad[1])
scatter!(ax[2], S[end], T[end]; markersize, label = "Initial upper density", color = σ_grad[2])
scatter!(ax[2], S[1], T[1]; markersize, label = "Initial lower density", color = σ_grad[3])
scatter!(ax[2], S_minmax[2], T_minmax[2]; markersize, label = "Maximum density", color = σ_grad[end])

Legend(fig[3, :], ax[2], orientation = :horizontal, nbanks = 2)
##
save("fig2.png", fig)
