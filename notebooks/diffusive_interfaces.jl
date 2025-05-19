### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 7d191973-92d1-4ed7-bcee-650f2185e36f
begin
	using Pkg
	Pkg.activate("..")
	using CairoMakie, GibbsSeaWater, PlutoUI, Statistics, JLD2
	using SpecialFunctions: erf
	using StaircaseShenanigans: total_density, CustomLinearEquationOfState
end

# ╔═╡ b51bd896-64f6-11ef-3b2b-3d362d0c3f8a
md"""
# Diffusive interfaces with a non-linear equation of state 

An initial value problem for a tracer (1D profile) can be posed as
```math
	\begin{aligned}
	\frac{\partial \Theta}{\partial t} &= \kappa_{\Theta}\frac{\partial^2 \Theta}{\partial z^2} \\
	\Theta(z, 0) &= \begin{cases} 0 & z < 0 \\ 1 & z > 0. \end{cases}
	\end{aligned}
```

This can be solved analytically to give the evolution of the tracer in time
```math
	\Theta(z, t) = \frac{1}{2}\left(1 + \mathrm{erf}\left(\frac{z}{\sqrt{4\kappa_{\Theta}t}}\right)\right).
```

The location of the step change and the magnitude of difference can be added so that the initial condition becomes
```math
	\Theta(z, 0) = \begin{cases} \Theta^{*} & z < -a \\ \Theta^{*} + \Delta\Theta & z > -a \end{cases}
```
which is the Heaviside function ``\Theta^{*} + \Delta\Theta H\left(z + a\right).``
This is the kind of step change profile we have looked already and it can be modelled analytically according to
```math
	\Theta(z, t) = \Theta^{*} + \frac{\Delta\Theta}{2}\left(1 + \mathrm{erf}\left(\frac{z + a}{\sqrt{4\kappa_{\Theta}t}}\right)\right)
```
with a similar expression for salinity.
So we set the lower layer temperature and salinity as ``\left(S^{*}, \Theta^{*}\right)`` then pick an interface and ``\left(\Delta S, \Delta\Theta\right)`` that the upper layer is and we can evolve salinity and temperature and see how density evolves due to mixing and different initial conditions.

I want to use this model to explore what can be said about a diffusive interface with a non-linear equation of state.
In particular when relating to a linear equation of state how do the fluxes, determined from the slope of the interface, change.
"""

# ╔═╡ 22c3e8cc-2d34-4fa4-9c2d-c807c6f53998
	erf_tracer_solution(z, Cₗ::Number, ΔC::Number, κ::Number, time, interface_location) =
    Cₗ + 0.5 * ΔC * (1 + erf((z - interface_location) / sqrt(4 * κ * time)))

# ╔═╡ a92c3315-ae34-447a-bc3e-7960bf2694d1
begin
	time = @bind time PlutoUI.Slider(0.00001:10000:500000.00001, default = 10000.00001)
	κₛ_var = @bind κₛ PlutoUI.Slider([1e-7, 1e-8, 1e-9], show_value = true, default = 1e-9)
	κₜ_var = @bind κₜ PlutoUI.Slider([1e-7, 1e-8, 1e-9], show_value = true, default = 1e-7)

	S_star, Θ_star = 34.7, 0.5
	ΔΘ = -2
	Θ_upper = Θ_star + ΔΘ
	S_stable, S_cab = 34.551, 34.58
	ΔSₛ = S_stable - S_star
	ΔS_c= S_cab - S_star
	select_upper_layer = @bind upper_layer Select(["Stable", "Cabbeling"])

	τ = κₛ / κₜ
	md"""
	# Non-linear equation of state 
	
	### Upper layer stability
	$(select_upper_layer)
	
	### Diffusivities
	Salt = $(κₛ_var) m²s⁻¹

	Temperature = $(κₜ_var) m²s⁻¹

	### Time
	t = $time
	"""
end

# ╔═╡ 8b0cb0d8-ad5b-46f8-8990-9438589976d8
let

	S_upper = upper_layer == "Stable" ? S_stable : S_cab
	ΔS = upper_layer == "Stable" ? ΔSₛ : ΔS_c
	
	z = range(-1, 0, length = 1400)
    interface_location = -0.5
	τ = κₛ / κₜ
	S = erf_tracer_solution.(z, S_star, ΔS, κₛ, time, interface_location)
	T = erf_tracer_solution.(z, Θ_star, ΔΘ, κₜ, time, interface_location)
	σ₀ = gsw_rho.(S, T, 0)

	fontsize = 22
	labelsize = 16
	fig = Figure(size = (900, 1000); fontsize)
	ax = [Axis(fig[1, i]) for i ∈ 1:2]
	lines!(ax[1], S, z; color = (:blue, 0.5), label = "Salinity")
	xlims!(ax[1], 34.52, 34.72)
	ax[1].xlabel = "Salinity (gkg⁻¹)"
	ax[1].xlabelcolor = :blue
	ax[1].xticklabelcolor = :blue
	ax[1].ylabel = "z (m)"
	axT = Axis(fig[1, 1];
	           xaxisposition = :top,
	           xticklabelcolor = :red,
	           xlabel = "Θ (°C)",
	           xlabelcolor = :red,
	           title = "Temperature and salinity profiles\nat t = $(round(time/60; digits = 4)) mins")
	lines!(axT, T, z; color = (:red, 0.5), label = "Temeperature")
	lines!(ax[2], σ₀, z; color = (:black, 0.5))
	ax[2].title = "Density at t = $(round(time/60; digits = 4)) mins."
	ax[2].xlabel = "σ₀ (kgm⁻³)"
	ax[2].xticklabelrotation = π/4
	axislegend(axT; labelsize)
	axislegend(ax[1], position = :lb; labelsize)
	hideydecorations!(ax[2], grid = false)

	linkyaxes!(ax[1], ax[2])
	colsize!(fig.layout, 1, Relative(3/5))
	fig

	ax2 = Axis(fig[2, :], title = "S-T evolution. τ = $(τ)", xlabel = "Absolute salinity (gkg⁻¹)", ylabel = "Conservative temperature (°C)")
	N = 2000
	S_range, Θ_range = range(minimum(S), maximum(S), length = N), range(minimum(T), maximum(T), length = N)
	S_grid, Θ_grid = ones(N) .* S_range', ones(N)' .* Θ_range
	ρ = gsw_rho.(S_grid, Θ_grid, 0)
	ρ_star = gsw_rho(S_star, Θ_star, 0)
	ρ_upper = gsw_rho(S_upper, Θ_upper, 0)

	scatter!(ax2, S, T, σ₀, markersize = 6, label = "S-T profile")
	σ₀_max, max_idx = findmax(σ₀)
	σ₀_min, min_idx = findmin(σ₀)

	S_minmax = [S[min_idx], S[max_idx]]
	T_minmax = [T[min_idx], T[max_idx]]

	scatter!(ax2, S_minmax[1], T_minmax[1], label = "Minimum density", color = :orange)
	scatter!(ax2, S_minmax[2], T_minmax[2], label = "Maximum density", color = :magenta)

	ρ_min = gsw_rho(S[min_idx], T[min_idx], 0)
	ρ_max = gsw_rho(S[max_idx], T[max_idx], 0)
	ρ_diff = abs(z[max_idx] - z[min_idx])
	contour!(ax2, S_range, Θ_range, ρ'; levels = [ρ_min, ρ_upper, ρ_star, ρ_max], colormap = :dense, label = "Isopycnals")

	scatter!(ax[2], σ₀_min, z[min_idx], label = "Minimum density", color = :orange)
	scatter!(ax[2], σ₀_max, z[max_idx], label = "Maximum density", color = :magenta)

	lines!(ax2, S_minmax, T_minmax, color = :red, linestyle = :dash)

	axislegend(ax2, position = :lt)
	fig

	ΔS_minmax = diff(S_minmax)[1]
	ΔT_minmax = diff(T_minmax)[1]
	slope = ΔT_minmax / ΔS_minmax

	quarter_domain = round(Int, length(z) / 4)
	S̄ᵤ, S̄ₗ = mean(S[1:quarter_domain]), mean(S[3*quarter_domain:end])
	ΔS_interface = abs(S̄ᵤ - S̄ₗ)
	find_S_interface = findall(median(S)-ΔS_interface/8 .< S .< median(S) + ΔS_interface/8)
	S_mat = [ones(length(find_S_interface)) z[find_S_interface]]
	S_linfit = S_mat \ S[find_S_interface]
	hₛ = abs(ΔS_interface / S_linfit[2])

	T̄ᵤ, T̄ₗ = mean(T[1:quarter_domain]), mean(T[3*quarter_domain:end])
	ΔΘ_interface = abs(T̄ᵤ - T̄ₗ)
	find_Θ_interface = findall(median(T)-ΔΘ_interface/8 .< T .< median(T)+ΔΘ_interface/8)
	Θ_mat = [ones(length(find_Θ_interface)) z[find_Θ_interface]]
	Θ_linfit = Θ_mat \ T[find_Θ_interface]
	hₜ = abs(ΔΘ_interface / Θ_linfit[2])

	scatter!(ax[1], median(S)-ΔS_interface/8, z[find_S_interface][end], color = :blue)
	scatter!(ax[1], median(S)+ΔS_interface/8, z[find_S_interface][1], color = :blue)
	scatter!(axT, median(T)-ΔΘ_interface/8, z[find_Θ_interface][end], color = :red)
	scatter!(axT, median(T)+ΔΘ_interface/8, z[find_Θ_interface][1], color = :red)
	
	Δρ = gsw_rho(S_upper, Θ_upper, 0) - gsw_rho(S_star, Θ_star, 0)

	md"""
	$(fig)

	``\Delta\rho_{\mathrm{lower}} =`` $(σ₀_max - ρ_star) (maximum density - density of lower layer)

	``\Delta\rho_{\mathrm{upper}} =`` $(σ₀_min - ρ_upper) (minimum density - density of upper layer)

	Asymmetry in the density differences!
	"""
end

# ╔═╡ 5a8f995b-f45d-47f0-bcdf-8df553c21d17
md"""
# Linear equation of state
"""

# ╔═╡ 1e626657-a23c-42b5-bbae-056c6a61948b
let

	S_upper = upper_layer == "Stable" ? S_stable : S_cab
	ΔS = upper_layer == "Stable" ? ΔSₛ : ΔS_c
	
	z = range(-1, 0, length = 1400)
	eos = CustomLinearEquationOfState(-0.5, 34.6)
	eos_vec = fill(eos, length(z))
    interface_location = -0.5
	S = erf_tracer_solution.(z, S_star, ΔS, κₛ, time, interface_location)
	T = erf_tracer_solution.(z, Θ_star, ΔΘ, κₜ, time, interface_location)
	τ = κₛ / κₜ
	σ₀ = total_density.(T, S, 0, eos_vec)

	fontsize = 22
	labelsize = 16
	fig = Figure(size = (900, 1000); fontsize)
	ax = [Axis(fig[1, i]) for i ∈ 1:2]
	lines!(ax[1], S, z; color = (:blue, 0.5), label = "Salinity")
	xlims!(ax[1], 34.52, 34.72)
	ax[1].xlabel = "Salinity (gkg⁻¹)"
	ax[1].xlabelcolor = :blue
	ax[1].xticklabelcolor = :blue
	ax[1].ylabel = "z (m)"
	axT = Axis(fig[1, 1];
	           xaxisposition = :top,
	           xticklabelcolor = :red,
	           xlabel = "Θ (°C)",
	           xlabelcolor = :red,
	           title = "Temperature and salinity profiles\nat t = $(round(time/60; digits = 4)) mins")
	lines!(axT, T, z; color = (:red, 0.5), label = "Temeperature")
	lines!(ax[2], σ₀, z; color = (:black, 0.5))
	ax[2].title = "Density at t = $(round(time/60; digits = 4)) mins."
	ax[2].xlabel = "σ₀ (kgm⁻³)"
	ax[2].xticklabelrotation = π/4
	axislegend(axT; labelsize)
	axislegend(ax[1], position = :lb; labelsize)
	hideydecorations!(ax[2], grid = false)

	linkyaxes!(ax[1], ax[2])
	colsize!(fig.layout, 1, Relative(3/5))
	fig

	ax2 = Axis(fig[2, :], title = "S-T evolution. τ = $(τ)", xlabel = "Absolute salinity (gkg⁻¹)", ylabel = "Conservative temperature (°C)")
	N = 2000
	S_range, Θ_range = range(minimum(S), maximum(S), length = N), range(minimum(T), maximum(T), length = N)
	S_grid, Θ_grid = ones(N) .* S_range', ones(N)' .* Θ_range
	ρ = total_density.(Θ_grid, S_grid, 0, fill(eos, size(S_grid)))
	ρ_star = total_density(Θ_star, S_star, 0, eos)
	ρ_upper = total_density(Θ_upper, S_upper, 0, eos)

	scatter!(ax2, S, T, σ₀, markersize = 6, label = "S-T profile")
	σ₀_max, max_idx = findmax(σ₀)
	σ₀_min, min_idx = findmin(σ₀)

	S_minmax = [S[min_idx], S[max_idx]]
	T_minmax = [T[min_idx], T[max_idx]]

	scatter!(ax2, S_minmax[1], T_minmax[1], label = "Minimum density", color = :orange)
	scatter!(ax2, S_minmax[2], T_minmax[2], label = "Maximum density", color = :magenta)

	ρ_min = total_density(T[min_idx], S[min_idx], 0, eos)
	ρ_max = total_density(T[max_idx], S[max_idx], 0, eos)
	ρ_diff = abs(z[max_idx] - z[min_idx])
	contour!(ax2, S_range, Θ_range, ρ'; levels = [ρ_min, ρ_upper, ρ_star, ρ_max], colormap = :dense, label = "Isopycnals")

	scatter!(ax[2], σ₀_min, z[min_idx], label = "Minimum density", color = :orange)
	scatter!(ax[2], σ₀_max, z[max_idx], label = "Maximum density", color = :magenta)

	lines!(ax2, S_minmax, T_minmax, color = :red, linestyle = :dash)

	axislegend(ax2, position = :lt)
	fig

	ΔS_minmax = diff(S_minmax)[1]
	ΔT_minmax = diff(T_minmax)[1]
	slope = ΔT_minmax / ΔS_minmax

	quarter_domain = round(Int, length(z) / 4)
	S̄ᵤ, S̄ₗ = mean(S[1:quarter_domain]), mean(S[3*quarter_domain:end])
	ΔS_interface = abs(S̄ᵤ - S̄ₗ)
	find_S_interface = findall(median(S)-ΔS_interface/8 .< S .< median(S) + ΔS_interface/8)
	S_mat = [ones(length(find_S_interface)) z[find_S_interface]]
	S_linfit = S_mat \ S[find_S_interface]
	hₛ = abs(ΔS_interface / S_linfit[2])

	T̄ᵤ, T̄ₗ = mean(T[1:quarter_domain]), mean(T[3*quarter_domain:end])
	ΔΘ_interface = abs(T̄ᵤ - T̄ₗ)
	find_Θ_interface = findall(median(T)-ΔΘ_interface/8 .< T .< median(T)+ΔΘ_interface/8)
	Θ_mat = [ones(length(find_Θ_interface)) z[find_Θ_interface]]
	Θ_linfit = Θ_mat \ T[find_Θ_interface]
	hₜ = abs(ΔΘ_interface / Θ_linfit[2])

	scatter!(ax[1], median(S)-ΔS_interface/8, z[find_S_interface][end], color = :blue)
	scatter!(ax[1], median(S)+ΔS_interface/8, z[find_S_interface][1], color = :blue)
	scatter!(axT, median(T)-ΔΘ_interface/8, z[find_Θ_interface][end], color = :red)
	scatter!(axT, median(T)+ΔΘ_interface/8, z[find_Θ_interface][1], color = :red)

	md"""
	$(fig)

	``\Delta\rho_{\mathrm{lower}} =`` $(σ₀_max - ρ_star) (maximum density - density of lower layer)

	``\Delta\rho_{\mathrm{upper}} =`` $(σ₀_min - ρ_upper) (minimum density - density of upper layer)

	In linear case as expected they are even!
	"""
end

# ╔═╡ d8e4ad73-b296-462b-a161-9666ba6322b6
md"""
# Asymmetry of density difference dependence on ``\tau``
"""

# ╔═╡ 48eb84eb-42bf-4938-9d9e-0c8794fa5a5c
begin
	S_upper = upper_layer == "Stable" ? S_stable : S_cab
	ΔS = upper_layer == "Stable" ? ΔSₛ : ΔS_c
	
	z = range(-1, 0, length = 1400)
	eos = CustomLinearEquationOfState(-0.5, 34.6)
	eos_vec = fill(eos, length(z))
    interface_location = -0.5
	κₛ_range = range(1e-9, κₜ, length = 100)
	τ_range = κₛ_range ./ κₜ
	σ₀_nonlinear_max = similar(κₛ_range)
	σ₀_nonlinear_min = similar(κₛ_range)
	σ₀_linear_max = similar(κₛ_range)
	σ₀_linear_min = similar(κₛ_range)
	t = 10000 # seconds, this value is not stricly important as maximum density does not increase
	for i ∈ eachindex(κₛ_range)
		S = erf_tracer_solution.(z, S_star, ΔS, κₛ_range[i], t, interface_location)
		T = erf_tracer_solution.(z, Θ_star, ΔΘ, κₜ, t, interface_location)
		σ₀_nonlinear = gsw_rho.(S, T, 0)
		σ₀_nonlinear_max[i] = maximum(σ₀_nonlinear)
		σ₀_nonlinear_min[i] = minimum(σ₀_nonlinear)
		
		σ₀_linear = total_density.(T, S, 0, eos_vec)
		σ₀_linear_max[i] = maximum(σ₀_linear)
		σ₀_linear_min[i] = minimum(σ₀_linear)
	end
	
	Δσ_lower_nonlinear = σ₀_nonlinear_max .- gsw_rho(S_star, Θ_star, 0)
	Δσ_upper_nonlinear = abs.(σ₀_nonlinear_min .- gsw_rho(S_upper, Θ_upper, 0))
	Δσ_nonlinear = Δσ_upper_nonlinear ./ Δσ_lower_nonlinear

	Δσ_lower_linear = σ₀_linear_max .- total_density(Θ_star, S_star, 0, eos)
	Δσ_upper_linear = abs.(σ₀_linear_min .- total_density(Θ_upper, S_upper, 0, eos))
	Δσ_linear = Δσ_upper_linear ./ Δσ_lower_linear
	
	fig, ax = lines(τ_range, Δσ_nonlinear, label = "nonlinear eos")
	lines!(ax, τ_range, Δσ_linear, label = "linear eos")
	ax.xlabel = "τ"
	ax.ylabel = "Δσ_upper/Δσ_lower"
	axislegend(ax, position = :rc)
	fig
end

# ╔═╡ 9695b41b-4dd0-429b-a3c5-c37c2a6dbaf8
md"""
This plot shows:

- asymmetry in density difference increases (lower density difference is the denomiator) as the diffusivity ratio decreases
- asymmetry becomes more apparent as the once less than around 0.1 or so
- almost linear dependence
- could I confirm this relationship using a model run with various $\tau$ values? I think could get away with five here but would you just show diffusion or include turbulence?
"""

# ╔═╡ 7ad7693e-6f18-474a-89e8-b2d433aea261
TableOfContents()

# ╔═╡ Cell order:
# ╟─7d191973-92d1-4ed7-bcee-650f2185e36f
# ╟─b51bd896-64f6-11ef-3b2b-3d362d0c3f8a
# ╟─22c3e8cc-2d34-4fa4-9c2d-c807c6f53998
# ╟─a92c3315-ae34-447a-bc3e-7960bf2694d1
# ╟─8b0cb0d8-ad5b-46f8-8990-9438589976d8
# ╟─5a8f995b-f45d-47f0-bcdf-8df553c21d17
# ╟─1e626657-a23c-42b5-bbae-056c6a61948b
# ╟─d8e4ad73-b296-462b-a161-9666ba6322b6
# ╟─48eb84eb-42bf-4938-9d9e-0c8794fa5a5c
# ╟─9695b41b-4dd0-429b-a3c5-c37c2a6dbaf8
# ╟─7ad7693e-6f18-474a-89e8-b2d433aea261
