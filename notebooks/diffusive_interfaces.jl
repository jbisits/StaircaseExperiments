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
	using StaircaseShenanigans: ρ, RoquetEquationOfState, CustomLinearEquationOfState
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

	lines!(ax2, S_minmax, T_minmax, color = :red, linestyle = :dash, label = "Interface (?)")

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
	
	r = hₜ / hₛ
	Fₜ = κₜ * (ΔΘ_interface) / hₜ
	Fₛ = κₛ * (ΔS_interface) / hₛ
	F_ratio = Fₛ / Fₜ
	
	Δρ = gsw_rho(S_upper, Θ_upper, 0) - gsw_rho(S_star, Θ_star, 0)
	md"""
	$(fig)

	Initial ``Δρ = `` $(round(Δρ, digits = 5)).

	Slope connecting points of minimum and maximum density on the ``S - \Theta`` figure is $(round(slope[1], digits = 3)).
	The interface thickness can be estimated as ``h_{s} = `` $(round(hₛ, digits = 3)) and ``h_{t} = `` $(round(hₜ, digits = 3)) which then give the interface thickness ratio ``r = h_{t} / h_{s} = `` $(round(r, digits = 3)) where .
	[Carpenter (2012)](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/simulations-of-a-doublediffusive-interface-in-the-diffusive-convection-regime/63D2ECE2AA41439E01A01F9A0D76F2E2) set an initial interface thickness so I am not quite sure how I can compare these numbers to their work plus this has non-linear equation of state.

	**Currently I am still a little unsure about the ``h_{s}`` and ``h_{t}`` calculations.**

	## Molecular flux

	Carpenter (2012) also defin the *molecular flux*, in equation (5.3), as
	```math
	F_{φ} = κ_{φ}\frac{Δφ}{h_{φ}}.
	```
	This assumes a purely molecular flux through the interface.
	Computing these ``F_{Θ} = `` $(Fₜ) and ``F_{S} = `` $(Fₛ).
	The ratio of these fluxes is then ``R_{F} = F_{S} / F_{Θ} = `` $(F_ratio).

	At this stage might be good to run some of these calculations for the non-linear eos run I have and see what I can find.
	"""
end

# ╔═╡ feb1de5f-8dfd-478d-89d9-ff668c55b892
begin
	animation_path = "../../CabbelingExperiments/data/animations/Double diffusion/Cabbeling_DD_600min_densityratio100"
	tracers = joinpath(animation_path, "tracers.mp4")
	LocalResource(tracers)
end

# ╔═╡ 056df154-3413-4459-85c2-a1aeaeed9bdb
begin
	density = joinpath(animation_path, "density.mp4")
	LocalResource(density)
end

# ╔═╡ 28033f6f-2107-4b8d-a2fe-78e4eacd8d14
md"""
## Model output

I have a simulation with ``τ = 0.01`` that I ran as part of my cabbeling work.
I have attempted to compute the interface thickness using the methods above and the horizontally averaged flux and diffusivity of both temperature and salinity.
This was a first pass and the averaging was because Gadi was going to be offline for a few days so I wanted to get something done.

The tracers and density evolution for this simulation are also shown.

### Interface thickness
"""

# ╔═╡ 78a2038d-2da0-43e9-9839-cdce03002927
begin
	computed_output = "single_interface_fluxes_nonlinear.jld2"
	co = load(computed_output)
	replace!(co["ha_κₜ"], Inf => 0)
	replace!(co["ha_κₜ"], -Inf => 0)
	replace!(co["ha_κₛ"], Inf => 0)
	replace!(co["ha_κₛ"], -Inf => 0)
	reverse!(co["ha_κₛ"], dims = 1)
	reverse!(co["ha_κₜ"], dims = 1)
	co["ha_κₜ"] *= 100 # 100 = 1 / area where area = horizontal area 0.1^2
	co["ha_κₛ"] *= 100
	reverse!(co["ha_Fₛ"], dims = 1)
	reverse!(co["ha_Fₜ"], dims = 1)
	replace!(co["hₜ"], Inf => 0)
	replace!(co["hₛ"], Inf => 0)
	nothing
end

# ╔═╡ 120a212f-f5c1-4520-802d-28ad737ebdee
@bind interface_plot_window PlutoUI.Slider(3:length(co["hₜ"]), default = length(co["hₜ"]))

# ╔═╡ 106ec4cd-a244-4457-b6b9-3b70a4b6db3e
let
	fig = Figure(size = (500, 500))
	ax = Axis(fig[1, 1], xlabel = "time (s)", ylabel = "Interface thickness (m, log10)", title = "Interface thickness")
	lines!(ax, log10.(co["hₜ"][2:interface_plot_window]), label = "hₜ")
	# hlines!(ax, log10.(10), label = "thickness = 10cm", linestyle = :dash, color = :black)
	axislegend(ax, position = :lt)
	fig
end

# ╔═╡ f2a9ff46-033c-4d12-94f8-3f12b4894a11
md"""
Clearly there are some issues here; mainly that these numbers look much too large.
As the simulation goes on the criteria of searching over ``\text{median}(φ) - Δφ/8 < φ < \text{median}(φ) + Δφ/8`` might need to be made more narrow.
There are also large fluctuations due to the turbulence that would likely not be present in lower density ratio flows so I also wonder if horizontal averaging first might be better.
This was done quite quickly by me so it could be that it just needs to be properly debugged and then will get more accurate results.

I also ended up only saving the temperature interface thickness so cannot get the thickness ratio unfortunately.

### Horizontally averaged fluxes
"""

# ╔═╡ 054c0968-fbc4-4ec2-b9a1-fd7c005f49cc
begin
	zoom_var = @bind zoom Select(["false", "true"])
	md"""
	Zoom on middle of domain $(zoom_var)
	"""
end

# ╔═╡ 214a62b5-d46f-4636-bc26-3d5a4372d0d9
begin
	z = range(-1, 0, length = 1399)
	zrange = zoom == "false" ? eachindex(z) : 680:800
	nothing
end

# ╔═╡ f2163623-d32f-4203-b1ce-9f87478fc9b3
let
	fig = Figure(size = (800, 500))
	axS = Axis(fig[1, 1], ylabel = "z (m)", title = "Horizontally averaged salinity flux")
	hidexdecorations!(axS)
	hmS = heatmap!(axS, 1:660, z[zrange], co["ha_Fₛ"][zrange, :]', colormap = :haline)
	Colorbar(fig[1, 2], hmS, label = "Fₛ")
	axT = Axis(fig[2, 1], xlabel = "time (min)", ylabel = "z (m)", title = "Horizontally averaged temperature flux")
	hmT = heatmap!(axT, 1:660, z[zrange], co["ha_Fₜ"][zrange, :]', colormap = :thermal)
	Colorbar(fig[2, 2], hmT, label = "Fₜ")
	fig
end

# ╔═╡ 88ad24ee-072c-4927-a21a-8159c80d75b0
let
	timestamps = 1:600
	int_idx = findfirst(z .> -0.4685)
	fig, ax = lines(timestamps[400:end], co["ha_Fₛ"][int_idx, 400:end], label = "Salinity flux")
	lines!(ax, timestamps[400:end], co["ha_Fₜ"][int_idx, 400:end], label = "Temperature flux")
	ax.xlabel = "time (mins)"
	axislegend(ax)
	fig
end

# ╔═╡ 3face522-479b-4f87-a332-b858877094bc
md"""
A figure like this for the horizontally averaged temperature and salinity might also be useful.
"""

# ╔═╡ e5079cbd-7878-4e71-96e1-71f7ceccfb29
md"""
### Horizontally averaged effective diffusivity

I think these might be out by a factor of 100 (i.e. need to multiply by surface area) but will check.

From the figures below it looks like the effective diffusivity could be used as an alternate measure for interface thickness.
Need to get the Carpenter (2012) method working first so can compare but in the hovmoller figure below there is a clear region (around the blue line) where there is lower, approximately molecular diffusivity values which appear to wider (especially in the temperature) as time rolls on.
"""

# ╔═╡ c521d03e-929c-493e-ba8b-db51983a2c2a
let
	fig = Figure(size = (800, 500))
	axS = Axis(fig[1, 1], ylabel = "z (m)", title = "Horizontally averaged effective diffusivity", subtitle = "salinity")
	hidexdecorations!(axS)
	z = range(-1, 0, length = 1399)
	hmS = heatmap!(axS, 1:660, z[zrange], log10.(abs.(co["ha_κₛ"][zrange, :]')), colormap = :haline)
	Colorbar(fig[1, 2], hmS, label = "κₛ")
	axT = Axis(fig[2, 1], xlabel = "time (min)", ylabel = "z (m)", title = "Horizontally averaged effective diffusivity", subtitle = "temperature")
	hmT = heatmap!(axT, 1:660, z[zrange], log10.(abs.(co["ha_κₜ"][zrange, :]')), colormap = :thermal)
	hlines!(axS, -0.4685, linestyle = :dot, color = :yellow)
	hlines!(axT, -0.4685, linestyle = :dot, color = :yellow)
	Colorbar(fig[2, 2], hmT, label = "κₜ")
	fig
end

# ╔═╡ d7aa3853-6c71-46da-983a-2c1eb0c817a7
let
	find = findfirst(z .< -0.4685)
	Θinterface = log10.(abs.(co["ha_κₜ"][find, :]))
	Θinterface_mean = mean(Θinterface[.!isnan.(Θinterface)])
	replace!(Θinterface, -Inf => NaN)
	Θinterface_mean = mean(Θinterface[.!isnan.(Θinterface)])
	Sinterface = log10.(abs.(co["ha_κₛ"][find, :]))
	replace!(Sinterface, -Inf => NaN)
	Sinterface_mean = mean(Sinterface[.!isnan.(Sinterface)])
	fig, ax = lines(Θinterface, label = "κₜ, time mean = $(round(Θinterface_mean, digits = 2))")
	ax.xlabel = "time (min)"
	ax.ylabel = "Effective diffusivity (log10)"
	ax.title = "Effective diffusivity at (or close to) interface"
	lines!(ax, Sinterface, label = "κₛ, time mean = $(round(Sinterface_mean, digits = 2))")
	axislegend(ax)
	fig
end

# ╔═╡ 7f417e03-72ef-47f5-923a-b932de687a37
@bind eff_diff_time_window PlutoUI.Slider(eachindex(co["ha_κₛ"][1, :]))

# ╔═╡ 7b798e2e-4370-43fc-a3d7-26587aeec349
let
	fig, ax = lines(vec(log10.(abs.(mean(co["ha_κₜ"][:, eff_diff_time_window:end], dims = 2)))), z, label = "temperature")
	lines!(ax, vec(log10.(abs.(mean(co["ha_κₛ"][:, eff_diff_time_window:end], dims = 2)))), z, label = "salinity")
	hlines!(ax, -0.4685, linestyle = :dot, color = :black)
	ax.title = "Time mean effective diffusivity at each depth level"
	ax.subtitle = "window = $(eff_diff_time_window) - 600min"
	ax.xlabel = "Effective diffusivity (log10)"
	ax.ylabel = "z (m)"
	axislegend(ax, position = :rc)
	fig
end

# ╔═╡ 7c657d45-cd7a-4106-b23b-77dccf5c982f
md"""
# Linear equation of state
"""

# ╔═╡ 308bc1ca-673d-48a9-8938-58c64c3d66a1
time_linear = @bind time_linear PlutoUI.Slider(0.00001:10000:500000.00001)

# ╔═╡ ac2e0369-626f-40f4-a405-6941ff60f6d2
let

	S_upper = upper_layer == "Stable" ? 34.5341 : S_cab # S_stable : S_cab
	ΔS = upper_layer == "Stable" ? S_upper - S_star : ΔS_c
	
	z = range(-1, 0, length = 800)
    interface_location = -0.5
	S = erf_tracer_solution.(z, S_star, ΔS, κₛ, time_linear, interface_location)
	T = erf_tracer_solution.(z, Θ_star, ΔΘ, κₜ, time_linear, interface_location)
	eos = CustomLinearEquationOfState(0.0, 34.6)
	σ₀ = ρ.(T, S, 0, fill(eos, length(S)))

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
	           title = "Temperature and salinity profiles\nat t = $(round(time_linear/60; digits = 4)) mins")
	lines!(axT, T, z; color = (:red, 0.5), label = "Temeperature")
	lines!(ax[2], σ₀, z; color = (:black, 0.5))
	ax[2].title = "Density at t = $(round(time_linear/60; digits = 4)) mins."
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
	ρ_ = ρ.(Θ_grid, S_grid, 0, fill(eos, (N, N)))
	ρ_star = ρ(Θ_star, S_star, 0, eos)
	ρ_upper = ρ(Θ_upper, S_upper, 0, eos)

	scatter!(ax2, S, T, σ₀, markersize = 6, label = "S-T profile")
	σ₀_max, max_idx = findmax(σ₀)
	σ₀_min, min_idx = findmin(σ₀)

	S_minmax = [S[min_idx], S[max_idx]]
	T_minmax = [T[min_idx], T[max_idx]]

	scatter!(ax2, S_minmax[1], T_minmax[1], label = "Minimum density", color = :orange)
	scatter!(ax2, S_minmax[2], T_minmax[2], label = "Maximum density", color = :magenta)

	ρ_min = ρ(T[min_idx], S[min_idx], 0, eos)
	ρ_max = ρ(T[max_idx], S[max_idx], 0, eos)
	contour!(ax2, S_range, Θ_range, ρ_'; levels = [ρ_min, ρ_upper, ρ_star, ρ_max], colormap = :dense, label = "Isopycnals")

	scatter!(ax[2], σ₀_min, z[min_idx], label = "Minimum density", color = :orange)
	scatter!(ax[2], σ₀_max, z[max_idx], label = "Maximum density", color = :magenta)

	lines!(ax2, S_minmax, T_minmax, color = :red, linestyle = :dash, label = "Interface (?)")

	axislegend(ax2, position = :lt)
	fig

	ΔS_minmax = diff(S_minmax)[1]
	ΔT_minmax = diff(T_minmax)[1]
	slope = ΔT_minmax / ΔS_minmax

	Δρ = ρ(Θ_upper, S_upper, 0, eos) - ρ(Θ_star, S_star, 0, eos)
	md"""
	$(fig)

	Initial ``Δρ = `` $(round(Δρ, digits = 5)).
	I have adjusted the salinity in the upper layer to give the same ``Δρ`` initially.
	
	Slope connecting points of minimum and maximum density on the ``S - \Theta`` figure is $(round(slope[1], digits = 3)).
	"""
end

# ╔═╡ b1a1dd8b-756a-4eff-9607-5408c65dae4e
md"""
## Model output
$(LocalResource("../single_interface/density.mp4"))

$(LocalResource("../single_interface/tracers.mp4"))
"""

# ╔═╡ 7ad7693e-6f18-474a-89e8-b2d433aea261
TableOfContents()

# ╔═╡ Cell order:
# ╟─7d191973-92d1-4ed7-bcee-650f2185e36f
# ╟─b51bd896-64f6-11ef-3b2b-3d362d0c3f8a
# ╟─22c3e8cc-2d34-4fa4-9c2d-c807c6f53998
# ╟─a92c3315-ae34-447a-bc3e-7960bf2694d1
# ╟─8b0cb0d8-ad5b-46f8-8990-9438589976d8
# ╟─feb1de5f-8dfd-478d-89d9-ff668c55b892
# ╟─056df154-3413-4459-85c2-a1aeaeed9bdb
# ╟─28033f6f-2107-4b8d-a2fe-78e4eacd8d14
# ╟─78a2038d-2da0-43e9-9839-cdce03002927
# ╟─120a212f-f5c1-4520-802d-28ad737ebdee
# ╟─106ec4cd-a244-4457-b6b9-3b70a4b6db3e
# ╟─f2a9ff46-033c-4d12-94f8-3f12b4894a11
# ╟─054c0968-fbc4-4ec2-b9a1-fd7c005f49cc
# ╟─214a62b5-d46f-4636-bc26-3d5a4372d0d9
# ╟─f2163623-d32f-4203-b1ce-9f87478fc9b3
# ╠═88ad24ee-072c-4927-a21a-8159c80d75b0
# ╟─3face522-479b-4f87-a332-b858877094bc
# ╟─e5079cbd-7878-4e71-96e1-71f7ceccfb29
# ╟─c521d03e-929c-493e-ba8b-db51983a2c2a
# ╠═d7aa3853-6c71-46da-983a-2c1eb0c817a7
# ╟─7f417e03-72ef-47f5-923a-b932de687a37
# ╟─7b798e2e-4370-43fc-a3d7-26587aeec349
# ╟─7c657d45-cd7a-4106-b23b-77dccf5c982f
# ╟─308bc1ca-673d-48a9-8938-58c64c3d66a1
# ╟─ac2e0369-626f-40f4-a405-6941ff60f6d2
# ╟─b1a1dd8b-756a-4eff-9607-5408c65dae4e
# ╟─7ad7693e-6f18-474a-89e8-b2d433aea261
