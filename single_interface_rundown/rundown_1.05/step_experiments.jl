### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 010ecdc3-51d6-41a6-9bc5-6efbba0723a6
begin
	using Pkg
	Pkg.activate("../..")
	using JLD2, CairoMakie, PlutoUI, Dates, Statistics, TimeSeries
end

# ╔═╡ 6301138c-fa0b-11ef-0f3b-39dac35db063
begin
initial_state = @bind is Select(["step", "tanh", "diff_initial_h"])
md"""
# Single interface experiments

This notebook contains diagnostics computed from two single interface experiments.
Both experiments have the same initial salinity and temperature within each layer and the two layers meet at a sharp interface.
One has a linear equation of state the other uses 55 term non-linear eos.
Initialy velocity noise ``\mathcal{O}(10^{-2})`` is used to kick off the simulation and they are then run for 18 hours.

The diagnostics, from [McDougall (1981)](https://www.sciencedirect.com/science/article/abs/pii/0079661181900021) and [Carpenter el al (2012)](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/simulations-of-a-doublediffusive-interface-in-the-diffusive-convection-regime/63D2ECE2AA41439E01A01F9A0D76F2E2), are:
- salinity, temperature and density in each layer
- entrainment parameter (computed from conservative temperature)
- ``\Delta S``, ``\Delta T`` and ``R_{\rho}`` computed from top and bottom quarter of domain for upper and lower layer (respectively)
- salinity and temperature flux across interface
- interface height
- interface thickness

I have two versions: one with an initial step change an one with `tanh` profiles that change at different rates over the interface to setup the interface.
Can choose which one to view: $initial_state
"""
end

# ╔═╡ dc0440f3-b915-4974-9b3f-76b34e20b8e0
begin
	linear_animation_path = joinpath(@__DIR__, is*"/lineareos")
	nonlinear_animation_path = joinpath(@__DIR__, is*"/nonlineareos")
	dims = load(joinpath(is, is*"_diagnostics.jld2"), "dims")
	# leos = load(is*"_diagnostics.jld2", "lineareos")
	# nleos = load(joinpath(is, is*"_diagnostics.jld2"), "nonlineareos")
	@info "$is initial condition output loaded"
end

# ╔═╡ 38b12f7f-4d53-4302-a4ce-7e8c07d49ca9
begin
	eos_select = @bind eos Select(["Nonlinear", "Linear"])
	md"""
	# Equation of state
	
	I have the data saved for both linear and non linear equation of state. To save writing all the same code twice select the eos: $(eos_select)
	"""
end

# ╔═╡ e177c879-b7d0-4328-b5ad-776f8c64e050
begin
	# 	expt_data = eos == "Linear" ? leos : nleos
	# animation_path = eos == "Linear" ? linear_animation_path : nonlinear_animation_path
			expt_data = load(joinpath(is, is*"_diagnostics.jld2"))
	# animation_path = eos == "Linear" ? linear_animation_path : nonlinear_animation_path
	md"""
	## Animations
	
	### Density
	$(LocalResource(is*"/density.mp4"))
	
	### Salinity and temperature
	$(LocalResource(is*"/tracers.mp4"))
	"""
end

# ╔═╡ 07089057-5b2f-40e5-a485-0eeac1e9b348
md"""
## Diagnostics

### Salinity, temperature and density in each layer
"""

# ╔═╡ c2dce901-8578-448c-8c6e-ec7bb3e6d71b
let
	fig = Figure(size = (600, 1000))
	axT = Axis(fig[1, 1], xlabel = "Θₗ", ylabel = "Θᵤ")
	lines!(axT, expt_data["Tₗ_Tᵤ_ts"][:, 1], expt_data["Tₗ_Tᵤ_ts"][:, 2])
	scatter!(axT, expt_data["Tₗ_Tᵤ_ts"][1, 1], expt_data["Tₗ_Tᵤ_ts"][1, 2], label = "Initial values", color = :orange)
	axislegend(axT)
	axS = Axis(fig[2, 1], xlabel = "Sₗ", ylabel = "Sᵤ")
	lines!(axS, expt_data["Sₗ_Sᵤ_ts"][:, 1], expt_data["Sₗ_Sᵤ_ts"][:, 2])
	scatter!(axS, expt_data["Sₗ_Sᵤ_ts"][1, 1], expt_data["Sₗ_Sᵤ_ts"][1, 2], label = "Initial values", color = :orange)
	axρ = Axis(fig[3, 1], xlabel = "ρₗ", ylabel = "ρᵤ")
	lines!(axρ, expt_data["ρₗ_ρᵤ_ts"][:, 1], expt_data["ρₗ_ρᵤ_ts"][:, 2])
	scatter!(axρ, expt_data["ρₗ_ρᵤ_ts"][1, 1], expt_data["ρₗ_ρᵤ_ts"][1, 2], label = "Initial values", color = :orange)
	fig
end

# ╔═╡ 5581d9c6-197b-4901-bb4e-e515ef249836
begin
	window_slider =	@bind window PlutoUI.Slider(1:600, show_value = true)
md"""
### Entrainment

Top panel is the entrainment every minute, the bottom panel takes a windowed average with the window adjusted by this slider: $(window_slider)
"""
end

# ╔═╡ a6403686-dc8a-480d-9d77-82a0562e4665
let
	# mins = round.(Int64, dims["time"] ./ 60)
	# = is == "step" ? 1 : is == "tanh" ? 8 : 1
	# time_hours = 1
	# timestamps = Time(0, 1, 0):Minute(1):Time(time_hours, 0, 0) 
	R_ρ_interp = 0.5 * (expt_data["R_ρ"][1:end-1] .+ expt_data["R_ρ"][2:end])
	# ta = TimeArray((;timestamps, Rᵨ = R_ρ_interp, Ẽ = expt_data["Ẽ"]), timestamp = :timestamps)
	# ta_mean = moving(mean, ta, window)
	
	fig = Figure(size = (500, 500))
	ax1 = Axis(fig[1, 1], xlabel = "R_ρ", ylabel = "Ẽ")
	lines!(ax1, R_ρ_interp, expt_data["Ẽ"])
	# vlines!(ax1, 1.6, color = :red, linestyle = :dash)

	# ax2 = Axis(fig[2, 1], xlabel = "R_ρ", ylabel = "Ẽ")
	# lines!(ax2, values(ta_mean), label = "Averaging window = $window mins")
	# vlines!(ax2, 1.6, color = :red, linestyle = :dash)
	# axislegend(ax2)
	fig
end

# ╔═╡ bbdef33d-6493-4f95-ba92-92d08e75c69a
begin
md"""
### Salinity and temperature flux
Salinity flux looks wrong when using whole domain.
Also have horizontal average flux of salintiy and temperature.
"""
end

# ╔═╡ c852e2d3-f489-4b98-a183-dae4c3594947
let
	R_ρ_interp = 0.5 * (expt_data["R_ρ"][1:end-1] .+ expt_data["R_ρ"][2:end])
	fig = Figure(size = (800, 800))
	axT = Axis(fig[1, 1], ylabel = "T flux")
	lines!(axT, R_ρ_interp, expt_data["T_flux"][2, :])
	axS = Axis(fig[2, 1], ylabel = "S flux")
	lines!(axS, R_ρ_interp, expt_data["S_flux"][2, :])
	R_f = expt_data["S_flux"][2, :] ./ expt_data["T_flux"][2, :]
	axf = Axis(fig[3, 1], xlabel = "Rᵨ",ylabel = "R_f")
	lines!(axf, R_ρ_interp, R_f)
	fig
end

# ╔═╡ baa37deb-ea5a-4437-afa0-0e30a755df9c
let
	R_ρ_interp = 0.5 * (expt_data["R_ρ"][1:end-1] .+ expt_data["R_ρ"][2:end])

	fig = Figure(size = (800, 800))
	axT = Axis(fig[1, 1], ylabel = "T flux")
	T_interface_idx = expt_data["ha_T_interface_idx"]
	T_flux_interface = [expt_data["ha_T_flux"][idx, i] for (i, idx) ∈ enumerate(T_interface_idx)]
	lines!(axT, R_ρ_interp, T_flux_interface)
	
	axS = Axis(fig[2, 1], ylabel = "S flux")
	S_interface_idx = expt_data["ha_S_interface_idx"]
	S_flux_interface = [expt_data["ha_S_flux"][idx, i] for (i, idx) ∈ enumerate(S_interface_idx)]
	lines!(axS, R_ρ_interp, S_flux_interface)
	
	R_f = S_flux_interface ./ T_flux_interface
	axf = Axis(fig[3, 1], xlabel = "Rᵨ",ylabel = "R_f")
	lines!(axf, R_ρ_interp, R_f)

	a, b = [R_ρ_interp.^0 R_ρ_interp] \ R_f
	lines!(axf, R_ρ_interp, a .+ b .* R_ρ_interp, label = "Linear fit, slope = $(b)")
	axislegend(axf, position = :rb)
	
	Label(fig[0, 1], "Horizontally averaged flux through interface", tellwidth = false, font = :bold)
	fig
end

# ╔═╡ 0a9d245d-8285-4a22-9edf-178d9e85addb
md"""
### Interface thickness, height and change in concentation between layers
"""

# ╔═╡ 9a8041ad-6b12-4ef5-9f2f-44189de067f9
let
	timestamps = dims["time"][2:end] ./ 60
	fig = Figure(size = (800, 800))
	axhₜ = Axis(fig[1, 1], ylabel = "hₜ (cm)")
	lines!(axhₜ, timestamps, 100 * expt_data["hₜ"])
	axhₛ = Axis(fig[2, 1], ylabel = "hₛ (cm)")
	lines!(axhₛ, timestamps, 100 * expt_data["hₛ"])
	axr = Axis(fig[3, 1], xlabel = "time (mins)", ylabel = "r = hₜ / hₛ (cm)")
	lines!(axr, timestamps, expt_data["r"])
	fig
end

# ╔═╡ 3e422d6d-912f-4119-a290-648dbe036dde
let
	mins = dims["time"] ./ 60
	fig = Figure(size = (500, 500))
	ax = Axis(fig[1, 1], xlabel = "time (mins)")
	lines!(ax, mins, expt_data["ΔS"] ./ expt_data["ΔS"][1], color = :blue, label = "ΔS / ΔS₀")
	lines!(ax, mins, expt_data["ΔT"] ./ expt_data["ΔT"][1], color = :red, label = "ΔT / ΔT₀")
	ylims!(ax, 0, 1.1)
	ax2 = Axis(fig[1, 1], yaxisposition = :right,
			    yticklabelcolor = :dodgerblue,
			    rightspinecolor = :dodgerblue,
				ylabelcolor = :dodgerblue,
			    ytickcolor = :dodgerblue,
				ylabel = "Rᵨ")
	lines!(ax2, mins, expt_data["R_ρ"], color = :dodgerblue, label = "Rᵨ")
	axislegend(ax, merge = true, position = :rc)
	fig
end

# ╔═╡ d0148931-4198-4bb3-893c-a9b73e1ec7a9
let
	fig = Figure(size = (500, 500))
	ax = Axis(fig[1, 1], xlabel = "Rᵨ", ylabel = "z✶")
	lines!(ax, expt_data["R_ρ"][2:end], dims["z✶"][expt_data["S_interface_idx"]], label = "salinity")
	lines!(ax, expt_data["R_ρ"][2:end], dims["z✶"][expt_data["T_interface_idx"]], label = "temperature", linestyle = :dot)
	# ylims!(ax, 0.49, 0.54)
	vlines!(ax, 1.6, color = :red, linestyle = :dash)
	axislegend(ax, position = :rb)
	fig
end

# ╔═╡ 6ce43b6e-c3fa-408f-8702-900eaeb17bf5
md"""
# Length scales

## Batchelor length
"""

# ╔═╡ 4538f159-01d9-45fd-9fa5-d7463c506a77
begin
	η = 1e3 * expt_data["η"][2:end]
	Ba = 1e3 * expt_data["Ba"][2:end]
	Ba_fig, Ba_ax = lines(log10.(Ba), label = "Ba")
	lines!(Ba_ax, log10.(η), label = "η")
	Ba_ax.xlabel = "time (mins)"
	Ba_ax.ylabel = "Length (log10(mm))"
	min_Ba = minimum(Ba)
	min_η = minimum(η)
	Ba_ax.title = "η minimum = $(min_η)mm,  Ba minimum = $(min_Ba)mm"
	axislegend(Ba_ax)
	Ba_fig
end

# ╔═╡ d9422085-e838-44a1-91be-b81458dc3013
begin
	Lx = Ly = 0.07
	Nx = 35
	Δx = diff(dims["x_caa"])[1]
	Lz = 0.5
	Nz = 250
	Δz = diff(dims["z_aaf"])[1]
	md"""
	Above can see a figure of the Batchelor scale with minimal length $(round(min_Ba, digits = 2))mm.
	To achieve this in the domain size I have this would need resolution of $(min_Ba * 1e-3) everywhere so around 4e-4.
	This simulation was run with:
	- Δx = $(Δx)
	- Δz = $(Δz).

	So with what I have done I am an order of magnitiude away.
	But I have the leeway that people use ``Δ < 2.5 Ba`` as the upper limit so what I need to resolve is $(round(min_Ba, digits = 2) * 2.5), so $(round(min_Ba, digits = 2) * 2.5 * 1e-3)m.
	
	If I set ``N_{x} = N_{y} = 100`` and ``N_{z} = 700``, for ``L_{x} = L_{y} = `` $(Lx) and ``L_{z} = `` $(Lz), I should be able to get to DNS resolution with the 2.5Ba argument provided that I get closed energy budget.
	This resolution is *less* than what I ran with the cabbeling DNS so it should be possible --- just that so far this simulation has required a significantly smaller timestep.

	One other option is to look at a further reduction to the diffusivities then ramp up and provided everything is scaled correctly it should still work but will check with supervisors what they think first.
	"""
end

# ╔═╡ 3c0e1dfd-e4ba-448f-8475-ada056c8b5fe
md"""
## Energy budget

I am not quite sure why this is not closed -- could be resolution?
*Update:* does not appear to be resolution as this is more or less at Batchelor scale.
I would say it looks like spikes in ``\varepsilon`` except that this is smooth.
I have also computed the buoyancy flux in two ways and they are equal.
"""

# ╔═╡ c576c4cd-1101-46e2-b6fa-b574f0b13dfe
let
	Δt = diff(dims["time"])
	ek = expt_data["∫Eₖ"]
	dₜek = diff(ek) ./ Δt
	ε = 0.5 * (expt_data["∫ε"][1:end-1] .+ expt_data["∫ε"][2:end])
	∫wb = 0.5 * (expt_data["∫wb"][1:end-1] .+ expt_data["∫wb"][2:end])
	∫gρw = 0.5 * (expt_data["∫gρw"][1:end-1] .+ expt_data["∫gρw"][2:end])
	RHS = ∫wb .- ε
	fig, ax = lines(eachindex(Δt)[2:end], dₜek[2:end], label = "dₜek")
	lines!(ax, eachindex(Δt)[2:end], RHS[2:end], label = "∫wb - ε")
	ax.title = "Energy  budget"
	ax.xlabel = "time (minutes)"
	ax.ylabel = "Watts / ρ₀"
	axislegend(ax, position = :rb)

	ax2 = Axis(fig[2, 1], title = "Absolute error")
	abs_err = abs.(dₜek .- RHS)
	lines!(ax2, abs_err[2:end])
	fig
end

# ╔═╡ f200b8e0-2b14-4270-963b-6bb1b154d550
let
	fig, ax = lines(expt_data["∫wb"], label = "wb")
	lines!(ax, -expt_data["∫gρw"], label = "∫gρw (post processing)", linestyle = :dash)
	axislegend(ax)
	fig
end

# ╔═╡ 50e87efc-a49c-4ffd-bfbd-cd5dfad40639
md"""
# Setup for boundary conditions

To get a model that does not run down it seems like best option is including some kind of boundary conditions.
Ideally we would use periodic or jump boundary conditions but nonlinear eos causes problems with this.
"""

# ╔═╡ ee9c0edb-477b-4cc0-8c57-36845a90bbaf
@bind Rᵨ_val PlutoUI.Slider(round.(expt_data["R_ρ"][:], digits =3), show_value = true, default=1.415)

# ╔═╡ 68a0a47e-e919-4d9d-b1a5-090d69bf633e
begin
	find_Rᵨ = findfirst(expt_data["R_ρ"] .> Rᵨ_val)
	Shaflux = expt_data["ha_S_flux"][expt_data["ha_S_interface_idx"][find_Rᵨ], find_Rᵨ]
	Thaflux = expt_data["ha_T_flux"][expt_data["ha_T_interface_idx"][find_Rᵨ], find_Rᵨ]
	md"""
	Horizontally averaged salinity flux through interface for ``R_{\rho} = `` $(Rᵨ_val) is ``J_{S} = `` $(round(Shaflux, digits = 10)) with 
	Horizontally averaged temperature flux through interface for ``R_{\rho} = `` $(Rᵨ_val) is ``J_{T} = `` $(round(Thaflux, digits = 7)).
	"""
end

# ╔═╡ 963fa274-2d8f-47fd-b227-4d7b3275d7ad
TableOfContents()

# ╔═╡ Cell order:
# ╟─6301138c-fa0b-11ef-0f3b-39dac35db063
# ╟─010ecdc3-51d6-41a6-9bc5-6efbba0723a6
# ╟─dc0440f3-b915-4974-9b3f-76b34e20b8e0
# ╟─38b12f7f-4d53-4302-a4ce-7e8c07d49ca9
# ╟─e177c879-b7d0-4328-b5ad-776f8c64e050
# ╟─07089057-5b2f-40e5-a485-0eeac1e9b348
# ╟─c2dce901-8578-448c-8c6e-ec7bb3e6d71b
# ╟─5581d9c6-197b-4901-bb4e-e515ef249836
# ╟─a6403686-dc8a-480d-9d77-82a0562e4665
# ╟─bbdef33d-6493-4f95-ba92-92d08e75c69a
# ╟─c852e2d3-f489-4b98-a183-dae4c3594947
# ╟─baa37deb-ea5a-4437-afa0-0e30a755df9c
# ╟─0a9d245d-8285-4a22-9edf-178d9e85addb
# ╟─9a8041ad-6b12-4ef5-9f2f-44189de067f9
# ╟─3e422d6d-912f-4119-a290-648dbe036dde
# ╟─d0148931-4198-4bb3-893c-a9b73e1ec7a9
# ╟─6ce43b6e-c3fa-408f-8702-900eaeb17bf5
# ╟─4538f159-01d9-45fd-9fa5-d7463c506a77
# ╟─d9422085-e838-44a1-91be-b81458dc3013
# ╟─3c0e1dfd-e4ba-448f-8475-ada056c8b5fe
# ╟─c576c4cd-1101-46e2-b6fa-b574f0b13dfe
# ╟─f200b8e0-2b14-4270-963b-6bb1b154d550
# ╟─50e87efc-a49c-4ffd-bfbd-cd5dfad40639
# ╟─ee9c0edb-477b-4cc0-8c57-36845a90bbaf
# ╟─68a0a47e-e919-4d9d-b1a5-090d69bf633e
# ╟─963fa274-2d8f-47fd-b227-4d7b3275d7ad
