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
choose_expt = @bind expt Select(["1e-2", "1e-8", "no noise"])
md"""
# McDougall expt 25 October 1979

This notebook shows the diagnostics listed below for the experiment ran on 25th October in [McDougall (1981)](https://www.sciencedirect.com/science/article/abs/pii/0079661181900021).
Perhaps it is because of the smaller domain size (only horizontally) or the initial noise is too intense but we have not got the best match between these two experiments, especially the salinity in the upper vs lower layer.
Currently I am not sure why this is but density and temperature match well.
The entrainment shows a similar patter to the lab experiment but is smaller (though I need to check I identified the correct expt) and the **salinity and temperature flux just look wrong.**
I thought I had used this function successfully but clearly I need to revisit what it is doing.
It does produce something sensible for the temperature in my other experiments though.

The diagnostics, from [McDougall (1981)](https://www.sciencedirect.com/science/article/abs/pii/0079661181900021) and [Carpenter el al (2012)](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/simulations-of-a-doublediffusive-interface-in-the-diffusive-convection-regime/63D2ECE2AA41439E01A01F9A0D76F2E2), are:
- salinity, temperature and density in each layer
- entrainment parameter (computed from conservative temperature)
- ``\Delta S``, ``\Delta T`` and ``R_{\rho}`` computed from top and bottom quarter of domain for upper and lower layer (respectively)
- salinity and temperature flux across interface
- interface height
- interface thickness.

I have run the same experiment twice with different noise magnitude because I am having issues with the flux calculation and I think it is because the turbulence is so vigorous in part due to the random noise.
Noise magnitude $(choose_expt).
"""
end

# ╔═╡ dc0440f3-b915-4974-9b3f-76b34e20b8e0
begin
	dims = load("diagnostics.jld2", "dims")
	expt_data = expt == "1e-2" ? load("diagnostics.jld2", "noise") : expt == "1e-8" ? load("diagnostics.jld2", "lessnoise") : load("diagnostics.jld2", "nonoise")
	animation_path = @__DIR__
	@info "Output loaded"
end

# ╔═╡ e177c879-b7d0-4328-b5ad-776f8c64e050
begin
	md"""
	## Animations

	### Density
	$(LocalResource(joinpath(animation_path, "density.mp4")))

	### Salinity and temperature
	$(LocalResource(joinpath(animation_path, "tracers.mp4")))
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
	lines!(axS, expt_data["Sₗ_Sᵤ_ts"][:, 1], expt_data["Sₗ_Sᵤ_ts"][:, 2], color = expt_data["R_ρ"])
	scatter!(axS, expt_data["Sₗ_Sᵤ_ts"][1, 1], expt_data["Sₗ_Sᵤ_ts"][1, 2], label = "Initial values")
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
	sim_hours = expt == "1e-2" ? 8 : expt == "1e-8" ? 4 : 1
	timestamps = Time(0, 0, 0):Minute(1):Time(sim_hours-1, 59, 0)
	R_ρ_interp = 0.5 * (expt_data["R_ρ"][1:end-1] .+ expt_data["R_ρ"][2:end])
	ta = TimeArray((;timestamps, Rᵨ = R_ρ_interp, Ẽ = expt_data["Ẽ"]), timestamp = :timestamps)
	ta_mean = moving(mean, ta, window)

	fig = Figure(size = (600, 800))
	ax1 = Axis(fig[1, 1], xlabel = "R_ρ", ylabel = "Ẽ")
	lines!(ax1, R_ρ_interp, expt_data["Ẽ"])
	vlines!(ax1, 1.6, color = :red, linestyle = :dash)

	ax2 = Axis(fig[2, 1], xlabel = "R_ρ", ylabel = "Ẽ")
	lines!(ax2, values(ta_mean), label = "Averaging window = $window mins")
	vlines!(ax2, 1.6, color = :red, linestyle = :dash)
	axislegend(ax2)
	fig
end

# ╔═╡ bbdef33d-6493-4f95-ba92-92d08e75c69a
md"""
### Salinity and temperature flux
Both of these look wrong now.
"""

# ╔═╡ c852e2d3-f489-4b98-a183-dae4c3594947
let
	R_ρ_interp = 0.5 * (expt_data["R_ρ"][1:end-1] .+ expt_data["R_ρ"][2:end])
	fig = Figure(size = (800, 800))
	axT = Axis(fig[1, 1], ylabel = "T flux")
	lines!(axT, R_ρ_interp, expt_data["T_flux"][3, :])
	axS = Axis(fig[2, 1], ylabel = "S flux")
	lines!(axS, R_ρ_interp, expt_data["S_flux"][3, :])
	R_f = expt_data["S_flux"][3, :] ./ expt_data["T_flux"][3, :]
	axf = Axis(fig[3, 1], xlabel = "Rᵨ",ylabel = "R_f")
	lines!(axf, R_ρ_interp, R_f)
	fig
end

# ╔═╡ fb112b31-b30a-4a37-8593-46506114b32f
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
	# lines!(axf, R_ρ_interp[100:end], 0.01/2 .* R_ρ_interp[100:end])
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
	timestamps = eachindex(expt_data["R_ρ"][2:end])
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
	mins = eachindex(expt_data["R_ρ"])
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
	lines!(ax, expt_data["R_ρ"][2:end], dims["z✶"][expt_data["T_interface_idx"]], label = "salinity")
	lines!(ax, expt_data["R_ρ"][2:end], dims["z✶"][expt_data["T_interface_idx"]], label = "temperature", linestyle = :dot)
	vlines!(ax, 1.6, color = :red, linestyle = :dash)
	axislegend(ax, position = :rb)
	fig
end

# ╔═╡ 963fa274-2d8f-47fd-b227-4d7b3275d7ad
TableOfContents()

# ╔═╡ Cell order:
# ╟─6301138c-fa0b-11ef-0f3b-39dac35db063
# ╟─010ecdc3-51d6-41a6-9bc5-6efbba0723a6
# ╟─dc0440f3-b915-4974-9b3f-76b34e20b8e0
# ╟─e177c879-b7d0-4328-b5ad-776f8c64e050
# ╟─07089057-5b2f-40e5-a485-0eeac1e9b348
# ╟─c2dce901-8578-448c-8c6e-ec7bb3e6d71b
# ╟─5581d9c6-197b-4901-bb4e-e515ef249836
# ╟─a6403686-dc8a-480d-9d77-82a0562e4665
# ╟─bbdef33d-6493-4f95-ba92-92d08e75c69a
# ╟─c852e2d3-f489-4b98-a183-dae4c3594947
# ╟─fb112b31-b30a-4a37-8593-46506114b32f
# ╟─0a9d245d-8285-4a22-9edf-178d9e85addb
# ╟─9a8041ad-6b12-4ef5-9f2f-44189de067f9
# ╟─3e422d6d-912f-4119-a290-648dbe036dde
# ╟─d0148931-4198-4bb3-893c-a9b73e1ec7a9
# ╟─963fa274-2d8f-47fd-b227-4d7b3275d7ad
