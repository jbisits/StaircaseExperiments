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

# ╔═╡ 010ecdc3-51d6-41a6-9bc5-6efbba0723a6
begin
	using Pkg
	Pkg.activate("../..")
	using JLD2, CairoMakie, PlutoUI, Dates, Statistics, TimeSeries
	using SpecialFunctions: erf
end

# ╔═╡ 6301138c-fa0b-11ef-0f3b-39dac35db063
begin
# initial_state = @bind is Select(["step", "tanh", "diff_initial_h"])
is = "step"
md"""
# Single interface experiments

This notebook contains diagnostics computed from two single interface experiments.
Both experiments have the same initial salinity and temperature within each layer and the two layers meet at a sharp interface.
One has a linear equation of state the other uses 55 term non-linear eos.

The diagnostics, from [McDougall (1981)](https://www.sciencedirect.com/science/article/abs/pii/0079661181900021) and [Carpenter el al (2012)](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/simulations-of-a-doublediffusive-interface-in-the-diffusive-convection-regime/63D2ECE2AA41439E01A01F9A0D76F2E2), are:
- salinity, temperature and density in each layer
- entrainment parameter (computed from conservative temperature)
- ``\Delta S``, ``\Delta T`` and ``R_{\rho}`` computed from top and bottom quarter of domain for upper and lower layer (respectively)
- salinity and temperature flux across interface
- interface height
- interface thickness
"""
end

# ╔═╡ 68d31cca-3f29-4402-ac79-8deaef98ef50
begin
	eos_select = @bind eos Select(["higher_res_nonlinear", "higher_res_linear", "R_rho_1.4_nonlinear", "largerdiffrationonlinear"])
	md"""
	# Equation of state
	
	I have the data saved for both linear and non linear equation of state. 
	To save writing all the same code twice select the eos: $(eos_select)
	"""
end

# ╔═╡ 087d2583-ee90-437a-97ec-0ab607337e30
begin
	expt_eos = eos*"eos"
	output_path = joinpath(@__DIR__, is, expt_eos)
	expt_data = load(joinpath(output_path, is*"_diagnostics.jld2"))
	dims = load(joinpath(output_path, is*"_diagnostics.jld2"), "dims")
	@info "$is initial condition with $(eos) eos output loaded"
end

# ╔═╡ e177c879-b7d0-4328-b5ad-776f8c64e050
begin
	md"""
	## Animations
	
	### Density
	$(LocalResource(joinpath(output_path, "density_Nsquared.mp4")))
	
	### Salinity and temperature
	$(LocalResource(joinpath(output_path,"tracers.mp4")))
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
	# R_ρ_interp = 0.5 * (expt_data["R_ρ"][1:end-1] .+ expt_data["R_ρ"][2:end])
	R_ρ_interp = dims["time"][1:end-1]
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

# ╔═╡ c3f03eaf-0c45-477a-ba2b-c411be6d07c8
begin
	R_ρ_interp = 0.5 * (expt_data["R_ρ"][1:end-1] .+ expt_data["R_ρ"][2:end])
	
	T_interface_idx = expt_data["T_ha_interface_idx"]
	T_flux_interface = [expt_data["T_ha_flux"][idx, i] for (i, idx) ∈ enumerate(T_interface_idx)]
	a_T, b_T = [R_ρ_interp[12:end].^0 R_ρ_interp[12:end]] \ T_flux_interface[12:end]
	
	S_interface_idx = expt_data["S_ha_interface_idx"]
	S_flux_interface = [expt_data["S_ha_flux"][idx, i] for (i, idx) ∈ enumerate(S_interface_idx)]
	lfit_T_flux = a_T .+ b_T .* R_ρ_interp[12:end]
	a_S, b_S = [R_ρ_interp[12:end].^0 R_ρ_interp[12:end]] \ S_flux_interface[12:end]
	lfit_S_flux = a_S .+ b_S .* R_ρ_interp[12:end]
	
	nothing
end

# ╔═╡ baa37deb-ea5a-4437-afa0-0e30a755df9c
let
	# R_ρ_interp = dims["time"][1:end-1] ./ 60
	fig = Figure(size = (800, 800))
	axT = Axis(fig[1, 1], ylabel = "T flux")
	lines!(axT, R_ρ_interp, T_flux_interface)

	lines!(axT, R_ρ_interp[12:end], lfit_T_flux, label = "Linear fit, slope = $(b_T)")
	axislegend(axT, position = :rb)
	
	axS = Axis(fig[2, 1], ylabel = "S flux")
	lines!(axS, R_ρ_interp, S_flux_interface)

	lines!(axS, R_ρ_interp[12:end], lfit_S_flux, label = "Linear fit, slope = $(b_S)")
	axislegend(axS, position = :rb)
	
	R_f = S_flux_interface ./ T_flux_interface
	axf = Axis(fig[3, 1], xlabel = "Rᵨ",ylabel = "R_f")
	lines!(axf, R_ρ_interp, R_f)

	a, b = [R_ρ_interp[12:end].^0 R_ρ_interp[12:end]] \ R_f[12:end]
	lines!(axf, R_ρ_interp[12:end], a .+ b .* R_ρ_interp[12:end], label = "Linear fit, slope = $(b)")
	lines!(axf, R_ρ_interp[12:end], lfit_S_flux ./ lfit_T_flux, linestyle = :dash, label = "Ratio of linear fits")
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
	lines!(ax, mins[2:end], expt_data["ΔS"] ./ expt_data["ΔS"][1], color = :blue, label = "ΔS / ΔS₀")
	lines!(ax, mins[2:end], expt_data["ΔT"] ./ expt_data["ΔT"][1], color = :red, label = "ΔT / ΔT₀")
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
	# vlines!(ax, 1.6, color = :red, linestyle = :dash)
	axislegend(ax, position = :rb)
	fig
end

# ╔═╡ 26776c6e-2864-4b0f-8ddd-d3fb16aa8779
begin
	figinterface = Figure(size = (500, 500))
	axinterface = Axis(figinterface[1, 1], xlabel = "Rᵨ", ylabel = "z")
	
	findS = [findfirst(expt_data["S_ha"][:, i] .≤ 0.5 * (34.58 + 34.7)) for i in eachindex(expt_data["S_ha"][1, 2:end])]
	S_interface = dims["z_aac"][findS]
	lines!(axinterface, expt_data["R_ρ"][2:end], S_interface, label = "salinity")
	
	findT = [findfirst(expt_data["T_ha"][:, i] .≤ 0.5 * (-1.5 + 0.5)) for i in eachindex(expt_data["T_ha"][1, 2:end])]
	T_interface = dims["z_aac"][findT]
	lines!(axinterface, expt_data["R_ρ"][2:end], T_interface, label = "temperature", linestyle = :dot)
	
	ylims!(axinterface, -0.25, -0.24)
	# vlines!(ax, 1.6, color = :red, linestyle = :dash)
	axislegend(axinterface, position = :rb)
	figinterface
end

# ╔═╡ 82964820-7220-4e21-b263-20754b6a3a33
md"""
# Salinity-temperature space

**Note:** the below will only run with later experiments where the horizontal average salinity and temperature profiles are saved.

This is trying to capture asymmetry in density difference and compare it to theory (most of which can be found in the `diffusive_interfaces` Pluto notebook).
"""

# ╔═╡ 2536f46d-bbe8-4f85-a16a-82afef16fef5
let
	fig = Figure(size = (1000, 1000))
	t, z = dims["time"] / 60, dims["z_aac"]
	axS = Axis(fig[1, 1], title = "Hovmoller for saliniy and temperature (ha profile)", ylabel = "z (m)")
	hmS = heatmap!(axS, t, z, expt_data["S_ha"]', colormap = :haline)
	Colorbar(fig[1, 2], hmS)
	axT = Axis(fig[2, 1], xlabel = "time (mins)", ylabel = "z (m)")
	hmT = heatmap!(axT, t, z, expt_data["T_ha"]', colormap = :thermal)
	Colorbar(fig[2, 2], hmT)
	fig
end

# ╔═╡ 3130c878-fda2-4d11-a658-748d6a15b2b8
let
	fig = Figure(size = (1000, 1000))
	t, z = dims["time"] / 60, dims["z_aac"]
	axS = Axis(fig[1, 1], title = "Hovmoller for anomaly saliniy and temperature (ha profile)", ylabel = "z (m)")
	hmS = heatmap!(axS, t, z, (expt_data["S_ha"] .- expt_data["S_ha"][:, 1])', colormap = :curl, colorrange = (-0.09, 0.09))
	Colorbar(fig[1, 2], hmS)
	axT = Axis(fig[2, 1], xlabel = "time (mins)", ylabel = "z (m)")
	hmT = heatmap!(axT, t, z, (expt_data["T_ha"] .- expt_data["T_ha"][:, 1])', colormap = :delta, colorrange = (-1.5, 1.5))
	Colorbar(fig[2, 2], hmT)
	fig
end

# ╔═╡ ca6991b4-ac76-435e-bcff-82103b6abdc7
erf_tracer_solution(z, Cₗ::Number, ΔC::Number, κ::Number, time, interface_location) =
    Cₗ + 0.5 * ΔC * (1 + erf((z - interface_location) / sqrt(4 * κ * time)))

# ╔═╡ b1fd4a61-796e-4b39-b175-971311c9c62a
@bind t PlutoUI.Slider(eachindex(dims["time"]))

# ╔═╡ 6e7d38ab-2824-474b-b90c-75a3a5a05e57
md"""
This is the rundown case here but can see the not well mixed layers here I think with the less vertical S-T relationship at the start and end of the curves.
I can still look at the ratio between lower and upper layer density differences by taking maximum, minimum and the top and bottom values.
From theory, for linear equation of state this should be around 1 and for non-linear equation of state this is around 0.6.
"""

# ╔═╡ e972f242-3b58-4581-87ae-437533b9fba1
begin
	σ_ha = expt_data["σ_ha"]
	R_Δσ = Array{Float64}(undef, length(σ_ha[1, :]))
	for (i, col) ∈ enumerate(eachcol(σ_ha))
		σ_min = minimum(col)
		σ_top = col[end]
		Δσ_upper = abs(σ_min - σ_top)
		σ_max = maximum(col)
		σ_bottom = col[1]
		Δσ_lower = abs(σ_max - σ_bottom)
		R_Δσ[i] = Δσ_upper / Δσ_lower
	end
end

# ╔═╡ 8711fdef-0da7-46bf-aa82-2f32b0590f7b
let
	Sₜ, Tₜ = expt_data["S_ha"][:, t], expt_data["T_ha"][:, t]
	fig, ax = lines(Sₜ, Tₜ, label = "Model output")
	Slims = extrema(expt_data["S_ha"][:, 1]) .+ [-0.01, 0.01]
	Tlims = extrema(expt_data["T_ha"][:, 1]) .+ [-0.1, 0.1]
	xlims!(Slims...)
	ylims!(Tlims...)
	ax.xlabel = "Salinity (gkg⁻¹)"
	ax.ylabel = "Temperature (°C)"
	ax.title = "Ha S and T profiles at time t = $(dims["time"][t] / 60)min"
	ax.subtitle = "R_Δσ = $(round(R_Δσ[t], digits = 1))"

	κₛ, κₜ = expt_data["attrib/κₛ (m²s⁻¹)"], expt_data["attrib/κₜ (m²s⁻¹)"]
	z = dims["z_aac"]
	ΔS = Sₜ[end] - Sₜ[1]
	ΔT = Tₜ[end] - Tₜ[1]
	id = expt_data["attrib/interface_depth"]
	interfaceS = vcat(id, S_interface)
	interfaceT = vcat(id, T_interface)
	S = erf_tracer_solution.(z, Sₜ[1], ΔS, κₛ, t, interfaceT[t])
	T = erf_tracer_solution.(z, Tₜ[1], ΔT, κₜ, t, interfaceT[t])
	lines!(ax, S, T, color = :orange, label = "Theoretical model")

	axislegend(ax, position = :lt)

	σₜ = expt_data["σ_ha"][:, t] .- mean(expt_data["σ_ha"][:, t])
	ax2 = Axis(fig[1, 2], 
			   title = "Ha σ′ profile at time t = $(dims["time"][t] / 60)min", 
			   subtitle = "Anomaly from mean",
			   xlabel = "σ₀′", 
			   ylabel = "z (m)")
	lines!(ax2, σₜ, z)
	σ_lims = extrema(expt_data["σ_ha"])
	xlims!(ax2, (-0.05, 0.05))
	fig
end

# ╔═╡ d2e81b8b-4a1c-4330-8f2a-14a502390bcd
let
	fig, ax = lines(dims["time"] ./ 60, R_Δσ)
	ax.title = "Raito of density difference between upper and lower layer"
	ax.xlabel = "time (mins)"
	ax.ylabel = "R_Δσ"
	fig
end

# ╔═╡ 6ce43b6e-c3fa-408f-8702-900eaeb17bf5
md"""
# Length scales
"""

# ╔═╡ 4538f159-01d9-45fd-9fa5-d7463c506a77
begin
	η = 1e3 * expt_data["η"][3:end]
	Ba = 1e3 * expt_data["Ba"][3:end]
	Ba_fig, Ba_ax = lines(log10.(η), label = "η")
	lines!(Ba_ax, log10.(Ba), label = "Ba")
	Ba_ax.xlabel = "time (mins)"
	Ba_ax.ylabel = "Length (log10(mm))"
	min_Ba = minimum(Ba)
	min_η = minimum(η)
	Ba_ax.title = "η minimum = $(min_η)mm,  Ba minimum = $(min_Ba)mm"
	axislegend(Ba_ax, position = :rc)
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
	But I have the leeway that people use ``Δ < 2.5 Ba`` as the upper limit so what I need to resolve is $(round(min_Ba, digits = 2) * 2.5)mm, so $(round(min_Ba, digits = 2) * 2.5 * 1e-3)m.
	
	If I set ``N_{x} = N_{y} = 100`` and ``N_{z} = 700``, for ``L_{x} = L_{y} = `` $(Lx) and ``L_{z} = `` $(Lz), I should be able to get to DNS resolution with the 2.5Ba argument provided that I get closed energy budget.
	This resolution is *less* than what I ran with the cabbeling DNS so it should be possible --- just that so far this simulation has required a significantly smaller timestep.

	One other option is to look at a further reduction to the diffusivities then ramp up and provided everything is scaled correctly it should still work but will check with supervisors what they think first.
	"""
end

# ╔═╡ 3c0e1dfd-e4ba-448f-8475-ada056c8b5fe
md"""
# Energy budget

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
	# RHS =-∫gρw .- ε
	fig, ax = lines(eachindex(Δt)[1:end], dₜek[1:end], label = "dₜek")
	lines!(ax, eachindex(Δt)[1:end], RHS[1:end], label = "∫wb - ε")
	ax.title = "Energy  budget"
	ax.xlabel = "time (minutes)"
	ax.ylabel = "Watts / ρ₀"
	axislegend(ax, position = :rb)

	abs_err = abs.(dₜek .- RHS)[3:end]
	ax2 = Axis(fig[2, 1], title = "Absolute error, MAE = $(mean(abs_err))")
	lines!(ax2, abs_err[3:end])
	fig
end

# ╔═╡ f200b8e0-2b14-4270-963b-6bb1b154d550
let
	fig, ax = lines(log10.(abs.((expt_data["∫wb"][3:end]))), label = "wb")
	lines!(ax,  log10.(abs.(-expt_data["∫gρw"][3:end])), label = "∫gρw (post processing)", linestyle = :dash)
	lines!(ax, log10.(expt_data["∫ε"][3:end]), label = "∫ε", linestyle = :dot)
	ax.title = "Buoyancy flux and TKE dissipation (log10)"
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
	Shaflux = expt_data["S_ha_flux"][expt_data["S_ha_interface_idx"][find_Rᵨ], find_Rᵨ]
	Thaflux = expt_data["T_ha_flux"][expt_data["T_ha_interface_idx"][find_Rᵨ], find_Rᵨ]
	int_R_f = Shaflux / Thaflux
	Shaflux_fit = a_S + b_S * Rᵨ_val
	Thaflux_fit = a_T + b_T * Rᵨ_val
	int_R_f_fit = Shaflux_fit / Thaflux_fit
	md"""
	Horizontally averaged fluxes through interface for ``R_{\rho} = `` $(Rᵨ_val): 
	- ``J_{S} = `` $(round(Shaflux, digits = 10)); and
	- ``J_{T} = `` $(round(Thaflux, digits = 7)).
	with flux ratio: $(int_R_f).

	Fluxes through interface from linear fit:
	- ``J_{S} = `` $(round(Shaflux_fit, digits = 10)); and
	- ``J_{T} = `` $(round(Thaflux_fit, digits = 7)).
	with flux ratio: $(int_R_f_fit).
	"""
end

# ╔═╡ 963fa274-2d8f-47fd-b227-4d7b3275d7ad
TableOfContents()

# ╔═╡ Cell order:
# ╟─6301138c-fa0b-11ef-0f3b-39dac35db063
# ╟─010ecdc3-51d6-41a6-9bc5-6efbba0723a6
# ╟─68d31cca-3f29-4402-ac79-8deaef98ef50
# ╟─087d2583-ee90-437a-97ec-0ab607337e30
# ╟─e177c879-b7d0-4328-b5ad-776f8c64e050
# ╟─07089057-5b2f-40e5-a485-0eeac1e9b348
# ╟─c2dce901-8578-448c-8c6e-ec7bb3e6d71b
# ╟─5581d9c6-197b-4901-bb4e-e515ef249836
# ╟─a6403686-dc8a-480d-9d77-82a0562e4665
# ╟─bbdef33d-6493-4f95-ba92-92d08e75c69a
# ╟─c852e2d3-f489-4b98-a183-dae4c3594947
# ╟─c3f03eaf-0c45-477a-ba2b-c411be6d07c8
# ╟─baa37deb-ea5a-4437-afa0-0e30a755df9c
# ╟─0a9d245d-8285-4a22-9edf-178d9e85addb
# ╟─9a8041ad-6b12-4ef5-9f2f-44189de067f9
# ╟─3e422d6d-912f-4119-a290-648dbe036dde
# ╟─d0148931-4198-4bb3-893c-a9b73e1ec7a9
# ╟─26776c6e-2864-4b0f-8ddd-d3fb16aa8779
# ╟─82964820-7220-4e21-b263-20754b6a3a33
# ╟─2536f46d-bbe8-4f85-a16a-82afef16fef5
# ╟─3130c878-fda2-4d11-a658-748d6a15b2b8
# ╟─ca6991b4-ac76-435e-bcff-82103b6abdc7
# ╟─b1fd4a61-796e-4b39-b175-971311c9c62a
# ╟─8711fdef-0da7-46bf-aa82-2f32b0590f7b
# ╟─6e7d38ab-2824-474b-b90c-75a3a5a05e57
# ╟─e972f242-3b58-4581-87ae-437533b9fba1
# ╟─d2e81b8b-4a1c-4330-8f2a-14a502390bcd
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
