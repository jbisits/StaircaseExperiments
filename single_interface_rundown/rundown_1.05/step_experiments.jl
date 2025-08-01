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
	using JLD2, CairoMakie, PlutoUI, Dates, Statistics, TimeSeries, GibbsSeaWater
	using SpecialFunctions: erf
	using SeawaterPolynomials:  total_density, TEOS10EquationOfState, thermal_expansion, haline_contraction
	using StaircaseShenanigans: CustomLinearEquationOfState, compute_R_ρ
end

# ╔═╡ 6301138c-fa0b-11ef-0f3b-39dac35db063
begin
# initial_state = @bind is Select(["step", "tanh", "diff_initial_h"])
is = "step"
md"""
# Single interface rundown experiments

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

Since starting this work new parts/diagnostics have been added which are documented below.
"""
end

# ╔═╡ 68d31cca-3f29-4402-ac79-8deaef98ef50
begin
	eos_select = @bind eos Select(["dns_res_dT2_nonlinear", 
								   "dns_res_dT2_linear",
								   "dns_res_dT1_nonlinear", 
								   "dns_res_dT1_linear", 
								   "dns_res_dT05_nonlinear", 
								   "dns_res_dT05_linear",
								   "dns_res_dT05_j110_nonlinear",
								   "equaldiffs_nonlinear",
								   "equaldiffs_nu1_nonlinear",
								   "R_rho_1.03_deltatheta_0.5_nonlinear", 
								   "R_rho_1.05_deltatheta_1.0_nonlinear", 
								   "R_rho_1.05_deltatheta_0.5_linear", 
								   "R_rho_1.4_nonlinear", "R_rho_1.67_nonlinear", "
								   R_rho_1.76_nonlinear", "R_rho_2.2_nonlinear", 
								   "R_rho_2.2_linear", "Tau_0.01_nonlinear", 
								   "_nonlinear"])
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
	ρ₀ = gsw_rho(34.7, 0.5, 0)
	_idx = findlast('_', eos)+1
	eos_type = eos[_idx:end] .== "nonlinear" ? TEOS10EquationOfState(reference_density = ρ₀) :
	CustomLinearEquationOfState(-0.5, 34.64, reference_density = ρ₀)
	α_sign = ifelse(eos[_idx:end] .== "nonlinear", +, -)
	@info "$is initial condition with $(eos) eos output loaded"
end

# ╔═╡ 700b8f65-83af-4e13-96b4-7192cbff4e86
0.05^2

# ╔═╡ b3fd9c9c-24d0-4ed8-a952-aaed2581045d
expt_data["attrib/RaS"], expt_data["attrib/RaT"]

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
	R_ρ_interp = 0.5 * (expt_data["R_ρ"][1:end-1] .+ expt_data["R_ρ"][2:end])
	# R_ρ_interp = dims["time"][1:end-1]
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
	# R_ρ_interp = 1:180
	start_flux = 30
	T_interface_idx = expt_data["T_ha_interface_idx"]
	T_flux_interface = [expt_data["T_ha_flux"][idx, i] for (i, idx) ∈ enumerate(T_interface_idx)]
	a_T, b_T = [R_ρ_interp[start_flux:end].^0 R_ρ_interp[start_flux:end]] \ T_flux_interface[start_flux:end]
	
	S_interface_idx = expt_data["S_ha_interface_idx"]
	S_flux_interface = [expt_data["S_ha_flux"][idx, i] for (i, idx) ∈ enumerate(S_interface_idx)]
	lfit_T_flux = a_T .+ b_T .* R_ρ_interp[start_flux:end]
	a_S, b_S = [R_ρ_interp[start_flux:end].^0 R_ρ_interp[start_flux:end]] \ S_flux_interface[start_flux:end]
	lfit_S_flux = a_S .+ b_S .* R_ρ_interp[start_flux:end]
	
	nothing
end

# ╔═╡ baa37deb-ea5a-4437-afa0-0e30a755df9c
let
	# R_ρ_interp = dims["time"][1:end-1] ./ 60
	fig = Figure(size = (800, 800))
	axT = Axis(fig[1, 1], ylabel = "T flux")
	lines!(axT, R_ρ_interp, T_flux_interface)

	lines!(axT, R_ρ_interp[start_flux:end], lfit_T_flux, label = "Linear fit, slope = $(b_T)")
	axislegend(axT, position = :rb)
	
	axS = Axis(fig[2, 1], ylabel = "S flux")
	lines!(axS, R_ρ_interp, S_flux_interface)

	lines!(axS, R_ρ_interp[start_flux:end], lfit_S_flux, label = "Linear fit, slope = $(b_S)")
	axislegend(axS, position = :rb)
	
	R_f = S_flux_interface ./ T_flux_interface
	axf = Axis(fig[3, 1], xlabel = "Rᵨ",ylabel = "R_f")
	lines!(axf, R_ρ_interp, R_f)

	a, b = [R_ρ_interp[start_flux:end].^0 R_ρ_interp[start_flux:end]] \ R_f[start_flux:end]
	lines!(axf, R_ρ_interp[start_flux:end], a .+ b .* R_ρ_interp[start_flux:end], label = "Linear fit, slope = $(b)")
	lines!(axf, R_ρ_interp[start_flux:end], lfit_S_flux ./ lfit_T_flux, linestyle = :dash, label = "Ratio of linear fits")
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
	lines!(ax, mins[1:end], expt_data["ΔS"][:] ./ expt_data["ΔS"][1], color = :blue, label = "ΔS / ΔS₀")
	lines!(ax, mins[1:end], expt_data["ΔT"][:] ./ expt_data["ΔT"][1], color = :red, label = "ΔT / ΔT₀")
	ylims!(ax, 0, 1.1)
	ax2 = Axis(fig[1, 1], yaxisposition = :right,
			    yticklabelcolor = :dodgerblue,
			    rightspinecolor = :dodgerblue,
				ylabelcolor = :dodgerblue,
			    ytickcolor = :dodgerblue,
				ylabel = "Rᵨ")
	lines!(ax2, expt_data["R_ρ"], color = :dodgerblue, label = "Rᵨ")
	axislegend(ax, merge = true, position = :rc)
	fig
end

# ╔═╡ d0148931-4198-4bb3-893c-a9b73e1ec7a9
begin
	fig_interface = Figure(size = (500, 500))
	ax_interface = Axis(fig_interface[1, 1], xlabel = "Rᵨ", ylabel = "z✶")
	lines!(ax_interface, expt_data["R_ρ"][2:end], dims["z✶"][expt_data["S_interface_idx"]], label = "salinity")
	lines!(ax_interface, expt_data["R_ρ"][2:end], dims["z✶"][expt_data["T_interface_idx"]], label = "temperature")

	id = expt_data["attrib/interface_depth"]
	S_mid = 0.5 * (expt_data["S_ha"][1, 1] + expt_data["S_ha"][end, 1])
	findS = [findfirst(reverse(expt_data["S_ha"][:, i]) .≥ S_mid) for i in eachindex(expt_data["S_ha"][1, 2:end])]
	S_interface = vcat(id, dims["z_aac"][findS])
	lines!(ax_interface, expt_data["R_ρ"], abs.(S_interface), label = "salinity (ha profile)", linestyle  = :dash)

	T_mid = 0.5 * (expt_data["T_ha"][1, 1] + expt_data["T_ha"][end, 1])
	findT = [findfirst(reverse(expt_data["T_ha"][:, i]) .≥ T_mid) for i in eachindex(expt_data["T_ha"][1, 2:end])]
	T_interface = vcat(id, dims["z_aac"][findT])
	lines!(ax_interface, expt_data["R_ρ"], abs.(T_interface), label = "temperature (ha profile)", linestyle = :dash)
	# ylims!(ax_interface, 0.49, 0.54)
	# vlines!(ax_interface, 1.6, color = :red, linestyle = :dash)
	axislegend(ax_interface, position = :rb)
	fig_interface
end

# ╔═╡ 82964820-7220-4e21-b263-20754b6a3a33
md"""
# Salinity-temperature evolution

**Note:** the below will only run with later experiments where the horizontal average salinity and temperature profiles are saved.

This is trying to capture asymmetry in density difference and compare it to theory (most of which can be found in the `diffusive_interfaces` Pluto notebook).

## Hovmollers
"""

# ╔═╡ 2536f46d-bbe8-4f85-a16a-82afef16fef5
let
	fig = Figure(size = (1000, 1000))
	t, z = dims["time"] / 60, dims["z_aac"]
	axS = Axis(fig[1, 1], title = "Hovmoller for saliniy and temperature (ha profile)", ylabel = "z (m)")
	hmS = heatmap!(axS, t, z, expt_data["S_ha"]', colormap = :haline)
	Colorbar(fig[1, 2], hmS)
	axT = Axis(fig[2, 1], ylabel = "z (m)")
	hmT = heatmap!(axT, t, z, expt_data["T_ha"]', colormap = :thermal)
	Colorbar(fig[2, 2], hmT)
	axσ = Axis(fig[3, 1], xlabel = "time (mins)", ylabel = "z (m)")
	hmσ = heatmap!(axσ, t, z, expt_data["σ_ha"]', colormap = :dense)
	Colorbar(fig[3, 2], hmσ)
	fig
end

# ╔═╡ 3130c878-fda2-4d11-a658-748d6a15b2b8
let
	fig = Figure(size = (1000, 1000))
	t, z = dims["time"] / 60, dims["z_aac"]
	axS = Axis(fig[1, 1], title = "Hovmoller for anomaly saliniy and temperature (ha profile)", ylabel = "z (m)")
	hmS = heatmap!(axS, t, z, (expt_data["S_ha"] .- expt_data["S_ha"][:, 1])', colormap = :curl, colorrange = (-0.09, 0.09))
	Colorbar(fig[1, 2], hmS, label = "S′ (gkg⁻¹)")
	axT = Axis(fig[2, 1], ylabel = "z (m)")
	hmT = heatmap!(axT, t, z, (expt_data["T_ha"] .- expt_data["T_ha"][:, 1])', colormap = :delta, colorrange = (-1.5, 1.5))
	Colorbar(fig[2, 2], hmT, label = "Θ′ (°C)")
	axσ = Axis(fig[3, 1], xlabel = "time (mins)", ylabel = "z (m)")
	hmσ = heatmap!(axσ, t, z, (expt_data["σ_ha"] .- expt_data["σ_ha"][:, 1])', colormap = :diff, colorrange = (-0.025, 0.025))
	Colorbar(fig[3, 2], hmσ, label = "σ′ kgm⁻³")
	fig
	# save("S_and_T_l_hov.png", fig)
end

# ╔═╡ ca6991b4-ac76-435e-bcff-82103b6abdc7
erf_tracer_solution(z, Cₗ::Number, ΔC::Number, κ::Number, time, interface_location) =
    Cₗ + 0.5 * ΔC * (1 + erf((z - interface_location) / sqrt(4 * κ * time)))

# ╔═╡ 616f02b0-5958-4e43-b747-4fb277460b55
md"""
## ``S-\Theta`` space evolution
"""

# ╔═╡ 6a8cda4f-8d4a-456e-b48c-7df73f933347
begin
	findS_neg = [findfirst(expt_data["S_ha"][:, i] .≤ S_mid) for i in eachindex(expt_data["S_ha"][1, 2:end])]
	S_interface_neg = vcat(id, dims["z_aac"][findS_neg])
	findT_neg = [findfirst(expt_data["T_ha"][:, i] .≤ T_mid) for i in eachindex(expt_data["T_ha"][1, 2:end])]
	T_interface_neg = vcat(id, dims["z_aac"][findT_neg])
	σ_mid = 0.5 * (expt_data["σ_ha"][1, 1] + expt_data["σ_ha"][end, 1])
	findσ_neg = [findfirst(expt_data["σ_ha"][:, i] .≤ σ_mid) for i in eachindex(expt_data["σ_ha"][1, 2:end])]
	σ_interface_neg = vcat(id, dims["z_aac"][findσ_neg])
	nothing
end

# ╔═╡ 18edf87f-fc2b-4e32-969b-eba0e2a813c1
t_slider = @bind t PlutoUI.Slider(eachindex(dims["time"]))

# ╔═╡ b4eded6e-eb4e-44b6-a5ea-1a16b09fc924
let
	Sₜ, Tₜ = expt_data["S_ha"][:, t], expt_data["T_ha"][:, t]
	z = dims["z_aac"]
	Rᵨ = expt_data["R_ρ"][t]

	fig = Figure(size = (800, 600))
	ax = Axis(fig[1, 1])
	lines!(ax, Sₜ, z; color = :blue, label = "Salinity")
	scatter!(ax, [Sₜ[findS_neg[t]]], [S_interface_neg[t]], color = :blue, label = "S interface")
	xlims!(ax, 34.56, 34.72)
	ax.xlabel = "Salinity (gkg⁻¹)"
	ax.xlabelcolor = :blue
	ax.xticklabelcolor = :blue
	ax.ylabel = "z (m)"
	hlines!(ax, -0.25, linestyle = :dash, color = :black, label = "Initial interface height")
	axT = Axis(fig[1, 1];
	            xaxisposition = :top,
	            xticklabelcolor = :red,
	            xlabel = "Θ (°C)",
	            xlabelcolor = :red,
	            title = "Temperature and salinity profiles")
	lines!(axT, Tₜ, z; color = :red, label = "Temperature")
	scatter!(axT, [Tₜ[findT_neg[t]]], [T_interface_neg[t]], color = :red, label = "T interface")
	axislegend(axT)
	axislegend(ax, position = :lb)
	xlims!(axT, -1.6, 0.6)
	linkyaxes!(ax, axT)
	Slims = extrema(expt_data["S_ha"][:, 1]) .+ [-0.01, 0.01]
	Tlims = extrema(expt_data["T_ha"][:, 1]) .+ [-0.1, 0.1]
	xlims!(ax, Slims...)
	xlims!(axT, Tlims...)
	ylims!(ax, -0.27, -0.23)

	σₜ = expt_data["σ_ha"][:, t]
	
	ax2 = Axis(fig[1, 2], 
			   title = "Ha σ profile at time t = $(dims["time"][t] / 60)min", 
			   subtitle = "Anomaly from mean",
			   xlabel = "σ₀", 
			   ylabel = "z (m)")
	lines!(ax2, σₜ, z)
	scatter!(ax2, [σₜ[findS_neg[t]]], [S_interface_neg[t]], color = :blue)
	scatter!(ax2, [σₜ[findT_neg[t]]], [T_interface_neg[t]], color = :red)
	scatter!(ax2, [σₜ[findσ_neg[t]]], [σ_interface_neg[t]], color = :orange)
	hlines!(ax2, -0.25, linestyle = :dash, color = :black, label = "Initial interface height")
	σ_lims = extrema(σₜ) .+ [-0.01, 0.01]
	xlims!(ax2, σ_lims)
	ylims!(ax, -0.27, -0.23)
	fig
end

# ╔═╡ 96ac45e7-9927-460e-907f-1449e09263f3
md"""
**Where do you cross midpoint density (in height)**
"""

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
	Δσ_upper = similar(R_Δσ)
	Δσ_lower = similar(R_Δσ)
	σ_min = similar(R_Δσ)
	σ_max = similar(R_Δσ)
	for (i, col) ∈ enumerate(eachcol(σ_ha))
		σ_min[i] = minimum(col)
		σ_top = col[end]
		Δσ_upper[i] = abs(σ_min[i] - σ_top)
		σ_max[i] = maximum(col)
		σ_bottom = col[1]
		Δσ_lower[i] = abs(σ_max[i] - σ_bottom)
		R_Δσ[i] = Δσ_upper[i] / Δσ_lower[i]
	end
end

# ╔═╡ 8711fdef-0da7-46bf-aa82-2f32b0590f7b
let
	Sₜ, Tₜ = expt_data["S_ha"][:, t], expt_data["T_ha"][:, t]
	Rᵨ = expt_data["R_ρ"][t]
	
	κₛ, κₜ = expt_data["attrib/κₛ (m²s⁻¹)"], expt_data["attrib/κₜ (m²s⁻¹)"]
	z = dims["z_aac"][:]
	ΔS = Sₜ[end] - Sₜ[1]
	ΔT = Tₜ[end] - Tₜ[1]
	S = erf_tracer_solution.(z, Sₜ[1], ΔS, κₛ, t, T_interface[t])
	T = erf_tracer_solution.(z, Tₜ[1], ΔT, κₜ, t, T_interface[t])

	# S_mixing = range(expt_data["S_ha"][end, 1], expt_data["S_ha"][1, 1], length = 1000)
	# ΔS_mix = expt_data["S_ha"][1, 1] - expt_data["S_ha"][end, 1]
	# ΔT_mix = expt_data["T_ha"][1, 1] - expt_data["T_ha"][end, 1]
	# T_mixing = expt_data["T_ha"][end, 1] .+ (ΔT_mix / ΔS_mix) * (S_mixing .- expt_data["S_ha"][end, 1])
	# ρ_mixing = total_density.(T_mixing, S_mixing, fill(0, length(S_mixing)), fill(eos_type, length(S_mixing)))
	# ρ_mix_max = maximum(ρ_mixing)
	# Δρ′ = ρ_mix_max - expt_data["σ_ha"][1, 1] 
	
	fig = Figure(size = (600, 500))
	ax = Axis(fig[1, 1], 
			  xlabel = "Salinity (gkg⁻¹)", 
			  ylabel = "Temperature (°C)",
			  title = "Ha S and T profiles at time t = $(dims["time"][t] / 60)min",
			  subtitle = "Rᵨ = $(round(Rᵨ, digits = 1)), R_Δσ = $(round(R_Δσ[t], digits = 1))")
	lines!(ax, Sₜ, Tₜ, label = "Model output")
	scatter!(ax, Sₜ[findS_neg[t]], Tₜ[findS_neg[t]], color = :red, label = "ST interface")
	Slims = extrema(expt_data["S_ha"][:, 1]) .+ [-0.01, 0.01]
	Tlims = extrema(expt_data["T_ha"][:, 1]) .+ [-0.1, 0.1]
	xlims!(Slims...)
	ylims!(Tlims...)

	lines!(ax, S, T, color = :orange, label = "Theoretical model")


	N = 200
	S_range = range(extrema(Sₜ)..., length=N)
	T_range = range(extrema(Tₜ)..., length=N)
	eos_vec = fill(eos_type, (N, N))
	FT = eltype(Tₜ[1])
	α = α_sign(thermal_expansion(Tₜ[1], Sₜ[1], FT(0), eos_vec[1]))
	β = haline_contraction(Tₜ[1], Sₜ[1], FT(0), eos_vec[1])
	ρ_grid = total_density.(T_range' .* ones(length(T_range)), S_range .* ones(length(S_range))', fill(0, (N, N)), eos_vec)
	ρ_shallow = total_density(Tₜ[end], Sₜ[end], 0, eos_vec[1])
	ρ_deep = total_density(Tₜ[1], Sₜ[1], 0, eos_vec[1])
	ρ_max = maximum(expt_data["σ_ha"][:, t])
	ρ_min = minimum(expt_data["σ_ha"][:, t])
	T_tangent = Tₜ[1] .+ (β / α) * (S_range .- Sₜ[1])
	lines!(ax, S_range, T_tangent, color = :green, linestyle = :dot, label = "Tangent to density\nat deep water")
	contour!(ax, S_range, T_range, ρ_grid, levels = [ρ_shallow, ρ_deep, ρ_max, ρ_min], color = :grey, label = "Deep water isopycnal")

	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ 8f0d41f8-d287-4c0d-af00-5b53f92cff75
abs(σ_min[3] - σ_max[3])

# ╔═╡ d2e81b8b-4a1c-4330-8f2a-14a502390bcd
let
	fig, ax = lines(expt_data["R_ρ"], R_Δσ)
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

# ╔═╡ 34f7166e-99a5-45fa-8b9b-963758e0d9f1
begin
	κₛ, κₜ = expt_data["attrib/κₛ (m²s⁻¹)"], expt_data["attrib/κₜ (m²s⁻¹)"]
	Θ, Sₐ = expt_data["T_ha"], expt_data["S_ha"]
	
	α = thermal_expansion.(Θ, Sₐ, zeros(Float32, size(Θ)),  fill(eos_type, size(Θ)))
	α_interp = α_sign(0.5*(α[1:end-1, :] .+ α[2:end, :]))
	
	β = haline_contraction.(Θ, Sₐ, zeros(Float32, size(Θ)),  fill(eos_type, size(Θ)))
	β_interp = 0.5*(β[1:end-1, :] .+ β[2:end, :])
	
	dz = diff(dims["z_aac"])
	
	 ∂_zΘ = diff(Θ, dims = 1) ./ dz 
	∂_zSₐ = diff(Sₐ, dims = 1) ./ dz

	A = 0.05^2
	g = 9.81

	Θ_mol = A * g * κₜ * sum(α_interp .* ∂_zΘ .* dz, dims = 1)
	S_mol = A * g * κₛ * sum(β_interp .* ∂_zSₐ.* dz, dims = 1)
	diff_offset = vec(Θ_mol .+ S_mol)
	diff_offset_interp = 0.5*(diff_offset[1:end-1] .+ diff_offset[2:end])
end

# ╔═╡ c576c4cd-1101-46e2-b6fa-b574f0b13dfe
let
	Δt = diff(dims["time"])
	ek = expt_data["∫Eₖ"]
	dₜek = diff(ek) ./ Δt
	ε = 0.5 * (expt_data["∫ε"][1:end-1] .+ expt_data["∫ε"][2:end])
	∫wb = 0.5 * (expt_data["∫wb"][1:end-1] .+ expt_data["∫wb"][2:end])
	∫gρw = - 0.5 * (expt_data["∫gρw"][1:end-1] .+ expt_data["∫gρw"][2:end])
	RHS = ∫wb .- ε #.- diff_offset_interp/2
	# RHS = ∫wb .- ε
	# RHS = ∫gρw .- ε #.- diff_offset_interp
	RHS_per = abs.(RHS) ./ abs.(ε)
	fig, ax = lines(eachindex(Δt)[3:end], dₜek[3:end], label = "dₜek")
	lines!(ax, eachindex(Δt)[3:end], RHS[3:end], label = "∫wb - ε")
	ax.title = "Energy  budget"
	ax.xlabel = "time (minutes)"
	ax.ylabel = "Watts / ρ₀"
	axislegend(ax, position = :rb)

	abs_err = abs.(dₜek .- RHS)
	ax2 = Axis(fig[2, 1], title = "Absolute error, MAE = $(mean(abs_err))")
	lines!(ax2, abs_err)
	fig
	# mean(RHS_per)
end

# ╔═╡ f200b8e0-2b14-4270-963b-6bb1b154d550
let
	fig, ax = lines(log10.(abs.((expt_data["∫wb"][2:end]))), label = "wb")
	# lines!(ax,  log10.(abs.(-expt_data["∫gρw"][2:end])), label = "∫gρw (post processing)", linestyle = :dash)
	lines!(ax, log10.(expt_data["∫ε"][2:end]), label = "∫ε", linestyle = :dot)
	ax.title = "Buoyancy flux and TKE dissipation (log10)"
	axislegend(ax)
	fig
end

# ╔═╡ 2f79adc0-2b5d-455d-beb8-5907c38d6e8e
let
	Δt = diff(dims["time"])
	Ep, Eb = expt_data["Ep"], expt_data["Eb"]
	Ea = Ep .- Eb
	dₜEp = diff(Ep) ./ Δt
	dₜEb = diff(Eb) ./ Δt
	dₜEa = diff(Ea) ./ Δt
	∫gρw = - 0.5 * (expt_data["∫gρw"][1:end-1] .+ expt_data["∫gρw"][2:end])
	fig, ax = lines(dₜEa, label = "dₜEa")
	RHS_APE = -(dₜEb .+ ∫gρw)
	lines!(ax, RHS_APE, label = "dₜEb .+ ∫gρw")
	abs_err = mean(abs.(dₜEa .- RHS_APE))
	ax.title = "APE budget, MAE = $(round(abs_err, digits = 14))"
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

# ╔═╡ 9e6998c4-6cca-49d5-9fff-2c697296849b
md"""
# Further diagnostics

Connected to the asymmetry in density anomalies across the interface is the more turbulent activity in the lower layer.
Two diagnostics can help here:
- interface height; and
- potential energy budget.

## Interface migration

The interface height I am going to define as the initial median ``S`` and ``\Theta`` values.
I can then find this in both the sorted 3D fields and the horizontally averaged fields to track layer thickness.

## Potential energy budget

This will be computed in the same way as project two.
My thought here is to break up the ``E_{p}`` integral into either lower layer, interface and upper layer or just use the interface as defined by the interface migration 
```math
E_{p} = \int_{z < z_{*}}g\rho z \mathrm{d}V + \int_{z > z_{*}}g\rho z \mathrm{d}V.
```
This might be the best place to start.
Otherwise can use something like
```math
E_{p} = \int_{z < z(\rho_{max})}g\rho z \mathrm{d}V + \int_{z(\rho_{max}) < z < z(\rho_{min})}g\rho z \mathrm{d}V + \int_{z > z(\rho_{min})}g\rho z \mathrm{d}V.
```
which then sums to total PE.
For the BPE we would then just use the sorted 1D density field but integrate to the same limits as ``E_{p}`` above.

Ideally this would be done by returning a 1D profile that is the horizontally integrated PE but I do not think this will work for APE so just have to choose the heights to integrate to.
"""

# ╔═╡ 31bed7ce-49d4-4009-be36-efd6531c979d
let
	fig = Figure(size = (600, 600))
	ax = Axis(fig[1, 1], xlabel = "time (mins)", ylabel = "Potential energy")
	
	Ep₀, Eb₀ = expt_data["Ep"][1], expt_data["Eb"][1]
	Ep, Eb = expt_data["Ep"], expt_data["Eb"]
	# Ep, Eb = (expt_data["∫Ep"] .- Eb₀) ./ Eb₀, (expt_data["∫Eb"] .- Eb₀) / Eb₀
	Ea = Ep .- Eb

	Ep_comp = expt_data["Ep_lower"] .+ expt_data["Ep_upper"]
	Eb_comp = expt_data["Eb_lower"] .+ expt_data["Eb_upper"]
	# Ep_comp = (expt_data["∫Ep_lower"] .+ expt_data["∫Ep_upper"] .- Ep₀) ./ Ep₀
	# Eb_comp = (expt_data["∫Eb_lower"] .+ expt_data["∫Eb_upper"] .- Eb₀) ./ Eb₀
	
	lines!(ax, dims["time"] ./ 60, Ep, label = "Ep")
	lines!(ax, dims["time"] ./ 60, Eb, label = "Eb")
	# lines!(ax, dims["time"] ./ 60, Ea, label = "Ea")
	lines!(ax, dims["time"] ./ 60, Ep_comp, label = "Ep upper + lower", linestyle = :dash)
	lines!(ax, dims["time"] ./ 60, Eb_comp, label = "Eb upper + lower", linestyle = :dash)
	axislegend(ax)
	Δt = diff(dims["time"])
	dₜEp = diff(Ep) ./ Δt
	dₜEb = diff(Eb) ./ Δt
	dₜEa = diff(Ea) ./ Δt
	ax2 = Axis(fig[2, 1], xlabel = "time (mins)", ylabel = "Potential energy")
	lines!(ax2, dims["time"][2:end] ./ 60, dₜEp, label = "dₜEp")
	lines!(ax2, dims["time"][2:end] ./ 60, dₜEb, label = "dₜEb")
	lines!(ax2, dims["time"][2:end] ./ 60, dₜEa, label = "dₜEa")
	axislegend(ax2)
	fig
end

# ╔═╡ 72353d1c-855b-463d-9bdb-b33bafc426d2
let
	fig = Figure(size = (600, 500))
	ax = Axis(fig[1, 1], xlabel = "time (mins)", ylabel = "Potential energy")
	
	Ep₀_upper = expt_data["Ep_upper"][1]
	# Ep_upper, Eb_upper = (expt_data["∫Ep_upper"] .- Ep₀_upper) ./ Ep₀_upper, (expt_data["∫Eb_upper"] .- Ep₀_upper) / Ep₀_upper
	Ep_upper, Eb_upper = expt_data["Ep_upper"], expt_data["Eb_upper"]
	Ea_upper = Ep_upper .- Eb_upper
	
	lines!(ax, dims["time"] ./ 60, Ep_upper, label = "Ep_upper")
	lines!(ax, dims["time"] ./ 60, Eb_upper, label = "Eb_upper")
	# lines!(ax, dims["time"] ./ 60, Ea_upper, label = "Ea_upper")

	Ep₀_lower = expt_data["Ep_lower"][1]
	# Ep_lower, Eb_lower = (expt_data["∫Ep_lower"] .- Ep₀_lower) ./ Ep₀_upper, (expt_data["∫Eb_lower"] .- Ep₀_lower) / Ep₀_upper
	Ep_lower, Eb_lower = expt_data["Ep_lower"], expt_data["Eb_lower"]
	Ea_lower = Ep_lower .- Eb_lower
	
	# lines!(ax, dims["time"] ./ 60, Ep_lower, label = "Ep_lower", linestyle = :dash)
	# lines!(ax, dims["time"] ./ 60, Eb_lower, label = "Eb_lower", linestyle = :dash)
	# lines!(ax, dims["time"] ./ 60, Ea_lower, label = "Ea_lower", linestyle = :dash)
	
	axislegend(ax)
	fig
	# Δt = diff(dims["time"])
	# dₜEp = diff(Ep) ./ Δt
	# dₜEb = diff(Eb) ./ Δt
	# dₜEa = diff(Ea) ./ Δt
	# ax2 = Axis(fig[2, 1], xlabel = "time (mins)", ylabel = "Potential energy")
	# lines!(ax2, dims["time"][2:end] ./ 60, dₜEp, label = "dₜEp")
	# lines!(ax2, dims["time"][2:end] ./ 60, dₜEb, label = "dₜEb")
	# lines!(ax2, dims["time"][2:end] ./ 60, dₜEa, label = "dₜEa")
	# axislegend(ax2)
end

# ╔═╡ b5102aaf-5c42-4c91-8016-2e9bae2073d8
let
	fig = Figure(size = (600, 500))
	ax = Axis(fig[1, 1], xlabel = "time (mins)", ylabel = "Potential energy")
	
	Ep₀_upper = expt_data["Ep_upper"][1]
	# Ep_upper, Eb_upper = (expt_data["∫Ep_upper"] .- Ep₀_upper) ./ Ep₀_upper, (expt_data["∫Eb_upper"] .- Ep₀_upper) / Ep₀_upper
	Ep_upper, Eb_upper = expt_data["Ep_upper"], expt_data["Eb_upper"]
	Ea_upper = Ep_upper .- Eb_upper
	
	# lines!(ax, dims["time"] ./ 60, Ep_upper, label = "Ep_upper")
	# lines!(ax, dims["time"] ./ 60, Eb_upper, label = "Eb_upper")
	# lines!(ax, dims["time"] ./ 60, Ea_upper, label = "Ea_upper")

	Ep₀_lower = expt_data["Ep_lower"][1]
	# Ep_lower, Eb_lower = (expt_data["∫Ep_lower"] .- Ep₀_lower) ./ Ep₀_upper, (expt_data["∫Eb_lower"] .- Ep₀_lower) / Ep₀_upper
	Ep_lower, Eb_lower = expt_data["Ep_lower"], expt_data["Eb_lower"]
	Ea_lower = Ep_lower .- Eb_lower
	
	lines!(ax, dims["time"] ./ 60, Ep_lower, label = "Ep_lower", linestyle = :dash)
	lines!(ax, dims["time"] ./ 60, Eb_lower, label = "Eb_lower", linestyle = :dash)
	# lines!(ax, dims["time"] ./ 60, Ea_lower, label = "Ea_lower", linestyle = :dash)
	
	axislegend(ax)
	fig
	# Δt = diff(dims["time"])
	# dₜEp = diff(Ep) ./ Δt
	# dₜEb = diff(Eb) ./ Δt
	# dₜEa = diff(Ea) ./ Δt
	# ax2 = Axis(fig[2, 1], xlabel = "time (mins)", ylabel = "Potential energy")
	# lines!(ax2, dims["time"][2:end] ./ 60, dₜEp, label = "dₜEp")
	# lines!(ax2, dims["time"][2:end] ./ 60, dₜEb, label = "dₜEb")
	# lines!(ax2, dims["time"][2:end] ./ 60, dₜEa, label = "dₜEa")
	# axislegend(ax2)
end

# ╔═╡ e70a5e3c-99de-44e0-b302-857ae75de3bf
md"""
## Heat flux
"""

# ╔═╡ 185fce17-00a3-40a3-b02e-ee44a2b10a28
let
	fig, ax, hm = heatmap(expt_data["ha_wT"][:, 1:60]', colormap = :speed)
	Colorbar(fig[1, 2], hm)
	fig
end

# ╔═╡ d053f751-6b25-4516-b5ca-c18fa886df52
md"""
# McDougall 1981 experiments

Checking what happens to the density asymmetry in McDougall experiment.
"""

# ╔═╡ 3ab99fd4-0e04-4fb1-b6ee-b009db9c93ad
t_slider

# ╔═╡ 650e9373-c299-4518-b223-de83f969df13
let
	Sₚ = [1.01, 6.01]
	T = [15.5, 30.39]
	salinity = gsw_sa_from_sp.(Sₚ, 0, 149, -35)
	temperature = gsw_ct_from_t.(salinity, T, 0)

	κₛ, κₜ = 1e-9, 1e-7
	z = range(-0.5, 0, length = 2000)
	ΔS = salinity[1] - salinity[2]
	ΔT = temperature[1] - temperature[2]
	S = erf_tracer_solution.(z, salinity[2], ΔS, κₛ, t, -0.25)
	T = erf_tracer_solution.(z, temperature[2], ΔT, κₜ, t, -0.25)
	σ = total_density.(T, S, 0, fill(eos_type, length(z)))

	R_Δρ = abs(minimum(σ) - σ[end]) / abs(maximum(σ) - σ[1])

	Δρ = σ[1] - σ[end]
	slope = ΔT / ΔS
	S_mix = range(salinity..., step = 0.000001)
    T_mix = @. temperature[1] + slope * (S_mix - salinity[1])
	ρ_mix = total_density.(T_mix, S_mix, 0, fill(eos_type, length(S_mix)))
	ρ_max = maximum(ρ_mix)
	Δρ′ = ρ_max - σ[1]
	C = Δρ′ / Δρ
	δ = (4 + 8*C) * sqrt(C^2 + C)
	
	fig = Figure(size = (600, 500))
	ax = Axis(fig[1, 1], 
			  xlabel = "Salinity (gkg⁻¹)", 
			  ylabel = "Temperature (°C)",
			  title = "Ha S and T profiles at time t = $(dims["time"][t] / 60)min",
			  subtitle = "Δρ′ = $(round(Δρ′, digits = 2)), R_Δσ = $(round(R_Δρ, digits = 1))")
	# Slims = salinity .+ [-0.01, 0.01]
	# Tlims = temperature .+ [-0.1, 0.1]
	# xlims!(Slims...)
	# ylims!(Tlims...)

	lines!(ax, S, T, color = :orange, label = "Theoretical model")
	lines!(ax, S_mix, T_mix)

	N = 100
	S_range = range(extrema(S)..., length=N)
	T_range = range(extrema(T)..., length=N)
	eos_vec = fill(eos_type, (N, N))
	FT = eltype(T[1])
	α = α_sign(thermal_expansion(temperature[end], salinity[end], FT(0), eos_vec[1]))
	β = haline_contraction(temperature[end], salinity[end], FT(0), eos_vec[1])
	ρ_grid = total_density.(T_range' .* ones(length(T_range)), S_range .* ones(length(S_range))', fill(0, (N, N)), eos_vec)
	ρ_deep = total_density(temperature[end], salinity[end], 0, eos_vec[1])
	T_tangent = T[1] .+ (β / α) * (S_range .- S[1])
	lines!(ax, S_range, T_tangent, color = :green, linestyle = :dot, label = "Tangent to density at deep water")
	contour!(ax, S_range, T_range, ρ_grid, levels = [ρ_deep], color = :black, label = "Deep water isopycnal")
	
	Legend(fig[2, 1], ax)

	σₜ = σ .- mean(σ)
	ax2 = Axis(fig[1, 2], 
			   title = "Ha σ′ profile at time t = $(dims["time"][t] / 60)min", 
			   subtitle = "Anomaly from mean",
			   xlabel = "σ₀′", 
			   ylabel = "z (m)")
	hlines!(ax2, -0.25, linestyle = :dash, color = :red, label = "Initial interface height")
	lines!(ax2, σₜ, z)
	σ_lims = extrema(σₜ) .+ [-0.1, 0.1]
	xlims!(ax2, σ_lims)
	Legend(fig[2, 2], ax2)
	fig
end

# ╔═╡ 4f04677f-b9c5-4026-882e-750ce2a3f76a
md"""
# Other ``S - Θ`` initial conditions
"""

# ╔═╡ 5b237029-d8ce-4289-8a29-416ae9babf6e
let
	salinity = [34.663, 34.70]
	temperature = [0, 0.5]
	# salinity = [34.6, 34.70]
	# temperature = [-1.5, 0.5]
	

	κₛ, κₜ = 1e-8, 1e-7
	z = range(-0.5, 0, length = 2000)
	ΔS = salinity[1] - salinity[2]
	ΔT = temperature[1] - temperature[2]
	S = erf_tracer_solution.(z, salinity[2], ΔS, κₛ, t, -0.25)
	T = erf_tracer_solution.(z, temperature[2], ΔT, κₜ, t, -0.25)
	σ = total_density.(T, S, 0, fill(eos_type, length(z)))

	R_Δρ = abs(minimum(σ) - σ[end]) / abs(maximum(σ) - σ[1])
	Rᵨ = compute_R_ρ(salinity, temperature, -0.25, eos_type)

	Δρ = σ[1] - σ[end]
	slope = ΔT / ΔS
	S_mix = range(salinity..., step = 0.000001)
    T_mix = @. temperature[1] + slope * (S_mix - salinity[1])
	ρ_mix = total_density.(T_mix, S_mix, 0, fill(eos_type, length(S_mix)))
	ρ_max = maximum(ρ_mix)
	Δρ′ = ρ_max - σ[1]
	C = Δρ′ / Δρ
	δ = (4 + 8*C) * sqrt(C^2 + C)
	
	fig = Figure(size = (600, 500))
	ax = Axis(fig[1, 1], 
			  xlabel = "Salinity (gkg⁻¹)", 
			  ylabel = "Temperature (°C)",
			  title = "Ha S and T profiles at time t = $(dims["time"][t] / 60)min",
			  subtitle = "Rᵨ = $(Rᵨ)  Δρ′ = $(round(Δρ′, digits = 6)), R_Δσ = $(round(R_Δρ, digits = 1))")
	# Slims = salinity .+ [-0.01, 0.01]
	# Tlims = temperature .+ [-0.1, 0.1]
	# xlims!(Slims...)
	# ylims!(Tlims...)

	lines!(ax, S, T, color = :orange, label = "Theoretical model")


	N = 100
	S_range = range(extrema(S)..., length=N)
	T_range = range(extrema(T)..., length=N)
	eos_vec = fill(eos_type, (N, N))
	FT = eltype(T[1])
	α = α_sign(thermal_expansion(temperature[end], salinity[end], FT(0), eos_vec[1]))
	β = haline_contraction(temperature[end], salinity[end], FT(0), eos_vec[1])
	ρ_grid = total_density.(T_range' .* ones(length(T_range)), S_range .* ones(length(S_range))', fill(0, (N, N)), eos_vec)
	ρ_deep = total_density(temperature[end], salinity[end], 0, eos_vec[1])
	T_tangent = T[1] .+ (β / α) * (S_range .- S[1])
	lines!(ax, S_range, T_tangent, color = :green, linestyle = :dot, label = "Tangent to density at deep water")
	contour!(ax, S_range, T_range, ρ_grid, levels = [ρ_deep], color = :black, label = "Deep water isopycnal")
	
	Legend(fig[2, 1], ax)

	σₜ = σ .- mean(σ)
	ax2 = Axis(fig[1, 2], 
			   title = "Ha σ′ profile at time t = $(dims["time"][t] / 60)min", 
			   subtitle = "Anomaly from mean",
			   xlabel = "σ₀′", 
			   ylabel = "z (m)")
	hlines!(ax2, -0.25, linestyle = :dash, color = :red, label = "Initial interface height")
	lines!(ax2, σₜ, z)
	σ_lims = extrema(σₜ) .+ [-0.1, 0.1]
	xlims!(ax2, σ_lims)
	Legend(fig[2, 2], ax2)
	fig
	# σ[1], σ[end]
end

# ╔═╡ 963fa274-2d8f-47fd-b227-4d7b3275d7ad
TableOfContents()

# ╔═╡ Cell order:
# ╟─6301138c-fa0b-11ef-0f3b-39dac35db063
# ╟─010ecdc3-51d6-41a6-9bc5-6efbba0723a6
# ╟─68d31cca-3f29-4402-ac79-8deaef98ef50
# ╟─087d2583-ee90-437a-97ec-0ab607337e30
# ╠═700b8f65-83af-4e13-96b4-7192cbff4e86
# ╠═b3fd9c9c-24d0-4ed8-a952-aaed2581045d
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
# ╟─82964820-7220-4e21-b263-20754b6a3a33
# ╟─2536f46d-bbe8-4f85-a16a-82afef16fef5
# ╟─3130c878-fda2-4d11-a658-748d6a15b2b8
# ╟─ca6991b4-ac76-435e-bcff-82103b6abdc7
# ╟─616f02b0-5958-4e43-b747-4fb277460b55
# ╟─6a8cda4f-8d4a-456e-b48c-7df73f933347
# ╟─b4eded6e-eb4e-44b6-a5ea-1a16b09fc924
# ╟─18edf87f-fc2b-4e32-969b-eba0e2a813c1
# ╟─8711fdef-0da7-46bf-aa82-2f32b0590f7b
# ╟─96ac45e7-9927-460e-907f-1449e09263f3
# ╟─6e7d38ab-2824-474b-b90c-75a3a5a05e57
# ╠═8f0d41f8-d287-4c0d-af00-5b53f92cff75
# ╟─e972f242-3b58-4581-87ae-437533b9fba1
# ╟─d2e81b8b-4a1c-4330-8f2a-14a502390bcd
# ╟─6ce43b6e-c3fa-408f-8702-900eaeb17bf5
# ╟─4538f159-01d9-45fd-9fa5-d7463c506a77
# ╟─d9422085-e838-44a1-91be-b81458dc3013
# ╟─3c0e1dfd-e4ba-448f-8475-ada056c8b5fe
# ╠═34f7166e-99a5-45fa-8b9b-963758e0d9f1
# ╠═c576c4cd-1101-46e2-b6fa-b574f0b13dfe
# ╟─f200b8e0-2b14-4270-963b-6bb1b154d550
# ╟─2f79adc0-2b5d-455d-beb8-5907c38d6e8e
# ╟─50e87efc-a49c-4ffd-bfbd-cd5dfad40639
# ╟─ee9c0edb-477b-4cc0-8c57-36845a90bbaf
# ╟─68a0a47e-e919-4d9d-b1a5-090d69bf633e
# ╟─9e6998c4-6cca-49d5-9fff-2c697296849b
# ╠═31bed7ce-49d4-4009-be36-efd6531c979d
# ╟─72353d1c-855b-463d-9bdb-b33bafc426d2
# ╟─b5102aaf-5c42-4c91-8016-2e9bae2073d8
# ╟─e70a5e3c-99de-44e0-b302-857ae75de3bf
# ╟─185fce17-00a3-40a3-b02e-ee44a2b10a28
# ╟─d053f751-6b25-4516-b5ca-c18fa886df52
# ╟─3ab99fd4-0e04-4fb1-b6ee-b009db9c93ad
# ╟─650e9373-c299-4518-b223-de83f969df13
# ╟─4f04677f-b9c5-4026-882e-750ce2a3f76a
# ╟─5b237029-d8ce-4289-8a29-416ae9babf6e
# ╟─963fa274-2d8f-47fd-b227-4d7b3275d7ad
