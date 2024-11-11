### A Pluto.jl notebook ###
# v0.20.3

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

# ╔═╡ 13d8331a-5396-11ef-24d4-6f8f7f57b6aa
begin
	using Pkg
	Pkg.activate("..")
	using StaircaseShenanigans, CairoMakie, NCDatasets, PlutoUI, Statistics, GibbsSeaWater
end

# ╔═╡ a902348d-5912-49b8-b784-4d065b0640ff
md"""
# Experiments run locally to test package functionality
"""

# ╔═╡ a84a9de8-6b20-481a-b0ef-2ca114aabf73
md"""
## Reentrant tracer and velocity
"""

# ╔═╡ 7377c6a0-b85a-4ce7-b401-3c31af4b5c4b
begin
	reentrant_output = "../../StaircaseShenanigans/examples/output/lineareos_single_interface_210min"
	tracers_reentrant = joinpath(reentrant_output, "tracers.nc")
	velocities_reentrant = joinpath(reentrant_output, "velocities.nc")
	co_reentrant = joinpath(reentrant_output, "computed_output.nc")
	ds_reentrant = NCDataset(tracers_reentrant)
	t_reentrant = ds_reentrant[:time][:]
	x_reentrant = ds_reentrant[:xC][:]
	z_reentrant = ds_reentrant[:zC][:]
	S_reentrant = ds_reentrant[:S][:, :, :, :]
	T_reentrant = ds_reentrant[:T][:, :, :, :]
	close(ds_reentrant)
	ds_reentrant = NCDataset(velocities_reentrant)
	zF_reentrant = ds_reentrant[:zF][:]
	w_reentrant = ds_reentrant[:w][:, :, :, :]
	close(ds_reentrant)
	ds_reentrant = NCDataset(co_reentrant)
	σ_reentrant = ds_reentrant[:σ][:, :, :, :]
	R_ρ_reentrant = ds_reentrant[:R_ρ][:]
	close(ds_reentrant)
end

# ╔═╡ 729bc65c-75e6-417e-b8fc-9e724d2acbdf
let
	fig, ax = lines(t_reentrant ./ 60, R_ρ_reentrant)
	ax.xlabel = "time (min)"
	ax.ylabel = "R_ρ"
	ax.title = "Density ratio between upper and lower quarter of domain"
	fig
end

# ╔═╡ c4559d79-f5bb-44a5-969e-36e379d62772
md"""
## Triple periodic with background state
"""

# ╔═╡ ac36c5a3-63bb-4a55-a356-c3732e15d071
begin
	periodic_output = "../../StaircaseShenanigans/examples/output/lineareos_single_interface_200min"
	tracers_periodic = joinpath(periodic_output, "tracers.nc")
	velocities_periodic = joinpath(periodic_output, "velocities.nc")
	co_periodic = joinpath(periodic_output, "computed_output.nc")
	ds_periodic = NCDataset(tracers_periodic)
	t_periodic = ds_periodic[:time][:]
	x_periodic = ds_periodic[:xC][:]
	z_periodic = ds_periodic[:zC][:]
	S_periodic = ds_periodic[:S][:, :, :, :]
	T_periodic = ds_periodic[:T][:, :, :, :]
	close(ds_periodic)
	ds_periodic = NCDataset(velocities_periodic)
	w_periodic = ds_periodic[:w][:, :, :, :]
	close(ds_periodic)
	ds_periodic = NCDataset(co_periodic)
	σ_periodic = ds_periodic[:σ][:, :, :, :]
	R_ρ_periodic = ds_periodic[:R_ρ][:]
	close(ds_periodic)
end

# ╔═╡ ba2793d6-f001-47ae-ac2f-f7e57353c48f
time_slider = @bind slider_periodic PlutoUI.Slider(eachindex(t_periodic))

# ╔═╡ 7d8a8521-6872-4c25-b5f8-77b5553b3a2b
let
	fig = Figure(size = (1000, 500))
	climits = (-1.6, 0.6)
	axbottom = Axis(fig[1, 1], xlabel = "x (m)", ylabel = "y (m)", title = "bottom temperature")
	hm = heatmap!(axbottom, T_reentrant[:, :, 1, slider_periodic], colormap = :thermal, colorrange = climits)
	axtop = Axis(fig[1, 2], xlabel = "x (m)", ylabel = "y (m)", title = "top temperature")
	heatmap!(axtop, T_reentrant[:, :, 50, slider_periodic], colormap = :thermal, colorrange = climits)
	Colorbar(fig[1, 3], label = "Θ (°C)", limits = climits, colormap = :thermal)
	fig
end

# ╔═╡ 29b4c096-23c7-47d1-8979-a23d912bec1e
let
	σ_z_mean = mean(σ_reentrant[:, :, :, slider_periodic], dims = (1, 2)) |> vec
	# T_z_mean = mean(T_periodic[:, :, :, slider_periodic], dims = (1, 2)) |> vec
	# S_z_mean = mean(S_periodic[:, :, :, slider_periodic], dims = (1, 2)) |> vec
	# σ_z_mean = gsw_sigma0.(S_z_mean, T_z_mean)
	fig, ax = lines(σ_z_mean, z_reentrant)
	ax.title = "t = $(t_reentrant[slider_periodic] / 60)min"
	ax.xlabel = "σ₀ (kgm⁻³)"
	ax.ylabel = "z (m)"
	fig
end

# ╔═╡ 39e8244d-c819-4717-aadc-3067748d1901
let
	w_z_mean = mean(w_reentrant[:, :, :, slider_periodic], dims = (1, 2)) |> vec
	w_plot = w_reentrant[2, 3, :, slider_periodic]
	fig, ax = lines(w_plot, zF_reentrant)
	ax.title = "t = $(t_reentrant[slider_periodic] / 60)min"
	ax.xlabel = "w (ms⁻¹)"
	ax.ylabel = "z (m)"
	# xlims!(-1.9, 0.7)
	# ax2 = Axis(fig[1, 2], xlabel = "S (gkg⁻¹)")
	# lines!(S_z_mean, z_periodic)
	fig
end

# ╔═╡ fd7bb1a0-dcd4-4a6f-89af-59c49702caa0
let
	T_z_mean = mean(T_reentrant[:, :, :, slider_periodic], dims = (1, 2)) |> vec
	S_z_mean = mean(S_reentrant[:, :, :, slider_periodic], dims = (1, 2)) |> vec
	fig, ax = lines(T_z_mean, z_reentrant)
	ax.title = "t = $(t_reentrant[slider_periodic] / 60)min"
	ax.xlabel = "Θ (°C)"
	ax.ylabel = "z (m)"
	# xlims!(-1.9, 0.7)
	ax2 = Axis(fig[1, 2], xlabel = "S (gkg⁻¹)")
	lines!(S_z_mean, z_periodic)
	fig
end

# ╔═╡ da8feb0b-901c-4a8b-ba4b-df926a05d458
time_slider

# ╔═╡ 8d998667-1dfb-4ea9-8273-b385ada56744
let
	T_z_mean = mean(T_periodic[:, :, :, slider_periodic], dims = (1, 2)) |> vec
	S_z_mean = mean(S_periodic[:, :, :, slider_periodic], dims = (1, 2)) |> vec
	fig, ax = lines(T_z_mean, z_periodic)
	ax.title = "t = $(t_periodic[slider_periodic] / 60)min"
	ax.xlabel = "Θ (°C)"
	ax.ylabel = "z (m)"
	# xlims!(-1.9, 0.7)
	ax2 = Axis(fig[1, 2], xlabel = "S (gkg⁻¹)")
	lines!(S_z_mean, z_periodic)
	fig
end

# ╔═╡ 10905668-b9c8-488f-bc17-d97bd588cf72
let
	fig = Figure(size = (1000, 500))
	climits = (-1.6, 0.6)
	axbottom = Axis(fig[1, 1], xlabel = "x (m)", ylabel = "y (m)", title = "bottom temperature")
	hm = heatmap!(axbottom, T_periodic[:, :, 1, slider_periodic], colormap = :thermal, colorrange = climits)
	axtop = Axis(fig[1, 2], xlabel = "x (m)", ylabel = "y (m)", title = "top temperature")
	heatmap!(axtop, T_periodic[:, :, 50, slider_periodic], colormap = :thermal, colorrange = climits)
	Colorbar(fig[1, 3], label = "Θ (°C)", limits = climits, colormap = :thermal)
	fig
end

# ╔═╡ 43f4f97b-ed20-47a3-905c-f04ddb9c9e8c
abs(T_periodic[1, 1, 50, end] - T_periodic[1, 1, 1, end])

# ╔═╡ e3de0e70-c40d-4ea6-b6a1-7027f92b9f0a
T_periodic[1, 1, 1, end], T_periodic[1, 1, 50, end]

# ╔═╡ d3596125-8a10-4314-91b5-0d6f35a66ba6
let
	σ_z_mean = mean(σ_periodic[:, :, :, slider_periodic], dims = (1, 2)) |> vec
	# T_z_mean = mean(T_periodic[:, :, :, slider_periodic], dims = (1, 2)) |> vec
	# S_z_mean = mean(S_periodic[:, :, :, slider_periodic], dims = (1, 2)) |> vec
	# σ_z_mean = gsw_sigma0.(S_z_mean, T_z_mean)
	fig, ax = lines(σ_z_mean, z_periodic)
	ax.title = "t = $(t_periodic[slider_periodic] / 60)min"
	ax.xlabel = "σ₀ (kgm⁻³)"
	ax.ylabel = "z (m)"
	fig
end

# ╔═╡ fda231f0-d672-48a1-8bef-c8822883f2a2
let
	fig, ax = lines(t_periodic ./ 60, R_ρ_periodic)
	ax.xlabel = "time (min)"
	ax.ylabel = "R_ρ"
	ax.title = "Density ratio between upper and lower quarter of domain"
	fig
end

# ╔═╡ 2ab36986-263d-4863-8312-25b307a4b1e9
time_slider

# ╔═╡ 4c86065c-6c91-4958-89ed-ffea606c297d
let
	w_z_mean = mean(w_periodic[:, :, :, slider_periodic], dims = (1, 2)) |> vec
	w_plot = w_periodic[2, 3, :, slider_periodic]
	fig, ax = lines(w_plot, z_periodic)
	ax.title = "t = $(t_periodic[slider_periodic] / 60)min"
	ax.xlabel = "w (ms⁻¹)"
	ax.ylabel = "z (m)"
	# xlims!(-1.9, 0.7)
	# ax2 = Axis(fig[1, 2], xlabel = "S (gkg⁻¹)")
	# lines!(S_z_mean, z_periodic)
	fig
end

# ╔═╡ 28e796ff-7213-4252-939c-a1002a1b2e73
sum(w_periodic[:, :, 1, slider_periodic]), sum(w_periodic[:, :, 50, slider_periodic])

# ╔═╡ 2155191a-249f-4987-9ca5-80d92579108a
md"""
### Check tanh function
"""

# ╔═╡ efc56079-3d72-4eae-a6db-1ea74c74639b
let
	 sol = @. 0.5 + 0.5 * -2 * (1  + tanh(100(z_periodic - (-0.5))))
	lines(sol, z_periodic)
end

# ╔═╡ 2ce3fdf4-9788-455b-afe5-da1e13a70690
let
	sol = @. 34.58 - (34.7-34.58) * z_periodic / 1
	lines(sol, z_periodic)
end

# ╔═╡ Cell order:
# ╟─13d8331a-5396-11ef-24d4-6f8f7f57b6aa
# ╟─a902348d-5912-49b8-b784-4d065b0640ff
# ╟─a84a9de8-6b20-481a-b0ef-2ca114aabf73
# ╟─7377c6a0-b85a-4ce7-b401-3c31af4b5c4b
# ╟─ba2793d6-f001-47ae-ac2f-f7e57353c48f
# ╟─fd7bb1a0-dcd4-4a6f-89af-59c49702caa0
# ╟─7d8a8521-6872-4c25-b5f8-77b5553b3a2b
# ╟─29b4c096-23c7-47d1-8979-a23d912bec1e
# ╟─729bc65c-75e6-417e-b8fc-9e724d2acbdf
# ╟─39e8244d-c819-4717-aadc-3067748d1901
# ╟─c4559d79-f5bb-44a5-969e-36e379d62772
# ╟─ac36c5a3-63bb-4a55-a356-c3732e15d071
# ╟─da8feb0b-901c-4a8b-ba4b-df926a05d458
# ╟─8d998667-1dfb-4ea9-8273-b385ada56744
# ╟─10905668-b9c8-488f-bc17-d97bd588cf72
# ╟─43f4f97b-ed20-47a3-905c-f04ddb9c9e8c
# ╟─e3de0e70-c40d-4ea6-b6a1-7027f92b9f0a
# ╟─d3596125-8a10-4314-91b5-0d6f35a66ba6
# ╟─fda231f0-d672-48a1-8bef-c8822883f2a2
# ╟─2ab36986-263d-4863-8312-25b307a4b1e9
# ╟─4c86065c-6c91-4958-89ed-ffea606c297d
# ╠═28e796ff-7213-4252-939c-a1002a1b2e73
# ╟─2155191a-249f-4987-9ca5-80d92579108a
# ╟─efc56079-3d72-4eae-a6db-1ea74c74639b
# ╟─2ce3fdf4-9788-455b-afe5-da1e13a70690
