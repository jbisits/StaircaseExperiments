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
# Experiments run locally to validate triply perioidic domain
"""

# ╔═╡ a84a9de8-6b20-481a-b0ef-2ca114aabf73
md"""
# Tanh salt and temperature

## No background
"""

# ╔═╡ 7377c6a0-b85a-4ce7-b401-3c31af4b5c4b
begin
	no_background = "../../StaircaseShenanigans/validation/periodic_tanh_interface/lineareos_single_interface_240min"
	tracers_no_background = joinpath(no_background, "tracers.nc")
	velocities_no_background = joinpath(no_background, "velocities.nc")
	co_no_background = joinpath(no_background, "computed_output.nc")
	ds_no_background = NCDataset(tracers_no_background)
	t_no_background = ds_no_background[:time][:]
	x_no_background = ds_no_background[:xC][:]
	z_no_background = ds_no_background[:zC][:]
	S_no_background = ds_no_background[:S][:, :, :, :]
	T_no_background = ds_no_background[:T][:, :, :, :]
	close(ds_no_background)
	ds_no_background = NCDataset(velocities_no_background)
	zF_no_background = ds_no_background[:zF][:]
	w_no_background = ds_no_background[:w][:, :, :, :]
	close(ds_no_background)
	ds_no_background = NCDataset(co_no_background)
	σ_no_background = ds_no_background[:σ][:, :, :, :]
	R_ρ_no_background = ds_no_background[:R_ρ][:]
	close(ds_no_background)
end

# ╔═╡ 729bc65c-75e6-417e-b8fc-9e724d2acbdf
let
	fig, ax = lines(t_no_background ./ 60, R_ρ_no_background)
	ax.xlabel = "time (min)"
	ax.ylabel = "R_ρ"
	ax.title = "Density ratio between upper and lower quarter of domain"
	fig
end

# ╔═╡ c4559d79-f5bb-44a5-969e-36e379d62772
md"""
## With background state

Still `tanh` interface smoothing
"""

# ╔═╡ f8cd18d6-8ab1-490d-9b60-97a383ed8980
@bind background_state Select(["tanh", "linear"])

# ╔═╡ ac36c5a3-63bb-4a55-a356-c3732e15d071
begin
	with_background = "../../StaircaseShenanigans/validation/periodic_tanh_interface_$(background_state)_background/lineareos_single_interface_240min"
	tracers_with_background = joinpath(with_background, "tracers.nc")
	velocities_with_background = joinpath(with_background, "velocities.nc")
	co_with_background = joinpath(with_background, "computed_output.nc")
	ds_with_background = NCDataset(tracers_with_background)
	t_with_background = ds_with_background[:time][:]
	x_with_background = ds_with_background[:xC][:]
	z_with_background = ds_with_background[:zC][:]
	S_with_background = ds_with_background[:S][:, :, :, :]
	T_with_background = ds_with_background[:T][:, :, :, :]
	S′_with_background = ds_with_background[:S′][:, :, :, :]
	T′_with_background = ds_with_background[:T′][:, :, :, :]
	Sb_with_background = ds_with_background[:S_background][:, :, :]
	Tb_with_background = ds_with_background[:T_background][:, :, :]
	close(ds_with_background)
	ds_with_background = NCDataset(velocities_with_background)
	w_with_background = ds_with_background[:w][:, :, :, :]
	close(ds_with_background)
	ds_with_background = NCDataset(co_with_background)
	σ_with_background = ds_with_background[:σ][:, :, :, :]
	σb_with_background = ds_with_background[:σ_background][:, :, :]
	R_ρ_with_background = ds_with_background[:R_ρ][:]
	close(ds_with_background)
end

# ╔═╡ ba2793d6-f001-47ae-ac2f-f7e57353c48f
time_slider = @bind slider_with_background PlutoUI.Slider(eachindex(t_with_background))

# ╔═╡ 7d8a8521-6872-4c25-b5f8-77b5553b3a2b
let
	fig = Figure(size = (1000, 500))
	climits = (-1.6, 0.6)
	axbottom = Axis(fig[1, 1], xlabel = "x (m)", ylabel = "y (m)", title = "bottom temperature")
	hm = heatmap!(axbottom, T_no_background[:, :, 1, slider_with_background], colormap = :thermal, colorrange = climits)
	axtop = Axis(fig[1, 2], xlabel = "x (m)", ylabel = "y (m)", title = "top temperature")
	heatmap!(axtop, T_no_background[:, :, 50, slider_with_background], colormap = :thermal, colorrange = climits)
	Colorbar(fig[1, 3], label = "Θ (°C)", limits = climits, colormap = :thermal)
	fig
end

# ╔═╡ 29b4c096-23c7-47d1-8979-a23d912bec1e
let
	σ_z_mean = mean(σ_no_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec
	# T_z_mean = mean(T_with_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec
	# S_z_mean = mean(S_with_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec
	# σ_z_mean = gsw_sigma0.(S_z_mean, T_z_mean)
	fig, ax = lines(σ_z_mean, z_no_background)
	ax.title = "t = $(t_no_background[slider_with_background] / 60)min"
	ax.xlabel = "σ₀ (kgm⁻³)"
	ax.ylabel = "z (m)"
	fig
end

# ╔═╡ 39e8244d-c819-4717-aadc-3067748d1901
let
	w_z_mean = mean(w_no_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec
	w_plot = w_no_background[2, 3, :, slider_with_background]
	fig, ax = lines(w_plot, zF_no_background)
	ax.title = "t = $(t_no_background[slider_with_background] / 60)min"
	ax.xlabel = "w (ms⁻¹)"
	ax.ylabel = "z (m)"
	# xlims!(-1.9, 0.7)
	# ax2 = Axis(fig[1, 2], xlabel = "S (gkg⁻¹)")
	# lines!(S_z_mean, z_with_background)
	fig
end

# ╔═╡ fd7bb1a0-dcd4-4a6f-89af-59c49702caa0
let
	T_z_mean = mean(T_no_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec
	S_z_mean = mean(S_no_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec
	fig, ax = lines(T_z_mean, z_no_background)
	ax.title = "t = $(t_no_background[slider_with_background] / 60)min"
	ax.xlabel = "Θ (°C)"
	ax.ylabel = "z (m)"
	# xlims!(-1.9, 0.7)
	ax2 = Axis(fig[1, 2], xlabel = "S (gkg⁻¹)")
	lines!(S_z_mean, z_with_background)
	fig
end

# ╔═╡ da8feb0b-901c-4a8b-ba4b-df926a05d458
time_slider

# ╔═╡ 8d998667-1dfb-4ea9-8273-b385ada56744
let
	T_z_mean = mean(T_with_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec
	S_z_mean = mean(S_with_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec
	fig, ax = lines(T_z_mean, z_with_background)
	ax.title = "t = $(t_with_background[slider_with_background] / 60)min"
	ax.xlabel = "Θ (°C)"
	ax.ylabel = "z (m)"
	# xlims!(-1.9, 0.7)
	ax2 = Axis(fig[1, 2], xlabel = "S (gkg⁻¹)")
	lines!(S_z_mean, z_with_background)
	fig
end

# ╔═╡ 10905668-b9c8-488f-bc17-d97bd588cf72
let
	fig = Figure(size = (1000, 500))
	climits = (-1.6, 0.6)
	axbottom = Axis(fig[1, 1], xlabel = "x (m)", ylabel = "y (m)", title = "bottom temperature")
	hm = heatmap!(axbottom, T_with_background[:, :, 1, slider_with_background], colormap = :thermal, colorrange = climits)
	axtop = Axis(fig[1, 2], xlabel = "x (m)", ylabel = "y (m)", title = "top temperature")
	heatmap!(axtop, T_with_background[:, :, 50, slider_with_background], colormap = :thermal, colorrange = climits)
	Colorbar(fig[1, 3], label = "Θ (°C)", limits = climits, colormap = :thermal)
	fig
end

# ╔═╡ 43f4f97b-ed20-47a3-905c-f04ddb9c9e8c
abs(T_with_background[1, 1, 50, end] - T_with_background[1, 1, 1, end])

# ╔═╡ e3de0e70-c40d-4ea6-b6a1-7027f92b9f0a
T_with_background[1, 1, 1, end], T_with_background[1, 1, 50, end]

# ╔═╡ d3596125-8a10-4314-91b5-0d6f35a66ba6
let
	σ_z_mean = mean(σ_with_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec
	# T_z_mean = mean(T_with_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec
	# S_z_mean = mean(S_with_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec
	# σ_z_mean = gsw_sigma0.(S_z_mean, T_z_mean)
	fig, ax = lines(σ_z_mean, z_with_background)
	ax.title = "t = $(t_with_background[slider_with_background] / 60)min"
	ax.xlabel = "σ₀ (kgm⁻³)"
	ax.ylabel = "z (m)"
	fig
end

# ╔═╡ fda231f0-d672-48a1-8bef-c8822883f2a2
let
	fig, ax = lines(t_with_background ./ 60, R_ρ_with_background)
	ax.xlabel = "time (min)"
	ax.ylabel = "R_ρ"
	ax.title = "Density ratio between upper and lower quarter of domain"
	fig
end

# ╔═╡ 2ab36986-263d-4863-8312-25b307a4b1e9
time_slider

# ╔═╡ 4c86065c-6c91-4958-89ed-ffea606c297d
let
	w_z_mean = mean(w_with_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec
	w_plot = w_with_background[2, 3, :, slider_with_background]
	fig, ax = lines(w_plot, z_with_background)
	ax.title = "t = $(t_with_background[slider_with_background] / 60)min"
	ax.xlabel = "w (ms⁻¹)"
	ax.ylabel = "z (m)"
	# xlims!(-1.9, 0.7)
	# ax2 = Axis(fig[1, 2], xlabel = "S (gkg⁻¹)")
	# lines!(S_z_mean, z_with_background)
	fig
end

# ╔═╡ 28e796ff-7213-4252-939c-a1002a1b2e73
sum(w_with_background[:, :, 1, slider_with_background]), sum(w_with_background[:, :, 50, slider_with_background])

# ╔═╡ df5b471a-450b-46cb-ae84-a9abd8ab5a20
md"""
## Checking saved anomlay against Full - background

All looks good
"""

# ╔═╡ 99661d5a-03c9-47c3-ac4b-e6bc4d3bd3dc
let
	T_z_mean = mean(T′_with_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec
	S_z_mean = mean(S′_with_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec
	fig, ax = lines(T_z_mean, z_with_background)
	ax.title = "Saved anomaly t = $(t_with_background[slider_with_background] / 60)min"
	ax.xlabel = "Θ′ (°C)"
	ax.ylabel = "z (m)"
	# xlims!(-1.9, 0.7)
	ax2 = Axis(fig[1, 2], xlabel = "S′ (gkg⁻¹)")
	lines!(S_z_mean, z_with_background)
	fig
end

# ╔═╡ 6ff01c76-a46e-47e9-ab14-c35cf5894b2b
let
	T_z_mean = mean(T_with_background[:, :, :, slider_with_background] .- T′_with_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec
	S_z_mean = mean(S_with_background[:, :, :, slider_with_background] .- S′_with_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec
	fig, ax = lines(T_z_mean, z_with_background)
	ax.title = "Full - anomaly i.e. background t = $(t_with_background[slider_with_background] / 60)min"
	ax.xlabel = "Θ (°C)"
	ax.ylabel = "z (m)"
	# xlims!(-1.9, 0.7)
	ax2 = Axis(fig[1, 2], xlabel = "S (gkg⁻¹)")
	lines!(S_z_mean, z_with_background)
	fig
end

# ╔═╡ a903ec1f-7424-4f61-99f8-f072ee167650
time_slider

# ╔═╡ a48ec43b-395d-47d0-8e06-dbcceb77bf0b
let
	S_params = background_state == "tanh" ? (ΔC = 0.14, Cₗ = 34.7, Lz = 1.0, z_interface = -0.5, D = 100) : (ΔC = 0.14, Cᵤ = 34.56, Lz = 1.0)
	S_background = background_state == "tanh" ? tanh_background.(z_with_background, S_params...) : linear_background.(z_with_background, S_params...)
	T_params =  background_state == "tanh" ? (ΔC = 2.0, Cₗ = 0.5, Lz = 1.0, z_interface = -0.5, D = 100) : (ΔC = 2.0, Cᵤ = -1.5, Lz = 1.0)
	T_background =  background_state == "tanh" ? tanh_background.(z_with_background, T_params...) : linear_background.(z_with_background, T_params...)

	T_z_mean = mean(T_with_background[:, :, :, slider_with_background] .- T′_with_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec
	S_z_mean = mean(S_with_background[:, :, :, slider_with_background] .- S′_with_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec

	Tb_z_mean = mean(Tb_with_background[:, :, :], dims = (1, 2)) |> vec
	Sb_z_mean = mean(Sb_with_background[:, :, :], dims = (1, 2)) |> vec
	fig, ax = lines(T_z_mean, z_with_background, label = "Full - anomaly")
	ax.title = "Full - anomaly t = $(t_with_background[slider_with_background] / 60)min"
	ax.xlabel = "Θ (°C)"
	ax.ylabel = "z (m)"
	# xlims!(-1.9, 0.7)
	lines!(ax, T_background, z_with_background, label = "Computed from function", linestyle = :dash)
	lines!(ax, Tb_z_mean, z_with_background, label = "Saved initial field", linestyle = :dot)
	ax2 = Axis(fig[1, 2], xlabel = "S′ (gkg⁻¹)")
	lines!(ax2, S_z_mean, z_with_background, label = "Full - anomaly")
	lines!(ax2, S_background, z_with_background, label = "Computed from function", linestyle = :dash)
	lines!(ax2, Sb_z_mean, z_with_background, label = "Saved initial field", linestyle = :dot)
	axislegend(ax)
	axislegend(ax2)
	fig
end

# ╔═╡ 136e28d9-9c60-4b8a-8db7-646eb01e1827
let
	S_params = background_state == "tanh" ? (ΔC = 0.14, Cₗ = 34.7, Lz = 1.0, z_interface = -0.5, D = 100) : (ΔC = 0.14, Cᵤ = 34.56, Lz = 1.0)
	S_background = background_state == "tanh" ? tanh_background.(z_with_background, S_params...) : linear_background.(z_with_background, S_params...)
	T_params =  background_state == "tanh" ? (ΔC = 2.0, Cₗ = 0.5, Lz = 1.0, z_interface = -0.5, D = 100) : (ΔC = 2.0, Cᵤ = -1.5, Lz = 1.0)
	T_background =  background_state == "tanh" ? tanh_background.(z_with_background, T_params...) : linear_background.(z_with_background, T_params...)

	T_z_mean = mean(T_with_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec
	T_z_mean .-= T_background
	S_z_mean = mean(S_with_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec
	S_z_mean .-= S_background

	T′_z_mean = mean(T′_with_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec
	S′_z_mean = mean(S′_with_background[:, :, :, slider_with_background], dims = (1, 2)) |> vec

	fig, ax = lines(T_z_mean, z_with_background, label = "Full - background")
	ax.title = "Full - anomaly t = $(t_with_background[slider_with_background] / 60)min"
	ax.xlabel = "Θ (°C)"
	ax.ylabel = "z (m)"
	# xlims!(-1.9, 0.7)
	lines!(ax, T′_z_mean, z_with_background, label = "Saved anomaly", linestyle = :dash)
	ax2 = Axis(fig[1, 2], xlabel = "S′ (gkg⁻¹)")
	lines!(ax2, S_z_mean, z_with_background, label = "Full - background")
	lines!(ax2, S′_z_mean, z_with_background, label = "Saved anomaly", linestyle = :dash)
	axislegend(ax)
	axislegend(ax2)
	fig
end

# ╔═╡ 2155191a-249f-4987-9ca5-80d92579108a
md"""
### Check tanh function
"""

# ╔═╡ efc56079-3d72-4eae-a6db-1ea74c74639b
let
	 sol = @. 0.5 + 0.5 * -2 * (1  + tanh(100(z_with_background - (-0.5))))
	lines(sol, z_with_background)
end

# ╔═╡ 2ce3fdf4-9788-455b-afe5-da1e13a70690
let
	sol = @. 34.58 - (34.7-34.58) * z_with_background / 1
	lines(sol, z_with_background)
end

# ╔═╡ 344cec98-8c7f-49ed-9b2c-920527794000
TableOfContents()

# ╔═╡ f46c4a18-e3c2-4313-b9c7-cebeff9baef3
md"""
# Check the animations
"""

# ╔═╡ fb32d6b1-a6dd-4dc7-89e7-98c4db31e1a8
# ╠═╡ disabled = true
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	animate_density_anomaly(co_with_background, "σ", xslice = 2, yslice = 2)
	animate_density(co_with_background, "σ", xslice = 2, yslice = 2)
	animate_tracers_anomaly(tracers_with_background, xslice = 2, yslice = 2)
	animate_tracers(tracers_with_background, xslice = 2, yslice = 2)
end
  ╠═╡ =#

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
# ╟─f8cd18d6-8ab1-490d-9b60-97a383ed8980
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
# ╟─28e796ff-7213-4252-939c-a1002a1b2e73
# ╟─df5b471a-450b-46cb-ae84-a9abd8ab5a20
# ╟─99661d5a-03c9-47c3-ac4b-e6bc4d3bd3dc
# ╟─6ff01c76-a46e-47e9-ab14-c35cf5894b2b
# ╟─a903ec1f-7424-4f61-99f8-f072ee167650
# ╟─a48ec43b-395d-47d0-8e06-dbcceb77bf0b
# ╟─136e28d9-9c60-4b8a-8db7-646eb01e1827
# ╟─2155191a-249f-4987-9ca5-80d92579108a
# ╟─efc56079-3d72-4eae-a6db-1ea74c74639b
# ╟─2ce3fdf4-9788-455b-afe5-da1e13a70690
# ╠═344cec98-8c7f-49ed-9b2c-920527794000
# ╟─f46c4a18-e3c2-4313-b9c7-cebeff9baef3
# ╠═fb32d6b1-a6dd-4dc7-89e7-98c4db31e1a8
