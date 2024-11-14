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

# ╔═╡ 21b56e8c-d764-4a37-b770-93263e17a7e5
using Printf

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
	S′_periodic = ds_periodic[:S′][:, :, :, :]
	T′_periodic = ds_periodic[:T′][:, :, :, :]
	Sb_periodic = ds_periodic[:S_background][:, :, :]
	Tb_periodic = ds_periodic[:T_background][:, :, :]
	close(ds_periodic)
	ds_periodic = NCDataset(velocities_periodic)
	w_periodic = ds_periodic[:w][:, :, :, :]
	close(ds_periodic)
	ds_periodic = NCDataset(co_periodic)
	σ_periodic = ds_periodic[:σ][:, :, :, :]
	σb_periodic = ds_periodic[:σ_background][:, :, :]
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

# ╔═╡ df5b471a-450b-46cb-ae84-a9abd8ab5a20
md"""
## Checking saved anomlay against Full - background

All looks good
"""

# ╔═╡ 99661d5a-03c9-47c3-ac4b-e6bc4d3bd3dc
let
	T_z_mean = mean(T′_periodic[:, :, :, slider_periodic], dims = (1, 2)) |> vec
	S_z_mean = mean(S′_periodic[:, :, :, slider_periodic], dims = (1, 2)) |> vec
	fig, ax = lines(T_z_mean, z_periodic)
	ax.title = "Saved anomaly t = $(t_periodic[slider_periodic] / 60)min"
	ax.xlabel = "Θ′ (°C)"
	ax.ylabel = "z (m)"
	# xlims!(-1.9, 0.7)
	ax2 = Axis(fig[1, 2], xlabel = "S′ (gkg⁻¹)")
	lines!(S_z_mean, z_periodic)
	fig
end

# ╔═╡ 6ff01c76-a46e-47e9-ab14-c35cf5894b2b
let
	T_z_mean = mean(T_periodic[:, :, :, slider_periodic] .- T′_periodic[:, :, :, slider_periodic], dims = (1, 2)) |> vec
	S_z_mean = mean(S_periodic[:, :, :, slider_periodic] .- S′_periodic[:, :, :, slider_periodic], dims = (1, 2)) |> vec
	fig, ax = lines(T_z_mean, z_periodic)
	ax.title = "Full - anomaly i.e. background t = $(t_periodic[slider_periodic] / 60)min"
	ax.xlabel = "Θ (°C)"
	ax.ylabel = "z (m)"
	# xlims!(-1.9, 0.7)
	ax2 = Axis(fig[1, 2], xlabel = "S (gkg⁻¹)")
	lines!(S_z_mean, z_periodic)
	fig
end

# ╔═╡ a903ec1f-7424-4f61-99f8-f072ee167650
time_slider

# ╔═╡ a48ec43b-395d-47d0-8e06-dbcceb77bf0b
let
	S_params = (ΔC = 0.12000000000000455, Cᵤ = 34.58, Cₗ = 34.7, Lz = 1.0, z_interface = -0.5, D = 100)
	_tanh_background(z, ΔC, Cᵤ, Cₗ, Lz, z_interface, D) = Cₗ - 0.5 * ΔC * (1  + tanh(D * (z - z_interface) / Lz))
	S_background = _tanh_background.(z_periodic, S_params...)
	T_params = (ΔC = 2.0, Cᵤ = -1.5, Cₗ = 0.5, Lz = 1.0, z_interface = -0.5, D = 100)
	T_background = _tanh_background.(z_periodic, T_params...)
	
	T_z_mean = mean(T_periodic[:, :, :, slider_periodic] .- T′_periodic[:, :, :, slider_periodic], dims = (1, 2)) |> vec
	S_z_mean = mean(S_periodic[:, :, :, slider_periodic] .- S′_periodic[:, :, :, slider_periodic], dims = (1, 2)) |> vec

	Tb_z_mean = mean(Tb_periodic[:, :, :], dims = (1, 2)) |> vec
	Sb_z_mean = mean(Sb_periodic[:, :, :], dims = (1, 2)) |> vec	
	fig, ax = lines(T_z_mean, z_periodic, label = "Full - anomaly")
	ax.title = "Full - anomaly t = $(t_periodic[slider_periodic] / 60)min"
	ax.xlabel = "Θ (°C)"
	ax.ylabel = "z (m)"
	# xlims!(-1.9, 0.7)
	lines!(ax, T_background, z_periodic, label = "Computed from function", linestyle = :dash)
	lines!(ax, Tb_z_mean, z_periodic, label = "Saved initial field", linestyle = :dot)
	ax2 = Axis(fig[1, 2], xlabel = "S′ (gkg⁻¹)")
	lines!(ax2, S_z_mean, z_periodic, label = "Full - anomaly")
	lines!(ax2, S_background, z_periodic, label = "Computed from function", linestyle = :dash)
	lines!(ax2, Sb_z_mean, z_periodic, label = "Saved initial field", linestyle = :dot)
	axislegend(ax)
	axislegend(ax2)
	fig
end

# ╔═╡ 136e28d9-9c60-4b8a-8db7-646eb01e1827
let
	S_params = (ΔC = 0.12000000000000455, Cᵤ = 34.58, Cₗ = 34.7, Lz = 1.0, z_interface = -0.5, D = 100)
	_tanh_background(z, ΔC, Cᵤ, Cₗ, Lz, z_interface, D) = Cₗ - 0.5 * ΔC * (1  + tanh(D * (z - z_interface) / Lz))
	S_background = _tanh_background.(z_periodic, S_params...)
	T_params = (ΔC = 2.0, Cᵤ = -1.5, Cₗ = 0.5, Lz = 1.0, z_interface = -0.5, D = 100)
	T_background = _tanh_background.(z_periodic, T_params...)
	
	T_z_mean = mean(T_periodic[:, :, :, slider_periodic], dims = (1, 2)) |> vec 
	T_z_mean .-= T_background
	S_z_mean = mean(S_periodic[:, :, :, slider_periodic], dims = (1, 2)) |> vec 
	S_z_mean .-= S_background

	T′_z_mean = mean(T′_periodic[:, :, :, slider_periodic], dims = (1, 2)) |> vec
	S′_z_mean = mean(S′_periodic[:, :, :, slider_periodic], dims = (1, 2)) |> vec
	
	fig, ax = lines(T_z_mean, z_periodic, label = "Full - background")
	ax.title = "Full - anomaly t = $(t_periodic[slider_periodic] / 60)min"
	ax.xlabel = "Θ (°C)"
	ax.ylabel = "z (m)"
	# xlims!(-1.9, 0.7)
	lines!(ax, T′_z_mean, z_periodic, label = "Saved anomaly", linestyle = :dash)
	ax2 = Axis(fig[1, 2], xlabel = "S′ (gkg⁻¹)")
	lines!(ax2, S_z_mean, z_periodic, label = "Full - background")
	lines!(ax2, S′_z_mean, z_periodic, label = "Saved anomaly", linestyle = :dash)
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
	 sol = @. 0.5 + 0.5 * -2 * (1  + tanh(100(z_periodic - (-0.5))))
	lines(sol, z_periodic)
end

# ╔═╡ 2ce3fdf4-9788-455b-afe5-da1e13a70690
let
	sol = @. 34.58 - (34.7-34.58) * z_periodic / 1
	lines(sol, z_periodic)
end

# ╔═╡ f46c4a18-e3c2-4313-b9c7-cebeff9baef3
md"""
# Check the animations
"""

# ╔═╡ 2b05755e-a29c-4e0f-b873-6c73e7ea79c9
function animate_density_anomaly(computed_output::AbstractString, variable::AbstractString;
                               xslice = 3, yslice = 3)

    NCDataset(computed_output) do ds

        x = ds["xC"][:]
        z = ds["zC"][:]
        t = ds["time"][:]

        n = Observable(1)
        σ_backgroud = ds[variable*"_background"][:, yslice, :]
        σ = @lift ds[variable][:, yslice, :, $n] .- σ_backgroud
        σ_backgroud_profile =  ds[variable*"_background"][xslice, yslice, :]
        σ_profile = @lift ds[variable][xslice, yslice, :, $n] .- σ_backgroud_profile
        time_title = @lift @sprintf("t=%1.2f minutes", t[$n] / 60)

        fig = Figure(size = (1000, 500))
        ax = [Axis(fig[1, i], title = i == 1 ? time_title : "") for i ∈ 1:2]

        colorrange = extrema(ds[variable][:, :, :, end] .- ds[variable*"_background"][:, :, :])

        lines!(ax[1], σ_profile, z)
        ax[1].xlabel = "σ₀′ kgm⁻³"
        ax[1].ylabel = "z"
        ax[1].xaxisposition = :top
        ax[1].xticklabelrotation = π / 4
        xlims!(ax[1], colorrange)

        colormap = cgrad(:dense)[2:end-1]
        lowclip = cgrad(:dense)[1]
        highclip = cgrad(:dense)[end]
        hm = heatmap!(ax[2], x, z, σ; colorrange, colormap, lowclip, highclip)

        ax[2].xlabel = "x (m)"
        ax[2].ylabel = "z (m)"
        Colorbar(fig[1, 3], hm, label = "σ₀′ kgm⁻³")

        linkyaxes!(ax[1], ax[2])
        hideydecorations!(ax[2], ticks = false)

        frames = eachindex(t)
        record(fig, joinpath(pwd(), "density.mp4"),
            frames, framerate=8) do i
            msg = string("Plotting frame ", i, " of ", frames[end])
            print(msg * " \r")
            n[] = i
        end

    end

    return nothing
end

# ╔═╡ fb32d6b1-a6dd-4dc7-89e7-98c4db31e1a8
animate_density_anomaly(co_periodic, "σ", xslice = 2, yslice = 2)

# ╔═╡ 26be505b-2297-4215-81ad-0605513d4f67
function animate_tracers_anomaly(tracers::AbstractString; xslice = 3, yslice = 3)

    NCDataset(tracers) do ds

        x = ds["xC"][:]
        z = ds["zC"][:]
        t = ds["time"][:]

        n = Observable(1)
        S = @lift ds[:S′][:, yslice, :, $n]
        S_profile = @lift ds[:S′][xslice, yslice, :, $n]
        Θ = @lift ds[:T′][:, yslice, :, $n]
        Θ_profile = @lift ds[:T′][xslice, yslice, :, $n]
        time_title = @lift @sprintf("t=%1.2f minutes", t[$n] / 60)

        fig = Figure(size = (1000, 1000))
        ax = [Axis(fig[j, i], title = (i == 1 && j == 1) ? time_title : "") for i ∈ 1:2, j ∈ 1:2]

        # Salinity
        Srange = extrema(ds[:S′][:, :, :, end])
        lines!(ax[1], S_profile, z)
        ax[1].xlabel = "S′ gkg⁻¹"
        ax[1].ylabel = "z (m)"
        ax[1].xaxisposition = :top
        xlims!(ax[1], Srange)

        Scmap = cgrad(:haline)[2:end-1]
        Slow = cgrad(:haline)[1]
        Shigh = cgrad(:haline)[end]
        hmS = heatmap!(ax[2], x, z, S, colorrange = Srange, colormap = Scmap,
                        lowclip = Slow, highclip = Shigh)

        ax[2].xlabel = "x (m)"
        ax[2].ylabel = "z (m)"
        Colorbar(fig[1, 3], hmS, label = "S′ gkg⁻¹")

        linkyaxes!(ax[1], ax[2])
        hideydecorations!(ax[2], ticks = false)

        # Temperature
        Θrange = extrema(ds[:T′][:, :, :, end])
        lines!(ax[3], Θ_profile, z)
        ax[3].xlabel = "Θ°C"
        ax[3].ylabel = "z (m)"
        ax[3].xaxisposition = :top
        xlims!(ax[3], Θrange)

        Θcmap = cgrad(:thermal)[2:end-1]
        Θlow = cgrad(:thermal)[1]
        Θhigh = cgrad(:thermal)[end]
        hmΘ = heatmap!(ax[4], x, z, Θ, colorrange = Θrange, colormap = Θcmap,
                        lowclip = Θlow, highclip = Θhigh)

        ax[4].xlabel = "x (m)"
        ax[4].ylabel = "z (m)"
        Colorbar(fig[2, 3], hmΘ, label = "Θ°C")

        linkyaxes!(ax[3], ax[4])
        hideydecorations!(ax[4], ticks = false)

        frames = eachindex(t)
        record(fig, joinpath(pwd(), "tracers.mp4"),
            frames, framerate=8) do i
            msg = string("Plotting frame ", i, " of ", frames[end])
            print(msg * " \r")
            n[] = i
        end

    end

    return nothing
end

# ╔═╡ b5108b49-199d-49f2-b2bf-fb6761f823d6
animate_tracers_anomaly(tracers_periodic)

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
# ╟─f46c4a18-e3c2-4313-b9c7-cebeff9baef3
# ╠═21b56e8c-d764-4a37-b770-93263e17a7e5
# ╠═2b05755e-a29c-4e0f-b873-6c73e7ea79c9
# ╠═fb32d6b1-a6dd-4dc7-89e7-98c4db31e1a8
# ╠═26be505b-2297-4215-81ad-0605513d4f67
# ╠═b5108b49-199d-49f2-b2bf-fb6761f823d6
