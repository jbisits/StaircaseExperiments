### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 13d8331a-5396-11ef-24d4-6f8f7f57b6aa
begin
	using Pkg
	Pkg.activate("..")
	using StaircaseShenanigans, CairoMakie, NCDatasets, PlutoUI, Statistics
end

# ╔═╡ a902348d-5912-49b8-b784-4d065b0640ff
md"""
# Experiments run locally to test package functionality

## Restoring

Aiming to be able to restore the uppermost and lowermost stairs in the staircase.
"""

# ╔═╡ dcad7506-c25d-4db7-ad87-e7e2ca5d50c6
begin
	expt_output = "../../StaircaseShenanigans/examples/output/nonlineareos_step_change_50min"
	tracers = joinpath(expt_output, "tracers.nc")
	ds = NCDataset(tracers)
	t = ds[:time][:]
	x = ds[:xC][:]
	z = ds[:zC][:]
	S = ds[:S][:, :, :, :]
	T = ds[:T][:, :, :, :]
	close(ds)
end

# ╔═╡ 565b9dde-bf8a-4152-a3ef-12dd5e4fcf48
md"""
## No temperature restoring
"""

# ╔═╡ 6e072f82-d550-46de-a2a3-eced3fc73528
let
	fig = Figure(size = (1200, 1200))
	ax = [Axis(fig[i, j], xlabel = "Temperature", ylabel = "z (m)") for i ∈ 1:2, j ∈ 1:2]
	snapshots = [1, 2, 100, length(t)]
	for (i, t) ∈ enumerate(snapshots)
		lines!(ax[i], T[2, 2, :, t], z)
		ax[i].title = "t = $(t)"
	end
	linkxaxes!(ax[1], ax[2])
	linkxaxes!(ax[3], ax[4])
	hidexdecorations!(ax[1], grid = false)
	hidexdecorations!(ax[3], grid = false)
	hideydecorations!(ax[3], grid = false)
	hideydecorations!(ax[4], grid = false)
	fig
end

# ╔═╡ c4559d79-f5bb-44a5-969e-36e379d62772
md"""
## With temperature restoring
"""

# ╔═╡ ac36c5a3-63bb-4a55-a356-c3732e15d071
begin
	# resroting_output = "../../StaircaseShenanigans/examples/output/nonlineareos_step_change_50min"
	resroting_output = "../../StaircaseShenanigans/examples/output/lineareos_step_change_50min"
	tracers_restore = joinpath(resroting_output, "tracers.nc")
	ds_restore = NCDataset(tracers_restore)
	t_restore = ds_restore[:time][:]
	x_restore = ds_restore[:xC][:]
	z_restore = ds_restore[:zC][:]
	S_restore = ds_restore[:S][:, :, :, :]
	T_restore = ds_restore[:T][:, :, :, :]
	# σ_restore = ds_restore[:σ][:, :, :, :]
	close(ds_restore)
end

# ╔═╡ da8feb0b-901c-4a8b-ba4b-df926a05d458
@bind slider_restore PlutoUI.Slider(eachindex(t_restore))

# ╔═╡ 8d998667-1dfb-4ea9-8273-b385ada56744
let
	T_z_mean = mean(T_restore[:, :, :, slider_restore], dims = (1, 2)) |> vec
	T_plot = T_restore[2, 2, :, slider_restore]
	fig, ax = lines(T_plot, z_restore)
	ax.title = "t = $(t_restore[slider_restore] / 60)min"
	ax.xlabel = "Θ (°C)"
	ax.ylabel = "z (m)"
	# xlims!(-1.9, 0.7)
	hlines!(ax, -0.1, color = :grey, linestyle = :dash)
	hlines!(ax, -0.9, color = :grey, linestyle = :dash)
	hlines!(ax, -0.4, color = :grey, linestyle = :dash)
	hlines!(ax, -0.6, color = :grey, linestyle = :dash)
	fig
end

# ╔═╡ Cell order:
# ╟─13d8331a-5396-11ef-24d4-6f8f7f57b6aa
# ╟─a902348d-5912-49b8-b784-4d065b0640ff
# ╟─dcad7506-c25d-4db7-ad87-e7e2ca5d50c6
# ╟─565b9dde-bf8a-4152-a3ef-12dd5e4fcf48
# ╠═6e072f82-d550-46de-a2a3-eced3fc73528
# ╟─c4559d79-f5bb-44a5-969e-36e379d62772
# ╠═ac36c5a3-63bb-4a55-a356-c3732e15d071
# ╟─da8feb0b-901c-4a8b-ba4b-df926a05d458
# ╠═8d998667-1dfb-4ea9-8273-b385ada56744
