### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 13d8331a-5396-11ef-24d4-6f8f7f57b6aa
begin
	using Pkg
	Pkg.activate("..")
	using StaircaseShenanigans, CairoMakie, NCDatasets
end

# ╔═╡ a902348d-5912-49b8-b784-4d065b0640ff
md"""
# Experiments run locally to test package functionality

## Restoring

Aiming to be able to restore the uppermost and lowermost stairs in the staircase.
"""

# ╔═╡ dcad7506-c25d-4db7-ad87-e7e2ca5d50c6
begin
	expt_output = "../../StaircaseShenanigans/examples/output/nonlineareos_step_change_2min"
	tracers = joinpath(expt_output, "tracers.nc")
	ds = NCDataset(tracers)
	t = ds[:time][:]
	x = ds[:xC][:]
	z = ds[:zC][:]
	S = ds[:S][:, :, :, :]
	close(ds)
end

# ╔═╡ 6e072f82-d550-46de-a2a3-eced3fc73528
let
	fig = Figure(size = (1200, 1200))
	ax = [Axis(fig[i, j], xlabel = "Salinity", ylabel = "z (m)") for i ∈ 1:2, j ∈ 1:2]
	snapshots = [10, 11, 12, 13]
	for (i, t) ∈ enumerate(snapshots)
		lines!(ax[i], S[2, 2, :, t], z)
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

# ╔═╡ 7b8db315-a883-49bb-b2d2-199226dc21ff
S[2, 2, 1:10, :]

# ╔═╡ Cell order:
# ╟─13d8331a-5396-11ef-24d4-6f8f7f57b6aa
# ╟─a902348d-5912-49b8-b784-4d065b0640ff
# ╠═dcad7506-c25d-4db7-ad87-e7e2ca5d50c6
# ╠═6e072f82-d550-46de-a2a3-eced3fc73528
# ╠═7b8db315-a883-49bb-b2d2-199226dc21ff
