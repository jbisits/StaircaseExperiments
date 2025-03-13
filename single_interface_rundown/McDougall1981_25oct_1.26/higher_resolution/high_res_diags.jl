### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 7f0968c7-21f6-405c-b47d-bb757d5da951
begin
	using Pkg
	Pkg.activate("../../..")
	using CairoMakie, JLD2, Statistics, PlutoUI
end

# ╔═╡ d57c16f2-ffae-11ef-1517-e7d5dbde6cd4
md"""
# Higher resolution simulation

Trying to understand what is happening with the salinity in the McDougall expts.
"""

# ╔═╡ 8ebecf2f-ec22-46c1-9587-eced7e145b5e
begin
	output_file = "diagnostics_high_res.jld2"
	dims = load(output_file, "dims")
	diags = load(output_file, "high_res")
	@info "Output loaded as dims and diags"
end

# ╔═╡ ad224317-f40d-49e2-9cc6-58f24409acd8
md"""
## Upper vs lower layers
"""

# ╔═╡ f3fd0fbd-8eae-43bb-87ad-607ee2d63e91
let
	fig = Figure(size = (600, 1000))
	axT = Axis(fig[1, 1], xlabel = "Θₗ", ylabel = "Θᵤ")
	lines!(axT, diags["Tₗ_Tᵤ_ts"][:, 1], diags["Tₗ_Tᵤ_ts"][:, 2])
	scatter!(axT, diags["Tₗ_Tᵤ_ts"][1, 1], diags["Tₗ_Tᵤ_ts"][1, 2], label = "Initial values", color = :orange)
	axislegend(axT)
	axS = Axis(fig[2, 1], xlabel = "Sₗ", ylabel = "Sᵤ")
	lines!(axS, diags["Sₗ_Sᵤ_ts"][:, 1], diags["Sₗ_Sᵤ_ts"][:, 2])
	scatter!(axS, diags["Sₗ_Sᵤ_ts"][1, 1], diags["Sₗ_Sᵤ_ts"][1, 2], label = "Initial values", color = :orange)
	axρ = Axis(fig[3, 1], xlabel = "ρₗ", ylabel = "ρᵤ")
	lines!(axρ, diags["ρₗ_ρᵤ_ts"][:, 1], diags["ρₗ_ρᵤ_ts"][:, 2])
	scatter!(axρ, diags["ρₗ_ρᵤ_ts"][1, 1], diags["ρₗ_ρᵤ_ts"][1, 2], label = "Initial values", color = :orange)
	fig
end

# ╔═╡ c88d2894-0b29-44f0-9fd7-36d221c567a1
md"""
## Fluxes
"""

# ╔═╡ 57c905e0-0b52-49e3-8324-c42ba363ecb2
let
	R_ρ_interp = 0.5 * (diags["R_ρ"][1:end-1] .+ diags["R_ρ"][2:end])

	fig = Figure(size = (800, 800))
	axT = Axis(fig[1, 1], ylabel = "T flux")
	T_interface_idx = diags["ha_T_interface_idx"]
	T_flux_interface = [diags["ha_T_flux"][idx, i] for (i, idx) ∈ enumerate(T_interface_idx)]
	lines!(axT, R_ρ_interp, T_flux_interface)
	
	axS = Axis(fig[2, 1], ylabel = "S flux")
	S_interface_idx = diags["ha_S_interface_idx"]
	S_flux_interface = [diags["ha_S_flux"][idx, i] for (i, idx) ∈ enumerate(S_interface_idx)]
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

# ╔═╡ 35746847-3b7b-4e2f-89b3-5466ca4c507c
md"""
## Energetics

Currently these are all computed with [Oceanostics.jl](https://github.com/tomchor/Oceanostics.jl).
At this resolution not closing the energy budget.
"""

# ╔═╡ fbe7c939-f1d7-4eb0-a9a0-ad2af38bd18b
let
	t = dims["time"]
	Δt = diff(t)
	Ba, η = diags["Ba"], diags["η"]
	∫wb, ∫ε, ∫Eₖ = diags["∫wb"], diags["∫ε"], diags["∫Eₖ"]
	dₜ∫Eₖ = diff(∫Eₖ) ./ diff(t)
	∫wb_i, ∫ε_i = 0.5 * (∫wb[1:end-1] .+ ∫wb[2:end]), 0.5 * (∫ε[1:end-1] .+ ∫ε[2:end])
	fig = Figure(size = (800, 800))
	ax1 = Axis(fig[1, 1], ylabel = "Watts/ρ₀")
	lines!(ax1, t[1:end-1], dₜ∫Eₖ, label = "dₜ∫Eₖ")
	lines!(ax1, t[1:end-1], ∫wb_i - ∫ε_i, label = "∫wb_i - ∫ε_i")
	axislegend(ax1)
	hidexdecorations!(ax1, grid = false, ticks = false)
	ax2 = Axis(fig[2, 1], xlabel = "time (s)", ylabel = "length (m)")
	lines!(ax2, t[2:end], log10.(η[2:end]), label = "η")
	lines!(ax2, t[2:end], log10.(2.5*Ba[2:end]), label = "Ba")
	axislegend(ax2)
	fig
end

# ╔═╡ 17788df9-165b-4709-ac3c-e50827ca8f57
let
	Δx, Δz = 0.1 / 100, 0.5 / 1000
	Nx_required, Nz_required = 4000, 20000
	Δx_required, Δz_required = 0.1 / Nx_required, 0.5 / Nz_required
md"""
Batchelor scale ``≈ \mathcal{O}(10^{-5})`` so require ``\Delta_{\mathrm{spacing}} \approx 1\times 10^{-5}`` to resolve this length scale.
Simulations have claimed ``\Delta_{\mathrm{spacing}} < 2.5L_{\mathrm{Ba}}`` is sufficient (presumably closed energy budget).
So looking at ``\Delta_{\mathrm{spacing}} < 2.5\times 10^{-5}``.
This also requires a rethink of the maximum timestep as found in project 2.

The resolution for this simulation was ``\Delta x = \Delta y = `` $(Δx) and ``\Delta z = `` $(Δz) and it took around 10 hours for this hour to be simulated.
So two orders of magnitude to small horizontally and one order and a bit vertically.
Increasing ``Nx`` to $(Nx_required) and ``Nz`` to $(Nz_required) gives ``\Delta x = \Delta y = `` $(Δx_required) and ``\Delta z = `` $(Δz_required).
But this resolution is insane --- I am not sure if it will be able to run!

This is where choice of initial noise might make a difference or it necessitates starting with a larger ``R_{\rho}`` to avoid the intense turbulence.
But latter option is not possible because it is only at lower ``R_\rho`` we get non-linear effects.

Likely I will have to run one very high resolution simulation to at least confirm things.
Best case altering the random noise produces something more reasonable.

In conjunction with this we can reduce domain size then repeat one experiment with greater extent to show no differences.
"""
end

# ╔═╡ Cell order:
# ╟─d57c16f2-ffae-11ef-1517-e7d5dbde6cd4
# ╟─7f0968c7-21f6-405c-b47d-bb757d5da951
# ╟─8ebecf2f-ec22-46c1-9587-eced7e145b5e
# ╟─ad224317-f40d-49e2-9cc6-58f24409acd8
# ╟─f3fd0fbd-8eae-43bb-87ad-607ee2d63e91
# ╟─c88d2894-0b29-44f0-9fd7-36d221c567a1
# ╟─57c905e0-0b52-49e3-8324-c42ba363ecb2
# ╟─35746847-3b7b-4e2f-89b3-5466ca4c507c
# ╟─fbe7c939-f1d7-4eb0-a9a0-ad2af38bd18b
# ╟─17788df9-165b-4709-ac3c-e50827ca8f57
