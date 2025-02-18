### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ b8342884-a419-46ff-a2cf-3c7caac2f902
begin
	using Pkg
	Pkg.activate("..")
	using GibbsSeaWater, CairoMakie
	using StaircaseShenanigans: CustomLinearEquationOfState
	using SeawaterPolynomials: total_density
	using PlutoUI
end

# ╔═╡ 189277c0-ecdf-11ef-0b11-ad7fc2913d2d
md"""
# Diffusive interfaces part two

This is another notebook exploring diffusive interfaces mainly from an initial state to see how things evolve.
The motivation is that the simulations with non-linear eos don't produce an instability --- they just have a large density perturbation that grows.
Clearly the densification is from cabbeling but I cannot figure out why the system stays so stable.
"""

# ╔═╡ ff501fef-8e92-45f6-b26a-a9398c86957d
md"""
# `tanh` background
"""

# ╔═╡ f4c36c3e-e8f8-4d3e-872f-9a2a1a6fc3b0
begin
	S✶, Θ✶ = 34.7, 0.5
	Sᵤ, Θᵤ = 34.59, -1.5
	ΔS, ΔΘ = (S✶ - Sᵤ), (Θ✶ - Θᵤ)
	Sₘ, Θₘ = 0.5 * (S✶ + Sᵤ), 0.5 * (Θ✶ + Θᵤ)
	Lz = -1
	interface = Lz / 2
	z = range(Lz, 0, length = 1000)
	scale = 100
	initial_S = @. Sᵤ + (ΔS / 2) * (1 + tanh(100 * (z - interface)) / abs(Lz))
	reverse!(initial_S)
	initial_Θ = @. Θᵤ + (ΔΘ / 2) * (1 + tanh((100 / 3) * (z - interface)) / abs(Lz))
	reverse!(initial_Θ)
	intial_ρ = gsw_sigma0.(initial_S, initial_Θ)
	max_ρ, max_idx = findmax(intial_ρ)

	S_range = range(Sᵤ, S✶, length = 1000)
	mixing_line = @. Θ✶ + (ΔΘ / ΔS) * (S_range - S✶)
	ρ_mixing = gsw_sigma0.(S_range, mixing_line)
	max_ρ_mix, mix_idx = findmax(ρ_mixing)

	linear_eos = CustomLinearEquationOfState(Θₘ, Sₘ)
	ρ_linear = total_density.(initial_Θ, initial_S, 0, fill(linear_eos, length(initial_S))) .- 1024.694#1024.6885
end

# ╔═╡ 8667678b-a40b-4dfa-879a-0314616474fb
let
	fig = Figure(size = (1000,1000))
	ax1 = Axis(fig[1, 1],
			    xticklabelcolor = :blue,
			    bottomspinecolor = :blue,
			    xtickcolor = :blue,
				xlabelcolor = :blue,
				xlabel = "S (gkg⁻¹)",
				ylabel = "z (m)")
	ax2 = Axis(fig[1, 1],
				xaxisposition = :top,
			    xticklabelcolor = :red,
			    topspinecolor = :red,
			    xtickcolor = :red,
				xlabelcolor = :red,
				xlabel = "Θ (°C)")
	lines!(ax1, initial_S, z, color = :blue)
	lines!(ax2, initial_Θ, z, color = :red)

	ax3 = Axis(fig[1, 2],
				xlabel = "σ₀ (kgm⁻³)")
	lines!(ax3, intial_ρ, z, label = "Nonlinear eos")
	lines!(ax3, ρ_linear, z, label = "Linear eos")
	axislegend(ax3)
	ax4 = Axis(fig[2, :], xlabel = "S (gkg⁻¹)", ylabel = "Θ (°C)")
	lines!(ax4, initial_S, initial_Θ)
	scatter!(ax4, initial_S[max_idx], initial_Θ[max_idx], color = :orange)
	lines!(ax4, S_range, mixing_line, color = :purple, linestyle = :dash)
	scatter!(ax4, S_range[mix_idx], mixing_line[mix_idx], color = :red)
	fig
end

# ╔═╡ 67e66c9a-309b-43f3-8956-92404d77a3ef
md"""
# Linear background

Setting an appropriate linear background gradient is something I am not quite sure I am doing correctly.
One thing I could try is to make it so that the change between levels is some ``R_{\rho}'``, so that,
```math
\frac{S_{z}}{\Theta_{z}} = R_{\rho}'
```
but the units do not add up.
Discretely my thought is that if there is a uniform ``\Delta S`` and ``\Delta \Theta`` for each ``\Delta z``, then we could set this to satisfy an ``R_{\rho}`` criteria?
This might then actually get an instability because at present I just create a linear change between the top and bottom salinity and temperature.
"""

# ╔═╡ 67fc284f-71d1-494a-bdcf-efff768f0683
begin
	S_linear_bg = @. Sᵤ - ΔS * z / abs(Lz)
	Θ_linear_bg = @. Θᵤ - ΔΘ * z / abs(Lz)
	α, β = gsw_alpha.(S_linear_bg, Θ_linear_bg, 0), gsw_beta.(S_linear_bg, Θ_linear_bg, 0)
	α_interp = 0.5 * (α[1:end-1] .+ α[2:end])
	β_interp = 0.5 * (β[1:end-1] .+ β[2:end])
	ρ_linear_bg = gsw_rho.(S_linear_bg, Θ_linear_bg, 0) .- 1000
	ρ_linear_bg_lineareos = total_density.(Θ_linear_bg, S_linear_bg, 0, fill(linear_eos, length(S_linear_bg))) .- 1024.6942
	Δz = diff(z)
	S_z = diff(S_linear_bg) .* Δz
	Θ_z = diff(Θ_linear_bg) .* Δz
	(β_interp .* S_z) ./ (Θ_z .* α_interp)
end

# ╔═╡ f126f7a6-5e86-4f4c-a36f-10b9cf5228b4
let
	fig = Figure(size = (1000,1000))
	ax1 = Axis(fig[1, 1],
			    xticklabelcolor = :blue,
			    bottomspinecolor = :blue,
			    xtickcolor = :blue,
				xlabelcolor = :blue,
				xlabel = "S (gkg⁻¹)",
				ylabel = "z (m)")
	ax2 = Axis(fig[1, 1],
				xaxisposition = :top,
			    xticklabelcolor = :red,
			    topspinecolor = :red,
			    xtickcolor = :red,
				xlabelcolor = :red,
				xlabel = "Θ (°C)")
	lines!(ax1, S_linear_bg, z, color = :blue)
	lines!(ax2, Θ_linear_bg, z, color = :red, linestyle = :dash)

	ax3 = Axis(fig[1, 2],
				xlabel = "σ₀ (kgm⁻³)")
	lines!(ax3, ρ_linear_bg, z, label = "Nonlinear eos")
	lines!(ax3, ρ_linear_bg_lineareos, z, label = "Linear eos")
	axislegend(ax3)
	ax4 = Axis(fig[2, :], xlabel = "S (gkg⁻¹)", ylabel = "Θ (°C)")
	lines!(ax4,  S_linear_bg, Θ_linear_bg)
	# scatter!(ax4, initial_S[max_idx], initial_Θ[max_idx], color = :orange)
	# lines!(ax4, S_range, mixing_line, color = :purple, linestyle = :dash)
	# scatter!(ax4, S_range[mix_idx], mixing_line[mix_idx], color = :red)
	fig
end

# ╔═╡ Cell order:
# ╟─189277c0-ecdf-11ef-0b11-ad7fc2913d2d
# ╟─b8342884-a419-46ff-a2cf-3c7caac2f902
# ╟─ff501fef-8e92-45f6-b26a-a9398c86957d
# ╟─f4c36c3e-e8f8-4d3e-872f-9a2a1a6fc3b0
# ╟─8667678b-a40b-4dfa-879a-0314616474fb
# ╟─67e66c9a-309b-43f3-8956-92404d77a3ef
# ╟─67fc284f-71d1-494a-bdcf-efff768f0683
# ╟─f126f7a6-5e86-4f4c-a36f-10b9cf5228b4
