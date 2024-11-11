using NCDatasets, JLD2, StatsBase

"""
    function mean_interface_thickness(φ)
Calculate the `mean_interface_thickness` of a single interface system. The `mean_interface_thickness`
uses the method of [Carpenter et al. (2012)](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/simulations-of-a-doublediffusive-interface-in-the-diffusive-convection-regime/63D2ECE2AA41439E01A01F9A0D76F2E2)
where a line is fitted in the region median(φ)-Δφ/8 < φ < median(φ) + Δφ/8 at each horizontal
location. Then the `mean` is taken to give the interface thickness of the tracer `φ`.
"""
function mean_interface_thickness(φ, z)

    Lx, Ly, Lz = size(φ)
    quarter_domain = round(Int, Lz / 4)
    φᵤ, φₗ = mean(φ[:, :, 1:quarter_domain]), mean(φ[:, :, 3*quarter_domain:end])
    Δφ = abs(φᵤ - φₗ)
    h_φ = Array{Float64}(undef, Lx, Ly)

    for i ∈ 1:Lx

        for j ∈ 1:Ly

            φ_profile = vec(φ[i, j, :])
            find_interface = findall(median(φ_profile) - Δφ / 8 .<  φ_profile .< median(φ_profile) + Δφ / 8)
            linfit = [ones(length(find_interface)) z[find_interface]] \ φ_profile[find_interface]
            h_φ[i, j] = abs(Δφ / linfit[2])

        end

    end

    return mean(h_φ)
end
"""
    function φ_flux(φ)
Calculate the flux for `φ` by sorting the `reshape`d 1D profile, then cumulatively
summing the result. **Note:** only appropriate where total amount of tracer within single
interface model is conserved.
"""
function φ_flux(φ, Δz, Δt)

    sort!(φ, dims = 1)
    ∫φ = cumsum(φ .* Δz, dims = 1)
    F_φ = vec(diff(∫φ  .* Δz, dims = 2) ./ Δt)

    return F_φ
end
"""
    function ha_flux_effective_κ(φ)
Calculate the horizontally averaged flux and effective diffusivity.
"""
function ha_flux_effective_κ(φ, Δz, Δt)

    ha_φ = [reshape(mean(φ[:, :, :, 1], dims = (1, 2)), :) reshape(mean(φ[:, :, :, 2], dims = (1, 2)), :)]
    ha_φ_flux = φ_flux(ha_φ, Δz, Δt)
    ∂φ_∂z = diff(ha_φ, dims = 1) ./ Δz
    ha_κ = (0.5 * (ha_φ_flux[1:end-1] .+ ha_φ_flux[2:end])) ./ (0.5 * (∂φ_∂z[:, 1] .+ ∂φ_∂z[:, 2]))

    return ha_φ_flux, ha_κ
end
