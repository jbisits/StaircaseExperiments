## I already have some model output of a single interface with double diffusive effects.
## Here I use this output to try and calculate fluxes through the interface, along with
## interface thickness

using NCDatasets, JLD2, StatsBase

data_path = "/g/data/e14/jb2381/CabbelingExperiments/outputs_doublediffusion/cabbeling_stepchange_nothing_600min"
tracers = joinpath(data_path, "tracers.nc")
fluxes_file = joinpath(@__DIR__, "single_interface_fluxes.jld2")

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
            find_interface = (median(φ_profile)) - Δφ / 8 .<  φ_profile .< (median(φ_profile) + Δφ / 8)
            linfit = [ones(length(find_interface)) z[find_interface]] \ φ_profile[find_interface]
            h_φ[i, j] = abs(Δφ / linfit[2])

        end

    end

    return mean(h_φ)
end
"""
    function φ_flux(φ)
Calculate the flux for `φ` by reshaping and sorting into a 1D profile, then cumulatively
summing the result. **Note:** only appropriate where total amount of tracer within single
interface model is conserved.
"""
function φ_flux(φ)

    return nothing
end
"""
    function ha_flux_effective_κ(φ)
Calculate the horizontally averaged flux and effective diffusivity.
"""
function ha_flux_effective_κ(φ, Δz, Δt)

    ha_φ = mean(φ, dims = (1, 2))
    sort!(ha_φ, dims = 1)
    ∫ha_φ = cumsum(ha_φ .* Δz, dims = 1)
    ha_φ_flux = vec(diff(∫ha_φ  .* Δz, dims = 2) ./ Δt)
    ∂φ_∂z = vec(diff(ha_φ, dims = 1) ./ Δz)
    ha_κ = 0.5 * (ha_φ_flux[1:end-1] .+ ha_φ_flux[2:end]) ./ ∂φ_∂z[:, 2]

    return ha_φ_flux, ha_κ
end

ds = NCDataset(tracers)

t = ds[:time][:]
z = ds[:zC][:]
Δz = diff(z)[1]
Δt = diff(t)

## Salinity
hₛ = similar(t)
ha_Fₛ = Array{Float64}(undef, length(z), length(t)-1)
ha_κₛ = Array{Float64}(undef, length(z)-1, length(t)-1)
for i ∈ eachindex(t)

    S = ds[:S][:, :, :, i]
    hₛ = mean_interface_thickness(S, z)
    if i >= 2
        S_ = S[:, :, :, i-1:i]
        F, κ = ha_flux_effective_κ(S_, Δz, Δt[i-1])
        ha_Fₛ[:, i-1] = F
        ha_κₛ[:, i-1] = κ
    end
end

## Temperature
hₜ = similar(t)
ha_Fₜ = Array{Float64}(undef, length(z), length(t)-1)
ha_κₜ = Array{Float64}(undef, length(z)-1, length(t)-1)
for i ∈ eachindex(t)

    T = T[:, :, :, i]
    hₜ = mean_interface_thickness(T, z)
    if i >= 2
        T_ = T[:, :, :, i-1:i]
        F, κ = ha_flux_effective_κ(T_, Δz, Δt[i-1])
        ha_Fₜ[:, i-1] = F
        ha_κₜ[:, i-1] = κ
    end

end

close(ds)

jldopen(fluxes_file, "w") do file
    file["hₜ"] = hₜ
    file["ha_Fₜ"] = ha_Fₜ
    file["ha_κₜ"] = ha_κₜ
    file["hₛ"] = hₜ
    file["ha_Fₛ"] = ha_Fₛ
    file["ha_κₛ"] = ha_κₛ
end
