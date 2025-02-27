using GibbsSeaWater, NCDatasets

"""
    function compute_Ẽ(co::AbstractString)
Compute the entrainment parameter using equation (25) from [McDougall (1981)](https://www.sciencedirect.com/science/article/pii/007966118190001X).
This is defined as
```math
Ẽ = \\frac{ρₗHₗdhₗ/dt + ρᵤHᵤdhᵤ/dt}{ρₗHₗdhₗ/dt - ρᵤHᵤdhᵤ/dt}.
```
The interface is found as the midpoint temperature and salinity of the two layers and the
heights Hᵤ and Hₗ are calculated from this but they take initial values of the domain divided
in two as that is how the experiments are set.
"""
function compute_Ẽ(co::AbstractString)

    ds = NCDataset(co)

    timestamps = ds[:time][:]
    z = ds[:zF][:]

    Δt = diff(timestamps)
    Ẽ = similar(Δt)
    σ = ds[:σ]

    for t ∈ eachindex(Δt)

    end

    close(ds)

    return Ẽ
end
