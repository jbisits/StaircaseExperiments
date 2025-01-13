using NCDatasets, JLD2

"""
    function φ_interface_flux!(flux_file::AbstractString, tracers::AbstractString, tracer::Symbol)
Calculate the flux through the diffusive interface for the `φ` tracer. The diffusive interface
is found in the sorted profile where φ > Δφ / 2 (with the previous two levels saved so it is
not a single value). The flux is then calculated as the change in φ content within the to
z✶(Δφ / 2) were z✶ is the backgrund z profile.
"""
function φ_interface_flux!(flux_file::AbstractString, tracers::AbstractString, tracer::Symbol)

    NCDataset(tracers) do ds

        φ = ds[tracer]
        Δφ₀ = φ[1, 1, 1, 1] - 0.5 * (φ[1, 1, 1, 1] - φ[1, 1, end, 1])
        timestamps = ds[:time][:]
        Δt = diff(timestamps)
        ΔV = diff(ds[:xC][1:2])[1] * diff(ds[:yC][1:2])[1] * diff(ds[:zC][1:2])[1]
        V = (1:length(reshape(φ[:, :, :, 1], :))) * ΔV
        SA = 0.07^2
        z✶ = V / SA
        Δz✶ = diff(z✶)[1]

        φ_interface_flux = Array{Float64}(undef, 3, length(Δt))
        interface_idx = Array{Int64}(undef, length(Δt))

        for i ∈ eachindex(Δt)

            φₜ = [reshape(φ[:, :, :, i], :) reshape(φ[:, :, :, i+1], :)]
            sort!(φₜ, dims = 1)
            ∫φdz✶ = cumsum(φₜ * Δz✶, dims = 1)
            dₜ∫φdz✶ = vec(diff(∫φdz✶, dims = 2) ./ Δt[i])

            φₜ_interp = 0.5 * vec(sum(φₜ, dims = 2))
            interface_idx[i] = ii = findfirst(φₜ_interp .> Δφ₀) - 1
            interface_idxs = [ii-1, ii, ii+1]

            φ_interface_flux[:, i] .= dₜ∫φdz✶[interface_idxs]

        end

        save_fluxes!(flux_file, φ_interface_flux, interface_idx, tracer)

    end

    return nothing
end
"""
    function save_fluxes!(flux_file, φ_interface_flux, interface_idx, tracer)
Save the flux through, and index of, the diffusive interface. **Note** the index is for the reshaped
and resorted vector.
"""
function save_fluxes!(flux_file, φ_interface_flux, interface_idx, tracer)

    if isfile(flux_file)
        jldopen(flux_file, "a+") do file
            file[string(tracer)*"_flux"] = φ_interface_flux
            file[string(tracer)*"_interface_idx"] = interface_idx
        end
    else
        jldopen(flux_file, "w") do file
            file[string(tracer)*"_flux"] = φ_interface_flux
            file[string(tracer)*"_interface_idx"] = interface_idx
        end
    end

    return nothing
end
