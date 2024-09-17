include("single_interface_diagnostics.jl")

data_path = "/g/data/e14/jb2381/StaircaseExperiments/single_interface/linear_simulation/lineareos_step_change_660min/tracers.nc"
tracers = joinpath(data_path, "tracers.nc")
fluxes_file = joinpath(@__DIR__, "single_interface_fluxes_linear.jld2")

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
    hₛ[i] = mean_interface_thickness(S, z)
    if i >= 2
        S_ = ds[:S][:, :, :, i-1:i]
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

    T = ds[:T][:, :, :, i]
    hₜ[i] = mean_interface_thickness(T, z)
    if i >= 2
        T_ = ds[:T][:, :, :, i-1:i]
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
    file["hₛ"] = hₛ
    file["ha_Fₛ"] = ha_Fₛ
    file["ha_κₛ"] = ha_κₛ
end
