using NCDatasets, JLD2, Statistics

output_path = joinpath(@__DIR__, "rundown_1.05/step/nonlineareos_single_interface_60min")
tracers = joinpath(output_path, "tracers.nc")
velocities = joinpath(output_path, "velocities.nc")
ps = joinpath(output_path, "plotting_snapshots.jld2")
snapshots = 1:25
jldopen(ps, "w") do file
    NCDataset(tracers) do ds

        t = ds["time"][:]
        for i ∈ snapshots
            file["S/xzslice_$(t[i])"] = ds[:S][:, 1, :, i]
            file["T/yzslice_$(t[i])"] = ds[:T][1, :, :, i]
        end

    end

    NCDataset(velocities) do ds

        t = ds["time"][:]
        for i ∈ snapshots
            file["w/w_zmean_$(t[i])"] = mean(ds[:w][:, :, :, i], dims = 3)
        end

    end
end
