using CairoMakie
using StaircaseShenanigans: animate_density, animate_density_anomaly,
                            animate_tracers, animate_tracers_anomaly

##
output_path = joinpath(@__DIR__, "lineareos_single_interface_240min_non_periodic/")
cd(output_path)
tracers = "tracers.nc"
co = "computed_output.nc"
xslice = 2
yslice = 2
animate_tracers(tracers; xslice, yslice)
animate_tracers_anomaly(tracers; xslice, yslice)
animate_density(co, "σ"; xslice, yslice)
animate_density_anomaly(co, "σ"; xslice, yslice)
