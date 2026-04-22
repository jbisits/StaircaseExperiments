include("plotting_env_and_load.jl")
## Figure five
# DNS flow evolution
files = (l_R_ρ_105_dT2_diagnostics, nl_R_ρ_105_dT2_diagnostics)
file = jldopen(files[2])
x = file["dims/x_caa"]
y = file["dims/y_aca"]
z = file["dims/z_aac"]
times = file["dims/time"]

x_xz = repeat(x, 1, length(z))
x_xy = repeat(x, 1, length(y))
y_xz = repeat(y, 1, length(z))
y_xy = repeat(y, 1, length(x))
z_xz = repeat(reshape(z, 1, length(z)), length(x), 1)
z_yz = repeat(reshape(z, 1, length(z)), length(x), 1)

x_xz_velocity = x[1] * ones(length(x), length(z))
x_xz_velocity_end = x[end] * ones(length(x), length(z))
y_xz_density = y[1] * ones(length(x), length(z))
y_xz_density_end = y[end] * ones(length(x), length(z))

z_xy_top = z[end] * ones(length(x), length(y))
close(file)

horizontal_ticks = LinRange(-0.025, 0.025, 5)
T_colorrange = (-1.5, 1.5)
S_colorrange = (-0.09, 0.09)
w_colorrnage =  (-0.0005, 0.0005)
fig = Figure(size = (1500, 1300))
snapshots = [[60.0 * 3, 60.0 * 6, 60 * 12, 60.0 * 24],
             [180.00000000000003, 60.0 * 6, 60 * 12, 60.0 * 24]]
panel_labels = ["(a)" "(b)" "(c)" "(d)"; "(e)" "(f)" "(g)" "(h)"]
for (j, path) ∈ enumerate(files)

    file = jldopen(path)

    initial_S = file["S"]["xzslice_0.0"]
    initial_T = file["T"]["yzslice_0.0"]
    slices = [(S_xz = file["S"]["xzslice_$(snapshot)"],
            T_yz = file["T"]["yzslice_$(snapshot)"],
            velocity_zmean = file["w"]["w_zmean_$(snapshot)"][:, :, 1]) for snapshot ∈ snapshots[j]]
    close(file)

    ax = [Axis3(fig[j, i],
                aspect=(1/3, 1/3, 1),
                titlefont = :regular,
                xlabel = "x (cm)",
                ylabel = "y (cm)",
                zlabel = "z (m)",
                xticks = (horizontal_ticks, string.(horizontal_ticks .* 100)),
                yticks = (horizontal_ticks, string.(horizontal_ticks .* 100)),
                xlabeloffset = 40,
                ylabeloffset = 40,
                zlabeloffset = 50,
                xlabelsize = 14,
                ylabelsize = 14,
                zlabelsize = 16,
                xticklabelsize = 12,
                yticklabelsize = 12,
                zticklabelsize = 16,
                zlabelrotation = π / 2,
                limits = ((x[1], x[end]), (y[1], y[end]), (z[1], z[end])),
                elevation = π / 6.5,
                azimuth = 1.25π,
                xspinesvisible = false,
                yspinesvisible = false,
                zspinesvisible = false,
                zgridvisible = false,
                protrusions = 20
                ) for i ∈ 1:4]

    for i ∈ eachindex(snapshots[j])
        sf_S = surface!(ax[i], x_xz, y_xz_density, z_xz; color = slices[i].S_xz .- initial_S,
                        colormap = :curl,
                        colorrange = S_colorrange)
        sf_T = surface!(ax[i], x_xz_velocity, y_xz, z_yz; color = slices[i].T_yz .- initial_T,
                        colormap = :delta, colorrange = T_colorrange,
                        backlight = 5f0, shading = FastShading)
        sf_w = surface!(ax[i], x, y, z_xy_top; color = slices[i].velocity_zmean, colormap = :balance,
                        colorrange = w_colorrnage)
        ax[i].title = j == 1 ? panel_labels[j, i] * " I, t = $(round(snapshots[j][i] / 60)) minutes" :
                               panel_labels[j, i] * " IV, t = $(round(snapshots[j][i] / 60)) minutes"
        if j == 2
            if i == 1
                Colorbar(fig[j+1, 1:2], sf_S, label = "S′ (gkg⁻¹)", vertical = false, flipaxis = false)
                Colorbar(fig[j+1, 3:4], sf_T, label = "Θ′ (°C)", vertical = false, flipaxis = false)
                Colorbar(fig[:, 5], sf_w, label = "w (ms⁻¹)")
            end
        end
        if i > 1
            hidezdecorations!(ax[i], ticks = false)
        end
    end
end
Label(fig[0, :], "Salinity and temperature evolution", fontsize = 22, font = :bold)
rowgap!(fig.layout, 2, Relative(0.05))
rowgap!(fig.layout, 3, Relative(0.05))
##
save("fig5.png", fig)
