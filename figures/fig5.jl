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

# horizontal_ticks = LinRange(-0.025, 0.025, 5)
horizontal_ticks = [-0.02, -0.01, 0, 0.01, 0.02]
st_horizontal_ticks = string.(round.(Int64, horizontal_ticks .* 100))
T_colorrange = (-1.5, 1.5)
S_colorrange = (-0.09, 0.09)
w_colorrnage = (-0.0005, 0.0005)
# fig = Figure(size=(800, 1800))
# snapshots = [[60.0 * 3, 60.0 * 6, 60 * 12, 60.0 * 24],
#     [180.00000000000003, 60.0 * 6, 60 * 12, 60.0 * 24]]
# panel_labels = ["(a)" "(b)" "(c)" "(d)"; "(e)" "(f)" "(g)" "(h)"]
# for j ∈ eachindex(files)
#
#     file = jldopen(files[j])
#
#     initial_S = file["S"]["xzslice_0.0"]
#     initial_T = file["T"]["yzslice_0.0"]
#     slices = [(S_xz=file["S"]["xzslice_$(snapshot)"],
#         T_yz=file["T"]["yzslice_$(snapshot)"],
#         velocity_zmean=file["w"]["w_zmean_$(snapshot)"][:, :, 1]) for snapshot ∈ snapshots[j]]
#     close(file)
#
#     ax = [Axis3(fig[k+(j-1)*3, i],
#         aspect=(1 / 3, 1 / 3, 1),
#         titlefont=:regular,
#         xlabel="x (cm)",
#         ylabel="y (cm)",
#         zlabel="z (m)",
#         xticks=(horizontal_ticks, string.(horizontal_ticks .* 100)),
#         yticks=(horizontal_ticks, string.(horizontal_ticks .* 100)),
#         xlabeloffset=40,
#         ylabeloffset=40,
#         zlabeloffset=60,
#         xlabelsize=18,
#         ylabelsize=18,
#         zlabelsize=18,
#         xticklabelsize=14,
#         yticklabelsize=14,
#         zticklabelsize=18,
#         zlabelrotation=π / 2,
#         limits=((x[1], x[end]), (y[1], y[end]), (z[1], z[end])),
#         elevation=π / 6.5,
#         azimuth=1.25π,
#         xspinesvisible=false,
#         yspinesvisible=false,
#         zspinesvisible=false,
#         zgridvisible=false,
#         protrusions=(0, 0, 7, 7)
#     ) for i ∈ 1:2, k ∈ 1:2]
#
#     if j == 2
#         for i ∈ 3:4
#             ax[i].protrusions = (0, 0, 40, 7)
#         end
#     end
#     for i ∈ eachindex(snapshots[j])
#         sf_S = surface!(ax[i], x_xz, y_xz_density, z_xz; color=slices[i].S_xz .- initial_S,
#             colormap=:curl,
#             colorrange=S_colorrange)
#         sf_T = surface!(ax[i], x_xz_velocity, y_xz, z_yz; color=slices[i].T_yz .- initial_T,
#             colormap=:delta, colorrange=T_colorrange,
#             backlight=5f0, shading=FastShading)
#         sf_w = surface!(ax[i], x, y, z_xy_top; color=slices[i].velocity_zmean, colormap=:balance,
#             colorrange=w_colorrnage)
#         ax[i].title = j == 1 ? panel_labels[j, i] * " t = $(round(snapshots[j][i] / 60)) minutes" :
#                       panel_labels[j, i] * " t = $(round(snapshots[j][i] / 60)) minutes"
#         if j == 2
#             if i == 1
#                 Colorbar(fig[0:1, 3], sf_T, label="Θ′ (°C)")
#                 Colorbar(fig[2:3, 3], sf_S, label="S′ (gkg⁻¹)")
#                 Colorbar(fig[4, 3], sf_w, label="w (ms⁻¹)")
#             end
#         end
#         if j == 1
#             hidexdecorations!(ax[i], ticks=false)
#             hideydecorations!(ax[i], ticks=false)
#         end
#         if j == 2 && i < 3
#             hidexdecorations!(ax[i], ticks=false)
#             hideydecorations!(ax[i], ticks=false)
#         end
#         if i % 2 == 0
#             hidezdecorations!(ax[i], ticks=false)
#         end
#     end
# end
# Label(fig[0, :], "Linear equation of state", fontsize=22, font=:bold)
# Label(fig[3, :], "Non-linear equation of state", fontsize=22, font=:bold)
# colgap!(fig.layout, 0)
fig = Figure(size=(1100, 1000))
snapshots = [[60.0 * 3, 60.0 * 6, 60 * 12, 60.0 * 24],
    [180.00000000000003, 60.0 * 6, 60 * 12, 60.0 * 24]]
for (j, path) ∈ enumerate(files)

    file = jldopen(path)

    initial_S = file["S"]["xzslice_0.0"]
    initial_T = file["T"]["yzslice_0.0"]
    slices = [(S_xz=file["S"]["xzslice_$(snapshot)"],
        T_yz=file["T"]["yzslice_$(snapshot)"],
        velocity_zmean=file["w"]["w_zmean_$(snapshot)"][:, :, 1]) for snapshot ∈ snapshots[j]]
    close(file)

    ax = [Axis3(fig[j+(j-1), i],
        aspect=(1, 1, 2.5),
        titlefont=:regular,
        xlabel="x (cm)",
        ylabel="y (cm)",
        zlabel="z (m)",
        xticks=(horizontal_ticks, st_horizontal_ticks),
        yticks=(horizontal_ticks, st_horizontal_ticks),
        xlabeloffset=40,
        ylabeloffset=40,
        zlabeloffset=60,
        xlabelsize=18,
        ylabelsize=18,
        zlabelsize=18,
        xticklabelsize=18,
        yticklabelsize=18,
        zticklabelsize=18,
        zlabelrotation=π / 2,
        limits=((x[1], x[end]), (y[1], y[end]), (z[1], z[end])),
        elevation=π / 6.5,
        azimuth=1.25π,
        xspinesvisible=false,
        yspinesvisible=false,
        zspinesvisible=false,
        zgridvisible=false,
        protrusions=j == 1 ? (8, 6, 1, 8) : (8, 6, 40, 8)
    ) for i ∈ 1:4]

    # j == 2 ? ax[1].protrusions = (50, 10, 40, 4) : (40, 10, 1, 8)
    ax[1].protrusions = (60, 6, 40, 4)
    for i ∈ eachindex(snapshots[j])
        sf_S = surface!(ax[i], x_xz, y_xz_density, z_xz; color=slices[i].S_xz .- initial_S,
            colormap=:curl,
            colorrange=S_colorrange)
        sf_T = surface!(ax[i], x_xz_velocity, y_xz, z_yz; color=slices[i].T_yz .- initial_T,
            colormap=:delta, colorrange=T_colorrange,
            backlight=5f0, shading=FastShading)
        sf_w = surface!(ax[i], x, y, z_xy_top; color=slices[i].velocity_zmean, colormap=:balance,
            colorrange=w_colorrnage)
        ax[i].title = panel_labels[j, i] * " t = $(round(snapshots[j][i] / 60)) minutes"
        if i == 1
            text!(ax[i], Point3f(y_xz_density[1, 1], 0, -0.125),
                text=L"$\Theta'$",
                color=:black,
                align=(:center, :center, :center),
                fontsize=22,
            )
            text!(ax[i], Point3f(0, x_xz_velocity[1, 1], -0.125),
                text=L"$S'$",
                color=:black,
                align=(:center, :center, :center),
                fontsize=22,
            )
            text!(ax[i], Point3f(0, 0, z_xy_top[1, 1]),
                text=L"$\langle w \rangle_{\!z}$",
                # rotation=3π / 4,
                color=:black,
                align=(:center, :center, :center),
                fontsize=2.5e-1,
                markerspace=:data
            )
        end
        if j == 2
            if i == 1
                Colorbar(fig[j+2, 1:2], sf_T, label=L"$\Theta'$ ($^{\circ}$C)", vertical=false, flipaxis=false, labelsize=22)
                Colorbar(fig[j+2, 3:4], sf_S, label=L"$S'$ (gkg$^{-1}$)", vertical=false, flipaxis=false, labelsize=22)
                Colorbar(fig[:, 5], sf_w, label=L"$\langle w \rangle_{z}$ (ms$^{-1}$)", labelsize=22)
            end
        end
        if i > 1
            hidexdecorations!(ax[i], ticks=false)
            hideydecorations!(ax[i], ticks=false)
            hidezdecorations!(ax[i], ticks=false)
        end
        # if j == 1 && i == 1
        #     hidexdecorations!(ax[i], ticks=false)
        #     hideydecorations!(ax[i], ticks=false)
        #     hidezdecorations!(ax[i], ticks=false)
        # end
    end
end
Label(fig[0, :], "Linear equation of state", fontsize=22, font=:bold)
Label(fig[2, :], "Non-linear equation of state", fontsize=22, font=:bold)
# rowgap!(fig.layout, 2, Relative(0.05))
# rowgap!(fig.layout, 3, Relative(0.05))
##
save("fig5.png", fig)
