using Plots

function save_binned_surface_plots(
    x,
    y,
    z,
    title,
    filename,
    n_plots = (3, 3),
    z_label_prefix = "z",
    n_digits = 5,
)
    n_z_partitions = prod(n_plots)
    z_min_global = min(z...)
    z_max_global = max(z...)
    Δz = (z_max_global - z_min_global) / n_z_partitions
    z_min = ntuple(i -> z_min_global + (i - 1) * Δz, n_z_partitions)
    z_max = ntuple(i -> z_min_global + (i) * Δz, n_z_partitions)
    p = []
    for i in 1:n_z_partitions
        mask = z_min[i] .<= z .<= z_max[i]
        x_i = x[mask]
        y_i = y[mask]
        sz_min = string(z_min[i])[1:min(n_digits, length(string(z_min[i])))]
        sz_max = string(z_max[i])[1:min(n_digits, length(string(z_max[i])))]
        p_i = plot(
            x_i,
            y_i,
            title = "$(title), in ($sz_min, $sz_max)",
            seriestype = :scatter,
            markersize = 5,
        )
        push!(p, p_i)
    end
    plot(p..., layout = n_plots, legend = false)
    savefig(filename)
end;
