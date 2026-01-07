import RootSolvers as RS
import Plots

import Thermodynamics as TD
import Thermodynamics.Parameters as TP
import ClimaParams as CP

FT = Float64

const param_set = TP.ThermodynamicsParameters(FT)

# Use joinpath and @__DIR__ to robustly locate the test file regardless of CWD
include(joinpath(@__DIR__, "../../test/TestedProfiles.jl"))
import .TestedProfiles

profiles = TestedProfiles.PhaseEquilProfiles(param_set, Array{FT});
(; ρ, q_tot) = profiles
T_true = profiles.T
prof_pts = (ρ, T_true, q_tot)

dims = (6, 6, 6);
ρ = range(min(ρ...), stop = max(ρ...), length = dims[1]);
T_true = range(min(T_true...), stop = max(T_true...), length = dims[2]);
q_tot = range(min(q_tot...), stop = max(q_tot...), length = dims[3]);

ρ_all = Array{FT}(undef, prod(dims));
T_true_all = Array{FT}(undef, prod(dims));
q_tot_all = Array{FT}(undef, prod(dims));

linear_indices = LinearIndices((1:dims[1], 1:dims[2], 1:dims[3]));

numerical_methods = (
    RS.SecantMethod,
    RS.NewtonsMethod,
    # RS.NewtonsMethodAD, # we need to relax diagonalization somewhere
    # RS.RegulaFalsiMethod, # hit assertion error, bounds need adjusted for 3D space
)

ts = Dict(
    NM => Array{Union{TD.ThermodynamicState, Nothing}}(undef, prod(dims))
    for NM in numerical_methods
)
ts_no_err = Dict(
    NM => Array{TD.ThermodynamicState}(undef, prod(dims)) for
    NM in numerical_methods
)

@inbounds for i in linear_indices.indices[1]
    @inbounds for j in linear_indices.indices[2]
        @inbounds for k in linear_indices.indices[3]
            n = linear_indices[i, j, k]
            ρ_all[n] = ρ[i]
            T_true_all[n] = T_true[j]
            q_tot_all[n] = q_tot[k]
            e_int = TD.internal_energy(
                param_set,
                T_true[j],
                TD.PhasePartition(q_tot[k]),
            )

            @inbounds for NM in numerical_methods
                # TD.set_error_on_non_convergence!(false)
                ts_no_err[NM][n] = TD.PhaseEquil_dev_only(
                    param_set,
                    ρ[i],
                    e_int,
                    q_tot[k];
                    sat_adjust_method = NM,
                    maxiter = 10,
                )
                # TD.set_error_on_non_convergence!(true)
                # @show n / prod(dims) * 100
                try
                    ts[NM][n] = TD.PhaseEquil_dev_only(
                        param_set,
                        ρ[i],
                        e_int,
                        q_tot[k];
                        sat_adjust_method = NM,
                        maxiter = 10,
                    )
                catch
                    ts[NM][n] = nothing
                end
            end
        end
    end
end

# folder = "sat_adjust_analysis"
folder = @__DIR__
mkpath(folder)

function save_binned_surface_plots(
    x,
    y,
    z,
    title,
    filename,
    n_plots = (3, 3),
    z_label_prefix = "z",
    n_digits = 5;
    xlims = (:auto, :auto),
    ylims = (:auto, :auto),
    label = label,
    ref_points = nothing,
)
    n_z_partitions = prod(n_plots)
    local z_min_global, z_max_global
    try
        z_min_global = min(z...)
        z_max_global = max(z...)
    catch
        z_min_global = 0
        z_max_global = 1
    end
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
        p_i = Plots.plot(
            x_i,
            y_i,
            title = "$(title), in ($sz_min, $sz_max)",
            seriestype = :scatter,
            markersize = 5,
            label = label,
            xlims = xlims,
            ylims = ylims,
        )
        if ref_points ≠ nothing
            ref_mask = z_min[i] .<= ref_points[3] .<= z_max[i]
            Plots.plot!(
                ref_points[1][ref_mask],
                ref_points[2][ref_mask],
                seriestype = :scatter,
                markersize = 5,
                label = "ref points",
            )
        end
        push!(p, p_i)
    end
    Plots.plot(p..., layout = n_plots, legend = false)
    Plots.savefig(filename)
end;
# Full 3D scatter plot
function plot3D(ts_no_err, ts, NM; converged)
    mask = converged ? ts .≠ nothing : ts .== nothing

    c_name = converged ? "converged" : "non_converged"
    label = converged ? "converged" : "non-converged"
    casename = converged ? "converged" : "non-converged"
    nm_name = nameof(NM)
    filename = "3DSpace_$(c_name)_$nm_name.svg"

    ρ_mask = ρ_all[mask]
    q_tot_mask = q_tot_all[mask]
    T_mask = T_true_all[mask]
    pts = (ρ_mask, T_mask, q_tot_mask)
    Plots.plot(
        pts...,
        color = "blue",
        seriestype = :scatter,
        markersize = 7,
        label = casename,
    )
    Plots.plot!(
        prof_pts...,
        color = "red",
        seriestype = :scatter,
        markersize = 7,
        label = "tested thermo profiles",
    )
    Plots.plot!(
        xlabel = "Density",
        ylabel = "Temperature",
        zlabel = "Total specific humidity",
        title = "3D input to PhaseEquil",
        xlims = (min(ρ_all...), max(ρ_all...)),
        ylims = (min(T_true_all...), max(T_true_all...)),
        zlims = (min(q_tot_all...), max(q_tot_all...)),
        camera = (50, 50),
        # camera = (50,70),
    )
    Plots.savefig(joinpath(folder, filename))
end

# 2D binned scatter plots
function plot2D_slices(ts_no_err, ts, NM; converged)
    mask = converged ? ts .≠ nothing : ts .== nothing
    ρ_mask = ρ_all[mask]
    T_mask = T_true_all[mask]
    q_tot_mask = q_tot_all[mask]
    c_name = converged ? "converged" : "non_converged"
    label = converged ? "converged" : "non-converged"
    short_name = converged ? "C" : "NC"
    nm_name = nameof(NM)
    filename = "2DSlice_$(c_name)_$nm_name.svg"
    filename = joinpath(folder, filename)
    save_binned_surface_plots(
        ρ_mask,
        T_mask,
        q_tot_mask,
        short_name,
        filename;
        xlims = (min(ρ_all...), max(ρ_all...)),
        ylims = (min(T_true_all...), max(T_true_all...)),
        label = label,
        ref_points = prof_pts,
    )
end

# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "nul"

for NM in numerical_methods
    plot3D(ts_no_err[NM], ts[NM], NM; converged = false)
    plot3D(ts_no_err[NM], ts[NM], NM; converged = true)
    plot2D_slices(ts_no_err[NM], ts[NM], NM; converged = true)
    plot2D_slices(ts_no_err[NM], ts[NM], NM; converged = false)
end

convergence_percent = Dict()
for NM in numerical_methods
    convergence_percent[NM] = count(ts[NM] .≠ nothing) / length(ts[NM])
end
println("Convergence percentages:")
for (k, v) in convergence_percent
    println("$k = $v")
end
