using Thermodynamics
using Thermodynamics.TemperatureProfiles
using Thermodynamics.TestedProfiles
using UnPack
using CLIMAParameters
using CLIMAParameters.Planet
using Plots
import Thermodynamics
Thermodynamics.print_warning() = false

struct EarthParameterSet <: AbstractEarthParameterSet end;
const param_set = EarthParameterSet();
FT = Float64;
thermo_dir = dirname(dirname(pathof(Thermodynamics)));
include(joinpath(thermo_dir, "docs", "plot_helpers.jl"));
profiles = TestedProfiles.PhaseEquilProfiles(param_set, Array{FT});
@unpack ρ, e_int, q_tot = profiles

dims = (10, 10, 10);
ρ = range(min(ρ...), stop = max(ρ...), length = dims[1]);
e_int = range(min(e_int...), stop = max(e_int...), length = dims[2]);
q_tot = range(min(q_tot...), stop = max(q_tot...), length = dims[3]);

ρ_all = Array{FT}(undef, prod(dims));
e_int_all = Array{FT}(undef, prod(dims));
q_tot_all = Array{FT}(undef, prod(dims));

linear_indices = LinearIndices((1:dims[1], 1:dims[2], 1:dims[3]));
TS = Array{Union{ThermodynamicState, Nothing}}(undef, prod(dims));
TS_no_err = Array{ThermodynamicState}(undef, prod(dims));

@inbounds for i in linear_indices.indices[1]
    @inbounds for j in linear_indices.indices[2]
        @inbounds for k in linear_indices.indices[3]
            p = linear_indices[i, j, k]
            ρ_all[p] = ρ[i]
            e_int_all[p] = e_int[j]
            q_tot_all[p] = q_tot[k]

            Thermodynamics.error_on_non_convergence() = false
            TS_no_err[p] = PhaseEquil(param_set, e_int[j], ρ[i], q_tot[k])
            Thermodynamics.error_on_non_convergence() = true
            # @show p/prod(linear_indices.indices)*100
            try
                TS[p] = PhaseEquil(param_set, e_int[j], ρ[i], q_tot[k])
            catch
                TS[p] = nothing
            end
        end
    end
end

# Full 3D scatter plot
function save_masked_plot3D(TS_no_err, mask, title, filename)
    ρ_mask = ρ_all[mask]
    T_mask = air_temperature.(TS_no_err[mask])
    q_tot_mask = q_tot_all[mask]
    pts = (ρ_mask, T_mask, q_tot_mask)
    Plots.plot(pts..., seriestype = :scatter, markersize = 7)
    plot!(
        xlabel = "Density",
        ylabel = "Temperature",
        zlabel = "Total specific humidity",
        title = "$title",
    )
    savefig(filename)
end;

save_masked_plot3D(TS_no_err, TS .== nothing, "NC", "Scatter3DNonConverged.svg");
save_masked_plot3D(TS_no_err, TS .!= nothing, "C", "Scatter3DConverged.svg");

# 2D binned scatter plots
function save_masked_plot2D_slices(TS_no_err, mask, title, filename)
    ρ_mask = ρ_all[mask]
    T_mask = air_temperature.(TS_no_err[mask])
    q_tot_mask = q_tot_all[mask]
    save_binned_surface_plots(ρ_mask, T_mask, q_tot_mask, title, filename)
end;

save_masked_plot2D_slices(
    TS_no_err,
    TS .== nothing,
    "NC",
    "Slices2DNonConverged.svg",
);
save_masked_plot2D_slices(
    TS_no_err,
    TS .!= nothing,
    "C",
    "Slices2DConverged.svg",
);

@warn "Note that the temperature axis for the non-converged
plot is not necessarily accurate, since the temperatures are
the result of a non-converged saturation adjustment"
