#=
This file is for updating the NCDataset database that stores
input values to the thermodynamic state constructors
which have caused convergence issues. Updating this database
allows us to optimize the convergence rate of the thermodynamics
constructor for a variety of realistic input values.

The database stores problematic input combinations that cause
saturation adjustment algorithms to fail or converge slowly.
This helps developers identify and fix convergence issues.
=#

using DelimitedFiles

"""
    compress_targz(file)

Platform-independent file compression using tar and gzip.
Currently only supports Unix-like systems (Linux, macOS).

# Arguments
- `file`: Output file path for the compressed archive

# Throws
- `ErrorException`: On Windows systems (not yet implemented)
"""
function compress_targz(file)
    if Sys.iswindows()
        error("Windows compression not yet implemented")
    else
        run(`tar -zcvf $file $(readdir())`)
    end
end

# Configuration flags for database management
const archive = true    # Whether to create compressed archive
const clean = true      # Whether to clean existing database files
const append = false    # Whether to append to existing database
const create = true     # Whether to create new database files

# Validate configuration - prevent conflicting operations
@assert !(append && create) # cannot create and append
@assert !(append && clean)  # cannot clean and append

# Define paths for database management
folder = joinpath(@__DIR__, "MTConstructorData")
output_file = joinpath(@__DIR__, "MTConstructorDataZipped.tar.gz")
mkpath(folder)

# Define constructors and their required parameters
# This maps constructor names to their input parameter names
constructors = Dict("PhaseEquil" => (:ρ, :e_int, :q_tot))

"""
    get_nc(k)

Get the NetCDF file path for a given constructor key.

# Arguments
- `k`: Constructor key (e.g., "PhaseEquil")

# Returns
- Path to the corresponding NetCDF file
"""
get_nc(k) = joinpath(folder, "test_data_$(k).nc")

"""
    get_data_to_append(k)

Get the CSV file path for data to append to a constructor database.

# Arguments
- `k`: Constructor key (e.g., "PhaseEquil")

# Returns
- Path to the corresponding CSV file
"""
get_data_to_append(k) = joinpath(@__DIR__, "test_data_$(k).csv")

# Clean existing database files if requested
if clean
    for k in keys(constructors)
        rm(get_nc(k); force = true)
    end
end

# Process data for each constructor
FT = Float64
if create || append
    for (k, v) in constructors
        if isfile(get_data_to_append(k))
            # Open NetCDF file in appropriate mode
            if append
                ds = Dataset(get_nc(k), "a")
            elseif create
                ds = Dataset(get_nc(k), "c")
            end

            # Read CSV data and convert to NetCDF format
            data_all = readdlm(get_data_to_append(k), ',')

            # Store each parameter as a separate variable
            for (i, _v) in enumerate(v)
                vardata = Array{FT}(data_all[2:end, i])  # Skip header row
                s = string(_v)
                defVar(ds, s, vardata, ("datapoint",))
            end

            close(ds)
        end
    end
end

# Create compressed archive if requested
if archive
    cd(folder) do
        compress_targz(output_file)
    end
end
