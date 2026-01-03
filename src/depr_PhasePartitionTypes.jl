
"""
    PhasePartitionTypes

Defines the PhasePartition struct for representing mass fractions of moist air.
"""

import DocStringExtensions

# Exports
export PhasePartition, promote_phase_partition

"""
    PhasePartition

Represents the mass fractions of the moist air mixture (the partitioning of water substance 
between vapor, liquid, and ice phases).

The total specific humidity `q_tot` represents the total water content, while `q_liq` and
`q_ice` represent the liquid and ice specific humidities, respectively. The vapor specific
humidity is computed as `q_vap = q_tot - q_liq - q_ice`.

# Constructors

    PhasePartition(q_tot::Real[, q_liq::Real[, q_ice::Real]])
    PhasePartition(param_set::APS, ts::ThermodynamicState)

See also [`PhasePartition_equil`](@ref)

# Fields

$(DocStringExtensions.FIELDS)
"""
struct PhasePartition{FT <: Real}
    "total specific humidity"
    tot::FT
    "liquid water specific humidity (default: `0`)"
    liq::FT
    "ice specific humidity (default: `0`)"
    ice::FT
    @inline function PhasePartition(tot::FT, liq::FT, ice::FT) where {FT}
        q_tot_safe = max(tot, 0)
        q_liq_safe = max(liq, 0)
        q_ice_safe = max(ice, 0)
        return new{FT}(q_tot_safe, q_liq_safe, q_ice_safe)
    end
end

PhasePartition(tot, liq, ice) = PhasePartition(promote(tot, liq, ice)...)

@inline Base.zero(::Type{PhasePartition{FT}}) where {FT} =
    PhasePartition(FT(0), FT(0), FT(0))

@inline PhasePartition(q_tot::FT, q_liq::FT) where {FT <: Real} =
    PhasePartition(q_tot, q_liq, zero(FT))
@inline PhasePartition(q_tot::FT) where {FT <: Real} =
    PhasePartition(q_tot, zero(FT), zero(FT))

Base.convert(::Type{PhasePartition{FT}}, q_pt::PhasePartition) where {FT} =
    PhasePartition(FT(q_pt.tot), FT(q_pt.liq), FT(q_pt.ice))

function promote_phase_partition(x, q_pt::PhasePartition)
    (x′, tot, liq, ice) = promote(x, q_pt.tot, q_pt.liq, q_pt.ice)
    return (x′, PhasePartition(tot, liq, ice))
end
function promote_phase_partition(x1, x2, q_pt::PhasePartition)
    (x1′, x2′, tot, liq, ice) = promote(x1, x2, q_pt.tot, q_pt.liq, q_pt.ice)
    return (x1′, x2′, PhasePartition(tot, liq, ice))
end
