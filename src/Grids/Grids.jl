module Grids

export AbstractGrid
export RectilinearGrid
export LatLonGrid
export LoadLatLonGrid
export Periodic, Bounded

using Adapt

using ..Architectures
using ..Architectures: array_type

"""
    AbstractGrid{FT, TX, TY, TZ}
Abstract type for grids with elements of type `FT` and topology `{TX, TY, TZ}`.
"""
abstract type AbstractGrid{FT, TX, TY, YZ} end

"""
    AbstractTopology
Abstract type for grid topologies.
"""
abstract type AbstractTopology end

"""
    Periodic
Grid topology for periodic dimensions.
"""
struct Periodic <: AbstractTopology end

"""
    Bounded
Grid topology for bounded dimensions.
"""
struct Bounded <: AbstractTopology end

import Base: show

include("rectilinear_grid.jl")
include("lat_lon_grid.jl")
include("utils.jl")
include("areas_volumes_spacings.jl")

end
