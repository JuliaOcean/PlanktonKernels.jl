module Fields
using KernelAbstractions
using Random
using Serialization: deserialize
using ..Architectures: Architecture, array_type, device
using ..Grids
using ..Grids: Ax, Ay, Az, volume

export vel_copy!, copy_interior!
export Field, BoundaryConditions, interior, init_tracers, initialize_tracer!, generate_tracers
export set_bc!

mutable struct BoundaryConditions
    west::Union{Nothing, Number, AbstractArray}
    east::Union{Nothing, Number, AbstractArray}
    north::Union{Nothing, Number, AbstractArray}
    south::Union{Nothing, Number, AbstractArray}
    top::Union{Nothing, Number, AbstractArray}
    bottom::Union{Nothing, Number, AbstractArray}
end

struct Field{FT}
    data::AbstractArray{FT,3}
    bc::BoundaryConditions
end

include("halo_regions.jl")
include("copy_fields.jl")
include("boundary_conditions.jl")
include("tracer_fields.jl")
include("apply_bcs.jl")

end
