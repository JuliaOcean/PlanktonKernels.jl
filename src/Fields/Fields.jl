module Fields
using KernelAbstractions
using Random
using Adapt
using Serialization: deserialize
using ..Architectures: Architecture, array_type, device
using ..Grids
using ..Grids: Ax, Ay, Az, volume

export vel_copy!, copy_interior!
export Field, BoundaryConditions, interior, init_tracers, initialize_tracer!, generate_tracers
export set_bc!

mutable struct BoundaryCondition
    c::Union{Nothing, Number, AbstractArray}  # input concentration or flux boundary value （mmol/m³)
    v::Union{Number, AbstractArray}           # input velocity boundary value (m/s)
end

mutable struct BoundaryConditions
    west::Union{Nothing, BoundaryCondition}
    east::Union{Nothing, BoundaryCondition}
    north::Union{Nothing, BoundaryCondition}
    south::Union{Nothing, BoundaryCondition}
    top::Union{Nothing, BoundaryCondition}
    bottom::Union{Nothing, BoundaryCondition}
end

struct Field{FT}
    data::AbstractArray{FT,3}
    bc::BoundaryConditions
end

include("halo_regions.jl")
include("copy_fields.jl")
include("boundary_conditions.jl")
include("tracer_fields.jl")
include("utils.jl")
include("apply_bcs.jl")

end
