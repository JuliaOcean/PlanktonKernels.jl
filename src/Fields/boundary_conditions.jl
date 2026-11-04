#####
##### generate boundary conditions for each nutrient field
#####

# generate the default boundary conditions for one field
function default_bcs()
    bcs = BoundaryConditions(nothing, nothing, nothing, nothing, nothing, nothing)
    return bcs
end

"Set a flux boundary value on a field, converting it to the field precision and architecture."
function set_bc!(field::Field{FT}, arch::Architecture; pos::Symbol, 
                 bc_c::Union{Nothing,Number,AbstractArray}, bc_v::Union{Number,AbstractArray}) where FT
    pos in fieldnames(BoundaryConditions) || throw(ArgumentError("Unknown boundary: $pos"))
    c = bc_c === nothing ? nothing : FT.(bc_c)
    v = FT.(bc_v)
    c isa AbstractArray && (c = array_type(arch)(c))
    v isa AbstractArray && (v = array_type(arch)(v))
    setproperty!(field.bc, pos, BoundaryCondition(c, v))
    return nothing
end

# Read one boundary value without allocating an array inside a kernel.
@inline boundary_value(x::Number, i, j, t) = x
@inline boundary_value(x::AbstractArray{T,2}, i, j, t) where T = x[i,j]
@inline boundary_value(x::AbstractArray{T,3}, i, j, t) where T = x[i,j,t]

# get boundary condition at each grid point
@inline function getbc(c, v, tr, pos::Symbol, i, j, t, g::AbstractGrid)
    if c === nothing
        if pos == :west
            concentration = tr[1+g.Hx, i+g.Hy, j+g.Hz]
        elseif pos == :east
            concentration = tr[g.Nx+g.Hx, i+g.Hy, j+g.Hz]
        elseif pos == :south
            concentration = tr[i+g.Hx, 1+g.Hy, j+g.Hz]
        elseif pos == :north
            concentration = tr[i+g.Hx, g.Ny+g.Hy, j+g.Hz]
        elseif pos == :top
            concentration = tr[i+g.Hx, j+g.Hy, 1+g.Hz]
        elseif pos == :bottom
            concentration = tr[i+g.Hx, j+g.Hy, g.Nz+g.Hz]
        end
    else
        concentration = boundary_value(c, i, j, t)
    end
    return boundary_value(v, i, j, t) * concentration
end

# Prescribed concentration and velocity at one boundary point.
@inline function getbc(c, v, i, j, t)
    return boundary_value(v, i, j, t) * boundary_value(c, i, j, t)
end

# validate boundary conditions, check if the grid information is compatible with nutrient field
function validate_bc(bc::Union{Nothing, BoundaryCondition}, bc_size, nΔT)
    bc === nothing && return nothing

    # Validate inputs independently, without allocating a flux array.
    # c === nothing selects the local tracer concentration for dilution.
    for (name, value) in ((:c, bc.c), (:v, bc.v))
        value === nothing && continue
        value isa Number && continue
        ndims(value) in (2, 3) ||
            throw(ArgumentError("BC $name arrays must have two or three dimensions"))
        expected = ndims(value) == 2 ? bc_size : (bc_size..., nΔT)
        size(value) == expected ||
            throw(ArgumentError("BC $name: grid/time mismatch, size must equal $expected, got $(size(value))."))
    end
    return nothing
end

function validate_bcs(nut, grid, nΔT)
    for name in keys(nut)
        bc_west    = nut[name].bc.west
        bc_east    = nut[name].bc.east
        bc_south   = nut[name].bc.south
        bc_north   = nut[name].bc.north
        bc_bottom  = nut[name].bc.bottom
        bc_top     = nut[name].bc.top

        bc_size_x = (grid.Ny, grid.Nz)
        bc_size_y = (grid.Nx, grid.Nz)
        bc_size_z = (grid.Nx, grid.Ny)
        
        validate_bc(bc_west, bc_size_x, nΔT)
        validate_bc(bc_east, bc_size_x, nΔT)
        validate_bc(bc_south, bc_size_y, nΔT)
        validate_bc(bc_north, bc_size_y, nΔT)
        validate_bc(bc_bottom, bc_size_z, nΔT)
        validate_bc(bc_top, bc_size_z, nΔT)
    end
end





