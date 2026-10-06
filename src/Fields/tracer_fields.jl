"""
    Field(arch::Architecture, grid::AbstractGrid, FT::DataType; bcs = default_bcs())
Construct a `Field` on `grid` with data and boundary conditions on architecture `arch`
with DataType `FT`.
"""
function Field(arch::Architecture, grid::AbstractGrid, FT::DataType; bcs = default_bcs())
    total_size = (grid.Nx+grid.Hx*2, grid.Ny+grid.Hy*2, grid.Nz+grid.Hz*2)
    data = zeros(FT, total_size) |> array_type(arch)
    return Field{FT}(data,bcs)
end

@inline interior(c, grid) = c[grid.Hx+1:grid.Hx+grid.Nx, grid.Hy+1:grid.Hy+grid.Ny, grid.Hz+1:grid.Hz+grid.Nz]
@inline interior(field::Field, grid::AbstractGrid) = interior(field.data, grid)

function zero_fields!(a)
    for tr in keys(a)
        @inbounds a[tr].data .= 0.0f0
    end
end

"""
    init_tracers(arch, grid, names::Tuple, FT=Float32)
Allocate independent fields for the caller-supplied tracer names.
"""
function init_tracers(arch, g, names::Tuple, FT=Float32)
    all(n -> n isa Symbol, names) || throw(ArgumentError("Tracer names must be Symbols"))
    length(unique(names)) == length(names) || throw(ArgumentError("Tracer names must be unique"))
    return NamedTuple{names}(Tuple(Field(arch, g, FT) for _ in names))
end

function initialize_tracer!(tr, g::AbstractGrid, source::Union{String, NamedTuple}, arch::Architecture, FT::DataType=Float32)
    if typeof(source) <: String
        tmp = deserialize(source)
        all(x -> x >= 0, tmp) || throw(ArgumentError("TRAC_INIT: initial concentrations must be nonnegative"))
        tmp = tmp |> array_type(arch)
        if size(tmp) == (g.Nx, g.Ny, g.Nz)
            @views @. tr.data[g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz] = FT.(tmp[:,:,:])
        else
            throw(ArgumentError("TRAC_INIT:  grid mismatch"))
        end
    elseif typeof(source) <: NamedTuple
        if source.init < 0.0
            throw(ArgumentError("TRAC_INIT:  The initial condition should be none-negetive."))
        end
        lower = FT(1.0 - source.rand_noise)
        upper = FT(1.0 + source.rand_noise)
        tr.data .= fill(FT(source.init), size(tr.data)) .* rand(lower:1.0f-4:upper, size(tr.data)) |> array_type(arch)
    else
        throw(ArgumentError("TRAC_INIT:  source should be a String or NamedTuple"))
    end
    return nothing
end

"""
    generate_tracers(arch, grid, tracer_names, source, FT)
Set up initial tracer fields according to `grid`.

Arguments
=================
- `arch`: `CPU()` or `GPU()`. The computer architecture used to time-step `model`.
- `grid`: The resolution and discrete geometry on which nutrient fields are solved.
- `tracer_names`: A tuple of symbols representing the tracer names.
- `source`: A `NamedTuple` mapping each tracer to `(init=value, rand_noise=noise)` or a serialized array path.
- `FT`: Floating point data type. Default: `Float32`.
"""
function generate_tracers(arch::Architecture, grid::AbstractGrid, tracer_names,
                            source::NamedTuple, FT::DataType=Float32)
    tracers = init_tracers(arch, grid, tracer_names, FT)

    for name in tracer_names
        if !haskey(source, name)
            throw(ArgumentError("TRAC_INIT: tracer not found $(name)"))
        else
            initialize_tracer!(tracers[name], grid, source[name], arch, FT)
        end

        @views @. tracers[name].data *= grid.landmask
    end

    fill_halo_tracer!(tracers,grid)

    return tracers
end
