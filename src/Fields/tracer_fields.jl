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
    tracers_init(arch, grid, names::Tuple, FT=Float32)
Allocate independent fields for the caller-supplied tracer names.
"""
function tracers_init(arch, g, names::Tuple, FT=Float32)
    all(n -> n isa Symbol, names) || throw(ArgumentError("Tracer names must be Symbols"))
    length(unique(names)) == length(names) || throw(ArgumentError("Tracer names must be unique"))
    return NamedTuple{names}(Tuple(Field(arch, g, FT) for _ in names))
end

"""
    generate_tracers(arch, grid, tracer_names, source, FT)
Set up initial tracer fields according to `grid`.

Arguments
=================
- `arch`: `CPU()` or `GPU()`. The computer architecture used to time-step `model`.
- `grid`: The resolution and discrete geometry on which nutrient fields are solved.
- `tracer_names`: A tuple of symbols representing the tracer names.
- `source`: A `NamedTuple` with `initial_condition` and `rand_noise` NamedTuples,
            or a `Dict` mapping tracer names to serialized arrays. Names are inferred
            from the source. File arrays must match the interior grid size.
- `FT`: Floating point data type. Default: `Float32`.
"""
function generate_tracers(arch::Architecture, grid::AbstractGrid, tracer_names,
                            source::Union{Dict,NamedTuple}, FT::DataType=Float32)
    total_size = (grid.Nx+grid.Hx*2, grid.Ny+grid.Hy*2, grid.Nz+grid.Hz*2)
    tracers = tracers_init(arch, grid, tracer_names, FT)
    pathkeys = collect(keys(source))

    if typeof(source) <: NamedTuple
        pathkeys = collect(keys(source.initial_condition))
    end

    for name in tracer_names
        if length(findall(x->x==name, pathkeys)) == 0
            throw(ArgumentError("TRAC_INIT: tracer not found $(name)"))
        else
            if typeof(source) <: Dict # file paths
                tmp = deserialize(source[name]) |> array_type(arch)
                if sum(tmp .< 0.0) > 0
                    throw(ArgumentError("TRAC_INIT: The initial condition should be none-negetive."))
                end

                if size(tmp) == (grid.Nx, grid.Ny, grid.Nz)
                    @views @. tracers[name].data[grid.Hx+1:grid.Hx+grid.Nx, grid.Hy+1:grid.Hy+grid.Ny, grid.Hz+1:grid.Hz+grid.Nz] = FT.(tmp[:,:,:])
                else
                    throw(ArgumentError("TRAC_INIT:  grid mismatch"))
                end
            elseif typeof(source) <: NamedTuple # value
                if source.initial_condition[name] < 0.0
                    throw(ArgumentError("TRAC_INIT:  The initial condition should be none-negetive."))
                end
                lower = FT(1.0 - source.rand_noise[name])
                upper = FT(1.0 + source.rand_noise[name])
                tracers[name].data .= fill(FT(source.initial_condition[name]),total_size) .* rand(lower:1.0f-4:upper, total_size) |> array_type(arch)
            end
        end

        @views @. tracers[name].data *= grid.landmask
    end

    fill_halo_tracer!(tracers,grid)

    return tracers
end
