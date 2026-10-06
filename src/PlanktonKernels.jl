module PlanktonKernels

include("Architectures.jl")
include("Units.jl")
include("Grids/Grids.jl")
include("Fields/Fields.jl")
include("Transport/Transport.jl")
include("Biogeochemistry/Biogeochemistry.jl")

using .Architectures
using .Units
using .Grids
using .Fields
using .Transport
using .Biogeochemistry

export
    # Architectures
    Architecture, CPU, GPU,

    # Grids
    AbstractGrid, RectilinearGrid, LatLonGrid, LoadLatLonGrid, Periodic, Bounded,

    # Fields
    Field, BoundaryConditions, interior, init_tracers, initialize_tracer!, generate_tracers, set_bc!,

    # Transport
    tracer_advection!, tracer_diffusion!, tracer_sinking!,

    # Biogeochemistry
    bgc_tracer_names, bgc_tracer_default_init, bgc_params_default, update_bgc_params, generate_bgc_tracers, bgc_tracer_update!,

    # Units
    second, minute, hour, meter, kilometer, seconds, minutes, hours, meters, kilometers, KiB, MiB, GiB, TiB

end
