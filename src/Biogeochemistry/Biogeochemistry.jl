module Biogeochemistry
using KernelAbstractions
using ..Architectures: Architecture, device
using ..Grids: AbstractGrid, volume
import ..Fields
using ..Fields: zero_fields!, apply_bcs!, fill_halo_tracer!
using ..Transport: tracer_advection!, tracer_diffusion!, tracer_sinking!

export bgc_tracer_names, bgc_tracer_default_init, bgc_params_default, update_bgc_params
export generate_bgc_tracers, bgc_tracer_update!

include("bgc_tracer_fields.jl")
include("bgc_params.jl")
include("bgc_tracer_forcings.jl")
include("bgc_tracer_update.jl")

end