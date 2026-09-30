module Transport

export tracer_advection!
export tracer_diffusion!
export tracer_sinking!

using KernelAbstractions

using ..Grids
using ..Grids: ΔxC, ΔyC, ΔzC, ΔxF, ΔyF, ΔzF, Ax, Ay, Az, volume
using ..Architectures: device, array_type, Architecture
using ..Fields: fill_halo_Gcs!, fill_halo_tracer!, fill_halo_flux_sink!

include("operaters.jl")
include("tracer_diffusion.jl")
include("DST3FL.jl")
include("multi_dim_adv.jl")
include("tracer_sinking.jl")

end

