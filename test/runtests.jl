using Test
using PlanktonKernels

@testset "PlanktonKernels" begin
    include("architecture_test.jl")
    include("unit_test.jl")
    include("grid_test.jl")
    include("field_test.jl")
    include("field_copy_test.jl")
    include("transport_test.jl")
    include("biogeochemistry_test.jl")
end
