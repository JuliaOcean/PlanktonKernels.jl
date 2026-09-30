using PlanktonKernels.Architectures
using PlanktonKernels.Architectures: array_type, rng_type, device, isfunctional, unsafe_free!
using Random: MersenneTwister
import KernelAbstractions

function test_cpu_architecture()
    @test CPU() isa Architecture
    @test GPU() isa Architecture
    @test device(CPU()) isa KernelAbstractions.CPU
    @test array_type(CPU()) === Array
    @test rng_type(CPU()) isa MersenneTwister
    @test isfunctional(CPU())
    @test unsafe_free!(zeros(2)) === nothing

    return nothing
end

@testset "Architectures" begin
    @testset "CPU and GPU types" begin
        test_cpu_architecture()
    end
end
