using PlanktonKernels.Units

function test_unit_constants()
    @test second === seconds === 1
    @test minute === minutes === 60seconds
    @test hour === hours === 60minutes
    @test meter === meters === 1.0f0
    @test kilometer === kilometers === 1000meters
    @test (KiB, MiB, GiB, TiB) == Tuple(Float32(1024)^i for i in 1:4)

    return nothing
end

@testset "Units" begin
    @testset "Conversion constants" begin
        test_unit_constants()
    end
end
