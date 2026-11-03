using PlanktonKernels.Fields: BoundaryCondition, validate_bc

@testset "Boundary validation" begin
    dims, steps = (3, 4), 2
    @test validate_bc(nothing, dims, steps) === nothing
    for FT in (Float32, Float64)
        scalar = FT(0.5)
        spatial = ones(FT, dims)
        temporal = ones(FT, dims..., steps)
        # Scalars, spatial arrays, and time-dependent arrays can be mixed.
        for c in (nothing, zero(FT), scalar, spatial, temporal),
            v in (zero(FT), scalar, spatial, temporal)
            @test validate_bc(BoundaryCondition(c, v), dims, steps) === nothing
        end
        invalid = (ones(FT, 3), ones(FT, 3, 4, 2, 1),
                   ones(FT, 2, 4), ones(FT, 3, 4, 3), ones(FT, 1, 4))
        for a in invalid
            @test_throws ArgumentError validate_bc(BoundaryCondition(nothing, a), dims, steps)
            @test_throws ArgumentError validate_bc(BoundaryCondition(a, scalar), dims, steps)
            @test_throws ArgumentError validate_bc(BoundaryCondition(scalar, a), dims, steps)
            # Broadcasting must not hide a malformed input's dimensions.
            @test_throws ArgumentError validate_bc(BoundaryCondition(a, temporal), dims, steps)
            @test_throws ArgumentError validate_bc(BoundaryCondition(temporal, a), dims, steps)
        end
    end
end
