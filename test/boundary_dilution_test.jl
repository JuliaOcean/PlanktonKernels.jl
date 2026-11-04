using PlanktonKernels.Architectures: CPU
using PlanktonKernels.Grids
using PlanktonKernels.Fields
using PlanktonKernels.Fields: BoundaryCondition, apply_bcs!, getbc

boundary_test_view(f, g) = @view f.data[g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz]

@testset "Boundary dilution" begin
    for FT in (Float32, Float64)
        @testset "$FT" begin
            @testset "Public boundary setter" begin
                g = RectilinearGrid(size=(3,4,2), x=(0,6), y=(0,12), z=(0,-8), FT=FT)
                c = init_tracers(CPU(), g, (:dye,), FT)
                # Concentration and velocity must create a usable boundary condition.
                set_bc!(c.dye, CPU(); pos=:east, bc_c=FT(4), bc_v=FT(0.5))
                @test c.dye.bc.east isa BoundaryCondition
                @test getbc(c.dye.bc.east.c, c.dye.bc.east.v, 1, 1, 1) == FT(2)
            end
            @testset "Local sampling does not mutate boundary configuration" begin
                g = RectilinearGrid(size=(3,4,2), x=(0,6), y=(0,12), z=(0,-8), FT=FT)
                c = init_tracers(CPU(), g, (:dye,), FT)
                c.dye.data .= FT(-999)
                boundary_test_view(c.dye,g) .= FT(8)
                c.dye.bc.east = BoundaryCondition(nothing, FT(0.5))
                @test getbc(c.dye.bc.east.c, c.dye.bc.east.v, c.dye.data, :east, 1, 1, 10, g) == FT(4)
                @test c.dye.bc.east.c === nothing
                boundary_test_view(c.dye,g) .= FT(12)
                @test getbc(c.dye.bc.east.c, c.dye.bc.east.v, c.dye.data, :east, 1, 1, 11, g) == FT(6)
            end
            @testset "Setter accepts local concentration mode" begin
                g = RectilinearGrid(size=(3,4,2), x=(0,6), y=(0,12), z=(0,-8), FT=FT)
                c = init_tracers(CPU(), g, (:dye,), FT)
                set_bc!(c.dye, CPU(); pos=:east, bc_c=nothing, bc_v=FT(0.5))
                @test c.dye.bc.east.c === nothing
            end
            @testset "Unset boundaries are no-ops" begin
                g = RectilinearGrid(size=(3,4,2), x=(0,6), y=(0,12), z=(0,-8), FT=FT)
                c = init_tracers(CPU(), g, (:dye,), FT)
                out = init_tracers(CPU(), g, (:dye,), FT)
                out.dye.data .= FT(7)
                apply_bcs!(out, c, g, 1, FT(0.5), CPU())
                @test all(==(FT(7)), out.dye.data)
            end
            # Test every face independently. Nonuniform interior and halo sentinels
            # distinguish boundary sampling from accidental halo or time indexing.
            for (face, axis, upper, sign, spacing) in
                ((:west,1,false,1,2), (:east,1,true,-1,2),
                 (:south,2,false,1,3), (:north,2,true,-1,3),
                 (:bottom,3,true,1,4), (:top,3,false,-1,4))
                @testset "$face local concentration" begin
                    g = RectilinearGrid(size=(3,4,2), x=(0,6), y=(0,12), z=(0,-8), FT=FT)
                    c = init_tracers(CPU(), g, (:dye,:salt), FT)
                    out = init_tracers(CPU(), g, keys(c), FT)
                    for (name, scale) in ((:dye,1), (:salt,2))
                        c[name].data .= FT(-999)
                        a = boundary_test_view(c[name], g)
                        for I in CartesianIndices(a)
                            a[I] = FT(scale * (I[1] + 10I[2] + 100I[3]))
                        end
                        # c=nothing selects the current local concentration.
                        # Choose outward velocity, so dilution decreases concentration.
                        setproperty!(c[name].bc, face, BoundaryCondition(nothing, FT(-sign*0.25)))
                        out[name].data .= FT(7)
                    end
                    before = deepcopy(c)
                    expected = deepcopy(out)
                    # Iteration deliberately exceeds Nz: time must not index tracer depth.
                    for step in 1:2
                        for name in keys(c)
                            a = boundary_test_view(expected[name], g)
                            b = boundary_test_view(c[name], g)
                            boundary = upper ? size(a, axis) : 1
                            for I in CartesianIndices(a)
                                if I[axis] == boundary
                                    a[I] -= FT(0.25) * b[I] * FT(0.5) / FT(spacing)
                                end
                            end
                        end
                        apply_bcs!(out, c, g, 10, FT(0.5), CPU())
                        for name in keys(c)
                            # Includes unchanged halo and non-boundary cells.
                            @test out[name].data ≈ expected[name].data
                            @test c[name].data == before[name].data
                        end
                    end
                end
            end
            @testset "Array boundary values through kernels" begin
                g = RectilinearGrid(size=(3,4,2), x=(0,6), y=(0,12), z=(0,-8), FT=FT)
                spatial = reshape(FT.(1:8), 4, 2)
                temporal = cat(spatial, FT(2) .* spatial; dims=3)
                for concentration in (nothing, zero(FT), FT(3), spatial, temporal),
                    velocity in (zero(FT), FT(0.25), spatial, temporal)
                    c = init_tracers(CPU(), g, (:dye,), FT)
                    out = init_tracers(CPU(), g, (:dye,), FT)
                    c.dye.data .= FT(5)
                    set_bc!(c.dye, CPU(); pos=:east, bc_c=concentration, bc_v=velocity)
                    expected = copy(out.dye.data)
                    for t in 1:2
                        for k in 1:g.Nz, j in 1:g.Ny
                            cv = concentration === nothing ? FT(5) :
                                 concentration isa Number ? concentration :
                                 ndims(concentration) == 2 ? concentration[j,k] : concentration[j,k,t]
                            vv = velocity isa Number ? velocity :
                                 ndims(velocity) == 2 ? velocity[j,k] : velocity[j,k,t]
                            expected[g.Hx+g.Nx,g.Hy+j,g.Hz+k] -= vv * cv * FT(0.5) / FT(2)
                        end
                        apply_bcs!(out, c, g, t, FT(0.5), CPU())
                        @test out.dye.data ≈ expected
                        @test all(==(FT(5)), c.dye.data)
                    end
                end
            end
            @testset "Array concentration and time-dependent velocity" begin
                g = RectilinearGrid(size=(3,4,2), x=(0,6), y=(0,12), z=(0,-8), FT=FT)
                c = init_tracers(CPU(), g, (:dye,), FT)
                concentrations = reshape(FT.(1:8), 4, 2)
                velocities = cat(fill(FT(0.25),4,2),fill(FT(0.5),4,2); dims=3)
                c.dye.bc.east = BoundaryCondition(concentrations, velocities)
                @test getbc(c.dye.bc.east.c, c.dye.bc.east.v, c.dye.data, :east, 3, 2, 2, g) == concentrations[3,2] * FT(0.5)
            end
        end
    end
end
