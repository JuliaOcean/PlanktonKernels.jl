using PlanktonKernels.Architectures: CPU
using PlanktonKernels.Grids
using PlanktonKernels.Transport
using PlanktonKernels.Fields
using PlanktonKernels.Fields: fill_halo_tracer!

function test_tracer_advection(FT)
    g = RectilinearGrid(size=(8,8,8), x=(0,8), y=(0,8), z=(0,-8),
                        topology=(Periodic,Periodic,Periodic), FT=FT)
    dims = (12,12,12)
    field(v) = (data=fill(FT(v), dims),)
    fields(v) = (dye=field(v), salt=field(2v))
    interior(a) = @view a[3:10,3:10,3:10]
    for axis in (:u, :v, :w), speed in (-0.2, 0.2)
        c, tmp, tendency = fields(3), fields(0), fields(0)
        vel = NamedTuple{(:u,:v,:w)}(Tuple(field(a == axis ? speed : 0) for a in (:u,:v,:w)))
        tracer_advection!(c, tmp, tendency, vel, g, FT(0.1), CPU())
        @test all(iszero, interior(tendency.dye.data))
        @test all(iszero, interior(tendency.salt.data))
        # A periodic pulse must move without creating or losing tracer mass.
        c.dye.data[5:6,5:6,5:6] .= 4
        fill_halo_tracer!(c, g)
        before = deepcopy(c)
        tracer_advection!(c, tmp, tendency, vel, g, FT(0.1), CPU())
        @test c == before
        @test any(!iszero, interior(tendency.dye.data))
        @test abs(sum(interior(tendency.dye.data))) < 100eps(FT)
        @test minimum(interior(tmp.dye.data)) >= 3 - 10eps(FT)
        @test maximum(interior(tmp.dye.data)) <= 4 + 10eps(FT)
    end

    return nothing
end

function test_tracer_diffusion(FT)
    g = RectilinearGrid(size=(8,8,8), x=(0,8), y=(0,8), z=(0,-8),
                        topology=(Periodic,Periodic,Periodic), FT=FT)
    dims = (12,12,12)
    field(v) = (data=fill(FT(v), dims),)
    fields(v) = (dye=field(v), salt=field(2v))
    interior(a) = @view a[3:10,3:10,3:10]
    c, tendency = fields(3), fields(0)
    tracer_diffusion!(tendency, CPU(), g, c, FT(0.1), FT(0.1), FT(0.1), FT(0.1))
    @test all(iszero, interior(tendency.dye.data))
    c.dye.data[6,6,6] += 1
    fill_halo_tracer!(c, g)
    tracer_diffusion!(tendency, CPU(), g, c, FT(0.1), FT(0.1), FT(0.1), FT(0.1))
    @test tendency.dye.data[6,6,6] ≈ FT(-0.06)
    @test tendency.dye.data[7,6,6] ≈ FT(0.01)
    @test abs(sum(interior(tendency.dye.data))) < 100eps(FT)

    return nothing
end

function test_tracer_sinking(FT)
    g = RectilinearGrid(size=(8,8,8), x=(0,8), y=(0,8), z=(0,-8),
                        topology=(Periodic,Periodic,Periodic), FT=FT)
    dims = (12,12,12)
    field(v) = (data=fill(FT(v), dims),)
    fields(v) = (dye=field(v), salt=field(2v))
    interior(a) = @view a[3:10,3:10,3:10]
    # Existing sinking API: closed surface/bottom redistribute a uniform tracer.
    bounded = RectilinearGrid(size=(8,8,8), x=(0,8), y=(0,8), z=(0,-8), FT=FT)
    c = (PFe_inorg=field(1), Dust=field(1), PFe_bio=field(1))
    out = (PFe_inorg=field(0), Dust=field(0), PFe_bio=field(0))
    tracer_sinking!(out, zeros(FT,dims), CPU(), bounded, c,
                    Dict("w_sink_inorg"=>FT(0.2), "w_sink_org"=>FT(0.1)), FT(0.1))
    for (name, rate) in ((:PFe_inorg,0.2),(:Dust,0.2),(:PFe_bio,0.1))
        @test out[name].data[3,3,3] ≈ FT(-rate*0.1)
        @test out[name].data[3,3,10] ≈ FT(rate*0.1)
        @test abs(sum(interior(out[name].data))) < 100eps(FT)
    end

    return nothing
end

@testset "Transport" begin
    @testset "Advection" begin
        for FT in (Float32, Float64)
            test_tracer_advection(FT)
        end
    end
    @testset "Diffusion" begin
        for FT in (Float32, Float64)
            test_tracer_diffusion(FT)
        end
    end
    @testset "Sinking" begin
        for FT in (Float32, Float64)
            test_tracer_sinking(FT)
        end
    end
end
