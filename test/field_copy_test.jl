using PlanktonKernels.Architectures: CPU
using PlanktonKernels.Grids
using PlanktonKernels.Fields: vel_copy!, copy_interior!, copy_interior_u!, copy_interior_v!, copy_interior_w!

function test_field_copying(FT, TX, TY, TZ)
    ns=(4,5,6); hs=(2,2,2); topology=(TX,TY,TZ)
    g=RectilinearGrid(size=ns,x=(0,4),y=(0,5),z=(0,-6),FT=FT,topology=topology)
    dims=ns .+ 2 .* hs
    source=reshape(collect(1:prod(ns)),ns)
    target=fill(FT(-1),dims)
    copy_interior!(target,source,g)
    expected=fill(FT(-1),dims)
    expected[3:6,3:7,3:8] .= source
    @test target == expected
    @test source == reshape(collect(1:prod(ns)),ns)
    vel=PlanktonKernels.Fields.init_tracers(CPU(),g,(:u,:v,:w),FT)
    inputs=ntuple(3) do axis
        sz=ntuple(d->ns[d]+(d==axis && topology[d]===Bounded),3)
        reshape(collect(1:prod(sz)),sz)
    end
    for (axis, copier) in enumerate((copy_interior_u!,copy_interior_v!,copy_interior_w!))
        dst=fill(FT(-1),dims)
        copier(dst,inputs[axis],g,topology[axis]())
        expected=fill(FT(-1),dims)
        ranges=ntuple(d->(hs[d]+1):(hs[d]+size(inputs[axis],d)),3)
        expected[ranges...] .= inputs[axis]
        @test dst == expected
    end
    vel_copy!(vel,inputs...,g)
    # Independently map every halo cell to its periodic or bounded source cell.
    for (axis, name) in enumerate((:u,:v,:w))
        for I in CartesianIndices(dims)
            J=ntuple(3) do d
                logical=I[d]-hs[d]
                topology[d]===Periodic ? mod1(logical,ns[d]) : clamp(logical,1,size(inputs[axis],d))
            end
            @test vel[name].data[I] == FT(inputs[axis][J...])
        end
    end

    return nothing
end

@testset "Field copying" begin
    @testset "Scalar and staggered velocity fields" begin
        for FT in (Float32, Float64), TX in (Periodic, Bounded), TY in (Periodic, Bounded), TZ in (Periodic, Bounded)
            test_field_copying(FT, TX, TY, TZ)
        end
    end
end
