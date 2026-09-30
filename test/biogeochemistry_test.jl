using PlanktonKernels.Architectures: CPU
using PlanktonKernels.Grids: RectilinearGrid
using PlanktonKernels.Fields: Field
import PlanktonKernels.Biogeochemistry as BGC
using PlanktonKernels.Fields: fill_halo_tracer!, apply_bcs!, zero_fields!
using PlanktonKernels.Grids: volume
using PlanktonKernels.Transport: tracer_advection!, tracer_diffusion!, tracer_sinking!

function test_biogeochemistry_integration()
    expected = (:DIC,:NH4,:NO3,:PO4,:DFe,:O2,:DOC,:DON,:DOP,:PFe_inorg,:POC,:PON,:POP,:PFe_bio,:Dust)
    @test BGC.bgc_tracer_names == expected
    defaults = BGC.bgc_tracer_init()
    @test Tuple(defaults.initial_condition) == (20.0,0.5,0.8,0.10,1e-6,160.0,10.0,0.1,0.05,0.0,0.0,0.0,0.0,0.0,0.0)
    @test all(==(0.1), defaults.rand_noise)
    g = RectilinearGrid(size=(4,4,4),x=(0,8),y=(0,8),z=(0,-8))
    source = (initial_condition=defaults.initial_condition,
              rand_noise=NamedTuple{expected}(ntuple(_->0.0,15)))
    @test_throws ArgumentError BGC.generate_bgc_tracers(CPU(),g,(initial_condition=(dye=1,),rand_noise=(dye=0,)))
    for FT in (Float32,Float64)
        c=BGC.generate_bgc_tracers(CPU(),g,source,FT)
        @test keys(c) == expected
        for n in expected
            @test c[n] isa Field{FT}
            @test all(==(FT(source.initial_condition[n])), c[n].data)
        end
        G=PlanktonKernels.Fields.tracers_init(CPU(),g,expected,FT); tmp=PlanktonKernels.Fields.tracers_init(CPU(),g,expected,FT)
        consume=PlanktonKernels.Fields.tracers_init(CPU(),g,expected,FT)
        vel=PlanktonKernels.Fields.tracers_init(CPU(),g,(:u,:v,:w),FT)
        flux=similar(c.DIC.data)
        params=BGC.bgc_params_default(FT)
        @test valtype(typeof(params)) === FT
        @test params["kDOC"] == FT(1/30/86400)
        overrides = Dict("κh"=>0.125, "kDOC"=>0.25)
        updated = BGC.update_bgc_params(overrides,FT)
        @test updated["κh"] === FT(0.125)
        @test updated["kDOC"] === FT(0.25)
        @test updated["Nit"] === params["Nit"]
        @test params["κh"] == 0
        @test overrides == Dict("κh"=>0.125,"kDOC"=>0.25)
        @test BGC.update_bgc_params(Dict(),FT) == params
        @test_throws ArgumentError BGC.update_bgc_params(Dict("unknown"=>1),FT)
        # Exercise all reaction pathways and verify elemental conservation.
        for n in expected; c[n].data .= FT(0.01); end
        BGC.bgc_tracer_forcing!(G,tmp,c,params,FT(60))
        @test G.DIC.data[3,3,3] ≈ FT(0.01)*params["kDOC"]*60
        @test G.NO3.data[3,3,3] ≈ FT(0.01)*params["Nit"]*60
        for group in ((:DIC,:DOC,:POC),(:NH4,:NO3,:DON,:PON),(:PO4,:DOP,:POP),(:DFe,:PFe_bio,:PFe_inorg))
            residual=sum(G[n].data[3,3,3] for n in group)
            @test abs(residual) <= 100eps(FT)
        end
        @test all(iszero,G.O2.data)
        @test all(iszero,G.Dust.data)
        # Full step must match composing Fields and Transport plus reaction increments.
        c=BGC.generate_bgc_tracers(CPU(),g,source,FT)
        c.DIC.data[4,4,4] += 1
        fill_halo_tracer!(c,g)
        vel.u.data .= FT(0.1)
        params["κh"]=FT(0.01)
        consume.DIC.data[4,4,4]=FT(8)
        PlanktonKernels.Fields.set_bc!(c.DIC,CPU();pos=:top,bc_value=FT(0.1))
        reference=deepcopy(c); dt=FT(0.5)
        tracer_advection!(reference,tmp,G,vel,g,dt,CPU())
        tracer_diffusion!(G,CPU(),g,reference,params["κh"],params["κh"],params["κv"],dt)
        tracer_sinking!(G,flux,CPU(),g,reference,params,dt)
        BGC.bgc_tracer_forcing!(G,tmp,reference,params,dt)
        apply_bcs!(G,reference,g,1,dt,CPU())
        for n in expected, k in 3:6, j in 3:6, i in 3:6
            reference[n].data[i,j,k] += G[n].data[i,j,k]+consume[n].data[i,j,k]/volume(i,j,k,g)
        end
        fill_halo_tracer!(reference,g)
        BGC.bgc_tracer_update!(c,G,tmp,flux,CPU(),g,params,vel,consume,dt,1)
        for n in expected; @test c[n].data ≈ reference[n].data; end
        # Zero transport/reactions leaves only per-cell coupling, applied once.
        params=Dict(k=>zero(v) for (k,v) in params)
        zero_fields!(vel); zero_fields!(consume)
        c=BGC.generate_bgc_tracers(CPU(),g,source,FT)
        before=copy(c.DIC.data); consume.DIC.data[4,4,4]=FT(8)
        BGC.bgc_tracer_update!(c,G,tmp,flux,CPU(),g,params,vel,consume,FT(9),1)
        @test c.DIC.data[4,4,4] == before[4,4,4]+1
        @test c.DIC.data[5,5,5] == before[5,5,5]
    end

    return nothing
end

@testset "Biogeochemistry" begin
    @testset "Defaults, reactions, and tracer updates" begin
        test_biogeochemistry_integration()
    end
end
