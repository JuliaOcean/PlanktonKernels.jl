using Serialization: serialize
using PlanktonKernels.Architectures: CPU
using PlanktonKernels.Grids
using PlanktonKernels.Transport
using PlanktonKernels.Fields
using PlanktonKernels.Fields: zero_fields!, init_tracers, validate_bc, validate_bcs, apply_bcs!
const field_test_names = (:DIC, :dye)
field_test_init() = (DIC=(init=20.0, rand_noise=0.1), dye=(init=1.0, rand_noise=0.0))

function test_fields()
    grid = RectilinearGrid(size = (4,6,2), x = (0,12), y = (0,12), z = (0,-8))

    tracers = init_tracers(CPU(), grid, field_test_names)

    @test Tuple(collect(keys(tracers))) == field_test_names
    @test tracers.DIC.data == zeros(8,10,6)
    @test interior(tracers.DIC.data, grid) == zeros(4,6,2) 
    
    tracers = generate_tracers(CPU(), grid, field_test_names, field_test_init(), Float32)
    @test maximum(tracers.DIC.data) < 23.0
    @test minimum(tracers.DIC.data) > 17.0

    zero_fields!(tracers)
    @test tracers.DIC.data == zeros(8,10,6)

    return nothing
end
function test_fill_halos()
    grid = RectilinearGrid(size = (4,6,2), x = (0,12), y = (0,12), z = (0,-8))
    tracers = generate_tracers(CPU(), grid, field_test_names, field_test_init(), Float32)
    Nx,Ny,Nz = grid.Nx, grid.Ny, grid.Nz
    Hx,Hy,Hz = grid.Hx, grid.Hy, grid.Hz

    PlanktonKernels.Fields.fill_halo_east!(tracers.DIC.data, Hx, Nx, Periodic())
    @test tracers.DIC.data[Nx+Hx+1:Nx+2*Hx, :, :] == tracers.DIC.data[1+Hx:2*Hx, :, :]

    PlanktonKernels.Fields.fill_halo_west!(tracers.DIC.data, Hx, Nx, Periodic())
    @test tracers.DIC.data[1:Hx, :, :] == tracers.DIC.data[Nx+1:Nx+Hx, :, :]

    PlanktonKernels.Fields.fill_halo_north!(tracers.DIC.data, Hy, Ny, Periodic())
    @test tracers.DIC.data[:, Ny+Hy+1:Ny+2*Hy, :] == tracers.DIC.data[:, 1+Hy:2*Hy, :]

    PlanktonKernels.Fields.fill_halo_south!(tracers.DIC.data, Hy, Ny, Periodic())
    @test tracers.DIC.data[:, 1:Hy, :] == tracers.DIC.data[:, Ny+1:Ny+Hy, :]

    PlanktonKernels.Fields.fill_halo_bottom!(tracers.DIC.data, Hz, Nz, Periodic())
    @test tracers.DIC.data[:, :, Nz+Hz+1:Nz+2*Hz] == tracers.DIC.data[:, :, 1+Hz:2*Hz]

    PlanktonKernels.Fields.fill_halo_top!(tracers.DIC.data, Hz, Nz, Periodic())
    @test tracers.DIC.data[:, :, 1:Hz] == tracers.DIC.data[:, :, Nz+1:Nz+Hz]

    PlanktonKernels.Fields.fill_halo_east!(tracers.DIC.data, Hx, Nx, Bounded())
    for i in 1:Hx
        @test tracers.DIC.data[Nx+Hx+i, :, :] == tracers.DIC.data[Nx+Hx, :, :]
    end

    PlanktonKernels.Fields.fill_halo_west!(tracers.DIC.data, Hx, Nx, Bounded())
    for i in 1:Hx
        @test tracers.DIC.data[i, :, :] == tracers.DIC.data[Hx+1, :, :]
    end

    PlanktonKernels.Fields.fill_halo_north!(tracers.DIC.data, Hy, Ny, Bounded())
    for i in 1:Hy
        @test tracers.DIC.data[:, Ny+Hy+i, :] == tracers.DIC.data[:, Ny+Hy, :]
    end

    PlanktonKernels.Fields.fill_halo_south!(tracers.DIC.data, Hy, Ny, Bounded())
    for i in 1:Hy
        @test tracers.DIC.data[:, i, :] == tracers.DIC.data[:, Hy+1, :]
    end

    PlanktonKernels.Fields.fill_halo_bottom!(tracers.DIC.data, Hz, Nz, Bounded())
    for i in 1:Hz
        @test tracers.DIC.data[:, :, Nz+Hz+i] == tracers.DIC.data[:, :, Nz+Hz]
    end

    PlanktonKernels.Fields.fill_halo_top!(tracers.DIC.data, Hz, Nz, Bounded())
    for i in 1:Hz
        @test tracers.DIC.data[:, :, i] == tracers.DIC.data[:, :, Hz+1]
    end

    PlanktonKernels.Fields.fill_halo_east_vel!(tracers.DIC.data, Hx, Nx, Bounded())
    for i in 1:Hx-1
        @test tracers.DIC.data[Nx+Hx+1+i, :, :] == tracers.DIC.data[Nx+Hx+1, :, :]
    end

    PlanktonKernels.Fields.fill_halo_north_vel!(tracers.DIC.data, Hy, Ny, Bounded())
    for i in 1:Hy-1
        @test tracers.DIC.data[:, Ny+Hx+1+i, :] == tracers.DIC.data[:, Ny+Hy+1, :]
    end

    PlanktonKernels.Fields.fill_halo_bottom_vel!(tracers.DIC.data, Hz, Nz, Bounded())
    for i in 1:Hz-1
        @test tracers.DIC.data[:, :, Nz+Hx+1+i] == tracers.DIC.data[:, :, Nz+Hz+1]
    end

    PlanktonKernels.Fields.fill_halo_east_Gc!(tracers.DIC.data, Hx, Nx, Bounded())
    for i in 1:Hx
        @test tracers.DIC.data[Nx+Hx+i, :, :] == zeros(10,6)
    end

    PlanktonKernels.Fields.fill_halo_north_Gc!(tracers.DIC.data, Hy, Ny, Bounded())
    for i in 1:Hy
        @test tracers.DIC.data[:, Ny+Hx+i, :] == zeros(8,6)
    end

    PlanktonKernels.Fields.fill_halo_bottom_Gc!(tracers.DIC.data, Hz, Nz, Bounded())
    for i in 1:Hz
        @test tracers.DIC.data[:, :, Nz+Hx+i] == zeros(8,10)
    end

    return nothing
end

function test_boundary_conditions()
    grid = RectilinearGrid(size = (4, 4, 4), x = (0,32), y = (0,32), z = (0,-32))
    tracers=init_tracers(CPU(), grid, field_test_names)
    FT = tracers.DIC.data |> eltype
    set_bc!(tracers.DIC,CPU(); pos = :west, bc_value = 0.1)
    @test tracers.DIC.bc.west == FT(0.1)
    set_bc!(tracers.DIC,CPU(); pos = :west, bc_value = ones(4,4))
    @test tracers.DIC.bc.west == ones(FT, 4,4)
    set_bc!(tracers.DIC,CPU(); pos = :west, bc_value = ones(4,4,10))
    @test tracers.DIC.bc.west == ones(FT, 4,4,10)

    Gcs = init_tracers(CPU(), grid, field_test_names)
    apply_bcs!(Gcs, tracers, grid, 10, 1, CPU())
    @test Gcs.DIC.data[3,3:6,3:6] == ones(FT, 4,4) ./ FT(8.0)

end

function test_generic_fields_and_fluxes()
    g = RectilinearGrid(size=(4,4,4), x=(0,8), y=(0,8), z=(0,-8))
    for FT in (Float32,Float64)
        c = init_tracers(CPU(), g, (:dye,:salt), FT)
        out = init_tracers(CPU(), g, keys(c), FT)
        @test eltype(c.dye.data) === FT
        @test c.dye.data !== c.salt.data
        @test c.dye.bc !== c.salt.bc
        @test interior(c.dye,g) == zeros(FT,4,4,4)
        @test all(isnothing, (getfield(c.dye.bc,f) for f in fieldnames(BoundaryConditions)))
        for (pos, index, sign) in ((:west,(3,4,4),1),(:east,(6,4,4),-1),
                                  (:south,(4,3,4),1),(:north,(4,6,4),-1),
                                  (:bottom,(4,4,6),1),(:top,(4,4,3),-1))
            for value in (2,fill(2,4,4),cat(fill(1,4,4),fill(2,4,4);dims=3))
                zero_fields!(out)
                set_bc!(c.dye,CPU();pos,bc_value=value)
                @test validate_bcs(c,g,2) === nothing
                apply_bcs!(out,c,g,2,FT(0.5),CPU())
                @test out.dye.data[index...] ≈ FT(sign*0.5)
                @test sum(interior(out.dye,g)) ≈ FT(sign*8)
                @test all(iszero,out.salt.data)
                setproperty!(c.dye.bc,pos,nothing)
            end
        end
        @test_throws ArgumentError set_bc!(c.dye,CPU();pos=:invalid,bc_value=1)
        @test_throws ArgumentError validate_bc(ones(3,4),(4,4),2)
        @test_throws ArgumentError validate_bc(ones(4,4,3),(4,4),2)
        @test_throws ArgumentError validate_bc(ones(4),(4,4),2)
        # Actual Field objects work with the migrated transport entry points.
        c.dye.data .= 1
        c.salt.data .= 2
        temp = init_tracers(CPU(),g,keys(c),FT)
        vel = init_tracers(CPU(),g,(:u,:v,:w),FT)
        zero_fields!(out)
        tracer_advection!(c,temp,out,vel,g,FT(0.1),CPU())
        tracer_diffusion!(out,CPU(),g,c,FT(0.1),FT(0.1),FT(0.1),FT(0.1))
        @test all(iszero,interior(out.dye,g))
        @test all(iszero,interior(out.salt,g))
    end
    @test_throws ArgumentError init_tracers(CPU(),g,(:dye,:dye))
    @test_throws ArgumentError generate_tracers(CPU(),g,(:dye,),(dye=(init=-1,rand_noise=0),))
    mktempdir() do dir
        path = joinpath(dir,"dye.bin")
        serialize(path,fill(3.0,4,4,4))
        c = generate_tracers(CPU(),g,(:dye,),(dye=path,),Float64)
        @test keys(c) == (:dye,)
        @test all(==(3.0),c.dye.data)
        serialize(path,zeros(2,2,2))
        @test_throws ArgumentError generate_tracers(CPU(),g,(:dye,),(dye=path,))
        serialize(path,fill(-1.0,4,4,4))
        @test_throws ArgumentError generate_tracers(CPU(),g,(:dye,),(dye=path,))
    end
    mask=ones(Bool,4,4,4); mask[2,2,2]=false
    masked=RectilinearGrid(size=(4,4,4), x=(0,8), y=(0,8), z=(0,-8),landmask=mask)
    c=generate_tracers(CPU(),masked,(:dye,),(dye=(init=2,rand_noise=0),))
    @test interior(c.dye,masked) == 2 .* mask

    return nothing
end

function test_explicit_tracer_selection()
    grid = RectilinearGrid(size=(4,4,4), x=(0,4), y=(0,4), z=(0,-4))
    source = (a=(init=1.0,rand_noise=0.0), b=(init=2.0,rand_noise=0.0), c=(init=3.0,rand_noise=0.0))
    fields = generate_tracers(CPU(),grid,(:c,:a),source,Float64)
    @test keys(fields) == (:c,:a)
    @test all(==(3.0),fields.c.data)
    @test all(==(1.0),fields.a.data)
    @test_throws ArgumentError generate_tracers(CPU(),grid,(:missing,),source)
    mktempdir() do dir
        path=joinpath(dir,"field.bin")
        serialize(path,fill(2.0,4,4,4))
        fields=generate_tracers(CPU(),grid,(:b,:a),(a=path,b=path,ignored="unused"))
        @test keys(fields) == (:b,:a)
        @test all(==(2.0f0),fields.b.data)
        @test_throws ArgumentError generate_tracers(CPU(),grid,(:missing,),(a=path,))
    end

    return nothing
end

@testset "Fields" begin
    @testset "Allocation and initialization" begin
        test_fields()
    end
    @testset "Halo regions" begin
        test_fill_halos()
    end
    @testset "Boundary conditions" begin
        test_boundary_conditions()
    end
    @testset "Generic fields and fluxes" begin
        test_generic_fields_and_fluxes()
    end
    @testset "Tracer selection" begin
        test_explicit_tracer_selection()
    end
end
