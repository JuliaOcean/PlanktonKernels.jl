using PlanktonKernels.Architectures: CPU
import Adapt
using PlanktonKernels.Grids
using PlanktonKernels.Grids: short_show, replace_grid_storage, ΔxC, ΔyC, ΔzC, ΔxF, ΔyF, ΔzF, Ax, Ay, Az, volume



function test_rectilinear_grid()
    grid = RectilinearGrid(size = (4,6,2), x = (0,12), y = (0,12), z = (0,-8), halo = (2,2,2))

    @test grid.Nx == 4
    @test grid.Ny == 6
    @test grid.Nz == 2

    @test grid.Δx == 3.0
    @test grid.Δy == 2.0
    @test grid.dzC == [4.0,4.0,4.0,4.0,4.0,4.0]
    @test grid.dzF == [4.0,4.0,4.0,4.0,4.0,4.0]

    @test length(grid.xC) == 4+2*2
    @test length(grid.yC) == 6+2*2
    @test length(grid.zC) == 2+2*2
    @test length(grid.xF) == 4+2*2
    @test length(grid.yF) == 6+2*2
    @test length(grid.zF) == 2+2*2

    return nothing
end

function test_vertically_stretched_rectilinear_grid()
    zfs = [0.0, -1.0, -3.0, -5.0, -10.0]
    grid = RectilinearGrid(size = (4,6,4), x = (0,12), y = (0,12), z = zfs, halo = (2,2,2))

    @test grid.Nx == 4
    @test grid.Ny == 6
    @test grid.Nz == 4

    @test grid.Δx == 3.0
    @test grid.Δy == 2.0
    @test grid.dzC == [1.0,1.0,1.5,2.0,3.5,5.0,5.0,5.0]
    @test grid.dzF == [1.0,1.0,1.0,2.0,2.0,5.0,5.0,5.0]

    @test grid.zF == [2.0,1.0,-0.0,-1.0,-3.0,-5.0,-10.0,-15.0]
    @test grid.zC == [1.5,0.5,-0.5,-2.0,-4.0,-7.5,-12.5,-17.5]

    @test length(grid.xC) == 4+2*2
    @test length(grid.yC) == 6+2*2
    @test length(grid.zC) == 4+2*2
    @test length(grid.xF) == 4+2*2
    @test length(grid.yF) == 6+2*2
    @test length(grid.zF) == 4+2*2

    return nothing
end

function test_rectilinear_areas_volumes()
    zfs = [0.0, -1.0, -3.0, -5.0, -10.0]
    grid = RectilinearGrid(size = (4,6,4), x = (0,12), y = (0,12), z = zfs, halo = (2,2,2))
    @test ΔxC(1,1,1,grid) == grid.Δx
    @test ΔyC(1,1,1,grid) == grid.Δy
    @test ΔzC(1,1,1,grid) == grid.dzC[1]
    @test ΔxF(1,1,1,grid) == grid.Δx
    @test ΔyF(1,1,1,grid) == grid.Δy
    @test ΔzF(1,1,1,grid) == grid.dzF[1]
    @test Ax(1,1,1,grid) == grid.Δy*grid.dzF[1]
    @test Ay(1,1,1,grid) == grid.Δx*grid.dzF[1]
    @test Az(1,1,1,grid) == grid.Δx*grid.Δy
    @test volume(1,1,1,grid) == grid.Δx*grid.Δy*grid.dzF[1]

    return nothing
end

function test_lat_lon_grid()
    grid = LatLonGrid(size = (360,160,10), lat = (-80,80), lon = (-180,180), z = (0,-20))
    @test grid.Nx == 360
    @test grid.Ny == 160
    @test grid.Nz == 10

    @test length(grid.xC) == 360+2*2
    @test length(grid.yC) == 160+2*2
    @test length(grid.zC) == 10+2*2
    @test length(grid.xF) == 360+2*2
    @test length(grid.yF) == 160+2*2
    @test length(grid.zF) == 10+2*2
end

function test_load_lat_lon_grid()
    dims = (4, 4, 2)
    grid_info = (RF=[0, -2, -6], RC=[-1, -4], DRF=[2, 4], DRC=[1, 3],
                 DXC=fill(3, 4, 4), DXG=fill(3, 4, 4),
                 DYC=fill(5, 4, 4), DYG=fill(5, 4, 4), RAC=fill(15, 4, 4),
                 hFacW=fill(0.5, dims), hFacS=fill(0.25, dims), hFacC=fill(0.75, dims))
    for lon in ((-180, 180), (-20, 20))
        grid = LoadLatLonGrid(; grid_info, size=dims, lat=(-20,20), lon)
        @test (grid.Nx, grid.Ny, grid.Nz) == dims
        @test grid.dzF == [2, 2, 2, 4, 4, 4]
        @test grid.dzC == [2, 2, 2, 3, 3, 3]
        @test all(grid.Ax[:, :, 3] .== 5)
        @test all(grid.Ay[:, :, 3] .== 1.5)
        @test all(grid.Vol[:, :, 3] .== 22.5)
        @test all(grid.Vol[:, :, 4] .== 45)
    end

    return nothing
end

function test_lat_lon_areas_volumes()
    grid = LatLonGrid(size = (360,160,10), lat = (-80,80), lon = (-180,180), z = (0,-20))
    @test ΔxC(1,1,1,grid) == grid.dxC[1,1]
    @test ΔyC(1,1,1,grid) == grid.dyC[1,1]
    @test ΔzC(1,1,1,grid) == grid.dzC[1]
    @test ΔxF(1,1,1,grid) == grid.dxF[1,1]
    @test ΔyF(1,1,1,grid) == grid.dyF[1,1]
    @test ΔzF(1,1,1,grid) == grid.dzF[1]
    @test Ax(1,1,1,grid) == grid.Ax[1,1,1]
    @test Ay(1,1,1,grid) == grid.Ay[1,1,1]
    @test Az(1,1,1,grid) == grid.Az[1,1]
    @test volume(1,1,1,grid) == grid.Vol[1,1,1]

    return nothing
end

function test_grid_storage_and_masks()
    mask = ones(Bool, 4, 4, 2)
    mask[1, 1, 1] = false
    for FT in (Float32, Float64), T in (Periodic, Bounded)
        grids = (RectilinearGrid(size=(4,4,2), x=(0,12), y=(0,20), z=(0,-6), FT=FT,
                                 topology=(T,T,Bounded), landmask=mask),
                 LatLonGrid(size=(4,4,2), lat=(-20,20), lon=(-20,20), z=(0,-6), FT=FT))
        for g in grids
            @test eltype(g.xF) === FT
            for converted in (replace_grid_storage(CPU(), g), Adapt.adapt(Array, g))
                @test typeof(converted) === typeof(g)
                for f in fieldnames(typeof(g))
                    @test getfield(converted, f) == getfield(g, f)
                end
            end
            @test !isempty(short_show(g))
        end
        g = first(grids)
        @test g.landmask[3:6,3:6,3:4] == mask
        @test g.landmask[1:2,:,:] == (T === Periodic ? g.landmask[5:6,:,:] : repeat(g.landmask[3:3,:,:],2,1,1))
        @test g.landmask[:,1:2,:] == (T === Periodic ? g.landmask[:,5:6,:] : repeat(g.landmask[:,3:3,:],1,2,1))
    end
    @test_throws ArgumentError RectilinearGrid(size=(4,4,2), x=(0,1), y=(0,1), z=(0,-1), landmask=zeros(4,4,2))
    @test_throws ArgumentError RectilinearGrid(size=(4,4,2), x=(0,1), y=(0,1), z=(0,-1), landmask=ones(Bool,2,2,2))

    return nothing
end

@testset "Grids" begin
    @testset "Rectilinear Grid" begin
        test_rectilinear_grid()
        test_vertically_stretched_rectilinear_grid()
        test_rectilinear_areas_volumes()
        grid = RectilinearGrid(size = (4,6,2), x = (0,12), y = (0,12), z = (0,-8), halo = (2,2,2))
        @test try
            show(grid); println()
            true
        catch err
            println("error in show(::RectilinearGrid)")
            println(sprint(showerror, err))
            false
        end
        @test grid isa RectilinearGrid
    end
    @testset "Latitude Longitude Grid" begin
        test_lat_lon_grid()
        test_load_lat_lon_grid()
        test_lat_lon_areas_volumes()
        grid = LatLonGrid(size = (360,160,10), lat = (-80,80), lon = (-180,180), z = (0,-20))
        @test try
            show(grid); println()
            true
        catch err
            println("error in show(::LatLonGrid)")
            println(sprint(showerror, err))
            false
        end
        @test grid isa LatLonGrid
    end
    @testset "Storage and masks" begin
        test_grid_storage_and_masks()
    end
end
