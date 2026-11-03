#####
##### apply flux boundary conditions to tracer tendency `Gc`
##### 

# positive west/south/bottom bc values are associated with an *increase* in `Gc` near the boundary.
@kernel function apply_west_bcs_kernel!(Gc, grid, tr, c, v, iter, ΔT)
    j, k = @index(Global, NTuple)
    jj = j + grid.Hy
    kk = k + grid.Hz
    @inbounds Gc[1+grid.Hx, jj, kk] += getbc(c, v, tr, :west, j, k, iter, grid) * ΔT * Ax(1+grid.Hx, jj, kk, grid) / volume(1+grid.Hx, jj, kk, grid)
end

@kernel function apply_south_bcs_kernel!(Gc, grid, tr, c, v, iter, ΔT)
    i, k = @index(Global, NTuple)
    ii = i + grid.Hx
    kk = k + grid.Hz
    @inbounds Gc[ii, 1+grid.Hy, kk] += getbc(c, v, tr, :south, i, k, iter, grid) * ΔT * Ay(ii, 1+grid.Hy, kk, grid) / volume(ii, 1+grid.Hy, kk, grid)
end

@kernel function apply_bottom_bcs_kernel!(Gc, grid, tr, c, v, iter, ΔT)
    i, j = @index(Global, NTuple)
    ii = i + grid.Hx
    jj = j + grid.Hy
    @inbounds Gc[ii, jj, grid.Nz+grid.Hz] += getbc(c, v, tr, :bottom, i, j, iter, grid) * ΔT * Az(ii, jj, grid.Nz+grid.Hz, grid) / volume(ii, jj, grid.Nz+grid.Hz, grid)
end

# positive east/north/top bc values are associated with a *decrease* in `Gc` near the boundary.
@kernel function apply_east_bcs_kernel!(Gc, grid, tr, c, v, iter, ΔT)
    j, k = @index(Global, NTuple)
    jj = j + grid.Hy
    kk = k + grid.Hz
    @inbounds Gc[grid.Nx+grid.Hx, jj, kk] -= getbc(c, v, tr, :east, j, k, iter, grid) * ΔT * Ax(grid.Nx+grid.Hx, jj, kk, grid) / volume(grid.Nx+grid.Hx, jj, kk, grid)
end

@kernel function apply_north_bcs_kernel!(Gc, grid, tr, c, v, iter, ΔT)
    i, k = @index(Global, NTuple)
    ii = i + grid.Hx
    kk = k + grid.Hz
    @inbounds Gc[ii, grid.Ny+grid.Hy, kk] -= getbc(c, v, tr, :north, i, k, iter, grid) * ΔT * Ay(ii, grid.Ny+grid.Hy, kk, grid) / volume(ii, grid.Ny+grid.Hy, kk, grid)
end

@kernel function apply_top_bcs_kernel!(Gc, grid, tr, c, v, iter, ΔT)
    i, j = @index(Global, NTuple)
    ii = i + grid.Hx
    jj = j + grid.Hy
    @inbounds Gc[ii, jj, 1+grid.Hz] -= getbc(c, v, tr, :top, i, j, iter, grid) * ΔT * Az(ii, jj, 1+grid.Hz, grid) / volume(ii, jj, 1+grid.Hz, grid)
end

##### apply flux boundary conditions in x direction (west and east)
function apply_west_bcs!(Gc, grid, tr, bc, iter, ΔT, arch)
    kernel! = apply_west_bcs_kernel!(device(arch), (16,16), (grid.Ny, grid.Nz))
    kernel!(Gc, grid, tr, bc.c, bc.v, iter, ΔT)
    return nothing
end
function apply_east_bcs!(Gc, grid, tr, bc, iter, ΔT, arch)
    kernel! = apply_east_bcs_kernel!(device(arch), (16,16), (grid.Ny, grid.Nz))
    kernel!(Gc, grid, tr, bc.c, bc.v, iter, ΔT)
    return nothing
end

##### apply flux boundary conditions in y direction (south and north)
function apply_south_bcs!(Gc, grid, tr, bc, iter, ΔT, arch)
    kernel! = apply_south_bcs_kernel!(device(arch), (16,16), (grid.Nx, grid.Nz))
    kernel!(Gc, grid, tr, bc.c, bc.v, iter, ΔT)
    return nothing
end
function apply_north_bcs!(Gc, grid, tr, bc, iter, ΔT, arch)
    kernel! = apply_north_bcs_kernel!(device(arch), (16,16), (grid.Nx, grid.Nz))
    kernel!(Gc, grid, tr, bc.c, bc.v, iter, ΔT)
    return nothing
end

##### apply flux boundary conditions in z direction (bottom and top)
function apply_bottom_bcs!(Gc, grid, tr, bc, iter, ΔT, arch)
    kernel! = apply_bottom_bcs_kernel!(device(arch), (16,16), (grid.Nx, grid.Ny))
    kernel!(Gc, grid, tr, bc.c, bc.v, iter, ΔT)
    return nothing
end
function apply_top_bcs!(Gc, grid, tr, bc, iter, ΔT, arch)
    kernel! = apply_top_bcs_kernel!(device(arch), (16,16), (grid.Nx, grid.Ny))
    kernel!(Gc, grid, tr, bc.c, bc.v, iter, ΔT)
    return nothing
end


##### apply flux boundary conditions in all three directions, six faces.
function apply_bcs!(Gcs, tracers, grid, iter, ΔT, arch)
    for name in keys(tracers)
          apply_west_bcs!(Gcs[name].data, grid, tracers[name].data, tracers[name].bc.west, iter, ΔT, arch)
          apply_east_bcs!(Gcs[name].data, grid, tracers[name].data, tracers[name].bc.east, iter, ΔT, arch)
         apply_south_bcs!(Gcs[name].data, grid, tracers[name].data, tracers[name].bc.south, iter, ΔT, arch)
         apply_north_bcs!(Gcs[name].data, grid, tracers[name].data, tracers[name].bc.north, iter, ΔT, arch)
        apply_bottom_bcs!(Gcs[name].data, grid, tracers[name].data, tracers[name].bc.bottom, iter, ΔT, arch)
           apply_top_bcs!(Gcs[name].data, grid, tracers[name].data, tracers[name].bc.top, iter, ΔT, arch)
    end

    return nothing
end

# Avoid some computation / memory accesses when no flux boundary conditions will be applied
  apply_west_bcs!(Gc, grid, tr, ::Nothing, iter, ΔT, arch) = nothing
  apply_east_bcs!(Gc, grid, tr, ::Nothing, iter, ΔT, arch) = nothing
 apply_north_bcs!(Gc, grid, tr, ::Nothing, iter, ΔT, arch) = nothing
 apply_south_bcs!(Gc, grid, tr, ::Nothing, iter, ΔT, arch) = nothing
apply_bottom_bcs!(Gc, grid, tr, ::Nothing, iter, ΔT, arch) = nothing
   apply_top_bcs!(Gc, grid, tr, ::Nothing, iter, ΔT, arch) = nothing
