"""
    vel_copy!(vel::NamedTuple, u, v, w, g::AbstractGrid{FT, TX, TY, TZ}) where {FT, TX, TY, TZ}
Copy external velocity arrays into the fields `vel.u`, `vel.v`, and `vel.w`,
then fill velocity halos. For bounded dimensions, the corresponding component
includes the extra boundary face (Nx+1, Ny+1, or Nz+1).
"""
function vel_copy!(vel::NamedTuple, u, v, w, g::AbstractGrid{FT, TX, TY, TZ}) where {FT, TX, TY, TZ}
    copy_interior_u!(vel.u.data, u, g, TX())
    copy_interior_v!(vel.v.data, v, g, TY())
    copy_interior_w!(vel.w.data, w, g, TZ())

    fill_halo_vel!(vel, g)
end

"""
    copy_interior!(destination, source, grid)
Copy a cell-centered array into the interior of a halo-sized destination array.
Halo cells are left unchanged.
"""
@inline function copy_interior!(c, t, g::AbstractGrid)
    copyto!(view(c, g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz), t)
end
@inline function copy_interior_u!(c, t, g::AbstractGrid, ::Periodic)
    copyto!(view(c, g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz), t)
end
@inline function copy_interior_v!(c, t, g::AbstractGrid, ::Periodic)
    copyto!(view(c, g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz), t)
end
@inline function copy_interior_w!(c, t, g::AbstractGrid, ::Periodic)
    copyto!(view(c, g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz), t)
end

@inline function copy_interior_u!(c, t, g::AbstractGrid, ::Bounded)
    copyto!(view(c, g.Hx+1:g.Hx+g.Nx+1, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz), t)
end
@inline function copy_interior_v!(c, t, g::AbstractGrid, ::Bounded)
    copyto!(view(c, g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny+1, g.Hz+1:g.Hz+g.Nz), t)
end
@inline function copy_interior_w!(c, t, g::AbstractGrid, ::Bounded)
    copyto!(view(c, g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz+1), t)
end

