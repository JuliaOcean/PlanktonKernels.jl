##### calculate the sinking flux for tracer c
@inline function calc_sinking_flux!(flux_sink, grid::AbstractGrid, w::Number, c)
    flux_sink .= w .* c
    return nothing
end

##### calculate the sinking tendency for tracer c
@kernel function calc_sinking_kernel!(Gc, flux_sink, grid::AbstractGrid, ΔT)
    i, j, k = @index(Global, NTuple)
    ii = i + grid.Hx
    jj = j + grid.Hy
    kk = k + grid.Hz
    @inbounds Gc[ii, jj, kk] = Gc[ii, jj, kk] + (flux_sink[ii, jj, kk-1] - flux_sink[ii, jj, kk]) / grid.dzC[k] * ΔT
end

function calc_sinking!(Gc, flux_sink, grid::AbstractGrid, ΔT, arch::Architecture)
    kernel! = calc_sinking_kernel!(device(arch), (16,16), (grid.Nx, grid.Ny, grid.Nz))
    kernel!(Gc, flux_sink, grid, ΔT)
    return nothing
end

function tracer_sinking!(Gcs, flux_sink, arch::Architecture, g::AbstractGrid, tracers, w_sink::Number, ΔT, tr_names::Tuple)
    
    for name in tr_names
        calc_sinking_flux!(flux_sink, g, w_sink, tracers[name].data)
        fill_halo_flux_sink!(flux_sink, g)
        calc_sinking!(Gcs[name].data, flux_sink, g, ΔT, arch)
    end
    return nothing
end