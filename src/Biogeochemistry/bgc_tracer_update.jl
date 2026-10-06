#####
##### update tracer fields
#####

##### apply tendency to tracer field
@kernel function apply_tendency_kernel!(tracer, Gc, consume, g::AbstractGrid)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds tracer[ii, jj, kk] += Gc[ii, jj, kk] + consume[ii, jj, kk]
end
function apply_tendency!(tracers, Gcs, consume, g::AbstractGrid, arch::Architecture)
    kernel! = apply_tendency_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    for name in bgc_tracer_names
        kernel!(tracers[name].data, Gcs[name].data, consume[name].data, g)
    end
    return nothing
end
    
"""
    tracer_update!(tracers, Gcs, bgc_tracer_tmp, flux_sink, arch, grid, params, vel, consume, ΔT, iter)
Advance the standard tracers with transport, reactions, boundary fluxes, and external
coupling. Gcs and bgc_tracer_tmp are independent standard-tracer workspaces; flux_sink
is a halo-sized array. consume contains concentration increments already integrated over the
time step; callers convert particle amounts to concentrations before this call. Fill input tracer/velocity halos first and
validate boundary arrays with Fields.validate_bcs. Sinking requires bounded z.
"""
function bgc_tracer_update!(tracers, Gcs, bgc_tracer_tmp, flux_sink, arch::Architecture, g::AbstractGrid, params, vel, consume, ΔT, iter)
    ##### compute advection tendency
    tracer_advection!(tracers, bgc_tracer_tmp, Gcs, vel, g, ΔT, arch)

    ##### compute tracer diffusion,for each time step
    tracer_diffusion!(Gcs, arch, g, tracers, params["κh"], params["κh"], params["κv"], ΔT)

    ##### compute tracer sinking,for each time step
    tracer_sinking!(Gcs, flux_sink, arch, g, tracers, params, ΔT)

    ##### compute biogeochemical forcings of tracers,for each time step
    zero_fields!(bgc_tracer_tmp)
    bgc_tracer_forcing!(Gcs, bgc_tracer_tmp, tracers, params, ΔT)

    ##### apply boundary conditions
    apply_bcs!(Gcs, tracers, g, iter, ΔT, arch)

    ##### apply diffusion and forcing tendency
    apply_tendency!(tracers, Gcs, consume, g, arch)

    fill_halo_tracer!(tracers, g)
end
