# PlanktonKernels.jl

**PlanktonKernels.jl** provides shared computational infrastructure for plankton ecosystem modeling in Julia. It is designed as a lightweight, reusable foundation for different plankton modeling frameworks, including both individual-based and continuum approaches.

The package contains common numerical and modeling components such as computational architectures, spatial grids, fields, boundary conditions, tracer transport, and biogeochemical utilities. By separating these general-purpose components from model-specific implementations, PlanktonKernels.jl allows different plankton models to share the same numerical foundation while retaining their own representations of plankton populations and physiological processes.

PlanktonKernels.jl is developed as the common backend for **PlanktonIndividuals.jl**, an individual-based plankton model, and **PlanktonContinuum.jl**, a continuum-based plankton modeling framework. This shared infrastructure facilitates consistent implementations of grids, environmental tracers, transport processes, and biogeochemical coupling across different modeling approaches, making it easier to develop, compare, and extend plankton ecosystem models.
## Biogeochemistry

`Biogeochemistry` owns the standard tracer names, initialization defaults, and
reaction formulas. It uses `Fields` for allocation and boundary conditions and
`Transport` for advection, diffusion, and sinking. The 15 tracers are DIC, NH4,
NO3, PO4, DFe, O2, DOC, DON, DOP, PFe_inorg, POC, PON, POP, PFe_bio, and Dust.

```julia
import PlanktonKernels.Biogeochemistry as BGC
tracers = BGC.generate_tracers(CPU(), grid) # original concentrations and random noise
Gcs = BGC.tracers_init(CPU(), grid)
temp = BGC.tracers_init(CPU(), grid)
consume = BGC.tracers_init(CPU(), grid) # external per-cell amounts per step
vel = PlanktonKernels.Fields.tracers_init(CPU(), grid, (:u, :v, :w))
flux_sink = similar(tracers.DIC.data)
params = BGC.bgc_params_default()
BGC.tracer_update!(tracers, Gcs, temp, flux_sink, CPU(), grid, params, vel, consume, 1.0f0, 1)
```

Use qualified `BGC.generate_tracers` / `BGC.tracers_init` for the standard tracer
set; the generic Fields functions keep their existing APIs. Reaction and tracer
transport parameter defaults match PlanktonIndividuals; light attenuation and
individual-particle parameters are not included. The source update order and
reaction formulas are preserved. `consume` is an amount increment, not a rate,
and is not reset by the update. Default sinking uses bounded vertical boundaries;
the previously documented stretched-grid sinking limitation still applies.

## Copying external fields

`Fields.copy_interior!(destination, source, grid)` copies cell-centered data
into a halo-sized array without changing its halos.
`Fields.vel_copy!(vel, u, v, w, grid)` copies external velocity arrays into a
NamedTuple of fields and fills their halos. For a bounded x, y, or z dimension,
the corresponding u, v, or w input includes one extra face along that dimension.
The lower-level `copy_interior_u!`, `copy_interior_v!`, and `copy_interior_w!`
accept destination/source arrays, a grid, and a `Periodic()` or `Bounded()`
topology instance; they do not fill halos.
