# PlanktonKernels.jl

[![CI](https://github.com/JuliaOcean/PlanktonKernels.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaOcean/PlanktonKernels.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/JuliaOcean/PlanktonKernels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaOcean/PlanktonKernels.jl)

**PlanktonKernels.jl** provides shared computational infrastructure for plankton ecosystem modeling in Julia. It is designed as a lightweight, reusable foundation for different plankton modeling frameworks, including both individual-based and continuum approaches.

The package contains common numerical and modeling components such as computational architectures, spatial grids, fields, boundary conditions, tracer transport, and biogeochemical utilities. By separating these general-purpose components from model-specific implementations, PlanktonKernels.jl allows different plankton models to share the same numerical foundation while retaining their own representations of plankton populations and physiological processes.

PlanktonKernels.jl is developed as the common backend for **PlanktonIndividuals.jl**, an individual-based plankton model, and **PlanktonContinuum.jl**, a continuum-based plankton modeling framework. This shared infrastructure facilitates consistent implementations of grids, environmental tracers, transport processes, and biogeochemical coupling across different modeling approaches, making it easier to develop, compare, and extend plankton ecosystem models.

## Installation

Requires Julia 1.11 or newer. While registration is pending, install from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/JuliaOcean/PlanktonKernels.jl")
```

Once the package is available in General, use `Pkg.add("PlanktonKernels")`.

## Quick start

Create a CPU grid with 8 × 8 × 4 cells over a 100 m × 100 m domain, extending
20 m below the surface. Horizontal boundaries are periodic; vertical boundaries
are bounded. Depth coordinates decrease from the surface.

```julia
using PlanktonKernels
using PlanktonKernels.Fields: copy_interior!

arch = CPU()
grid = RectilinearGrid(
    size = (8, 8, 4),
    x = (0, 100), y = (0, 100), z = (0, -20),
    topology = (Periodic, Periodic, Bounded),
)

# Allocate a scalar field and initialize its physical cells.
concentration = Field(arch, grid, Float32)
copy_interior!(concentration.data, fill(2.0f0, 8, 8, 4), grid)

values = interior(concentration, grid)
@assert size(values) == (8, 8, 4)
@assert all(values .== 2.0f0)
```

Fields also allocate halo cells used by numerical transport schemes.
`interior(field, grid)` returns the physical cells without halos as a copy;
use `copy_interior!` to write an array into a field.

### Biogeochemical tracers

Using the same grid, initialize the standard tracer set, including dissolved
inorganic carbon (`DIC`), nitrate (`NO3`), phosphate (`PO4`), and oxygen (`O2`):

```julia
initial = bgc_tracer_init()
# Disable the default random perturbations for reproducible uniform fields.
source = (
    initial_condition = initial.initial_condition,
    rand_noise = map(_ -> 0.0, initial.rand_noise),
)
tracers = generate_bgc_tracers(arch, grid, source, Float32)

@assert keys(tracers) == bgc_tracer_names
@assert all(interior(tracers.NO3, grid) .== 0.8f0)
```

`source.initial_condition` sets each tracer's starting concentration, and
`source.rand_noise` sets its relative perturbation amplitude. Custom tracer sets
can be created with `generate_tracers`. These examples initialize fields;
advancing them requires velocities, a timestep, and transport or biogeochemical
updates supplied by the calling model.

## Package components

- `Architectures`: CPU storage and optional CUDA/Metal GPU extensions.
- `Grids`: rectilinear and latitude–longitude grids, topology, and geometry.
- `Fields`: scalar fields, boundary conditions, halo operations, and field copying.
- `Transport`: tracer advection, diffusion, and sinking.
- `Biogeochemistry`: standard tracer initialization, parameters, and updates.
- `Units`: time, length, and storage-size constants.

Constructor documentation is available in Julia's help mode, for example
`?RectilinearGrid` and `?Field`. To run the package tests, use
`using Pkg; Pkg.test("PlanktonKernels")`.
