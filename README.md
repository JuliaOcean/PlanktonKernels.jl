# PlanktonKernels.jl

**PlanktonKernels.jl** provides shared computational infrastructure for plankton ecosystem modeling in Julia. It is designed as a lightweight, reusable foundation for different plankton modeling frameworks, including both individual-based and continuum approaches.

The package contains common numerical and modeling components such as computational architectures, spatial grids, fields, boundary conditions, tracer transport, and biogeochemical utilities. By separating these general-purpose components from model-specific implementations, PlanktonKernels.jl allows different plankton models to share the same numerical foundation while retaining their own representations of plankton populations and physiological processes.

PlanktonKernels.jl is developed as the common backend for **PlanktonIndividuals.jl**, an individual-based plankton model, and **PlanktonContinuum.jl**, a continuum-based plankton modeling framework. This shared infrastructure facilitates consistent implementations of grids, environmental tracers, transport processes, and biogeochemical coupling across different modeling approaches, making it easier to develop, compare, and extend plankton ecosystem models.
