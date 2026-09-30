module Architectures

export CPU, GPU, Architecture

using KernelAbstractions

using Random

"""
    Architecture
Abstract type for architectures supported by PlanktonKernels.
"""
abstract type Architecture end

"""
    CPU <: Architecture
Run PlanktonKernels on one CPU node.
"""
struct CPU <: Architecture end

"""
    GPU <: Architecture
Run PlanktonKernels on one CUDA GPU node.
"""
struct GPU <: Architecture end

##### CPU #####
device(::CPU) = KernelAbstractions.CPU()
array_type(::CPU) = Array
rng_type(::CPU) = MersenneTwister()
isfunctional(::CPU) = true
unsafe_free!(m::Array) = nothing

end
