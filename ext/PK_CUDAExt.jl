module PK_CUDAExt

using CUDA
import PlanktonKernels.Architectures: GPU, device, array_type, rng_type, isfunctional,  unsafe_free!

device(::GPU) = CUDABackend()
array_type(::GPU) = CuArray
rng_type(::GPU) = CUDA.default_rng()
isfunctional(::GPU) = CUDA.functional()
unsafe_free!(m::CuArray) = CUDA.unsafe_free!(m)

end
