module SL2CfoamGPU

using SL2Cfoam
using CUDA

export to_GPU, Base.only

# checks that GPU is available and working
__init__() = @assert CUDA.functional(true)

# allow extracting scalars from CUDA arrays
Base.only(a::CuArray{<:Number, 0}) = CUDA.@allowscalar @inbounds a[]

"Copies an object to GPU memory."
to_GPU(s::Union{SL2Cfoam.GenericVertex, CoherentState, Boosters}) = CuArray(s.a)

end # module GPU