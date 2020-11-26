module SL2Cfoam

using HalfIntegers
using LinearAlgebra
using Libdl

export Spin, dim, intertwiner_range,
       VerbosityOff, LowVerbosity, HighVerbosity,
       NormalAccuracy, HighAccuracy, VeryHighAccuracy,
       Vertex, Boosters, CoherentState,
       vertex_amplitude, vertex_compute, vertex_load,
       boosters_compute, boosters_load, b4_compute,
       LinkIn, LinkOut, coherentstate_compute,
       contract

const clib = "libsl2cfoam"

# Initialization of shared libraries.
function __init__()

    # check if sl2cfoam C library is loaded correctly
    dlopen(clib)

    nothing

end
  
###################################################################
# Spin functions.
###################################################################

const Spin = Union{Integer, HalfInteger}
const Cbool = Cuchar

"Dimension of spin representations."
dim(x) = 2x + 1
dim(x::HalfInteger) = twice(x) + 1

"Allowed intertwiner range for spins (j1, j2, j3, j4)
in the recoupling (j1, j2) -- (j3, j4)."
function intertwiner_range(j1::Spin, j2::Spin, j3::Spin, j4::Spin)

    if (2j1 + 2j2) % 2 != (2j3 + 2j4) % 2
        return (), 0
    end

    imin = max(abs(j1-j2), abs(j3-j4))
    imax = min(j1+j2, j3+j4)

    range = imax >= imin ? (imin, imax) : ()
    size = imax >= imin ? Int(imax - imin) + 1 : 0

    (range, size)

end

# Converts spins to two_spins for C library.
ctwo(j) = trunc(Cint, twice(j))

###################################################################
# C library configuration.
###################################################################  

# Checks C library initialization.
clib_initialized = false
function check_cinit()
    if !clib_initialized; error("libsl2cfoam not initialized") end
    nothing
end

"Verbosity level."
@enum Verbosity VerbosityOff LowVerbosity HighVerbosity

"Accuracy level."
@enum Accuracy NormalAccuracy HighAccuracy VeryHighAccuracy

"Specifies what to do at the end of a vertex computation."
VertexResult = @NamedTuple{ret::Bool, store::Bool, store_batches::Bool}
Base.cconvert(::Type{Cint}, vr::VertexResult) = Cint(vr[:ret] + 2 * vr[:store] + 4 * vr[:store_batches])

# C-struct of library configuration.
struct __C_config
    verbosity             :: Cint
    accuracy              :: Cint
    max_spin              :: Cint
    max_MB_mem_per_thread :: Culong
end

"Configuration for the library."
struct Config 
    verbosity             :: Verbosity
    accuracy              :: Accuracy
    max_spin              :: Spin
    max_MB_mem_per_thread :: Integer
end


"Initializes the library with given configuration."
function cinit(folder::String, Immirzi::Real, conf::Config)

    if !isdir(folder) throw(ArgumentError("error opening folder $folder")) end
    if Immirzi <= 0 throw(ArgumentError("Immirzi parameter must be strictly positive")) end

    cconf = __C_config(Int(conf.verbosity), Int(conf.accuracy), conf.max_spin, conf.max_MB_mem_per_thread)
    ccall((:sl2cfoam_init_conf, clib), Cvoid, (Cstring, Cdouble, Ref{__C_config}), folder, Immirzi, Ref(cconf))

    global clib_initialized = true
    nothing

end

"Clears the library."
function cclear()

    @ccall clib.sl2cfoam_free()::Cvoid

    global clib_initialized = false
    nothing

end

function set_verbosity(v::Verbosity)
    @ccall clib.sl2cfoam_set_verbosity(v::Cint)::Cvoid
end

function set_accuracy(a::Accuracy)
    @ccall clib.sl2cfoam_set_accuracy(a::Cint)::Cvoid
end

function set_Immirzi(Immirzi::Real)

    if Immirzi <= 0 throw(ArgumentError("Immirzi parameter must be strictly positive")) end
    @ccall clib.sl2cfoam_set_Immirzi(Immirzi::Cdouble)::Cvoid

end

function set_OMP(enable::Bool)

    @ccall clib.sl2cfoam_set_omp(enable::Cbool)::Cvoid

end

# sanity checks for arguments
check_spins(js, n) = if (ng = length(js)) != n; throw(ArgumentError("$n spins required, got $ng")) end


###################################################################
# C types mappings.
###################################################################  

# C-struct of SL2Cfoam tensor.
# TODO: read tag data
struct __C_tensor{N}
    num_keys :: Cuchar
    dims     :: NTuple{N, Cuint}
    strides  :: NTuple{N, Culong}
    dim      :: Culong
    d        :: Ptr{Cdouble}
    tag      :: Ptr{Cvoid}
end

# C-struct for vertex tensors.
__C_vertex_tensor = __C_tensor{5}

# C-struct for boosters tensors.
__C_boosters_tensor = __C_tensor{6}


###################################################################
# Vertex functions.
###################################################################  

"Vertex object. Contains the data array and the pointer
to the library tensor."
mutable struct Vertex

    a       :: Array{Float64, 5}
    cptr    :: Ptr{__C_vertex_tensor}

    function Vertex(a::Array, cptr)

        v = new(a, cptr)
        finalizer(v) do x
            x.cptr != C_NULL && ccall((:sl2cfoam_vertex_free, clib), Cvoid, (Ptr{__C_vertex_tensor},), x.cptr)
        end
        return v

    end

    function Vertex(cptr)

        if cptr == C_NULL; error("libsl2cfoam returned a NULL pointer") end
        ctens = unsafe_load(cptr)
        a = unsafe_wrap(Array, ctens.d, ctens.dims; own = false)
        return Vertex(a, cptr)

    end

end

"Computes a single vertex amplitude given 10 spins (j12, ...),
5 intertwiners (i1,...) and number of shells."
function vertex_amplitude(js, is, Dl)

    check_cinit()
    check_spins(js, 10)
    check_spins(is, 5)

    ccall((:sl2cfoam_vertex_amplitude, clib), Cdouble, 
          (Ref{Cint}, Ref{Cint}, Cint), ctwo.(js), ctwo.(is), Dl)

end

"Computes the vertex tensor given 10 spins (j12, ...) and number of shells.
It optionally takes an intertwiner range ((i1_min, i1_max), (i2_min, i2_max)...)
and a VertexResult object to specify what to to with the result."
function vertex_compute(js, Dl::Integer; result = VertexResult((true, true, false)))

    check_cinit()
    check_spins(js, 10)

    cptr = ccall((:sl2cfoam_vertex_fullrange, clib), Ptr{__C_vertex_tensor}, (Ref{Cint}, Cint, Cint), ctwo.(js), Dl, result)

    Vertex(cptr)

end

function vertex_compute(js, Dl::Integer, irange; result = VertexResult((true, true, false)))

    check_cinit()
    check_spins(js, 10)

    cptr = ccall((:sl2cfoam_vertex_range, clib), Ptr{__C_vertex_tensor}, 
                 (Ref{Cint}, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint),
                 ctwo.(js), 
                 ctwo(irange[1][1]), ctwo(irange[1][2]),
                 ctwo(irange[2][1]), ctwo(irange[2][2]),
                 ctwo(irange[3][1]), ctwo(irange[3][2]),
                 ctwo(irange[4][1]), ctwo(irange[4][2]),
                 ctwo(irange[5][1]), ctwo(irange[5][2]),
                 Dl, result)

    Vertex(cptr)

end

"Loads a vertex tensor given the Immirzi parameter,
10 spins, number of shells and optional intertwiners range."
function vertex_load(js, Dl::Integer)

    check_cinit()
    check_spins(js, 10)

    cptr = ccall((:sl2cfoam_vertex_fullrange_load, clib), Ptr{__C_vertex_tensor}, (Ref{Cint}, Cint), ctwo.(js), Dl)

    Vertex(cptr)

end

function vertex_load(js, Dl::Integer, irange)

    check_cinit()
    check_spins(js, 10)

    cptr = ccall((:sl2cfoam_vertex_range_load, clib), Ptr{__C_vertex_tensor}, 
                 (Ref{Cint}, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint),
                 ctwo.(js), 
                 ctwo(irange[1][1]), ctwo(irange[1][2]),
                 ctwo(irange[2][1]), ctwo(irange[2][2]),
                 ctwo(irange[3][1]), ctwo(irange[3][2]),
                 ctwo(irange[4][1]), ctwo(irange[4][2]),
                 ctwo(irange[5][1]), ctwo(irange[5][2]),
                 Dl)

    Vertex(cptr)

end

# The following function can be called without 
# initializing the library.

"Loads a vertex tensor computed with the C library
given the absolute path to the file."
function vertex_load(path::String)

    cptr = ccall((:sl2cfoam_vertex_load, clib), Ptr{__C_vertex_tensor}, (Cstring,), path)
    Vertex(cptr)
    
end


###################################################################
# Boosters functions.
###################################################################  

"Boosters object. Contains the data array and the pointer
to the library tensor."
mutable struct Boosters

    a       :: Array{Float64, 6}
    cptr    :: Ptr{__C_boosters_tensor}

    function Boosters(a::Array, cptr)

        v = new(a, cptr)
        finalizer(v) do x
            x.cptr != C_NULL && ccall((:sl2cfoam_boosters_free, clib), Cvoid, (Ptr{__C_boosters_tensor},), x.cptr)
        end
        return v

    end

    function Boosters(cptr)

        if cptr == C_NULL; error("libsl2cfoam returned a NULL pointer") end
        ctens = unsafe_load(cptr)
        a = unsafe_wrap(Array, ctens.d, ctens.dims; own = false)
        return Boosters(a, cptr)

    end

end

"Computes a boosters tensor, given the gauge-fixed index (1 to 4), 4 spins
and number of shells. Spins order must match the order of the symbol (anti-clockwise).
An optional store parameter sets if to store the tensor after computation."
function boosters_compute(gf, js, Dl::Integer; store = true)

    check_cinit()
    check_spins(js, 4)
    !(1 <= gf <= 4) && throw(ArgumentError("gauge-fixed index must be 1 to 4"))

    cptr = ccall((:sl2cfoam_boosters, clib), Ptr{__C_boosters_tensor}, (Cint, Cint, Cint, Cint, Cint, Cint, Cbool), 
                 gf, ctwo(js[1]), ctwo(js[2]), ctwo(js[3]), ctwo(js[4]), Dl, store)

    Boosters(cptr)

end

"Loads a computed tensor for the boosters given gauge-fixed index,
spins and number of shells."
function boosters_load(gf, js, Dl::Integer)

    check_cinit()
    check_spins(js, 4)
    !(1 <= gf <= 4) && throw(ArgumentError("gauge-fixed index must be 1 to 4"))


    cptr = ccall((:sl2cfoam_boosters_load, clib), Ptr{__C_boosters_tensor}, (Cint, Cint, Cint, Cint, Cint, Cint),
                 gf, ctwo(js[1]), ctwo(js[2]), ctwo(js[3]), ctwo(js[4]), Dl)

    Boosters(cptr)

end

"Computes a boosters tensor, given the gauge-fixed index (1 to 4), 4 spins
and number of shells. Spins order must match the order of the symbol (anti-clockwise).
An optional store parameter sets if to store the tensor after computation."
function b4_compute(js, ls)

    check_cinit()
    check_spins(js, 4)
    check_spins(ls, 4)

    for i in 1:4
        ls[i] < js[i] && throw(ArgumentError("spin l$i = $(ls[i]) must be greater or equal $(js[i])"))
    end 

    cptr = ccall((:sl2cfoam_b4, clib), Ptr{Cdouble}, (Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint), 
                 ctwo(js[1]), ctwo(js[2]), ctwo(js[3]), ctwo(js[4]),
                 ctwo(ls[1]), ctwo(ls[2]), ctwo(ls[3]), ctwo(ls[4]))

    # wrap matrix
    _, isize = intertwiner_range(js...)
    _, ksize = intertwiner_range(ls...)

    m = unsafe_wrap(Array, cptr, (isize, ksize); own = false)

    # copy and release C memory
    b4 = copy(m)
    ccall((:sl2cfoam_vector_free, clib), Cvoid, (Ptr{Cvoid},), cptr)

    b4

end


###################################################################
# Coherent states functions.
###################################################################  

"CoherentState object. Contains the data array and the pointer
to the library array."
mutable struct CoherentState

    a     :: Vector{ComplexF64}
    cptr  :: Ptr{Cdouble}

    function CoherentState(a::Array, cptr)

        v = new(a, cptr)
        finalizer(v) do x
            x.cptr != C_NULL && ccall((:sl2cfoam_vector_free, clib), Cvoid, (Ptr{Cvoid},), x.cptr)
        end
        return v

    end

    function CoherentState(size::Integer, cptr)

        if cptr == C_NULL; error("libsl2cfoam returned a NULL pointer") end
        a = unsafe_wrap(Array, cptr, size; own = false)
        return CoherentState(a, cptr)

    end

end

"Orientation of the links for a coherent state."
@enum CsLink LinkOut=+1 LinkIn=-1

"Computes a vector of coefficients of (normalized) Livine-Speziale coherent states.
Angles are passed in a 4x2 matrix (4 normals x (theta, phi))."
function coherentstate_compute(js, angles, inout::NTuple{4, CsLink})

    check_cinit()
    check_spins(js, 4)

    cptr = ccall((:sl2cfoam_coherentstate_fullrange, clib), Ptr{ComplexF64}, 
                 (Ref{Cint}, Ref{Cdouble}, Ref{NTuple{4,Int32}}), ctwo.(js), Float64.(transpose(angles)), Cint.(inout))

    _, size = intertwiner_range(js...) 

    CoherentState(size, cptr)

end

function coherentstate_compute(js, angles, inout::NTuple{4, CsLink}, irange)

    check_cinit()
    check_spins(js, 4)

    cptr = ccall((:sl2cfoam_coherentstate_range, clib), Ptr{ComplexF64}, 
                 (Ref{Cint}, Cint, Cint, Ref{Cdouble}, Ref{NTuple{4,Int32}}),
                 ctwo.(js), ctwo(irange[1]), ctwo(irange[2]), Float64.(transpose(angles)), Cint.(inout))

    size = Int(irange[2] - irange[1]) + 1

    CoherentState(size, cptr)

end


###################################################################
# Functions for contracting vertices and coherent states.
###################################################################

"Contracts a multidimensional array with a vector
over dimension with stride 1 (leftmost index)."
function contract(a::Array, v::Vector)

    if size(a,1) != length(v)
        throw(ArgumentError("first dimension of the array does not match vector length"))
    end

    av = reshape(a, (length(v), prod(size(a)[2:end])))

    reshape(transpose(av) * v, size(a)[2:end])

end


"Contracts an array with many vectors."
function contract(a::Array, vs::Vararg{Vector, N}) where N

    dima = length(size(a))

    if N < 1 || N > dima
        throw(ArgumentError("1 to $dima vectors required"))
    end

    cn = contract(a, vs[1])

    for v in vs[2:end]
        cn = contract(cn, v)
    end

    if isempty(size(cn)) return cn[1] end

    cn
    
end

"Contracts a vertex with 1 to 5 coherent states vectors
starting from leftmost index (i5 -> i4 -> i3 ...)."
function contract(v::Vertex, css::Vararg{CoherentState, N}) where N

    if N < 1 || N > 5
        throw(ArgumentError("1 to 5 coherent states required"))
    end

    # first contraction is between REAL vertex amplitudes
    # and COMPLEX coherent state coefficients
    re_cn = contract(v.a, real.(css[1].a))
    im_cn = contract(v.a, imag.(css[1].a))

    cn = re_cn .+ im .* im_cn

    for cs in css[2:end]
        cn = contract(cn, cs.a)
    end

    if isempty(size(cn)) return cn[1] end

    cn
    
end


end # module
