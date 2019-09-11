module FHE

export encrypt, decrypt, keygen

abstract type SHEShemeParams end

using Random

# For now
using BitIntegers
Base.widen(::Type{Int128}) = Int256

# TODO: CSPRNG here
keygen(params::SHEShemeParams) = keygen(Random.GLOBAL_RNG, params)

function encrypt end
function decrypt end

include("poly.jl")
include("Karney.jl")
include("NTT.jl")
using .NTT
include("utils.jl")
include("cryptparams.jl")
include("bgv.jl")
include("bfv.jl")
include("nemo.jl")

end

