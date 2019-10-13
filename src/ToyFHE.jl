module ToyFHE

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
include("signedmod.jl")
include("pow2_cyc_rings.jl")
using .NTT
include("utils.jl")
include("cryptparams.jl")
include("bgv.jl")
include("bfv.jl")
include("ckks.jl")
include("nemo.jl")
include("crt.jl")

if VERSION < v"1.4.0-DEV.208"
    include("div_hacks.jl")
end

end

