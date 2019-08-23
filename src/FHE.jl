module FHE

export encrypt, decrypt, keygen

abstract SHEShemeParams end

using Random

# TODO: CSPRNG here
keygen(params::SHEShemeParams) = keygen(Random.GLOBAL_RNG, params)

function encrypt end
function decrypt end

include("Karney.jl")
include("NTT.jl")
using .NTT
include("utils.jl")
include("bgv.jl")
include("cryptparams.jl")


end

