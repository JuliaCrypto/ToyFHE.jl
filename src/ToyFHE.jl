module ToyFHE

export encrypt, decrypt, keygen, keyswitch, PolyCRTEncoding
export BFVParams, BGVParams, CKKSParams
export PubKey, PrivKey, EvalKey, KeyPair

using Random
using Distributions
using Nemo
using Hecke
using GaloisFields
using Mods
using BitIntegers

include("poly.jl")
include("signedmod.jl")
include("pow2_cyc_rings.jl")
using .NTT
export NegacyclicRing
include("utils.jl")
using .Utils: @fields_as_locals
include("rlwe_she.jl")
include("cryptparams.jl")
include("bgv.jl")
include("bfv.jl")
include("ckks.jl")
include("nemo.jl")
include("encoding.jl")
include("polycrtencoding.jl")
include("crt.jl")

import .Utils: plaintext_space
export plaintext_space

if VERSION < v"1.4.0-DEV.208"
    include("div_hacks.jl")
end

end

