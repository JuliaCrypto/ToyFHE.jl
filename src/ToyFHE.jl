module ToyFHE

export encrypt, decrypt, keygen, keyswitch, PolyCRTEncoding, CKKSEncoding
export BFVParams, BGVParams, CKKSParams
export PubKey, PrivKey, EvalKey, EvalMultKey, GaloisKey, GaloisKeys, KeyPair, CipherText
export modswitch, modswitch_drop
export ℛ_cipher, ℛ_key, ℛ_plain

export coefftype

using Random
using Distributions
using Nemo
using Hecke
using GaloisFields
using Mods
using BitIntegers
using StructArrays

include("poly.jl")
include("signedmod.jl")
include("pow2_cyc_rings.jl")
using .NTT
export NegacyclicRing
include("utils.jl")
using .Utils: @fields_as_locals
include("rlwe_she.jl")
include("cryptparams.jl")
include("crt.jl")
include("bgv.jl")
include("bfv.jl")
include("ckks.jl")
include("ckksencoding.jl")
include("nemo.jl")
include("encoding.jl")
include("polycrtencoding.jl")
include("insecuredebug.jl")
include("modulusraising.jl")

import .Utils: plaintext_space
export plaintext_space

if VERSION < v"1.4.0-DEV.208"
    include("div_hacks.jl")
end

end

