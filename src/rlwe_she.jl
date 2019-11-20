# This file contains common definitions for the three RLWE based SHE
# schemes we support (BGV, BFV, CKKS). See the cryptographic background section
# in the manual for an overview.

################################################################################
#                        SHEShemeParams
################################################################################

abstract type SHEShemeParams end

# These four functions (π, π⁻¹, 𝒩, 𝒢) determine the differences between the
# three schemes we support.
function π end
function π⁻¹ end
function 𝒩 end
function 𝒢 end

"""
    ℛ_cipher(params::SHEShemeParams)

Given the parameters of the SHE scheme, return the ciphertext ring.
"""
function ℛ_cipher end
const ciphertext_space = ℛ_cipher

ℛ_key(p::SHEShemeParams) = ℛ_cipher(p)

"""
    ℛ_plain(params::SHEShemeParams)

Given the parameters of the SHE scheme, return the plaintext ring.
"""
function ℛ_plain end
const plaintext_space = ℛ_plain

relin_window(params::SHEShemeParams) = params.relin_window

# These optional function change multiplication. By default they do nothing
mul_expand(params::SHEShemeParams, x) = x
mul_contract(params::SHEShemeParams, x) = x

# Some utilities
Base.show(io::IO, params::SHEShemeParams) = print(io, scheme_name(params), " parameters")
scheme_name(params::SHEShemeParams) = scheme_name(typeof(params))
Broadcast.broadcastable(params::SHEShemeParams) = Ref(params)


################################################################################
#                        PassthroughParams
################################################################################

abstract type PassthroughParams{P} <: SHEShemeParams; end
parent_params(p::PassthroughParams) = p.params
ℛ_plain(p::PassthroughParams) = ℛ_plain(parent_params(p))
ℛ_cipher(p::PassthroughParams) = ℛ_cipher(parent_params(p))
ℛ_key(p::PassthroughParams) = ℛ_key(parent_params(p))
relin_window(p::PassthroughParams) = relin_window(parent_params(p))
π⁻¹(p::PassthroughParams, plaintext) = π⁻¹(parent_params(p), plaintext)
π(p::PassthroughParams, b) = π(parent_params(p), b)
𝒩(p::PassthroughParams) = 𝒩(parent_params(p))
𝒢(p::PassthroughParams) = 𝒢(parent_params(p))

################################################################################
#                        The various FHE key types
################################################################################

struct PrivKey{P <: SHEShemeParams}
    params::P
    secret
end
Base.show(io::IO, kp::PrivKey{P}) where {P} = print(io, scheme_name(P), " private key")

struct KeyComponent
    mask
    masked
end

struct PubKey{P <: SHEShemeParams}
    params::P
    # One can think of a public key as a keyswitching key from 0 to the secret
    # secret key of the scheme. Of course, for 0, we have no noise growth
    # concerns, so we can get away with a single KeyComponent
    key::KeyComponent
end
Base.show(io::IO, kp::PubKey{P}) where {P} = print(io, scheme_name(P), " public key")

struct KeySwitchKey{P <: SHEShemeParams}
    params::P
    key::Vector{KeyComponent}
end
Base.show(io::IO, kp::KeySwitchKey{P}) where {P} = print(io, scheme_name(P), " key-switching key")

abstract type EvalKey end
struct EvalMultKey{P <: SHEShemeParams} <: EvalKey
    key::KeySwitchKey{P}
end
Base.show(io::IO, gk::EvalMultKey{P}) where {P} = print(io, scheme_name(P), " multiplication key")

struct GaloisKey{P <: SHEShemeParams} <: EvalKey
    galois_element::Int
    key::KeySwitchKey{P}
end
Base.show(io::IO, gk::GaloisKey{P}) where {P} = print(io, scheme_name(P), " galois key (element ", gk.galois_element, ")")

struct GaloisKeys{P <: SHEShemeParams}
    # Sorted collection of galois keys for various step sizes
    keys::Vector{GaloisKey{P}}
end
Base.show(io::IO, gk::GaloisKeys{P}) where {P} = print(io, scheme_name(P), " galois keys (elements ", join(map(k->k.galois_element, gk.keys), ", ", ")"))

struct KeyPair{P <: SHEShemeParams}
    priv::PrivKey{P}
    pub::PubKey{P}
end
Base.show(io::IO, kp::KeyPair{P}) where {P} = print(io, scheme_name(P), " key pair")

Base.broadcastable(k::Union{PrivKey, PubKey, EvalKey, KeyPair}) = Ref(k)

################################################################################
#                        CipherText
################################################################################

"""
    CipherText{P <: SHEShemeParams, Enc, T, N}

The CipherText for an RLWE-based FHE scheme (the scheme is determined by the
parameter `p`). Optionally, a corresponding plaintext encoding `Enc` may be
specified that will be applied on decryption (`Any` to obtain the raw
decryption result in the scheme's native plaintext space).
"""
struct CipherText{Plain, P <: SHEShemeParams, T, N}
    params::P
    cs::NTuple{N, T}
end
CipherText(params::P, cs::NTuple{N,T}) where {P <: SHEShemeParams,T,N} =
    CipherText{Any,P,T,N}(params, cs)
CipherText{Plain}(params::P, cs::NTuple{N,T}) where {Plain, P <: SHEShemeParams,T,N} =
    CipherText{Plain,P,T,N}(params, cs)

Base.length(c::CipherText) = length(c.cs)
Base.getindex(c::CipherText, i::Integer) = c.cs[i]
Base.lastindex(c::CipherText) = length(c)
Base.broadcastable(c::CipherText) = Ref(c)
function Base.show(io::IO, kp::CipherText{Enc, P, <:Any, N}) where {P, Enc, N}
    print(io, scheme_name(P), " ciphertext (length ", N)
    Enc != Any && print(io, ", encoding $Enc")
    print(io, ")")
end

################################################################################
#                        Key generation
################################################################################

function keygen(rng::AbstractRNG, params::SHEShemeParams)
    𝒰 = RingSampler(ℛ_key(params), DiscreteUniform(coefftype(ℛ_key(params))))

    mask = rand(rng, 𝒰)
    secret = rand(rng, 𝒢(params))
    error = rand(rng, 𝒩(params))

    masked = -(mask*secret + error)

    KeyPair(
        PrivKey(params, secret),
        PubKey(params, KeyComponent(mask, masked)))
end

# TODO: CSPRNG here
keygen(params::SHEShemeParams) = keygen(Random.GLOBAL_RNG, params)

################################################################################
#                        encryption/decryption
################################################################################

function encrypt(rng::AbstractRNG, pk::PubKey, ::Zero)
    @fields_as_locals pk::PubKey
    @fields_as_locals key::KeyComponent

    u = rand(rng, 𝒢(params))
    e₁, e₂ = rand(rng, 𝒩(params), 2)

    c₁ = masked*u + e₁
    c₂ = mask*u + e₂
    CipherText{Zero}(params, (c₁, c₂))
end

function encrypt(rng::AbstractRNG, key::PubKey, plaintext)
    c = encrypt(rng, key, Zero())
    c += oftype(c.cs[1], π⁻¹(key.params, plaintext))

    EncT = typeof(plaintext)
    EncT <: NTT.RingElement && (EncT = Any)
    return CipherText{EncT}(key.params, c.cs)
end
encrypt(rng::AbstractRNG, kp::KeyPair, plaintext) = encrypt(rng, kp.pub, plaintext)
encrypt(key::KeyPair, plaintext) = encrypt(Random.GLOBAL_RNG, key, plaintext)

function decrypt(key::PrivKey, c::CipherText{T}) where T
    @fields_as_locals key::PrivKey

    while eltype(secret) != eltype(c[1])
        secret = modswitch_drop(secret)
    end

    b = c[1]
    spow = secret

    for i = 2:length(c)
        b += spow*c[i]
        spow *= secret
    end

    dec = π(params, b)
    T === Any ? dec : T(dec)
end
decrypt(key::KeyPair, plaintext) = decrypt(key.priv, plaintext)

################################################################################
#                            error handling
################################################################################

struct UsageError
    msg::AbstractString
end

################################################################################
#                        Homomorphic arithmetic
################################################################################

for f in (:+, :-)
    @eval function $f(c1::CipherText{P, Enc, T, N1}, c2::CipherText{P, Enc, T,N2}) where {P, Enc, T, N1, N2}
        if c1.params !== c2.params
            throw(UsageError("Attempting to add ciphertexts with differing parameters"))
        end
        CipherText{P, Enc, T, max(N1, N2)}(c1.params, tuple((
            i > length(c1) ? c2[i] :
            i > length(c2) ? c1[i] :
            $f(c1[i], c2[i]) for i in 1:max(N1, N2))...))
    end
end

function +(c::CipherText{Enc, P, T, N}, b::T) where {T, Enc, P, N}
    CipherText{Enc, P, T, N}(c.params, (c.cs[1] + b, c.cs[2:end]...))
end

function enc_mul(c1, c2)
    if c1.params !== c2.params
        throw(UsageError("Attempting to multiply ciphertexts with differing parameters"))
    end
    params = c1.params

    (c1, c2) = mul_expand.(params, (c1, c2))

    c = [zero(c1[1]) for i = 1:(length(c1) + length(c2) - 1)]
    for i = 1:length(c1), j = 1:length(c2)
        c[i+j-1] += c1[i] * c2[j]
    end

    c = mul_contract(params, c)
    (c...,)
end

function *(c1::CipherText{P, Enc, T}, c2::CipherText{P, Enc, T}) where {P, Enc, T}
    CipherText(c1.params, enc_mul(c1, c2))
end

################################################################################
#                        Key switching
################################################################################

function make_eval_key(rng::AbstractRNG, (old, new)::Pair{<:Any, <:PrivKey})
    @fields_as_locals new::PrivKey

    ℛ = ring(old)

    𝒰 = RingSampler(ℛ, DiscreteUniform(coefftype(ℛ)))
    𝒩gen = 𝒩(params)

    if relin_window(params) != 0
        nwindows = ndigits(modulus(coefftype(ℛ)), base=2^relin_window(params))
        evala = [old * coefftype(ℛ)(2)^(i*relin_window(params)) for i = 0:nwindows-1]
        evalb = eltype(evala)[]
    else
        # CRT basis decomposition
        evala = [typeof(old)(OffsetArray(map(x->eltype(ℛ)(CRTResidual(x)), arr), axes(old)), nothing) for arr in StructArrays.fieldarrays(parent(NTT.coeffs_primal(old)))]
        evalb = eltype(evala)[]
    end

    for i = 1:length(evala)
        mask = rand(rng, 𝒰)
        e = rand(rng, 𝒩gen)
        push!(evalb, mask)
        evala[i] -= mask*new.secret + e
    end
    KeySwitchKey(new.params, map(KeyComponent, evalb, evala))
end
keygen(rng::AbstractRNG, ::Type{EvalMultKey}, priv::PrivKey) = EvalMultKey(make_eval_key(rng, priv.secret^2=>priv))
function keygen(rng::AbstractRNG, ::Type{GaloisKey}, priv::PrivKey; galois_element = nothing, steps=nothing)
    @assert galois_element === nothing || steps === nothing && !(galois_element === nothing && steps === nothing)
    if galois_element === nothing
        @assert steps !== nothing
        galois_element = mod(steps > 0 ? 2degree(priv.secret)-3^steps : 3^-steps, 2degree(priv.secret))
    else
        @assert steps === nothing
    end
    GaloisKey(galois_element, make_eval_key(rng, NTT.apply_galois_element(priv.secret, galois_element)=>priv))
end
keygen(T::Type{<:EvalKey}, priv::PrivKey; kwargs...) = keygen(Random.GLOBAL_RNG, T, priv; kwargs...)

keyswitch_expand(ek, re) = re
keyswitch_contract(ek, c) = c
function keyswitch(ek::KeySwitchKey, c::CipherText{Enc}) where {Enc}
    @fields_as_locals ek::KeySwitchKey

    @assert length(c.cs) in (2,3)
    ℛ = ring(c.cs[1])

    any_mask = ek.key[1].mask

    c1 = keyswitch_expand(ek, c[1])
    c2 = length(c) == 2 ? zero(any_mask) : c[2]

    if relin_window(params) == 0
        # CRT basis decomposition
        cendcoeffs = NTT.coeffs_primal(c[end])
        ps = [typeof(any_mask)(eltype(any_mask).(convert.(Integer, SignedMod.(OffsetArray(arr, axes(cendcoeffs))))), nothing) for arr in StructArrays.fieldarrays(parent(cendcoeffs))]
    else
        # digit-wise decomposition
        cendcoeffs = NTT.coeffs_primal(c[end])
        nwindows = ndigits(modulus(coefftype(ℛ)), base=2^relin_window(params))
        ds = Any[digits(convert(Integer, x), base=2^relin_window(params), pad=nwindows) for x in cendcoeffs]
        ps = map(1:nwindows) do i
            typeof(any_mask)([coefftype(ℛ)(ds[j][i]) for j in eachindex(cendcoeffs)], nothing)
        end
    end

    for i in eachindex(ps)
        c2 += ek.key[i].mask * ps[i]
        c1 += ek.key[i].masked * ps[i]
    end

    CipherText{Enc}(ek.params, (keyswitch_contract(ek, c1), keyswitch_contract(ek, c2)))
end
keyswitch(k::GaloisKey, c::CipherText) = keyswitch(k.key, c)
keyswitch(k::EvalMultKey, c::CipherText) = keyswitch(k.key, c)

################################################################################
#     Rotations
################################################################################

function NTT.apply_galois_element(c::CipherText{Enc}, galois_element) where {Enc}
    CipherText{Enc}(c.params, map(re->NTT.apply_galois_element(re, galois_element), c.cs))
end

rotate(gk::GaloisKey, c::CipherText) = keyswitch(gk, NTT.apply_galois_element(c, gk.galois_element))

################################################################################
#     Modulus switching
################################################################################

function modswitch(c::CipherText{Any}, new_modulus)
    error("Not implemented")
end

################################################################################
#     Plaintext spaces for a given ciphertext space an plain modulus
################################################################################

function plaintext_space(r::ResRing, p)
    ℤp = ResidueRing(Nemo.ZZ, p)
    ℤpx = PolynomialRing(ℤp, "x")[1]
    ResidueRing(ℤpx, Nemo.lift(Nemo.ZZ["x"][1], modulus(r)))
end

function plaintext_space(r::NegacyclicRing, p)
    coefft = Primes.isprime(p) ? GaloisField(p) :
        p == 256 ? UInt8 :
        Mod(p)
    if Primes.isprime(p) && p > 2degree(modulus(r))
        # TODO: Also needs to check here if the prime admits 2n-th roots of
        # unities.
        NegacyclicRing{coefft, degree(modulus(r))}(
            GaloisFields.minimal_primitive_root(coefft, 2degree(modulus(r))))
    else
        NegacyclicRing{coefft, degree(modulus(r))}()
    end
end
