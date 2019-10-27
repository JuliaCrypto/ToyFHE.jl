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

"""
    ℛ_plain(params::SHEShemeParams)

Given the parameters of the SHE scheme, return the plaintext ring.
"""
function ℛ_plain end
const plaintext_space = ℛ_plain

# These optional function change multiplication. By default they do nothing
mul_expand(params::SHEShemeParams, x) = x
mul_contract(params::SHEShemeParams, x) = x

# Some utilities
Base.show(io::IO, params::SHEShemeParams) = print(io, scheme_name(params), " parameters")
scheme_name(params::SHEShemeParams) = scheme_name(typeof(params))
Broadcast.broadcastable(params::SHEShemeParams) = Ref(params)

################################################################################
#                        The various FHE key types
################################################################################

struct PrivKey{P <: SHEShemeParams}
    params::P
    secret
end
Base.show(io::IO, kp::PrivKey{P}) where {P} = print(io, scheme_name(P), " private key")

struct PubKey{P <: SHEShemeParams}
    params::P
    mask
    masked
end
Base.show(io::IO, kp::PubKey{P}) where {P} = print(io, scheme_name(P), " public key")

struct EvalKey{P <: SHEShemeParams}
    params::P
    masked::Vector
    mask::Vector
end
Base.show(io::IO, kp::EvalKey{P}) where {P} = print(io, scheme_name(P), " evaluation key")

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
function Base.show(io::IO, kp::CipherText{Enc, P, <:Any, N}) where {P, Enc, N}
    print(io, scheme_name(P), " ciphertext (length ", N)
    Enc != Any && print(io, ", encoding $Enc")
    print(io, ")")
end

################################################################################
#                        Key generation
################################################################################

function keygen(rng::AbstractRNG, params::SHEShemeParams)
    𝒰 = RingSampler(ℛ_cipher(params), DiscreteUniform(coefftype(ℛ_cipher(params))))

    mask = rand(rng, 𝒰)
    secret = rand(rng, 𝒢(params))
    error = rand(rng, 𝒩(params))

    masked = -(mask*secret + error)

    KeyPair(
        PrivKey(params, secret),
        PubKey(params, mask, masked))
end

# TODO: CSPRNG here
keygen(params::SHEShemeParams) = keygen(Random.GLOBAL_RNG, params)

################################################################################
#                        encryption/decryption
################################################################################

function encrypt(rng::AbstractRNG, key::PubKey, plaintext)
    @fields_as_locals key::PubKey

    u = rand(rng, 𝒢(params))
    e₁, e₂ = rand(rng, 𝒩(params), 2)

    c₁ = masked*u + e₁ + π⁻¹(params, plaintext)
    c₂ = mask*u + e₂

    return CipherText(params, (c₁, c₂))
end
encrypt(rng::AbstractRNG, kp::KeyPair, plaintext) = encrypt(rng, kp.pub, plaintext)
encrypt(key::KeyPair, plaintext) = encrypt(Random.GLOBAL_RNG, key, plaintext)

function decrypt(key::PrivKey, c::CipherText)
    @fields_as_locals key::PrivKey

    b = c[1]
    spow = secret

    for i = 2:length(c)
        b += spow*c[i]
        spow *= secret
    end

    π(params, b)
end
decrypt(key::KeyPair, plaintext) = decrypt(key.priv, plaintext)

################################################################################
#                        Homomorphic arithmetic
################################################################################

for f in (:+, :-)
    @eval function $f(c1::CipherText{T,N1}, c2::CipherText{T,N2}) where {T,N1,N2}
        throw(UsageError("Attempting to add ciphertexts with differing parameters"))
        CipherText((
            i > length(c1) ? c2[i] :
            i > length(c2) ? c1[i] :
            $f(c1[i], c2[i]) for i in max(N1, N2)))
    end
end

function *(c1::CipherText{P, Enc, T}, c2::CipherText{P, Enc, T}) where {P, Enc, T}
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

    CipherText(params, (c...,))
end

################################################################################
#                        Key switching
################################################################################

function make_eval_key(rng::AbstractRNG, ::Type{EvalKey}, (old, new)::Pair{<:Any, <:PrivKey})
    @fields_as_locals new::PrivKey

    𝒰 = RingSampler(ℛ_cipher(params), DiscreteUniform(coefftype(ℛ_cipher(params))))
    𝒩gen = 𝒩(params)

    nwindows = ndigits(modulus(coefftype(ℛ_cipher(params))), base=2^params.relin_window)
    evala = [old * coefftype(ℛ_cipher(params))(2)^(i*params.relin_window) for i = 0:nwindows-1]
    evalb = eltype(evala)[]

    for i = 1:length(evala)
        mask = rand(rng, 𝒰)
        e = rand(rng, 𝒩gen)
        push!(evalb, mask)
        evala[i] -= mask*new.secret + e
    end
    EvalKey(new.params, evala, evalb)
end
keygen(rng::AbstractRNG, ::Type{EvalKey}, priv::PrivKey) = make_eval_key(rng, EvalKey, priv.secret^2=>priv)
keygen(::Type{EvalKey}, priv::PrivKey) = keygen(Random.GLOBAL_RNG, EvalKey, priv)

function keyswitch(ek::EvalKey, c::CipherText)
    @fields_as_locals ek::EvalKey

    @assert length(c.cs) in (2,3)
    nwindows = ndigits(modulus(coefftype(ℛ_cipher(params))), base=2^params.relin_window)

    c1 = c[1]
    c2 = length(c) == 2 ? zero(c[2]) : c[2]

    cendcoeffs = NTT.coeffs_primal(c[end])
    ds = map(cendcoeffs) do x
        digits(x.n, base=2^params.relin_window, pad=nwindows)
    end
    ps = map(1:nwindows) do i
        ℛ_cipher(params)([coefftype(ℛ_cipher(params))(ds[j][i]) for j in eachindex(cendcoeffs)])
    end

    for i in eachindex(ek.masked)
        c2 += ek.mask[i] * ps[i]
        c1 += ek.masked[i] * ps[i]
    end

    CipherText(ek.params, (c1, c2))
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
