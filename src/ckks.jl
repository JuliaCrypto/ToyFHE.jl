module CKKS

    using Random
    using Distributions
    using GaloisFields
    using ..NTT
    using ..CryptParameters
    using Primes
    using BitIntegers
    using Nemo
    using AbstractAlgebra
    using Mods

    import GaloisFields: PrimeField
    import ..Utils: @fields_as_locals, plaintext_space
    import ..ToyFHE: SHEShemeParams, RingSampler, modulus, degree
    export BFVParams

    import ToyFHE: keygen, encrypt, decrypt, coefftype
    import Base: +, *, -

    export CKKSParams, FixedRational

    """
        FixedRational{T<:Integer, den} <: Real
    Rational number type, with numerator of type T and fixed denominator `den`.
    """
    struct FixedRational{T, denom}
        x::T
        Base.reinterpret(::Type{FixedRational{T, denom}}, x::T) where {T,denom} =
            new{T, denom}(x)
        function FixedRational{T, denom}(x::Real) where {T, denom}
            n = round(BigInt, big(x)*denom)
            if n < 0
                n = modulus(T) + n
            end
            new{T, denom}(convert(T, n))
        end
    end
    function Base.convert(::Type{Float64}, fr::FixedRational{T, denom}) where {T, denom}
        n = fr.x.n
        if n > div(modulus(T), 2)
            n = n - modulus(T)
        end
        n/denom
    end
    function Base.show(io::IO, fr::FixedRational{<:Any, denom}) where {denom}
        print(io, convert(Float64, fr))
    end

    struct CKKSParams <: SHEShemeParams
        # The Cypertext ring over which operations are performed
        ℛ
        # The big ring used during multiplication
        ℛbig
        relin_window
        σ
    end

    struct PrivKey
        params::CKKSParams
        s
    end

    struct PubKey
        params::CKKSParams
        a
        b
    end

    struct EvalKey
        params::CKKSParams
        a
        b
    end

    struct CipherText{T, N}
        params::CKKSParams
        cs::NTuple{N, T}
    end
    Base.length(c::CipherText) = length(c.cs)
    Base.getindex(c::CipherText, i::Integer) = c.cs[i]
    Base.lastindex(c::CipherText) = length(c)

    struct KeyPair
        priv
        pub
    end
    Base.show(io::IO, kp::KeyPair) = print(io, "CKKS key pair")

    function keygen(rng, params::CKKSParams)
        @fields_as_locals params::CKKSParams

        dug = RingSampler(ℛ, DiscreteUniform(coefftype(ℛ)))
        dgg = RingSampler(ℛ, DiscreteNormal(0, σ))

        a = nntt_hint(rand(rng, dug))
        s = nntt_hint(rand(rng, dgg))

        e = nntt_hint(rand(rng, dgg))

        KeyPair(
            PrivKey(params, s),
            PubKey(params, a, -(a*s + e)))
    end

    function encrypt(rng::AbstractRNG, key::PubKey, plaintext)
        @fields_as_locals key::PubKey
        @fields_as_locals params::CKKSParams

        dgg = RingSampler(ℛ, DiscreteNormal(0, σ))

        u = nntt_hint(rand(rng, dgg))
        e₁ = nntt_hint(rand(rng, dgg))
        e₂ = nntt_hint(rand(rng, dgg))

        c₁ = b*u + e₁ + plaintext
        c₂ = a*u + e₂

        return CipherText(params, (c₁, c₂))
    end
    encrypt(rng::AbstractRNG, kp::KeyPair, plaintext) = encrypt(rng, kp.pub, plaintext)
    encrypt(key::KeyPair, plaintext) = encrypt(Random.GLOBAL_RNG, key, plaintext)

    function decrypt(key::PrivKey, c::CipherText)
        @fields_as_locals key::PrivKey
        @fields_as_locals params::CKKSParams

        b = nntt_hint(c[1])
        spow = s

        for i = 2:length(c)
            b += spow*nntt_hint(c[i])
            spow *= s
        end

        inntt_hint(b)
    end
    decrypt(key::KeyPair, plaintext) = decrypt(key.priv, plaintext)

end
