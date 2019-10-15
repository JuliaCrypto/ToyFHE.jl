module BGV

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
    import ..ToyFHE: SHEShemeParams, RingSampler, modulus, degree, SignedMod
    export BGVParams

    import ToyFHE: keygen, encrypt, decrypt, coefftype
    import Base: +, *, -


    struct BGVParams <: SHEShemeParams
        # The Cypertext ring over which operations are performed
        ℛ
        # The plaintext ring.
        ℛplain
        σ
    end
    BGVParams(ring, p::Integer, σ) =
        BGVParams(ring, plaintext_space(ring, p), σ)
    plaintext_space(params::BGVParams) = params.ℛplain

    struct PrivKey
        params::BGVParams
        s
    end

    struct PubKey
        params::BGVParams
        a
        b
    end

    struct EvalKey
        params::BGVParams
        a
        b
    end

    struct KeyPair
        priv
        pub
    end

    struct CipherText{T, N}
        params::BGVParams
        cs::NTuple{N, T}
    end
    Base.length(c::CipherText) = length(c.cs)
    Base.getindex(c::CipherText, i::Integer) = c.cs[i]
    Base.lastindex(c::CipherText) = length(c)

    function keygen(rng, params::BGVParams)
        @fields_as_locals params::BGVParams

        dug = RingSampler(ℛ, DiscreteUniform(coefftype(ℛ)))
        dgg = RingSampler(ℛ, DiscreteNormal(0, σ))

        a = rand(rng, dug)
        s = rand(rng, dgg)
        e = rand(rng, dgg)

        p = modulus(coefftype(ℛplain))
        b = a*s + e*p

        KeyPair(PrivKey(params, s), PubKey(params, a, b))
    end
    keygen(params::BGVParams) = keygen(Random.GLOBAL_RNG, params)

    function encrypt(rng::AbstractRNG, key::PubKey, plaintext)
        @fields_as_locals key::PubKey
        @fields_as_locals params::BGVParams

        dgg = RingSampler(ℛ, DiscreteNormal(0, σ))

        v = rand(rng, dgg)
        e0 = rand(rng, dgg)
        e1 = rand(rng, dgg)

        p = modulus(coefftype(ℛplain))
        c0 = b*v + p*e0 + oftype(v, plaintext)
        c1 = a*v + p*e1
        CipherText(params, (c0, c1))
    end
    encrypt(rng::AbstractRNG, kp::KeyPair, plaintext) = encrypt(rng, kp.pub, plaintext)
    encrypt(key, plaintext) = encrypt(Random.GLOBAL_RNG, key, plaintext)

    function decrypt(key::PrivKey, c::CipherText)
        @fields_as_locals key::PrivKey
        @fields_as_locals params::BGVParams

        b = c[1]
        spow = s
        for i = 2:length(c)
            b -= spow*c[i]
            spow *= s
        end

        ℛplain = plaintext_space(params)
        ℛplain(map(x->coefftype(ℛplain)(convert(Integer, mod(SignedMod(x), modulus(base_ring(ℛplain))))), NTT.coeffs_primal(b)))
    end
    decrypt(kp::KeyPair, plaintext) = decrypt(kp.priv, plaintext)

    function *(c1::CipherText, c2::CipherText)
        CipherText(c1.params,(
           c1[1] * c2[1],
           c1[1] * c2[2] + c1[2] * c2[1],
         -(c1[2] * c2[2])
        ))
    end

    function keyswitch(ek::EvalKey, c::CipherText)
        @fields_as_locals ek::EvalKey
        # TODO
    end
end
