module BFV

    using Random
    using Distributions
    using GaloisFields
    using ..Karney
    using ..NTT

    import GaloisFields: PrimeField
    import ..Utils: @fields_as_locals
    export BGVParams

    import FHE: keygen, encrypt, decrypt
    import Base: +, *

    struct BFVParams <: SHEShemeParams
        # The Cypertext ring over which operations are performed
        ℛ
        # The big ring used during multiplication
        ℛbig
        # The plain modulus. Plaintexts are elements mod p.
        p
        σ
    end

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
        cs::NTuple{N, T}
    end

    function keygen(rng, params::BFVParams)
        @fields_as_locals params::BFVParams

        dug = RingSampler{ℛ}(DiscreteUniform(eltype(ℛ)))
        dgg = RingSampler{ℛ}(DiscreteNormal{eltype(ℛ)}(0, σ))

        a = nntt(rand(rng, dug))
        s = nntt(rand(rng, dgg))

        e = nntt(rand(rng, dgg))

        KeyPair(
            PrivKey(params, s),
            PubKey(a, a*s - e))
    end

    function encrypt(rng::AbstractRNG, key::PubKey, plaintext)
        @fields_as_locals key::PubKey
        @fields_as_locals params::BFVParams

        dgg = RingSampler{ℛ}(DiscreteNormal{eltype(ℛ)}(0, σ))

        u = nntt(rand(rng, dgg))
        e₁ = nntt(rand(rang, dgg))
        e₂ = nntt(rand(rang, dgg))

        c₁ = b*u + e₁ + Δ * plaintext
        c₂ = a*u + e₂

        return CipherText((c₁, c₂))
    end

    for f in (:+, :-)
        @eval $f(c1::CipherText{T,N1}, c2::CipherText{T,N2}) where {T,N1,N2}
            CipherText((
                i > length(c1) ? c2[i] :
                i > length(c2) ? c1[i] :
                $f(c1[i], c2[i]) for i in max(N1, N2))
        end
    end

    function multround(e, a, b)
        div(e * a, b, round)
    end

    function *(c1::CipherText{T}, c2::CipherText{T}) where {T}
        @fields_as_locals c1.params::BFVParams

        modswitch(c) = nntt(switch(ℛbig, inntt(c)))
        c1 = map(modswitch, c1.cs)
        c2 = map(modswitch, c2.cs)

        c = [zero(typeof(c1)) for i = 1:(length(c1) + length(c2) - 1)]
        for i = 1:length(c1), j = 1:length(c2)
            c[i+j-1] += c1[i] * c2[j]
        end

        c = map(c) do e
            modswitch(ℛ, multround.(inntt(e), p, q))
        end
    end
end
