module BFV

    using Random
    using Distributions
    using GaloisFields
    using ..Karney
    using ..NTT
    using ..CryptParameters
    using Primes

    import GaloisFields: PrimeField
    import ..Utils: @fields_as_locals, fqmod
    import ..FHE: SHEShemeParams
    export BFVParams

    import FHE: keygen, encrypt, decrypt
    import Base: +, *, -

    struct BFVParams <: SHEShemeParams
        # The Cypertext ring over which operations are performed
        ℛ
        # The big ring used during multiplication
        ℛbig
        # The plain modulus. Plaintexts are elements mod p.
        p
        σ
        Δ
    end

    # Matches parameter generation in PALISADE
    function BFVParams(p, σ=8/√(2π), α=9, r=1; eval_mult_count = 0, security = HEStd_128_classic)
        @assert r >= 1
        Berr = σ*√(α)
        Bkey = Berr
        δ(n) = 2*√(n)
        Vnorm(n) = Berr * (1 + 2*δ(n)*Bkey)

        function nRLWE(q)
            if isa(security, StdSecurity)
                CryptParameters.std_ring_dim(HEStd_error, security, ceil(log2(q)))
            else
                # The security parameter is interpreted as the hermite factor as
                # in PALISADE.
                log2(q / σ) / (4 * log2(security));
            end
        end

        n = 512
        q = 0

        if eval_mult_count > 0
            w = 2^r
            ϵ₁(n) = 4 / δ(n)*Bkey
            C₁(n) = (1 + ϵ₁(n))*δ(n)^2*p*Bkey
            C₂(n, qPrev) =
                δ(n)^2*p*Bkey*(Bkey + p^2) +
                δ(n)*(floor(log2(qPrev) / r) + 1)*w*Berr
            qBFV(n, qPrev) =
                p^2 + 2p*(
                    C₁(n)^eval_mult_count * Vnorm(n) +
                    eval_mult_count*C₁(n)^(eval_mult_count-1)*C₂(n, qPrev))

            qPrev = 1e6
            q = qBFV(n, qPrev)
            qPrev = q

            while nRLWE(q) > n
                while nRLWE(q) > n
                    n *= 2
                    # TODO: So in original, but is this right?
                    # Shouldn't we set qPrev = q first.
                    q = qBFV(n, qPrev)
                    qPrev = q
                end

                q = qBFV(n, qPrev)

                while abs(q - qPrev) > 0.001q
                    qPrev = q
                    q = qBFV(n, qPrev)
                end
            end
        end

        qPrime = nextprime(Int128(2)^(ceil(Int, log2(q))+1) + 1, 1; interval=2n)
        qPrimeLarge = nextprime(Int128(2)^(2*ceil(Int, log2(q)) + ceil(Int, log2(p)) + 3) + 1, 1; interval=2n)

        Δ = div(qPrime, p)

        𝔽 = GaloisField(qPrime)
        ℛ = LWERing{𝔽, n}(GaloisFields.minimal_primitive_root(𝔽, 2n))
        𝔽big = GaloisField(qPrimeLarge)
        ℛbig = LWERing{𝔽big, n}(GaloisFields.minimal_primitive_root(𝔽big, 2n))

        BFVParams(ℛ, ℛbig, p, σ, Δ)
    end

    struct PrivKey
        params::BFVParams
        s
    end

    struct PubKey
        params::BFVParams
        a
        b
    end

    struct EvalKey
        params::BFVParams
        a
        b
    end

    struct KeyPair
        priv
        pub
    end

    struct CipherText{T, N}
        params::BFVParams
        cs::NTuple{N, T}
    end
    Base.length(c::CipherText) = length(c.cs)
    Base.getindex(c::CipherText, i::Integer) = c.cs[i]

    function keygen(rng, params::BFVParams)
        @fields_as_locals params::BFVParams

        dug = RingSampler{ℛ}(DiscreteUniform(eltype(ℛ)))
        dgg = RingSampler{ℛ}(DiscreteNormal{eltype(ℛ)}(0, σ))

        a = nntt(rand(rng, dug))
        s = nntt(rand(rng, dgg))

        e = nntt(rand(rng, dgg))

        KeyPair(
            PrivKey(params, s),
            PubKey(params, a, -(a*s + e)))
    end

    function encrypt(rng::AbstractRNG, key::PubKey, plaintext)
        @fields_as_locals key::PubKey
        @fields_as_locals params::BFVParams

        dgg = RingSampler{ℛ}(DiscreteNormal{eltype(ℛ)}(0, σ))

        u = nntt(rand(rng, dgg))
        e₁ = nntt(rand(rng, dgg))
        e₂ = nntt(rand(rng, dgg))

        c₁ = b*u + e₁ + Δ * plaintext
        c₂ = a*u + e₂

        return CipherText(params, (c₁, c₂))
    end
    encrypt(rng::AbstractRNG, kp::KeyPair, plaintext) = encrypt(rng, kp.pub, plaintext)
    encrypt(key::KeyPair, plaintext) = encrypt(Random.GLOBAL_RNG, key, plaintext)

    for f in (:+, :-)
        @eval function $f(c1::CipherText{T,N1}, c2::CipherText{T,N2}) where {T,N1,N2}
            CipherText((
                i > length(c1) ? c2[i] :
                i > length(c2) ? c1[i] :
                $f(c1[i], c2[i]) for i in max(N1, N2)))
        end
    end

    function multround(e::Integer, a::Integer, b::Integer)
        div(e * a, b, RoundNearestTiesAway)
    end

    function multround(e::PrimeField, a::Integer, b::Integer)
        q = GaloisFields.char(e)
        halfq = q >> 1
        if e.n > halfq
            return typeof(e)(q - multround(q - e.n, a, b))
        else
            return typeof(e)(multround(e.n, a, b))
        end
    end
    function multround(e::LWERingElement{ℛ}, a::Integer, b::Integer) where {ℛ}
        LWERingElement(ℛ)(FixedDegreePoly(map(e.p.p) do x
            multround(x, a, b)
        end))
    end

    divround(e::Integer, q::Integer) = div(e, q, RoundNearestTiesAway)
    function divround(e::PrimeField, d::Integer)
        q = GaloisFields.char(e)
        halfq = q >> 1
        if e.n > halfq
            return typeof(e)(q - divround(q - e.n, d))
        else
            return typeof(e)(divround(e.n, d))
        end
    end

    function switch(::Type{T}, e::S) where {T<:PrimeField, S<:PrimeField}
        q = GaloisFields.char(e)
        halfq = q >> 1
        diff = abs(char(T) - q)
        if (q < char(T))
            if e.n > halfq
                return T(e.n + diff)
            else
                return T(e.n)
            end
        else
            if e.n > halfq
                return T(e.n - diff)
            else
                return T(e.n)
            end
        end
    end

    function switch(ℛ::LWERing, e::LWERingElement)
        LWERingElement(ℛ)(FixedDegreePoly(map(e.p.p) do x
            switch(eltype(ℛ), x)
        end))
    end

    function *(c1::CipherText{T}, c2::CipherText{T}) where {T}
        @fields_as_locals c1.params::BFVParams

        modswitch(c) = nntt(switch(ℛbig, inntt(c)))
        c1 = map(modswitch, c1.cs)
        c2 = map(modswitch, c2.cs)

        c = [zero(typeof(c1[1])) for i = 1:(length(c1) + length(c2) - 1)]
        for i = 1:length(c1), j = 1:length(c2)
            c[i+j-1] += c1[i] * c2[j]
        end

        c = map(c) do e
            switch(ℛ, multround(inntt(e), p, char(eltype(ℛ))))
        end

        CipherText(c1.params, (c...,))
    end

    maybe_nntt(x::LWERingElement) = nntt(x)
    maybe_nntt(x::LWERingDualElement) = x

    function decrypt(key::PrivKey, c::CipherText)
        @fields_as_locals key::PrivKey
        @fields_as_locals params::BFVParams

        b = maybe_nntt(c[1])
        spow = s

        for i = 2:length(c)
            b += spow*maybe_nntt(c[i])
            spow *= s
        end

        b = inntt(b)
        FixedDegreePoly(map(x->UInt8(fqmod(divround(x, Δ), p)), b.p.p))
    end
    decrypt(key::KeyPair, plaintext) = decrypt(key.priv, plaintext)
end
