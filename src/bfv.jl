module BFV

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
    import ..Utils: @fields_as_locals, fqmod
    import ..ToyFHE: SHEShemeParams, RingSampler, modulus, degree
    export BFVParams

    import ToyFHE: keygen, encrypt, decrypt, coefftype
    import Base: +, *, -

    struct BFVParams <: SHEShemeParams
        # The Cypertext ring over which operations are performed
        ℛ
        # The big ring used during multiplication
        ℛbig
        # The plain ring.
        ℛplain
        relin_window
        σ
        Δ
    end

    plaintext_space(p::BFVParams) = p.ℛplain

    plaintext_space(r::ResRing, p) = PolynomialRing(ResidueRing(Nemo.ZZ, p), "x")[1]
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

    # Matches parameter generation in PALISADE
    function BFVParams(p, σ=8/√(2π), α=9, r=1; eval_mult_count = 0, security = HEStd_128_classic, relin_window=1)
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
        largebits = 2*ceil(Int, log2(q)) + ceil(Int, log2(p)) + 3
        Tlarge = largebits > 128 ? Int256 : Int128
        qLargeBig = nextprime(big(2)^largebits + 1, 1; interval=2n)
        qPrimeLarge = Tlarge(qLargeBig)

        Δ = div(qPrime, p)

        𝔽 = GaloisField(qPrime)
        ℛ = NegacyclicRing{𝔽, n}(GaloisFields.minimal_primitive_root(𝔽, 2n))
        𝔽big = GaloisField(qPrimeLarge)
        r = GaloisFields.minimal_primitive_root(𝔽big, 2n)
        ℛbig = NegacyclicRing{𝔽big, n}(r)

        BFVParams(ℛ, ℛbig, plaintext_space(ℛ, p), relin_window, σ, Δ)
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
    Base.show(io::IO, kp::KeyPair) = print(io, "BFV key pair")

    struct CipherText{T, N}
        params::BFVParams
        cs::NTuple{N, T}
    end
    Base.length(c::CipherText) = length(c.cs)
    Base.getindex(c::CipherText, i::Integer) = c.cs[i]
    Base.lastindex(c::CipherText) = length(c)

    nntt_hint(r) = r
    nntt_hint(r::NegacyclicRingElement) = nntt(r)
    inntt_hint(r) = r
    inntt_hint(r::NegacyclicRingDualElement) = inntt(r)

    function keygen(rng, params::BFVParams)
        @fields_as_locals params::BFVParams

        dug = RingSampler(ℛ, DiscreteUniform(coefftype(ℛ)))
        dgg = RingSampler(ℛ, DiscreteNormal(0, σ))

        a = nntt_hint(rand(rng, dug))
        s = nntt_hint(rand(rng, dgg))

        e = nntt_hint(rand(rng, dgg))

        KeyPair(
            PrivKey(params, s),
            PubKey(params, a, -(a*s + e)))
    end

    function make_eval_key(rng::AbstractRNG, ::Type{EvalKey}, (old, new)::Pair{<:Any, PrivKey})
        @fields_as_locals new::PrivKey
        @fields_as_locals params::BFVParams

        dug = RingSampler(ℛ, DiscreteUniform(coefftype(ℛ)))
        dgg = RingSampler(ℛ, DiscreteNormal(0, σ))

        nwindows = ndigits(modulus(coefftype(ℛ)), base=2^relin_window)
        evala = [old * coefftype(params.ℛ)(2)^(i*relin_window) for i = 0:nwindows-1]
        evalb = eltype(evala)[]

        for i = 1:length(evala)
            a = nntt_hint(rand(rng, dug))
            e = nntt_hint(rand(rng, dgg))
            push!(evalb, a)
            evala[i] -= a*new.s + e
        end
        EvalKey(new.params, evala, evalb)
    end
    keygen(rng::AbstractRNG, ::Type{EvalKey}, priv::PrivKey) = make_eval_key(rng, EvalKey, priv.s^2=>priv)
    keygen(::Type{EvalKey}, priv::PrivKey) = keygen(Random.GLOBAL_RNG, EvalKey, priv)

    function encrypt(rng::AbstractRNG, key::PubKey, plaintext)
        @fields_as_locals key::PubKey
        @fields_as_locals params::BFVParams

        dgg = RingSampler(ℛ, DiscreteNormal(0, σ))

        u = nntt_hint(rand(rng, dgg))
        e₁ = nntt_hint(rand(rng, dgg))
        e₂ = nntt_hint(rand(rng, dgg))

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
    multround(e::fmpz, a::Integer, b::Integer) = multround(BigInt(e), a, b)
    multround(e::fmpz, a::Integer, b::fmpz) = multround(BigInt(e), a, BigInt(b))

    function multround(e::Union{PrimeField, Nemo.nmod, AbstractAlgebra.Generic.Res{fmpz}}, a::Integer, b)
        q = modulus(e)
        halfq = q >> 1
        en = lift(e)
        if en > halfq
            return oftype(e, q - multround(q - en, a, b))
        else
            return oftype(e, multround(en, a, b))
        end
    end
    function multround(e, a::Integer, b)
        oftype(e, map(NTT.coeffs(e)) do x
            multround(x, a, b)
        end)
    end

    Nemo.modulus(e::PrimeField) = GaloisFields.char(e)
    Nemo.lift(e::PrimeField) = e.n
    Nemo.lift(e::Nemo.nmod) = lift(Nemo.ZZ, e)

    function fqmod(e::Union{PrimeField, Nemo.nmod, AbstractAlgebra.Generic.Res{fmpz}}, nq::Integer)
        q = modulus(e)
        halfq = q >> 1
        e = lift(e)
        if e > halfq
            mod(e - q, nq)
        else
            mod(e, nq)
        end
    end

    divround(e::Integer, q::Integer) = div(e, q, RoundNearestTiesAway)
    divround(e::fmpz, q::Integer) = divround(BigInt(e), q)
    function divround(e::Union{PrimeField, Nemo.nmod, AbstractAlgebra.Generic.Res{fmpz}}, d::Integer)
        q = modulus(e)
        halfq = q >> 1
        en = lift(e)
        if en > halfq
            return oftype(e, q - divround(q - en, d))
        else
            return oftype(e, divround(en, d))
        end
    end

    function switchel(T, e)
        q = modulus(e)
        halfq = q >> 1
        diff = modulus(T) > q ? modulus(T) - q : q - modulus(T)
        en = lift(e)
        if (q < modulus(T))
            if en > halfq
                return T(en + diff)
            else
                return T(en)
            end
        else
            if en > halfq
                return T(en - diff)
            else
                return T(en)
            end
        end
    end

    function switch(ℛ, e)
        ℛ(map(NTT.coeffs(e)) do x
            switchel(coefftype(ℛ), x)
        end)
    end

    function *(c1::CipherText{T}, c2::CipherText{T}) where {T}
        params = c1.params
        @fields_as_locals params::BFVParams

        modswitch(c) = nntt_hint(switch(ℛbig, inntt_hint(c)))
        c1 = map(modswitch, c1.cs)
        c2 = map(modswitch, c2.cs)

        c = [zero(c1[1]) for i = 1:(length(c1) + length(c2) - 1)]
        for i = 1:length(c1), j = 1:length(c2)
            c[i+j-1] += c1[i] * c2[j]
        end

        c = map(c) do e
            switch(ℛ, multround(inntt_hint(e), modulus(base_ring(ℛplain)), modulus(coefftype(ℛ))))
        end

        CipherText(params, (c...,))
    end

    function decrypt(key::PrivKey, c::CipherText)
        @fields_as_locals key::PrivKey
        @fields_as_locals params::BFVParams

        b = nntt_hint(c[1])
        spow = s

        for i = 2:length(c)
            b += spow*nntt_hint(c[i])
            spow *= s
        end

        b = inntt_hint(b)
        ℛplain = plaintext_space(params)
        ℛplain(map(x->coefftype(ℛplain)(fqmod(divround(x, Δ), modulus(base_ring(ℛplain)))), NTT.coeffs(b)))
    end
    decrypt(key::KeyPair, plaintext) = decrypt(key.priv, plaintext)

    function keyswitch(ek::EvalKey, c::CipherText)
        @fields_as_locals ek::EvalKey
        @fields_as_locals params::BFVParams
        @assert length(c.cs) in (2,3)
        nwindows = ndigits(modulus(coefftype(ℛ)), base=2^relin_window)

        c1 = nntt_hint(c[1])
        c2 = length(c) == 2 ? zero(c[2]) : nntt_hint(c[2])

        cendcoeffs = NTT.coeffs(inntt_hint(c[end]))
        ds = map(cendcoeffs) do x
            digits(x.n, base=2^params.relin_window, pad=nwindows)
        end
        ps = map(1:nwindows) do i
            nntt_hint(ℛ([coefftype(ℛ)(ds[j][i]) for j in eachindex(cendcoeffs)]))
        end

        for i in eachindex(a)
            c2 += b[i] * ps[i]
            c1 += a[i] * ps[i]
        end

        CipherText(ek.params, (c1, c2))
    end
end
