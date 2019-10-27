################################################################################
#                        BFV Scheme definition
################################################################################

struct BFVParams <: SHEShemeParams
    # The Cypertext ring over which operations are performed
    ℛ
    # The big ring used during multiplication
    ℛbig
    # The plaintext ring.
    ℛplain
    relin_window
    σ
    Δ
end
scheme_name(p::Type{BFVParams}) = "BFV"

ℛ_plain(p::BFVParams) = p.ℛplain
ℛ_cipher(p::BFVParams) = p.ℛ

function π⁻¹(params::BFVParams, plaintext)
    @fields_as_locals params::BFVParams
    params.Δ * oftype(zero(params.ℛ), ℛplain(plaintext))
end

function π(params::BFVParams, b)
    @fields_as_locals params::BFVParams
    ℛplain(map(x->coefftype(ℛplain)(convert(Integer, mod(divround(x, Δ), modulus(base_ring(ℛplain))))), NTT.coeffs_primal(b)))
end

𝒩(params::BFVParams) = RingSampler(params.ℛ, DiscreteNormal(0, params.σ))
𝒢(params::BFVParams) = RingSampler(params.ℛ, DiscreteNormal(0, params.σ))

mul_expand(params::BFVParams, c::CipherText) = map(c->switch(params.ℛbig, c), c.cs)
function mul_contract(params::BFVParams, c)
    @fields_as_locals params::BFVParams
    map(c) do e
        switch(ℛ, multround(e, modulus(base_ring(ℛplain)), modulus(coefftype(ℛ))))
    end
end

################################################################################
#                 BFV Noise Modeling / Parameter Generation
################################################################################

# Matches parameter generation in PALISADE
function BFVParams(p, σ=8/√(2pi), α=9, r=1; eval_mult_count = 0, security = CryptParameters.HEStd_128_classic, relin_window=1)
    @assert r >= 1
    Berr = σ*√(α)
    Bkey = Berr
    δ(n) = 2*√(n)
    Vnorm(n) = Berr * (1 + 2*δ(n)*Bkey)

    function nRLWE(q)
        if isa(security, CryptParameters.StdSecurity)
            CryptParameters.std_ring_dim(CryptParameters.HEStd_error, security, ceil(log2(q)))
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

    bits = ceil(Integer, log2(q))+1
    T = bits > 256 ? Int512 : bits > 128 ? Int256 : Int128
    qPrime = nextprime(T(2)^(ceil(Int, log2(q))+1) + 1, 1; interval=2n)
    largebits = 2*ceil(Int, log2(q)) + ceil(Int, log2(p)) + 3
    Tlarge = largebits > 256 ? Int512 : largebits > 128 ? Int256 : Int128
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


"""
Compute the *invariant noise budget*, defined by:

        -log2(2‖v‖) = log2(q) - log2(q‖v‖) - 1.

If this quantity is >0, the ciphertext is expected to decrypt correctly with
high probability.

This notion of noise was first introduced by the SEAL HE library. See [CLP19]
for details.

[CLP19] Anamaria Costache, Kim Laine, and Rachel Player
        "Homomorphic noise growth in practice: comparing BGV and FV"
        https://eprint.iacr.org/2019/493.pdf
"""
function invariant_noise_budget(pk::PrivKey{BFVParams}, c::CipherText{<:Any, BFVParams})
    @fields_as_locals pk::PrivKey
    @fields_as_locals params::BFVParams

    b = c[1]
    spow = secret

    for i = 2:length(c)
        b += spow*c[i]
        spow *= secret
    end

    ℛplain = plaintext_space(params)

    function birem(x)
        r = rem(x, Δ)
        if r > div(Δ, 2)
            return Δ - r
        else
            return r
        end
    end

    # -log2(2‖v‖) = log(q) - log(t) - 1 - max_i log2(Δ |v_i|)
    log2(modulus(coefftype(ℛ))) - log2(modulus(coefftype(ℛplain))) - 1 -
        maximum(log2(birem(c.n)) for c in NTT.coeffs_primal(b))
end
invariant_noise_budget(kp::KeyPair, c::CipherText) =
    invariant_noise_budget(kp.priv, c)

export invariant_noise_budget

################################################################################
#                 BFV Computational utilities
################################################################################

function multround(e::SignedMod, a::Integer, b::Integer)
    div(e * a, b, RoundNearestTiesAway)
end
function multround(e::BigInt, a::Integer, b::Integer)
    div(e * a, b, RoundNearestTiesAway)
end
multround(e, a::Integer, b::fmpz) = multround(e, a, BigInt(b))
multround(e::fmpz, a::Integer, b::Integer) = multround(BigInt(e), a, b)
multround(e::fmpz, a::Integer, b::fmpz) = multround(BigInt(e), a, BigInt(b))

function multround(e, a::Integer, b)
    oftype(e, broadcast(NTT.coeffs_primal(e)) do x
        if isa(x, AbstractAlgebra.Generic.Res{fmpz})
            multround(BigInt(Nemo.lift(x)), a, b)
        else
            multround(SignedMod(x), a, b).x
        end
    end)
end

Nemo.modulus(e::GaloisFields.PrimeField) = GaloisFields.char(e)
Nemo.lift(e::GaloisFields.PrimeField) = e.n
Nemo.lift(e::Nemo.nmod) = lift(Nemo.ZZ, e)

divround(e::Integer, q::Integer) = div(e, q, RoundNearestTiesAway)
divround(e::fmpz, q::Integer) = divround(BigInt(e), q)
function divround(e, d::Integer)
    div(SignedMod(e), d, RoundNearestTiesAway)
end

function switchel(T, e)
    q = modulus(e)
    halfq = q >> 1
    diff = modulus(T) > q ? modulus(T) - q : q - modulus(T)
    en = convert(Integer, e)
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
    ℛ(broadcast(NTT.coeffs_primal(e)) do x
        switchel(coefftype(ℛ), x)
    end)
end
