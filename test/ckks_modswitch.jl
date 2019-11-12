using ToyFHE
using OffsetArrays
using Test
using Primes
using GaloisFields

N = 2^5
q₀ = nextprime(Int128(2)^40 + 1, 1; interval=2N)
ps = nextprime(q₀ + 2N, 1; interval=2N)

ℛ = let CT = ToyFHE.CRTEncoded{2, Tuple{GaloisField.((q₀,ps))...}}
    ζ₂n = GaloisFields.minimal_primitive_root(CT, 2N)
    ToyFHE.NegacyclicRing{CT, N}(ζ₂n)
end

scale = 2^60
Tscale = FixedRational{coefftype(ℛ), scale}

plain = CKKSEncoding{Tscale}(zero(ℛ))
plain .= 2

Tscale2 = FixedRational{ToyFHE.drop_last(coefftype(ℛ)), scale/ps}
@test isapprox(CKKSEncoding{Tscale2}(ToyFHE.modswitch(convert(ToyFHE.NTT.RingElement, plain)))[0], 2.0; atol=10^-5)

# The same but with encryption noise
ckks_params = CKKSParams(ℛ, ℛ, 1, 3.2)
kp = keygen(ckks_params)

switched = let c = ToyFHE.modswitch(encrypt(kp, plain))
    pk′ = PrivKey(kp.priv.params, ToyFHE.modswitch_drop(kp.priv.secret))
    decrypt(pk′, c)
end
@test switched ≈ plain atol=10^-3
