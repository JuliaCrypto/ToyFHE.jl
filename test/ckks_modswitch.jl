using ToyFHE
using OffsetArrays
using Test
using Primes
using GaloisFields

N = 2^5
q₀ = nextprime(Int128(2)^40 + 1, 1; interval=2N)
q₁ = nextprime(q₀ + 2N, 1; interval=2N)
ps = nextprime(q₁ + 2N, 1; interval=2N)

ℛ = let CT = ToyFHE.CRTEncoded{3, Tuple{GaloisField.((q₀,q₁,ps))...}}
    ζ₂n = GaloisFields.minimal_primitive_root(CT, 2N)
    ToyFHE.NegacyclicRing{CT, N}(ζ₂n)
end

scale = 2^60
Tscale = FixedRational{scale}

plain = CKKSEncoding{Tscale}(zero(ℛ))
plain .= 2

Tscale2 = FixedRational{scale/ps}
@test isapprox(CKKSEncoding{Tscale2}(ToyFHE.modswitch(convert(ToyFHE.NTT.RingElement, plain)))[0], 2.0; atol=10^-5)

# The same but with encryption noise
ckks_params = CKKSParams(ℛ, 1, 3.2)
kp = keygen(ckks_params)

switched = let c = ToyFHE.modswitch(encrypt(kp, plain))
    decrypt(kp, c)
end
@test switched ≈ plain atol=10^-3
