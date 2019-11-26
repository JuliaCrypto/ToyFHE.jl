using ToyFHE
using OffsetArrays
using Test
using Primes
using GaloisFields
using LinearAlgebra
using Random


N = 2^5
q₀ = nextprime(Int128(2)^40 + 1, 1; interval=2N)
q₁ = nextprime(q₀ + 2N, 1; interval=2N)
ps = nextprime(q₁ + 2N, 1; interval=2N)

ℛ = let CT = ToyFHE.CRTEncoded{3, Tuple{GaloisField.((q₀,q₁,ps))...}}
    ζ₂n = GaloisFields.minimal_primitive_root(CT, 2N)
    ToyFHE.NegacyclicRing{CT, N}(ζ₂n)
end

ckks_params = ModulusRaised(CKKSParams(ℛ, 0, 3.2))
kp = keygen(ckks_params)

scale = 2^40
Tscale = FixedRational{scale}

plain = CKKSEncoding{Tscale}(zero(ℛ_cipher(ckks_params)))
plain .= OffsetArray(1:div(N,2), 0:div(N,2)-1)

c = encrypt(kp, plain)
@test decrypt(kp, keyswitch(ToyFHE.make_eval_key(Random.GLOBAL_RNG, kp.priv.secret=>kp.priv), c)) ≈ plain atol=10^-8
