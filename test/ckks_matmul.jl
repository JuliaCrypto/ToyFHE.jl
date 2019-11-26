using ToyFHE
using OffsetArrays
using Test
using Primes
using GaloisFields
using LinearAlgebra

N = 2^5
q₀ = nextprime(Int128(2)^40 + 1, 1; interval=2N)
q₁ = nextprime(q₀ + 2N, 1; interval=2N)
ps = nextprime(q₁ + 2N, 1; interval=2N)

ℛ = let CT = ToyFHE.CRTEncoded{3, Tuple{GaloisField.((q₀,q₁,ps))...}}
    ζ₂n = GaloisFields.minimal_primitive_root(CT, 2N)
    ToyFHE.NegacyclicRing{CT, N}(ζ₂n)
end

scale = 2^40
Tscale = FixedRational{scale}

plain = CKKSEncoding{Tscale}(zero(ℛ))
plain .= OffsetArray(1:div(N,2), 0:div(N,2)-1)

plain_matrix = ones(Float32, 4, 4)

ckks_params = CKKSParams(ℛ, 1, 3.2)
kp = keygen(ckks_params)

c = encrypt(kp, plain)

gk = keygen(GaloisKey, kp.priv; steps=4)
ek = keygen(EvalMultKey, kp.priv)

function encrypted_matmul(gk, weights, x::ToyFHE.CipherText)
    result = repeat(diag(weights), size(weights, 2)).*x
    rotated = x
    for k = 2:size(weights, 2)
        rotated = ToyFHE.rotate(gk, rotated)
        result += repeat(diag(circshift(weights, (0,(k-1)))), size(weights, 2)) .* rotated
    end
    result
end

@test reshape(collect(decrypt(kp, encrypted_matmul(gk, plain_matrix, c))), (4,4)) ≈ (plain_matrix*reshape(collect(plain), 4,4)')' atol=10^-5
