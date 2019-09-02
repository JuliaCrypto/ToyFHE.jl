using FHE
using FHE.BFV
using FHE.NTT
using GaloisFields
using Polynomials
using OffsetArrays
using Test
using Primes
using BitIntegers

Primes.isprime(x::Int256) = isprime(big(x))

params = BFVParams(
    65537, # plaintext modulus
    ; eval_mult_count = 1
)
kp = FHE.BFV.keygen(params)

F = GaloisField(65537)
ℛplain = LWERing{F, degree(params.ℛ)}(GaloisFields.minimal_primitive_root(F, 2degree(params.ℛ)))

plain = OffsetArray(zeros(F, degree(params.ℛ)), 0:degree(params.ℛ)-1)
plain[0] = 1
plain[1] = 1
encoded = inntt(RingCoeffs{ℛplain}(plain))

embedded = FHE.NTT.LWERingElement(
    FHE.NTT.RingCoeffs{params.ℛ}(map(encoded.coeffs) do x
        eltype(params.ℛ)(x.n)
    end)
)

plain2 = OffsetArray(zeros(F, degree(params.ℛ)), 0:degree(params.ℛ)-1)
fill!(plain2, F(10))
plain2[0] = 5
encoded2 = inntt(RingCoeffs{ℛplain}(plain2))

function switch_plain(x)
    FHE.NTT.RingCoeffs{ℛplain}(map(FHE.NTT.coeffs(x)) do x
        FHE.BFV.switch(eltype(ℛplain), x)
    end)
end

embedded2 = FHE.NTT.LWERingElement(
    FHE.NTT.RingCoeffs{params.ℛ}(map(encoded2.coeffs) do x
        eltype(params.ℛ)(x.n)
    end)
)

nntt(switch_plain(inntt(nntt(embedded)*nntt(embedded2))))

c1 = encrypt(kp, FHE.NTT.nntt(embedded))
c2 = encrypt(kp, FHE.NTT.nntt(embedded2))
y = c1*c2
xx = decrypt(kp, y)
nntt(RingCoeffs{ℛplain}(xx.p))
