using ToyFHE
using ToyFHE: degree
using ToyFHE.BFV
using ToyFHE.NTT
using GaloisFields
using OffsetArrays
using Test
using Primes
using BitIntegers

Primes.isprime(x::Int256) = isprime(big(x))

params = BFVParams(
    65537, # plaintext modulus
    ; eval_mult_count = 1
)
kp = ToyFHE.BFV.keygen(params)

F = GaloisField(65537)
ℛplain = ToyFHE.BFV.plaintext_space(params)

plain = OffsetArray(zeros(F, degree(params.ℛ)), 0:degree(params.ℛ)-1)
plain[0] = 1
plain[1] = 1
encoded = inntt(RingCoeffs{ℛplain}(plain))

embedded = ToyFHE.NTT.NegacyclicRingElement(
    ToyFHE.NTT.RingCoeffs{params.ℛ}(map(encoded.coeffs) do x
        eltype(params.ℛ)(x.n)
    end)
)

plain2 = OffsetArray(zeros(F, degree(params.ℛ)), 0:degree(params.ℛ)-1)
fill!(plain2, F(10))
plain2[0] = 5
encoded2 = inntt(RingCoeffs{ℛplain}(plain2))

function switch_plain(x)
    ToyFHE.NTT.RingCoeffs{ℛplain}(map(ToyFHE.NTT.coeffs(x)) do x
        ToyFHE.BFV.switch(eltype(ℛplain), x)
    end)
end

embedded2 = ToyFHE.NTT.NegacyclicRingElement(
    ToyFHE.NTT.RingCoeffs{params.ℛ}(map(encoded2.coeffs) do x
        eltype(params.ℛ)(x.n)
    end)
)

c1 = encrypt(kp, ToyFHE.NTT.nntt(embedded))
c2 = encrypt(kp, ToyFHE.NTT.nntt(embedded2))
y = c1*c2
xx = decrypt(kp, y)
let data = nntt(xx).data
    @test data[0] == 5
    @test data[1] == 10
    @test all(iszero, data[2:end])
end
