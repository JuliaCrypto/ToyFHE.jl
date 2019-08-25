using FHE
using FHE.BFV
using GaloisFields
using Polynomials
using OffsetArrays
using Test

params = BFVParams(
    256, # plaintext modulus
    ; eval_mult_count = 1
)
kp = FHE.BFV.keygen(params)

plain = OffsetArray(zeros(UInt8, degree(params.ℛ)), 0:degree(params.ℛ)-1)
plain[0] = 6
plain = FHE.NTT.FixedDegreePoly(plain)
encoded = FHE.NTT.LWERingElement(params.ℛ)(
    FHE.NTT.FixedDegreePoly(map(x->eltype(params.ℛ)(x), plain.p))
)
encoded = FHE.NTT.nntt(encoded)
#@test FHE.NTT.inntt(encoded)[0] == 6

c = encrypt(kp, encoded)
@test decrypt(kp, c).p[0] == 6

let y = c*c
    @test decrypt(kp, y) == 6^2
end

