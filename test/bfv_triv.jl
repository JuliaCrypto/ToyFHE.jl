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
encoded = FHE.NTT.NegacyclicRingElement(
    FHE.NTT.RingCoeffs{params.ℛ}(map(x->eltype(params.ℛ)(x), plain))
)
#@test FHE.NTT.inntt(encoded)[0] == 6

c = encrypt(kp, FHE.NTT.nntt(encoded))
@test decrypt(kp, c).p[0] == 6

let y = c*c
    @test decrypt(kp, y).p[0] == 0x24
end
