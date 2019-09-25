using ToyFHE.BGV
using GaloisFields
using Polynomials
using OffsetArrays
using Test

params = BGVParams(
    ToyFHE.CryptParameters.palisade_parameters[1024],
    256,
    8/√(2π)
)
kp = ToyFHE.BGV.keygen(params)

plain = OffsetArray(zeros(UInt8, degree(params.ℛ)), 0:degree(params.ℛ)-1)
plain[0] = 6
encoded = ToyFHE.NTT.NegacyclicRingElement(
    ToyFHE.NTT.RingCoeffs{params.ℛ}(map(x->eltype(params.ℛ)(x), plain))
)
#@test FHE.NTT.inntt(encoded)[0] == 6

c = encrypt(kp, ToyFHE.NTT.nntt(encoded))
@test decrypt(kp, c).p[0] == 6

let y = c*c
    @test decrypt(kp, y).p[0] == 0x24
end
