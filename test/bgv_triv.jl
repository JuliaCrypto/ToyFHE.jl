using ToyFHE
using ToyFHE.BGV
using GaloisFields
using OffsetArrays
using Test

params = BGVParams(
    ToyFHE.CryptParameters.palisade_parameters[4096],
    256,
    8/√(2π)
)
kp = ToyFHE.BGV.keygen(params)

plain = OffsetArray(zeros(UInt8, ToyFHE.degree(params.ℛ)), 0:ToyFHE.degree(params.ℛ)-1)
plain[0] = 6
encoded = ToyFHE.NTT.nntt(ToyFHE.NTT.NegacyclicRingElement(
    ToyFHE.NTT.RingCoeffs{params.ℛ}(map(x->eltype(params.ℛ)(x), plain))
))

c = encrypt(kp, encoded)
@test decrypt(kp, c).p[0] == 6

let y = c*c
    @test decrypt(kp, y).p[0] == 0x24
end
