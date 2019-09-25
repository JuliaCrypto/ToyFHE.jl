using FHE
using FHE.BGV
using GaloisFields
using Polynomials
using OffsetArrays
using Test

params = BGVParams(
    FHE.CryptParameters.palisade_parameters[1024],
    256,
    8/√(2π)
)
kp = FHE.BGV.keygen(params)

plain = OffsetArray(zeros(UInt8, degree(params.ℛ)), 0:degree(params.ℛ)-1)
plain[0] = 6
plain = FHE.NTT.FixedDegreePoly(plain)
encoded = FHE.NTT.NegacyclicRingElement(params.ℛ)(
    FHE.NTT.FixedDegreePoly(map(x->eltype(params.ℛ)(x), plain.p))
)
encoded = FHE.NTT.nntt(encoded)
@test FHE.NTT.inntt(encoded)[0] == 6

c = encrypt(kp, encoded)
@test decrypt(kp, c).p[0] == 6

let y = x*x
    @test decrypt(kp, y) == 6^2
end

using FHE
using FHE.BGV
using GaloisFields
using Polynomials
using OffsetArrays
using Test

params = BGVParams(
    FHE.CryptParameters.palisade_parameters[1024],
    256,
    8/√(2π)
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
