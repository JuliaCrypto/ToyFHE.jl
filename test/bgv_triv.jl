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

plain = zero(plaintext_space(params))
plain[0] = 6

c = encrypt(kp, plain)
@test decrypt(kp, c)[0] == 6

let y = c*c
    @test decrypt(kp, y)[0] == 0x24
end
