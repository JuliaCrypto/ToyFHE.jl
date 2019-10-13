using ToyFHE
using ToyFHE.BFV
using OffsetArrays
using Test

params = BFVParams(
    256, # plaintext modulus
    # We're only doing one multiplication, but let's give ourselves some room,
    # since this test is only testing the basic correctness of the implementation.
    # We need to validate the parameter selection elsewhere
    ; eval_mult_count = 2
)
kp = ToyFHE.BFV.keygen(params)

plain = OffsetArray(zeros(UInt8, ToyFHE.degree(params.ℛ)), 0:ToyFHE.degree(params.ℛ)-1)
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
