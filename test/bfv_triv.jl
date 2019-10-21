using ToyFHE
using OffsetArrays
using Test

params = BFVParams(
    53, # plaintext modulus
    # We're only doing one multiplication, but let's give ourselves some room,
    # since this test is only testing the basic correctness of the implementation.
    # We need to validate the parameter selection elsewhere
    ; eval_mult_count = 2
)
kp = keygen(params)

plain = zero(plaintext_space(params))
plain[0] = 6

c = encrypt(kp, plain)
@test decrypt(kp, c)[0] == 6

let y = c*c
    @test decrypt(kp, y)[0] == 0x24
end
