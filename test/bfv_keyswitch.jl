using ToyFHE
using ToyFHE.BFV
using OffsetArrays
using Test

params = BFVParams(
    7, # plaintext modulus
    ; eval_mult_count = 3
)

plain = zero(plaintext_space(params))
plain[0] = 2

# Generate a keypair and a corresponding relin key
kp1 = ToyFHE.BFV.keygen(params)
ek = ToyFHE.BFV.keygen(BFV.EvalKey, kp1.priv)

c1 = encrypt(kp1, plain)

@test decrypt(kp1, c1)[0] == 2
let c1squared = c1*c1
    @test decrypt(kp1, c1squared)[0] == 4
    cswitch = ToyFHE.BFV.keyswitch(ek, c1squared)
    @test length(cswitch.cs) == 2
    @test decrypt(kp1, cswitch)[0] == 4
    @test decrypt(kp1, cswitch*c1)[0] == 1
end
