using ToyFHE
using ToyFHE.BFV
using OffsetArrays
using Test

params = BFVParams(
    7, # plaintext modulus
    ; eval_mult_count = 3
)

plain = OffsetArray(zeros(UInt8, ToyFHE.degree(params.ℛ)), 0:ToyFHE.degree(params.ℛ)-1)
plain[0] = 2
encoded = ToyFHE.NTT.nntt(ToyFHE.NTT.NegacyclicRingElement(
    ToyFHE.NTT.RingCoeffs{params.ℛ}(map(x->eltype(params.ℛ)(x), plain))
))

# Generate a keypair and a corresponding relin key
kp1 = ToyFHE.BFV.keygen(params)
ek = ToyFHE.BFV.keygen(BFV.EvalKey, kp1.priv)

c1 = encrypt(kp1, encoded)

@test decrypt(kp1, c1).p[0] == 2
let c1squared = c1*c1
    @test decrypt(kp1, c1squared).p[0] == 4
    cswitch = ToyFHE.BFV.keyswitch(ek, c1squared)
    @test length(cswitch.cs) == 2
    @test decrypt(kp1, cswitch).p[0] == 4
    @test decrypt(kp1, cswitch*c1).p[0] == 1
end
