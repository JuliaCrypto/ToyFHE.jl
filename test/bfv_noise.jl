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

kp1 = ToyFHE.BFV.keygen(params)
ek = ToyFHE.BFV.keygen(BFV.EvalKey, kp1.priv)

c1 = encrypt(kp1, plain)
b1 = BFV.invariant_noise_budget(kp1.priv, c1)
c1squared = c1*c1
b2 = BFV.invariant_noise_budget(kp1.priv, c1squared)
@test b2 < b1
cswitch1 = ToyFHE.BFV.keyswitch(ek, c1squared)
bswitch1 = BFV.invariant_noise_budget(kp1.priv, cswitch1)
cswitchmul = cswitch1*c1
bswitchmul = BFV.invariant_noise_budget(kp1.priv, cswitchmul)
@test bswitchmul < bswitch1 < b1
cswitch2 = ToyFHE.BFV.keyswitch(ek, cswitchmul)
bswitch2 = BFV.invariant_noise_budget(kp1.priv, cswitch2)
cswitchmul2 = cswitch2*c1
bswitchmul2 = BFV.invariant_noise_budget(kp1.priv, cswitchmul2)
@test bswitchmul2 < bswitch2 < bswitch1

# This is mostly a test that our heuristics are doing a decent job of chosing
# parameters.
@test 1 < bswitchmul2 < 10
