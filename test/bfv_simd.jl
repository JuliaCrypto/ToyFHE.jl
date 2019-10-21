using ToyFHE
using ToyFHE: degree
using ToyFHE.NTT
using GaloisFields
using OffsetArrays
using Test
using Primes
using BitIntegers

params = BFVParams(
    65537, # plaintext modulus
    ; eval_mult_count = 1
)
kp = keygen(params)

plain = SlotEncoding(zero(plaintext_space(params)))
plain[0] = 1
plain[1] = 1

plain2 = SlotEncoding(zero(plaintext_space(params)))
plain2[:] .= 10
plain2[0] = 5

c1 = encrypt(kp, plain)
c2 = encrypt(kp, plain2)
y = c1*c2
let data = SlotEncoding(decrypt(kp, y))
    @test data[0] == 5
    @test data[1] == 10
    @test all(iszero, data[2:end])
end
