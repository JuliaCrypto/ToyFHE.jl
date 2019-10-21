using ToyFHE
using OffsetArrays
using Test
using GaloisFields
using Primes
using StructArrays

n = 2048
p1 = nextprime(Int128(2)^50 + 1, 1; interval=2n)
p2 = nextprime(p1+2n, interval=2n)

ℛ = let CT = ToyFHE.CRTEncoded{Int64, 2, Tuple{GaloisField(p1), GaloisField(p2)}}
    ζ₄₀₉₆ = GaloisFields.minimal_primitive_root(CT, 2n)
    ToyFHE.NegacyclicRing{CT, n}(ζ₄₀₉₆)
end

p3 = nextprime(p2+2n, 1; interval=2n)
p4 = nextprime(p3+2n, 1; interval=2n)
p5 = nextprime(p4+2n, 1; interval=2n)
p6 = nextprime(p5+2n, 1; interval=2n)

ℛbig = let CT = ToyFHE.CRTEncoded{Int64, 4, Tuple{GaloisField(p3), GaloisField(p4), GaloisField(p5), GaloisField(p6)}}
    ζ₄₀₉₆ = GaloisFields.minimal_primitive_root(CT, 2n)
    ToyFHE.NegacyclicRing{CT, n}(ζ₄₀₉₆)
end

ℛplain = plaintext_space(ℛ, 53)

params = BFVParams(
    ℛ,
    ℛbig,
    ℛplain,
    1,
    3.2,
    div(ToyFHE.NTT.modulus(ToyFHE.NTT.coefftype(ℛ)), ToyFHE.NTT.modulus(ToyFHE.NTT.coefftype(ℛplain)))
)
kp = keygen(params)

plain = zero(plaintext_space(params))
plain[0] = 6

c = encrypt(kp, plain)
@test decrypt(kp, c)[0] == 6

let y = c*c
    @test decrypt(kp, y)[0] == 0x24
end
