using Nemo
using Hecke
using ToyFHE
using Primes
using Test

const n_plaintext_slots = 3
const plaintext_modulus = 256

p_fact = Primes.factor(plaintext_modulus)
@assert length(p_fact) == 1
pbase = first(keys(p_fact))

const ℤ = Nemo.ZZ
const ℤx, x = PolynomialRing(ℤ, "x")
const ℤp = ResidueRing(ℤ, pbase)
const ℤpx, xp = PolynomialRing(ℤp, "xp")

# Find a prime cyclotomic that supports at least `n_plaintext_slots` slots
function find_cyclotomic(nslots)
    # TODO: What are the security considerations on ring dimensions
    for m in primes(2^20)
        poly = cyclotomic(m, x)
        fact = Nemo.factor(ℤpx(poly))
        if length(fact) > nslots
            return m
        end
    end
end

m = find_cyclotomic(n_plaintext_slots)
poly = cyclotomic(m, x)

# TODO: What are the security considerations for choosing the ciphertext modulus
q = nextprime(big(2)^51)
ℛ = ResidueRing(PolynomialRing(ResidueRing(ℤ, q), "x")[1], poly)
ℛbig = ResidueRing(PolynomialRing(ResidueRing(ℤ, nextprime(big(2)^111)), "x")[1], poly)

params = BFVParams(
    ℛ,
    ℛbig,
    plaintext_space(ℛ, plaintext_modulus),
    1,
    8/√(2π),
    div(q, plaintext_modulus)
)

kp = keygen(params)

plain = PolyCRTEncoding(zero(plaintext_space(params)))
plain[:] .= [1, 2, 3, 4, 5, 6]

c1 = encrypt(kp, plain)
let dec = PolyCRTEncoding(decrypt(kp, c1*c1))
    @test dec == [1, 4, 9, 16, 25, 36]
end
