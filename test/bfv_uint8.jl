using Nemo
using Hecke
using FHE
using FHE.BFV
using Primes

const n_plaintext_slots = 3
const plaintext_modulus = 256

p_fact = Primes.factor(plaintext_modulus)
@assert length(p_fact) == 1
pbase = first(keys(p_fact))

const ℤ = ZZ
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

function encode(factors, 𝔽p, data)
    f = prod(factors)
    ℤx = PolynomialRing(ZZ, "x")[1]
    crt = map(factors) do fᵢ
        # See page 3 of https://eprint.iacr.org/2011/133.pdf
        hᵢ = (f ÷ fᵢ)
        if parent(fᵢ) != 𝔽p
            fᵢl = lift(ℤx, fᵢ)
            hᵢl = lift(ℤx, hᵢ)
            fᵢp = 𝔽p(fᵢl)
            hᵢp = 𝔽p(hᵢl)
            s = invmod(hᵢp, fᵢp)
            t = (1-s*hᵢp) ÷ fᵢp
            ss, tt = lift(ℤx, s), lift(ℤx, t)
            pk = 1
            for i = 1:(Int(log(modulus(ℤpx), modulus(ℤplainx)))-1)
                # TODO: This is code is from Helib - it's some sort of hensel
                # lifting, but I don't really understand it. What is the
                # factorization we're actually lifting here?
                pk = pk *= modulus(ℤpx)
                d, r = divrem(1 - (ss*hᵢl + tt*fᵢl), ℤx(pk))
                @assert r == 0
                d = 𝔽p(d)
                s1 = (s * d) % fᵢp
                t1 = (d-s1*hᵢp) ÷ fᵢp
                ss += pk*lift(ℤx, s1)
                tt += pk*lift(ℤx, t1)
            end
            gᵢ = parent(fᵢ)(ss)
            @assert mod(gᵢ*hᵢ, fᵢ) == 1
        else
            gᵢ = invmod(hᵢ, fᵢ)
        end
        (hᵢ, gᵢ)
    end
    mapreduce(+, 1:length(factors)) do i
        (hᵢ, gᵢ) = crt[i]
        mulmod(data[i], gᵢ, factors[i]) * hᵢ
    end
end

decode(factors, encoded) = map(fᵢ->rem(encoded, fᵢ), factors)

const ℤplain = ResidueRing(ℤ, plaintext_modulus)
const ℤplainx, xplain = PolynomialRing(ℤplain, "x")

m = find_cyclotomic(n_plaintext_slots)
poly = cyclotomic(m, x)
factors = if plaintext_modulus != pbase
    map(ℤplainx, collect(keys(factor_mod_pk(poly, pbase, first(values(p_fact))))))
else
    collect(keys(Nemo.factor(ℤpx(poly)).fac))
end

# TODO: What are the security considerations for choosing the ciphertext modulus
q = nextprime(big(2)^51)
ℛ = ResidueRing(PolynomialRing(ResidueRing(ℤ, q), "x")[1], poly)
ℛbig = ResidueRing(PolynomialRing(ResidueRing(ℤ, nextprime(big(2)^111)), "x")[1], poly)

params = BFVParams(
    ℛ,
    ℛbig,
    BFV.plaintext_space(ℛ, plaintext_modulus),
    8/√(2π),
    div(q, plaintext_modulus)
)

kp = FHE.BFV.keygen(params)

ppoly = encode(factors, ℤpx, map(ℤplainx, [1, 2, 3, 4, 5, 6]))
epoly = params.ℛ(lift(PolynomialRing(ℤ,"x")[1], ppoly))

c1 = encrypt(kp, epoly)
let dec = decode(factors, decrypt(kp, c1*c1))
    @test dec == [1, 4, 9, 16, 25, 36]
end
