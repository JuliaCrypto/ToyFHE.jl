using Nemo
using Hecke
using ToyFHE
using Primes
using Test
using GaloisFields

const ℤ = ZZ
const ℤx, x = ZZ["x"]
const ℤpx, xp = ResidueRing(ZZ, 2)["xp"]
ℛ = ResidueRing(ℤpx, cyclotomic(7, x))

plain = PolyCRTEncoding(zero(ℛ))
β = GaloisFields.gen(eltype(plain))

plain[1] = β+1
plain[2] = β^2+1

p = convert(ℛ, plain)
enc = PolyCRTEncoding(p)
@test enc[1] == β+1
@test enc[2] == β^2+1
