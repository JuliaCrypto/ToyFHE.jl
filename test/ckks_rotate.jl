using ToyFHE
using OffsetArrays
using Test
using Primes
using GaloisFields
using Random

N = 2^4
q₀ = nextprime(Int128(2)^40 + 1, 1; interval=2N)
ps = nextprime(q₀ + 2N, 1; interval=2N)

ℛ = let CT = ToyFHE.CRTEncoded{2, Tuple{GaloisField.((q₀,ps))...}}
    ζ₂n = GaloisFields.minimal_primitive_root(CT, 2N)
    ToyFHE.NegacyclicRing{CT, N}(ζ₂n)
end

scale = 2^60
Tscale = FixedRational{scale}

plain = CKKSEncoding{Tscale}(zero(ℛ))
plain .= OffsetArray(1:div(N,2), 0:div(N,2)-1)
plain[0] += 1im

re = convert(ToyFHE.NTT.RingElement, plain)
@test CKKSEncoding{Tscale}(ToyFHE.NTT.apply_galois_element(re, 3)) ≈ circshift(plain, -1)

# The same, but encrypted
ckks_params = CKKSParams(ℛ, 1, 3.2)
kp = keygen(ckks_params)

function ToyFHE.NTT.apply_galois_element(c::CipherText{Enc}, galois_element) where {Enc}
    CipherText{Enc}(c.params, map(re->ToyFHE.NTT.apply_galois_element(re, galois_element), c.cs))
end

rt = let c = encrypt(kp, plain)
    c′ = ToyFHE.NTT.apply_galois_element(c, 3)
    ek = ToyFHE.make_eval_key(Random.GLOBAL_RNG, ToyFHE.NTT.apply_galois_element(kp.priv.secret, 3)=>kp.priv)
    c′ = ToyFHE.keyswitch(ek, c′)
    decrypt(kp, c′)
end
@test rt ≈ circshift(plain, -1)

let gk = ToyFHE.keygen(GaloisKey, kp.priv; steps=1)
    @test decrypt(kp, ToyFHE.rotate(gk, encrypt(kp, plain))) ≈ circshift(plain, 1)
end
