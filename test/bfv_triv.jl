using FHE
using FHE.BFV
using GaloisFields
using Polynomials
using OffsetArrays
using Test

params = BFVParams(
    256, # plaintext modulus
    ; eval_mult_count = 1
)
kp = FHE.BFV.keygen(params)


for i = 0:degree(params.ℛ)-1
    @show i
    plain = OffsetArray(zeros(UInt8, degree(params.ℛ)), 0:degree(params.ℛ)-1)
    plain[i] = 6
    encoded = FHE.NTT.LWERingElement(
        FHE.NTT.RingCoeffs{params.ℛ}(map(x->eltype(params.ℛ)(x), plain))
    )
    #@test FHE.NTT.inntt(encoded)[0] == 6

    c = encrypt(kp, FHE.NTT.nntt(encoded))
    @test decrypt(kp, c).p[i] == 6

    let y = c*c
        @test decrypt(kp, y).p[mod(2*i,degree(params.ℛ))] == 0x24
    end
end
