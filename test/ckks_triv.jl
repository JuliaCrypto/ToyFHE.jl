using ToyFHE
using OffsetArrays
using Test

# For now use the BFV parameter estimator and just copy it over.
# We need to figure out the right way to do this in CKKS
bfv_params = BFVParams(
    7, # plaintext modulus
    ; eval_mult_count = 3
)

params = CKKSParams(bfv_params.ℛ, bfv_params.ℛbig, bfv_params.relin_window,
    bfv_params.σ)

scale = 2^40
Tscale = FixedRational{scale}
Tscalesq = FixedRational{Int128(scale)^2}

plain = CKKSEncoding{Tscale}(zero(plaintext_space(params)))
plain .= OffsetArray(LinRange(0.0, 1.0, 2048), 0:2047)

re = convert(ToyFHE.NTT.RingElement, plain)

# TODO: Can we do better here by being more careful with the FFT?
@test real.(collect(CKKSEncoding{Tscalesq}(re*re))) ≈ LinRange(0.0, 1.0, 2048).^2 atol=10^-4

# First test the encoder in isolation without the extra noise from encryption

kp = keygen(params)
c = encrypt(kp, plain)

@test real.(collect(decrypt(kp, c))) ≈ LinRange(0.0, 1.0, 2048) atol=10^-4
@test real.(collect(decrypt(kp, c*c))) ≈ LinRange(0.0, 1.0, 2048).^2 atol=10^-4
