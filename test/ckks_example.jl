# This should match the SEAL CKKS example in
# https://github.com/microsoft/SEAL/blob/master/native/examples/4_ckks_basics.cpp

using ToyFHE
using ToyFHE.NTT
using ToyFHE.CKKS
using ToyFHE.BFV
using ToyFHE: coefftype
using OffsetArrays
using FFTW

# For now use the BFV parameter estimator and just copy it over.
# We need to figure out the right way to do this in CKKS
bfv_params = BFVParams(
    7, # plaintext modulus
    ; eval_mult_count = 3
)

params = CKKSParams(bfv_params.ℛ, bfv_params.ℛbig, bfv_params.relin_window,
    bfv_params.σ)

scale = 2^40
Tscale = FixedRational{coefftype(params.ℛ), scale}

points = LinRange(0.0, 1.0, 2048)
ipoints = OffsetArray(ifft([points; map(conj, reverse(points))]), 0:2*length(points)-1)
# Make it negacyclic
ipoints = ipoints .* [exp(2*k/(2*length(ipoints))*pi*im) for k in eachindex(ipoints)]

# Project to real coefficents
ipoints = map(ipoints) do p
    isapprox(imag(p), 0)
    real(p)
end

# Encode into the fixed fraction representation
encoded = map(Tscale, ipoints)

plaintext = params.ℛ(map(x->x.x, encoded))

kp = ToyFHE.CKKS.keygen(params)
c = encrypt(kp, ToyFHE.NTT.nntt(plaintext))

result = let dec = decrypt(kp, c)
    # Decode
    scaled = map(x->reinterpret(Tscale, x), ToyFHE.NTT.coeffs(dec))
    # Undo the root of unity premul
    multed = map(x->convert(Float64, x), scaled) .* [exp(-2*k/(2*length(scaled))*pi*im) for k in eachindex(scaled)]
    # FFT it and take only the non-conjugated coefficients
    fft(multed)[1:2048]
end
@show result
