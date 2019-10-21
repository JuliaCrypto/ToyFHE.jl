################################################################################
#                        BGV Scheme definition
################################################################################

struct BGVParams <: SHEShemeParams
    # The Cypertext ring over which operations are performed
    ℛ
    # The plaintext ring.
    ℛplain
    σ
end
scheme_name(p::Type{BGVParams}) = "BGV"
plaintext_modulus(p::BGVParams) = modulus(base_ring(p.ℛplain))

BGVParams(ring, p::Integer, σ) =
    BGVParams(ring, plaintext_space(ring, p), σ)

ℛ_plain(p::BGVParams) = p.ℛplain
ℛ_cipher(p::BGVParams) = p.ℛ

π⁻¹(params::BGVParams, plaintext) = convert(params.ℛ, params.ℛplain(plaintext))
function π(params::BGVParams, b)
    @fields_as_locals params::BGVParams
    ℛplain(map(x->coefftype(ℛplain)(convert(Integer, mod(SignedMod(x), plaintext_modulus(params)))), NTT.coeffs_primal(b)))
end

struct ShiftedDiscreteNormal
    p::Int
    dn::DiscreteNormal
end
Base.rand(rng::AbstractRNG, d::ShiftedDiscreteNormal) = d.p*rand(rng, d.dn)

𝒩(params::BGVParams) = RingSampler(params.ℛ, ShiftedDiscreteNormal(plaintext_modulus(params), DiscreteNormal(0, params.σ)))
𝒢(params::BGVParams) = RingSampler(params.ℛ, DiscreteNormal(0, params.σ))
