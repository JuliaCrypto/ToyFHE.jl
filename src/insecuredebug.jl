export InsecureDebug

struct InsecureDebug{P<:SHEShemeParams} <: PassthroughParams{P}
    params::P
end

scheme_name(::Type{InsecureDebug{P}}) where P = string("Insecure ", scheme_name(P))

struct ZeroSampler <: Random.Sampler{Zero}; end
Random.rand(rng::AbstractRNG, zs::ZeroSampler) = Zero()
Base.:+(a::NTT.RingElement, b::Zero) = a
Base.:+(a::Zero, b::NTT.RingElement) = a

𝒩(params::CKKSParams) = ZeroSampler()
