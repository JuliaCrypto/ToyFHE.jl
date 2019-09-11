using OffsetArrays

coefftype(ring::Type) = eltype(poly)

struct RingSampler{Ring} <: Random.Sampler{Ring}
    ring::Ring
    coeff_distribution::Any #DiscreteUnivariateDistribution
end
ring(r::RingSampler) = r.ring

function Random.rand(rng::Random.AbstractRNG, r::RingSampler)
    ℛ = ring(r)
    ℛ(OffsetArray(
        [rand(rng, r.coeff_distribution) for _ in 1:degree(modulus(ℛ))],
        0:degree(modulus(ℛ))-1))
end
