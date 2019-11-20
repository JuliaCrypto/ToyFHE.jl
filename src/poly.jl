using OffsetArrays

coefftype(ring::Type) = eltype(ring)
modulus(ring) = Nemo.modulus(ring)
degree(ring) = Nemo.degree(ring)

struct RingSampler{Ring} <: Random.Sampler{Ring}
    ring::Ring
    coeff_distribution::Any #DiscreteUnivariateDistribution
end
Random.gentype(r::RingSampler) = Any
ring(r::RingSampler) = r.ring

function sample_ring_array(rng, ℛ, coeff_distribution)
    [coefftype(ℛ)(rand(rng, coeff_distribution)) for _ in 1:degree(modulus(ℛ))]
end

function Random.rand(rng::Random.AbstractRNG, r::RingSampler)
    ℛ = ring(r)
    ℛ(OffsetArray(
        sample_ring_array(rng, ℛ, r.coeff_distribution),
        0:degree(modulus(ℛ))-1))
end

struct Zero; end
