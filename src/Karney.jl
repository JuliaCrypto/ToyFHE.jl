module Karney
    using Distributions
    using Random

    export DiscreteNormal, DiscreteGaussian

    # TODO: Algorithm H?
    struct BernoulliExpHalf
    end
    const exphalf = exp(-1/2)
    Random.rand(rng::AbstractRNG, ::BernoulliExpHalf) = rand(rng, Bernoulli(exphalf))


    struct HUntilFailure; end
    Random.rand(rng::AbstractRNG, ::HUntilFailure) = (n = 0; while rand(rng, BernoulliExpHalf()); n += 1; end; n)

    struct DiscreteNormal{S, R<:Real} <: DiscreteUnivariateDistribution
        T::S
        μ::R
        σ::R
    end
    const DiscreteGaussian = DiscreteNormal

    function DiscreteNormal(T::S, a::Real, b::Real) where {S}
        c = promote(a,b)
        DiscreteNormal{S,typeof(c[1])}(T, c...)
    end

    function _rand_karney(rng, d::DiscreteNormal)
        while true
            # Step D1
            k = rand(rng, HUntilFailure())

            # Step D2
            # TODO: Use Algorithm H instead?
            if !rand(rng, Bernoulli(exp(-k*(k+1)/2)))
                continue
            end

            # Step D3
            positive = rand(rng, Bool)

            # Step D4
            sμ = positive ? d.μ : -d.μ
            di₀ = d.σ * k + sμ
            i₀::Int64 = ceil(Int64, di₀)
            x₀ = i₀ - di₀
            j::Int64 = rand(rng, DiscreteUniform(0, ceil(Int64, d.σ) - 1))
            x = x₀ + j/d.σ

            # Step D5
            if (!isinteger(d.σ) && x >= 1)
                continue
            end

            # Step D6
            if k == 0 && x == 0 && !positive
                continue
            end

            # Step D7
            # TODO: Use Algorithm B instead?
            if !rand(rng, Bernoulli(exp(-1/2*x*(2k+x))))
                continue
            end

            # Step D8
            i::Int64 = i₀ + j
            positive || (i = -i)

            return d.T(i)
        end
    end

    Random.rand(rng::AbstractRNG, d::DiscreteNormal) = _rand_karney(rng, d)
    Random.rand(d::DiscreteNormal) = rand(Random.GLOBAL_RNG, d)
end
