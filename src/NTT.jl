module NTT

using GaloisFields
using Polynomials
using OffsetArrays
using Distributions
using Random
using AutoHashEquals

import Base: *, +, -

import GaloisFields: PrimeField

export LWERing, RingSampler, nntt, inntt, FixedDegreePoly,
    LWERingElement, LWERingDualElement

@auto_hash_equals struct FixedDegreePoly{N, T}
    p::OffsetVector{T}
end
function FixedDegreePoly(p::OffsetVector)
    @assert first(axes(p)[1]) == 0
    FixedDegreePoly{length(p), eltype(p)}(p)
end
Base.zero(::Type{FixedDegreePoly{N, T}}) where {N, T} =
    FixedDegreePoly(OffsetArray(zeros(T, N),0:N-1))
Polynomials.degree(p::FixedDegreePoly{N}) where {N} = N

"""
Represents the ring 𝔽q[x]/(xⁿ+1).
"""
struct LWERing{Field <: PrimeField, N}
    # 2N'th primitive root of unity in Field
    ψ::Field
end
Polynomials.degree(ℛ::LWERing{F,N}) where {F,N} = N
Base.eltype(ℛ::LWERing{F,N}) where {F,N} = F

"""
Represents an element of 𝔽q[x]/(xⁿ+1).
"""
@auto_hash_equals struct LWERingElement{ℛ #= ::LWERing{Field} =#, Field <: PrimeField,  N}
    p::FixedDegreePoly{N, Field}
end
LWERingElement{ℛ,Field,N}(coeffs::AbstractVector) where {ℛ,Field <: PrimeField,  N} = LWERingElement{ℛ,Field,N}(FixedDegreePoly(coeffs))
coeffs(e::LWERingElement) = e.p.p
LWERingElement(ℛ::LWERing) = LWERingElement{ℛ, eltype(ℛ), degree(ℛ)}
Base.zero(::Type{LWERingElement{ℛ,Field,N}}) where {ℛ,Field,N} =
    LWERingElement{ℛ,Field,N}(zero(FixedDegreePoly{N, Field}))

"""
Represents an ntt-dual element of 𝔽q[x]/(xⁿ+1).
"""
@auto_hash_equals struct LWERingDualElement{ ℛ #= ::LWERing{Field} =#, Field <: PrimeField}
    data::OffsetVector{Field}
end
Base.zero(::Type{LWERingDualElement{ℛ,Field}}) where {ℛ,Field} =
    LWERingDualElement{ℛ,Field}(OffsetArray(zeros(Field, degree(ℛ)),0:degree(params.ℛ)-1))
coeffs(e::LWERingDualElement) = e.data
LWERingDualElement(ℛ::LWERing) = LWERingDualElement{ℛ, eltype(ℛ)}

function *(a::LWERingDualElement{ℛ},
           b::LWERingDualElement{ℛ}) where {ℛ}
    LWERingDualElement(ℛ)(a.data .* b.data)
end
function *(a::LWERingDualElement{ℛ},
           b::Integer) where {ℛ}
    LWERingDualElement(ℛ)(a.data * b)
end
function *(a::Integer,
           b::LWERingDualElement{ℛ}) where {ℛ}
    LWERingDualElement(ℛ)(a * b.data)
end

for f in (:+, :-)
    for T in (LWERingElement, LWERingDualElement)
        @eval function $f(a::$T{ℛ},
                b::$T{ℛ}) where {ℛ}
            $T(ℛ)(map($f, coeffs(a), coeffs(b)))
        end
        @eval $f(a::$T{ℛ}) where {ℛ} = $T(ℛ)(map($f, coeffs(a)))
    end
end

function _ntt(ℛ::LWERing, v::AbstractVector)
    @assert first(axes(v)[1]) == 0
    ω = ℛ.ψ^2
    # TODO: Do this using the DFT algorithm
    [sum((v[j]*ω^(j*i)) for j in eachindex(v)) for i in eachindex(v)]
end

function _intt(ℛ::LWERing, v::AbstractVector)
    @assert first(axes(v)[1]) == 0
    ω = ℛ.ψ^2
    # TODO: Do this using the DFT algorithm
    [sum(v[j]*inv(ω)^(j*i) for j = eachindex(v)) for i in eachindex(v)]
end

"""
Perform a negacyclic on the coefficient vector of
p ∈ ℛ = 𝔽q[x]/(xⁿ+1).

Adopting the terminology from [1], we have

    nntt(p) = NTT(PowMulψ(a))

where ψ = √ω and ω is the n-th primitive (the n
here being the same n as in the definition of ℛ).

The additional PowMulψ turns the NTT negacyclic,
i.e. introducing the extra minus sign in the wrap
around required to achieve x^n = -1.

https://eprint.iacr.org/2015/382.pdf
"""
function nntt(p::LWERingElement{ℛ})::LWERingDualElement{ℛ} where {ℛ}
    ψ = ℛ.ψ
    ω = ψ^2
    c̃ = _ntt(ℛ, [x*ψ^i for (i,x) in pairs(p.p.p)])
    LWERingDualElement(ℛ)(c̃)
end

"""
Computes the inverse of nntt(p).
"""
function inntt(p̃::LWERingDualElement{ℛ})::LWERingElement{ℛ} where {ℛ}
    ψ = ℛ.ψ
    ψ⁻¹ = inv(ψ)
    pp = _intt(ℛ, p̃.data)
    n⁻¹ = inv(eltype(ℛ)(length(pp)))
    LWERingElement(ℛ)(FixedDegreePoly([x * n⁻¹ * ψ⁻¹^i for (i, x) in pairs(pp)]))
end

struct RingSampler{Ring} <: Random.Sampler{Ring}
    coeff_distribution::Any #DiscreteUnivariateDistribution
end

function Random.rand(rng::Random.AbstractRNG, r::RingSampler{ℛ}) where {ℛ}
    coeffs = OffsetArray(
        [rand(rng, r.coeff_distribution) for _ in 1:degree(ℛ)],
        0:degree(ℛ)-1)
    LWERingElement(ℛ)(FixedDegreePoly(coeffs))
end

end
