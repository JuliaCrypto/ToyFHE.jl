module NTT

using GaloisFields
using Polynomials
using OffsetArrays
using Distributions
using Random
using AutoHashEquals
using FourierTransforms
using LinearAlgebra

import Base: *, +, -

import GaloisFields: PrimeField

export LWERing, RingSampler, nntt, inntt, FixedDegreePoly,
    LWERingElement, LWERingDualElement, RingCoeffs, coeffs

@auto_hash_equals struct FixedDegreePoly{N, T <: AbstractVector}
    p::T
    function FixedDegreePoly{N,T}(p::T) where {N, T<:AbstractVector}
        @assert first(axes(p)[1]) == 0
        new{N,T}(p)
    end
end
FixedDegreePoly{N}(p::T) where {N, T<:AbstractVector} = FixedDegreePoly{N,T}(p)
function FixedDegreePoly(p::AbstractVector)
    @assert first(axes(p)[1]) == 0
    FixedDegreePoly{length(p), typeof(p)}(p)
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
    function LWERing{Field, N}(ψ::Field) where {Field, N}
        @assert ψ^2N == 1
        new{Field, N}(ψ)
    end
end
Polynomials.degree(ℛ::LWERing{F,N}) where {F,N} = N
Base.eltype(ℛ::LWERing{F,N}) where {F,N} = F

@auto_hash_equals struct RingCoeffs{ℛ, Field, T<:AbstractVector{Field}} <: AbstractVector{Field}
    coeffs::T
end
RingCoeffs{ℛ}(coeffs::T) where {ℛ,T} = RingCoeffs{ℛ, eltype(ℛ), T}(coeffs)
Base.axes(r::RingCoeffs) = axes(r.coeffs)
Base.size(r::RingCoeffs) = size(r.coeffs)
Base.getindex(r::RingCoeffs, idxs...) = getindex(r.coeffs, idxs...)

"""
Represents an element of 𝔽q[x]/(xⁿ+1).
"""
@auto_hash_equals struct LWERingElement{ℛ #= ::LWERing{Field} =#, Field <: PrimeField,  N}
    p::FixedDegreePoly{N, RingCoeffs{ℛ, Field, OffsetVector{Field, Vector{Field}}}}
end
LWERingElement{ℛ,Field,N}(coeffs::AbstractVector) where {ℛ,Field <: PrimeField,  N} = LWERingElement{ℛ,Field,N}(FixedDegreePoly(coeffs))
coeffs(e::LWERingElement) = e.p.p
LWERingElement(ℛ::LWERing) = LWERingElement{ℛ, eltype(ℛ), degree(ℛ)}
LWERingElement(coeffs::RingCoeffs{ℛ}) where {ℛ} = LWERingElement(ℛ)(coeffs)
Base.zero(::Type{LWERingElement{ℛ,Field,N}}) where {ℛ,Field,N} =
    LWERingElement{ℛ,Field,N}(zero(FixedDegreePoly{N, Field}))

"""
Represents an ntt-dual element of 𝔽q[x]/(xⁿ+1).
"""
@auto_hash_equals struct LWERingDualElement{ ℛ #= ::LWERing{Field} =#, Field <: PrimeField}
    data::RingCoeffs{ℛ, Field, OffsetVector{Field, Vector{Field}}}
end
Base.zero(::Type{LWERingDualElement{ℛ,Field}}) where {ℛ,Field} =
    LWERingDualElement(RingCoeffs{ℛ}(OffsetArray(zeros(Field, degree(ℛ)),0:degree(ℛ)-1)))
coeffs(e::LWERingDualElement) = e.data.coeffs
LWERingDualElement(ℛ::LWERing) = LWERingDualElement{ℛ, eltype(ℛ)}

function *(a::LWERingDualElement{ℛ},
           b::LWERingDualElement{ℛ}) where {ℛ}
    LWERingDualElement(RingCoeffs{ℛ}(coeffs(a) .* coeffs(b)))
end
function *(a::LWERingDualElement{ℛ},
           b::Integer) where {ℛ}
    LWERingDualElement(RingCoeffs{ℛ}(coeffs(a) * b))
end
function *(a::Integer,
           b::LWERingDualElement{ℛ}) where {ℛ}
    LWERingDualElement(RingCoeffs{ℛ}(a * coeffs(b)))
end

for f in (:+, :-)
    for T in (LWERingElement, LWERingDualElement)
        @eval function $f(a::$T{ℛ},
                b::$T{ℛ}) where {ℛ}
            $T(ℛ)(RingCoeffs{ℛ}(map($f, coeffs(a), coeffs(b))))
        end
        @eval $f(a::$T{ℛ}) where {ℛ} = $T(RingCoeffs{ℛ}(map($f, coeffs(a))))
    end
end

using FourierTransforms: NontwiddleKernelStep, TwiddleKernelStep, fftgen, CTPlan
@generated function FourierTransforms.applystep(ns::NontwiddleKernelStep{T,N,forward},
    vn::Integer,
    X::RingCoeffs{ℛ},
    x0::Integer, xs::Integer, xvs::Integer,
    Y::RingCoeffs{ℛ},
    y0::Integer, ys::Integer, yvs::Integer) where {T <: PrimeField, ℛ, N, forward}
    ω_n = ℛ.ψ^(2*(div(degree(ℛ), N)))
    forward || (ω_n = inv(ω_n))
    quote
        XV = X.coeffs
        YV = Y.coeffs
        @inbounds @simd for i in 0:vn-1
            $(fftgen(T, ω_n, forward, N,
                     j -> :(XV[(x0 + xvs*i) + xs*$j]),
                     k -> :(YV[(y0 + yvs*i) + ys*$k])))
        end
        Y
    end
end

@generated function FourierTransforms.applystep(ts::TwiddleKernelStep{T,N,forward},
    vn::Integer,
    X::RingCoeffs{ℛ},
    x0::Integer, xs::Integer, xvs::Integer,
    W::AbstractMatrix{T}) where {T <: PrimeField, ℛ, N, forward}
    ω_n = ℛ.ψ^(2*(div(degree(ℛ), N)))
    forward || (ω_n = inv(ω_n))
    quote
        XV = X.coeffs
        @inbounds @simd for i in 0:vn-1
            $(fftgen(T, ω_n, forward, N,
            j -> j == 0 ? :(XV[(x0 + xvs*i) + xs*$j]) :
            forward ? :(W[$j,i+1] * XV[(x0 + xvs*i) + xs*$j]) : :(inv(W[$j,i+1])*XV[(x0 + xvs*i) + xs*$j]),
            j -> :(XV[(x0 + xvs*i) + xs*$j])))
        end
        X
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

function LinearAlgebra.mul!(y::LWERingDualElement{ℛ}, p::CTPlan{T}, x::LWERingElement{ℛ}) where {T, ℛ}
    @assert p.n == length(y.data) == degree(ℛ)
    FourierTransforms.applystep(p, x, 0, 1, y, 0, 1, 1)
    return y
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

[1] https://eprint.iacr.org/2015/382.pdf
"""
function nntt(c::RingCoeffs{ℛ})::RingCoeffs{ℛ} where {ℛ}
    ψ = ℛ.ψ
    ω = ψ^2
    powmulp = RingCoeffs{ℛ}([x*ψ^i for (i,x) in pairs(c.coeffs)])
    ωpow(T, n, i) = (@assert T == eltype(ℛ); ω^(i*div(degree(ℛ),n)))
    c̃ = RingCoeffs{ℛ}(OffsetArray(Vector{eltype(ℛ)}(undef, degree(ℛ)), 0:degree(ℛ)-1))
    mul!(c̃, CTPlan(eltype(ℛ), true, degree(ℛ); ωpow=ωpow), powmulp)
    c̃
end

function nntt(p::LWERingElement{ℛ})::LWERingDualElement{ℛ} where {ℛ}
    LWERingDualElement(nntt(p.p.p))
end

"""
Computes the inverse of nntt(p).
"""
function inntt(c̃::RingCoeffs{ℛ})::RingCoeffs{ℛ} where {ℛ}
    ψ = ℛ.ψ
    ω = ψ^2
    ψ⁻¹ = inv(ψ)
    ω⁻¹ = ψ⁻¹^2
    ωpow(T, n, i) = (@assert T == eltype(ℛ); ω^(i*div(degree(ℛ),n)))
    c = RingCoeffs{ℛ}(OffsetArray(Vector{eltype(ℛ)}(undef, degree(ℛ)), 0:degree(ℛ)-1))
    mul!(c, CTPlan(eltype(ℛ), false, degree(ℛ); ωpow=ωpow), c̃)
    n⁻¹ = inv(eltype(ℛ)(degree(ℛ)))
    RingCoeffs{ℛ}([x * n⁻¹ * ψ⁻¹^i for (i, x) in pairs(c.coeffs)])
end

function inntt(p̃::LWERingDualElement{ℛ})::LWERingElement{ℛ} where {ℛ}
    LWERingElement(inntt(p̃.data))
end

struct RingSampler{Ring} <: Random.Sampler{Ring}
    coeff_distribution::Any #DiscreteUnivariateDistribution
end

function Random.rand(rng::Random.AbstractRNG, r::RingSampler{ℛ}) where {ℛ}
    coeffs = RingCoeffs{ℛ}(OffsetArray(
        [rand(rng, r.coeff_distribution) for _ in 1:degree(ℛ)],
        0:degree(ℛ)-1))
    LWERingElement(coeffs)
end

end
