module NTT

using GaloisFields
using Polynomials
using OffsetArrays
using Distributions
using Random
using AutoHashEquals
using FourierTransforms
using LinearAlgebra

import Base: *, +, -, ^

import GaloisFields: PrimeField
import ..ToyFHE: coefftype, modulus, degree
import Nemo: base_ring

export NegacyclicRing, RingSampler, nntt, inntt, FixedDegreePoly,
    NegacyclicRingElement, NegacyclicRingDualElement, RingCoeffs, coeffs,
    nntt_hint, inntt_hint

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
Base.getindex(p::FixedDegreePoly, args...) = getindex(p.p, args...)

"""
Represents the ring 𝔽q[x]/(xⁿ+1) with optional identified 2n-th primitive root of unity.
"""
struct NegacyclicRing{BaseRing, N}
    # 2N'th primitive root of unity in BaseRing (or zero if unused)
    ψ::BaseRing
    function NegacyclicRing{BaseRing, N}(ψ::BaseRing) where {BaseRing, N}
        @assert ψ^2N == 1
        new{BaseRing, N}(ψ)
    end
    function NegacyclicRing{BaseRing, N}() where {BaseRing, N}
        new{BaseRing, N}(zero(BaseRing))
    end
end
base_ring(R::NegacyclicRing{F}) where {F} = F
coefftype(::Type{NegacyclicRing{F}}) where {F} = F
coefftype(::NegacyclicRing{F}) where {F} = F

# We match the AbstractAlgebra interface here to support more general cycltomics
# elsewhere in the code.
struct PowTwoCyclotomic; n::Int; end
modulus(::NegacyclicRing{<:Any, N}) where {N} = (@assert ispow2(N); PowTwoCyclotomic(2N))
modulus(::Type{UInt8}) = 256
modulus(F::Type{<:PrimeField}) = char(F)
degree(ptc::PowTwoCyclotomic) = div(ptc.n, 2)

function NegacyclicRing{Field, N}(ψ::Integer) where {Field, N}
    @assert ψ < char(Field)
    NegacyclicRing{Field, N}(convert(Field, ψ))
end
degree(ℛ::NegacyclicRing{F,N}) where {F,N} = N
Base.eltype(ℛ::NegacyclicRing{F,N}) where {F,N} = F

@auto_hash_equals struct RingCoeffs{ℛ, Field, T<:AbstractVector{Field}} <: AbstractVector{Field}
    coeffs::T
end
RingCoeffs{ℛ}(coeffs::T) where {ℛ,T} = RingCoeffs{ℛ, eltype(ℛ), T}(coeffs)
Base.convert(::Type{RingCoeffs{ℛ,F,T}}, coeffs::T) where {ℛ,F,T} = RingCoeffs{ℛ,F,T}(coeffs)
Base.copy(r::RingCoeffs{ℛ,F,T}) where {ℛ,F,T} = RingCoeffs{ℛ,F,T}(copy(r.coeffs))
RingCoeffs{ℛ}(r::RingCoeffs{ℛ}) where {ℛ} = copy(r)
Base.axes(r::RingCoeffs) = axes(r.coeffs)
Base.size(r::RingCoeffs) = size(r.coeffs)
Base.getindex(r::RingCoeffs, idxs...) = getindex(r.coeffs, idxs...)

"""
Represents an element of 𝔽q[x]/(xⁿ+1).
"""
@auto_hash_equals struct NegacyclicRingElement{ℛ #= ::NegacyclicRing{Field} =#, Field,  N}
    p::FixedDegreePoly{N, RingCoeffs{ℛ, Field, OffsetVector{Field, Vector{Field}}}}
end
NegacyclicRingElement{ℛ,Field,N}(coeffs::AbstractVector) where {ℛ,Field,  N} =
    NegacyclicRingElement{ℛ,Field,N}(FixedDegreePoly(RingCoeffs{ℛ}(coeffs)))
Base.convert(::Type{NegacyclicRingElement{ℛ,Field,N}}, coeffs::OffsetVector{Field, Vector{Field}}) where {ℛ, Field, N} =
    NegacyclicRingElement{ℛ,Field,N}(FixedDegreePoly(RingCoeffs{ℛ}(coeffs)))
coeffs(e::NegacyclicRingElement) = e.p.p
NegacyclicRingElement(ℛ::NegacyclicRing) = NegacyclicRingElement{ℛ, eltype(ℛ), degree(modulus(ℛ))}
NegacyclicRingElement(coeffs::RingCoeffs{ℛ}) where {ℛ} = NegacyclicRingElement(ℛ)(coeffs)
Base.zero(::Type{NegacyclicRingElement{ℛ,Field,N}}) where {ℛ,Field,N} =
    NegacyclicRingElement{ℛ,Field,N}(zero(FixedDegreePoly{N, Field}))

(ℛ::NegacyclicRing)(coeffs) = NegacyclicRingElement(ℛ)(coeffs)

"""
Represents an ntt-dual element of 𝔽q[x]/(xⁿ+1).
"""
@auto_hash_equals struct NegacyclicRingDualElement{ ℛ #= ::NegacyclicRing{Field} =#, Field <: PrimeField}
    data::RingCoeffs{ℛ, Field, OffsetVector{Field, Vector{Field}}}
end
Base.zero(::Type{NegacyclicRingDualElement{ℛ,Field}}) where {ℛ,Field} =
    NegacyclicRingDualElement(RingCoeffs{ℛ}(OffsetArray(zeros(Field, degree(modulus(ℛ))),0:degree(modulus(ℛ))-1)))
Base.zero(d::NegacyclicRingDualElement) = zero(typeof(d))
coeffs(e::NegacyclicRingDualElement) = e.data.coeffs
NegacyclicRingDualElement(ℛ::NegacyclicRing) = NegacyclicRingDualElement{ℛ, eltype(ℛ)}

function *(a::NegacyclicRingDualElement{ℛ},
           b::NegacyclicRingDualElement{ℛ}) where {ℛ}
    NegacyclicRingDualElement(RingCoeffs{ℛ}(coeffs(a) .* coeffs(b)))
end
function *(a::NegacyclicRingDualElement{ℛ},
           b::Union{Integer, PrimeField}) where {ℛ}
    NegacyclicRingDualElement(RingCoeffs{ℛ}(coeffs(a) * b))
end
function *(a::Union{Integer, PrimeField},
           b::NegacyclicRingDualElement{ℛ}) where {ℛ}
    NegacyclicRingDualElement(RingCoeffs{ℛ}(a * coeffs(b)))
end
function ^(x::NegacyclicRingDualElement, n::Integer)
    @assert n >= 0
    Base.power_by_squaring(x,n)
end


for f in (:+, :-)
    for T in (NegacyclicRingElement, NegacyclicRingDualElement)
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
    ω_n = ℛ.ψ^(2*(div(degree(modulus(ℛ)), N)))
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
    ω_n = ℛ.ψ^(2*(div(degree(modulus(ℛ)), N)))
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

function _ntt(ℛ::NegacyclicRing, v::AbstractVector)
    @assert first(axes(v)[1]) == 0
    ω = ℛ.ψ^2
    # TODO: Do this using the DFT algorithm
    [sum((v[j]*ω^(j*i)) for j in eachindex(v)) for i in eachindex(v)]
end

function _intt(ℛ::NegacyclicRing, v::AbstractVector)
    @assert first(axes(v)[1]) == 0
    ω = ℛ.ψ^2
    # TODO: Do this using the DFT algorithm
    [sum(v[j]*inv(ω)^(j*i) for j = eachindex(v)) for i in eachindex(v)]
end

function LinearAlgebra.mul!(y::NegacyclicRingDualElement{ℛ}, p::CTPlan{T}, x::NegacyclicRingElement{ℛ}) where {T, ℛ}
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
    ωpow(T, n, i) = (@assert T == eltype(ℛ); ω^(i*div(degree(modulus(ℛ)),n)))
    c̃ = RingCoeffs{ℛ}(OffsetArray(Vector{eltype(ℛ)}(undef, degree(modulus(ℛ))), 0:degree(modulus(ℛ))-1))
    mul!(c̃, CTPlan(eltype(ℛ), true, degree(modulus(ℛ)); ωpow=ωpow), powmulp)
    c̃
end

function nntt(p::NegacyclicRingElement{ℛ})::NegacyclicRingDualElement{ℛ} where {ℛ}
    NegacyclicRingDualElement(nntt(p.p.p))
end

"""
Computes the inverse of nntt(p).
"""
function inntt(c̃::RingCoeffs{ℛ})::RingCoeffs{ℛ} where {ℛ}
    ψ = ℛ.ψ
    ω = ψ^2
    ψ⁻¹ = inv(ψ)
    ω⁻¹ = ψ⁻¹^2
    ωpow(T, n, i) = (@assert T == eltype(ℛ); ω^(i*div(degree(modulus(ℛ)),n)))
    c = RingCoeffs{ℛ}(OffsetArray(Vector{eltype(ℛ)}(undef, degree(modulus(ℛ))), 0:degree(modulus(ℛ))-1))
    mul!(c, CTPlan(eltype(ℛ), false, degree(modulus(ℛ)); ωpow=ωpow), c̃)
    n⁻¹ = inv(eltype(ℛ)(degree(modulus(ℛ))))
    RingCoeffs{ℛ}([x * n⁻¹ * ψ⁻¹^i for (i, x) in pairs(c.coeffs)])
end

function inntt(p̃::NegacyclicRingDualElement{ℛ})::NegacyclicRingElement{ℛ} where {ℛ}
    NegacyclicRingElement(inntt(p̃.data))
end

# Hints
nntt_hint(r) = r
nntt_hint(r::NegacyclicRingElement) = nntt(r)
inntt_hint(r) = r
inntt_hint(r::NegacyclicRingDualElement) = inntt(r)

end
