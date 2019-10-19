module NTT

using GaloisFields
using Polynomials
using OffsetArrays
using Distributions
using Random
using AutoHashEquals
using FourierTransforms
using LinearAlgebra
using StructArrays

import Base: *, +, -, ^

import GaloisFields: PrimeField
import ..ToyFHE: coefftype, modulus, degree
import Nemo: base_ring

export NegacyclicRing, RingSampler, nntt, inntt, RingCoeffs, coeffs,
    nntt_hint, inntt_hint

is_primitive_root(ψ, n) = ψ^n == 1

"""
Represents the ring 𝔽q[x]/(xⁿ+1) with optional identified 2n-th primitive root of unity.
"""
struct NegacyclicRing{BaseRing, N}
    # 2N'th primitive root of unity in BaseRing (or zero if unused)
    ψ::BaseRing
    function NegacyclicRing{BaseRing, N}(ψ::BaseRing) where {BaseRing, N}
        @assert is_primitive_root(ψ, 2N)
        new{BaseRing, N}(ψ)
    end
    function NegacyclicRing{BaseRing, N}(::Nothing) where {BaseRing, N}
        new{BaseRing, N}(zero(BaseRing))
    end
end
function NegacyclicRing{BaseRing, N}() where {BaseRing <: PrimeField, N}
    if gcd(char(BaseRing) - 1, 2N) == 2N
        NegacyclicRing{BaseRing, N}(GaloisFields.minimal_primitive_root(BaseRing, 2N))
    else
        NegacyclicRing{BaseRing, N}(nothing)
    end
end
function NegacyclicRing{BaseRing, N}() where {BaseRing, N}
    NegacyclicRing{BaseRing, N}(nothing)
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

function (ℛ::NegacyclicRing)(coeffs::OffsetVector)
    RingElement{ℛ}(coeffs, nothing)
end
(ℛ::NegacyclicRing)(coeffs) = convert(RingElement{ℛ}, coeffs)


function Base.zero(ℛ::NegacyclicRing)
    RingElement{ℛ}(OffsetArray(zeros(eltype(ℛ), degree(ℛ)), 0:degree(ℛ)-1), nothing)
end

"""
Represents an element of 𝔽q[x]/(xⁿ+1).

Also optionally caches its dual to efficiently perform multiplicative
operations.
"""
@auto_hash_equals mutable struct RingElement{ℛ #= ::NegacyclicRing{Field} =#, Field, Storage <: AbstractVector{Field}} <: AbstractVector{Field}
    primal::Union{Nothing, OffsetVector{Field, Storage}}
    dual::Union{Nothing, OffsetVector{Field, Storage}}
end

function Base.convert(::Type{<:RingElement{ℛ₁}}, r::RingElement{ℛ₂}) where {ℛ₁, ℛ₂}
    RingElement{ℛ₁}(map(c->convert(eltype(ℛ₁), c), coeffs_primal(r)), nothing)
end

function Base.convert(::Type{RingElement{ℛ,Field,Storage}}, primal::OffsetVector{Field, Storage}) where {ℛ, Field, Storage}
    RingElement{ℛ,Field,Storage}(primal, nothing)
end

function RingElement{ℛ}(primal::Union{Nothing, OffsetVector{Field, Storage}},
                        dual::Union{Nothing, OffsetVector{Field, Storage}}) where {ℛ, Field, Storage}
    @assert primal !== nothing || dual !== nothing
    RingElement{ℛ,Field,Storage}(primal, dual)
end
Base.axes(r::RingElement{ℛ}) where {ℛ} = (Base.IdentityUnitRange(0:degree(ℛ)-1),)
Base.size(r::RingElement{ℛ}) where {ℛ} = map(length, axes(r))
Base.zero(r::RingElement{ℛ}) where {ℛ} = RingElement{ℛ}(zero(r.primal), nothing)

function coeffs_primal(r::RingElement{ℛ}) where {ℛ}
    if r.primal === nothing
        @assert r.dual !== nothing
        r.primal = inntt(RingCoeffs{ℛ}(r.dual)).coeffs
    end
    return r.primal
end

function coeffs_dual(r::RingElement{ℛ}) where {ℛ}
    if r.dual === nothing
        @assert r.primal !== nothing
        r.dual = nntt(RingCoeffs{ℛ}(r.primal)).coeffs
    end
    return r.dual
end

Base.getindex(r::RingElement, idxs...) = getindex(coeffs_primal(r), idxs...)
function Base.setindex!(r::RingElement, v, idxs...)
    ret = setindex!(coeffs_primal(r), v, idxs...)
    r.dual = nothing
    ret
end

function ring_multiply(ℛ::NegacyclicRing{BaseRing}, a::RingElement, b::RingElement) where BaseRing
    # TODO: This should ideally use dispatch
    # https://github.com/JuliaLang/julia/issues/33387
    if iszero(ℛ.ψ)
        # Do the naive convolution. This is just for plaintexts and
        # testing.
        ca = coeffs_primal(a)
        cb = coeffs_primal(b)
        res = zero(ca)
        N = lastindex(ca)
        for (i, j) in Iterators.product(eachindex(ca), eachindex(cb))
            idx = i+j
            if idx <= N
                res[idx] += ca[i] * cb[j]
            else
                res[idx-N] -= ca[i] * cb[j]
            end
        end
        return RingElement{ℛ}(res, nothing)
    else
        return RingElement{ℛ}(nothing, coeffs_dual(a) .* coeffs_dual(b))
    end
end

function *(a::RingElement{ℛ}, b::RingElement{ℛ}) where {ℛ}
    ring_multiply(ℛ, a, b)
end

const RingScalar = Union{Integer, PrimeField}
function *(a::RingScalar, b::RingElement{ℛ}) where {ℛ}
    RingElement{ℛ}(b.primal === nothing ? nothing : a .* b.primal,
                    b.dual === nothing ? nothing : a .* b.dual)
end
*(a::RingElement{ℛ}, b::RingScalar) where {ℛ} = b*a

function -(a::RingElement{ℛ}) where {ℛ}
    RingElement{ℛ}(a.primal === nothing ? nothing : -a.primal,
                    a.dual === nothing ? nothing : -a.dual)
end

for f in (:+, :-)
    @eval function ($f)(a::RingElement{ℛ}, b::RingElement{ℛ}) where {ℛ}
        # If both have both primal and dual set, we sum both (we sort of
        # expect a higher order compiler to remove any unnecessary computation,
        # though julia doesn't currently do that). If only one of them is set in
        # each, we do that.
        new_primal = new_dual = nothing
        if a.primal !== nothing && b.primal !== nothing
            new_primal = broadcast($f, a.primal, b.primal)
        end
        if a.dual !== nothing && b.dual !== nothing
            new_dual = broadcast($f, a.dual, b.dual)
        end
        if new_primal === nothing && new_dual === nothing
            if a.primal === nothing
                new_primal = broadcast($f, coeffs_primal(a), b.primal)
            else
                new_primal = broadcast($f, a.primal, coeffs_primal(b))
            end
            if a.dual === nothing
                new_dual = broadcast($f, coeffs_dual(a), b.dual)
            else
                new_dual = broadcast($f, a.dual, coeffs_dual(b))
            end
        end
        RingElement{ℛ}(new_primal, new_dual)
    end
end

function ^(x::RingElement, n::Integer)
    @assert n >= 0
    Base.power_by_squaring(x,n)
end

@auto_hash_equals struct RingCoeffs{ℛ, Field, T<:AbstractVector{Field}} <: AbstractVector{Field}
    coeffs::T
end
RingCoeffs{ℛ}(coeffs::T) where {ℛ,T} = RingCoeffs{ℛ, eltype(ℛ), T}(coeffs)
Base.convert(::Type{RingCoeffs{ℛ,F,T}}, coeffs::T) where {ℛ,F,T} = RingCoeffs{ℛ,F,T}(coeffs)
Base.copy(r::RingCoeffs{ℛ,F,T}) where {ℛ,F,T} = RingCoeffs{ℛ,F,T}(copy(r.coeffs))
RingCoeffs{ℛ}(r::RingCoeffs{ℛ}) where {ℛ} = copy(r)
Base.axes(r::RingCoeffs) = axes(r.coeffs)
Base.size(r::RingCoeffs) = size(r.coeffs)

# Negacyclic Numbertheortic transform (i.e. FFT over finite fields)
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

end
