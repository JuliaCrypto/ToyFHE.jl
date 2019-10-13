using Hecke
using AbstractAlgebra
using GaloisFields
using GaloisFields: PrimeField
using Random
using StructArrays

struct CRTEncoded{T, N, M <: NTuple{N, PrimeField{T}}} <: Number
    c::M
end
Base.getproperty(c::CRTEncoded, i::Int64) = c.c[i]
NTT.modulus(::Type{CRTEncoded{T,N,mods}}) where {T,N,mods} = prod(x->BigInt(NTT.modulus(x)), fieldtypes(mods))
NTT.modulus(c::CRTEncoded) = NTT.modulus(typeof(c))
moduli(::Type{CRTEncoded{T,N,mods}}) where {T,N,mods} = mods

function Base.promote_rule(::Type{CRTEncoded{T, N, moduli}}, ::Type{<:Integer}) where {T, N, moduli}
    CRTEncoded{T, N, moduli}
end

function CRTEncoded{T, N, M}(x::Integer) where {T,N,M <: NTuple{N, PrimeField{T}}}
    CRTEncoded{T,N,M}(map(fieldtypes(M)) do T
        T(x)
    end)
end
CRTEncoded{T, N, M}(x::CRTEncoded{T, N, M}) where {T,N,M<: NTuple{N, PrimeField{T}}} = x

function AbstractAlgebra.crt(r1::Integer, m1::Integer, r2::Integer, m2::Integer)
    g, u, v = gcdx(m1, m2)
    @assert g==1
    m = m1*m2
    return mod(r1*v*m2 + r2*u*m1, m)
end

function Base.convert(::Type{Integer}, x::CRTEncoded{<:Any, <:Any, moduli}) where {moduli}
    # TODO: Do this ourselves?
    AbstractAlgebra.crt(collect(map(y->BigInt(convert(Integer, y)), x.c)), collect(map(x->BigInt(modulus(x)), fieldtypes(moduli))))
end

function Base.:+(a::CRTEncoded{T, N, moduli}, b::CRTEncoded{T, N, moduli}) where {T,N,moduli}
    CRTEncoded{T, N, moduli}(a.c .+ b.c)
end

function Base.:*(a::CRTEncoded{T, N, moduli}, b::CRTEncoded{T, N, moduli}) where {T,N,moduli}
    CRTEncoded{T, N, moduli}(a.c .* b.c)
end

function Base.:-(a::CRTEncoded{T, N, moduli}) where {T,N,moduli}
    CRTEncoded{T, N, moduli}(.-(a.c))
end

function NTT.is_primitive_root(ψ::CRTEncoded{<:Any, <:Any, moduli}, n) where {moduli}
    all(ψ.c) do r
        return NTT.is_primitive_root(r, n)
    end
end

function GaloisFields.minimal_primitive_root(::Type{CRTEncoded{T,N,moduli}}, n) where {T,N,moduli}
    CRTEncoded{T, N, moduli}(map(T->GaloisFields.minimal_primitive_root(T, n), fieldtypes(moduli)))
end

function Random.rand(rng::AbstractRNG, ::Random.SamplerType{CRTEncoded{T, N, moduli}}) where {T,N,moduli}
    CRTEncoded{T,N,moduli}(map(rand, fieldtypes(moduli)))
end

function StructArrays.staticschema(CT::Type{<:CRTEncoded{T,N}}) where {T,N}
    fieldtype(CT, 1)
end

function StructArrays.createinstance(::Type{T}, args...) where {T<:CRTEncoded}
    T(args)
end


# Integration with NTT
function NTT.nntt(p::NegacyclicRingElement{ℛ,<:Any,<:Any,<:StructArray{T}})::NegacyclicRingDualElement{ℛ} where {ℛ, T<:CRTEncoded}
    rcs = p.p.p
    oa = rcs.coeffs
    sa = oa.parent
    sa′ = StructArray{eltype(sa)}(tuple(map(enumerate(fieldarrays(sa))) do (i,a)
        oai = OffsetArray(a, axes(oa)...)
        ℛi = NegacyclicRing{eltype(oai), degree(ℛ)}(ℛ.ψ.c[i])
        nntt(RingCoeffs{ℛi}(oai)).coeffs.parent
    end...))
    NegacyclicRingDualElement(RingCoeffs{ℛ}(OffsetArray(sa′, axes(oa)...)))
end

function NTT.inntt(p::NegacyclicRingDualElement{ℛ,<:Any,<:StructArray{T}})::NegacyclicRingElement{ℛ} where {ℛ, T<:CRTEncoded}
    rcs = p.data
    oa = rcs.coeffs
    sa = oa.parent
    sa′ = StructArray{eltype(sa)}(tuple(map(enumerate(fieldarrays(sa))) do (i,a)
        oai = OffsetArray(a, axes(oa)...)
        ℛi = NegacyclicRing{eltype(oai), degree(ℛ)}(ℛ.ψ.c[i])
        inntt(RingCoeffs{ℛi}(oai)).coeffs.parent
    end...))
    NegacyclicRingElement(RingCoeffs{ℛ}(OffsetArray(sa′, axes(oa)...)))
end

function sample_ring_array(rng::Random.AbstractRNG, ℛ::NegacyclicRing{<:CRTEncoded{T, N, moduli}}, coeff_distribution) where {T,N,moduli}
    StructArray(coefftype(ℛ)(rand(rng, coeff_distribution)) for _ in 1:degree(modulus(ℛ)))
end
