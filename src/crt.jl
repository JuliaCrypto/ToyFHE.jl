using Hecke
using AbstractAlgebra
using GaloisFields
using GaloisFields: PrimeField
using Random
using StructArrays

struct CRTEncoded{N, M <: NTuple{N, PrimeField}} <: Number
    c::M
end
Base.getproperty(c::CRTEncoded, i::Int64) = c.c[i]
NTT.modulus(::Type{CRTEncoded{N,mods}}) where {N,mods} = prod(x->BigInt(NTT.modulus(x)), fieldtypes(mods))
NTT.modulus(c::CRTEncoded) = NTT.modulus(typeof(c))
moduli(::Type{CRTEncoded{N,mods}}) where {N,mods} = mods

function Base.promote_rule(::Type{CRTEncoded{N, moduli}}, ::Type{<:Integer}) where {N, moduli}
    CRTEncoded{N, moduli}
end

function CRTEncoded{N, M}(x::Integer) where {N,M <: NTuple{N, PrimeField}}
    CRTEncoded{N,M}(map(fieldtypes(M)) do T
        T(x)
    end)
end
CRTEncoded{N, M}(x::CRTEncoded{N, M}) where {T,N,M<: NTuple{N, PrimeField}} = x

function AbstractAlgebra.crt(r1::Integer, m1::Integer, r2::Integer, m2::Integer)
    g, u, v = gcdx(m1, m2)
    @assert g==1
    m = m1*m2
    return mod(r1*v*m2 + r2*u*m1, m)
end

function Base.convert(::Type{Integer}, x::CRTEncoded{<:Any, moduli}) where {moduli}
    # TODO: Do this ourselves?
    AbstractAlgebra.crt(collect(map(y->BigInt(convert(Integer, y)), x.c)), collect(map(x->BigInt(modulus(x)), fieldtypes(moduli))))
end

function Base.convert(T::Type{<:CRTEncoded}, x::PrimeField)
    x = convert(Integer, x)
    x > modulus(T) && throw(InexactError(:convert, T, x))
    T(x)
end

function Base.:+(a::CRTEncoded{N, moduli}, b::CRTEncoded{N, moduli}) where {T,N,moduli}
    CRTEncoded{N, moduli}(a.c .+ b.c)
end

function Base.:*(a::CRTEncoded{N, moduli}, b::CRTEncoded{N, moduli}) where {N,moduli}
    CRTEncoded{N, moduli}(a.c .* b.c)
end

function Base.:-(a::CRTEncoded{N, moduli}) where {N,moduli}
    CRTEncoded{N, moduli}(.-(a.c))
end

function NTT.is_primitive_root(ψ::CRTEncoded{<:Any, moduli}, n) where {moduli}
    all(ψ.c) do r
        return NTT.is_primitive_root(r, n)
    end
end

function GaloisFields.minimal_primitive_root(::Type{CRTEncoded{N,moduli}}, n) where {T,N,moduli}
    CRTEncoded{N, moduli}(map(T->GaloisFields.minimal_primitive_root(T, n), fieldtypes(moduli)))
end

function Random.rand(rng::AbstractRNG, ::Random.SamplerType{CRTEncoded{N, moduli}}) where {T,N,moduli}
    CRTEncoded{N,moduli}(map(rand, fieldtypes(moduli)))
end

function StructArrays.staticschema(CT::Type{<:CRTEncoded})
    fieldtype(CT, 1)
end

function StructArrays.createinstance(::Type{T}, args...) where {T<:CRTEncoded}
    T(args)
end


# Integration with NTT
function NTT.nntt(rcs::RingCoeffs{ℛ,T,OffsetVector{T, S}})::RingCoeffs{ℛ} where {ℛ, T<:CRTEncoded, S<:StructArray{T}}
    oa = rcs.coeffs
    sa = oa.parent
    sa′ = StructArray{eltype(sa)}(tuple(map(enumerate(fieldarrays(sa))) do (i,a)
        oai = OffsetArray(a, axes(oa)...)
        ℛi = NegacyclicRing{eltype(oai), degree(ℛ)}(ℛ.ψ.c[i])
        nntt(RingCoeffs{ℛi}(oai)).coeffs.parent
    end...))
    RingCoeffs{ℛ}(OffsetArray(sa′, axes(oa)...))
end

function NTT.inntt(rcs::RingCoeffs{ℛ,T,OffsetVector{T, S}})::RingCoeffs{ℛ} where {ℛ, T<:CRTEncoded, S<:StructArray{T}}
    oa = rcs.coeffs
    sa = oa.parent
    sa′ = StructArray{eltype(sa)}(tuple(map(enumerate(fieldarrays(sa))) do (i,a)
        oai = OffsetArray(a, axes(oa)...)
        ℛi = NegacyclicRing{eltype(oai), degree(ℛ)}(ℛ.ψ.c[i])
        inntt(RingCoeffs{ℛi}(oai)).coeffs.parent
    end...))
    RingCoeffs{ℛ}(OffsetArray(sa′, axes(oa)...))
end

function NTT.nntt(rcs::RingCoeffs{ℛ,T,OffsetVector{T, S}})::RingCoeffs{ℛ} where {ℛ, T<:CRTEncoded, S<:AbstractArray{T}}
    error("NNTT methods currently only implemented for CRTEncoded in StructArray format")
end

function NTT.inntt(rcs::RingCoeffs{ℛ,T,OffsetVector{T, S}})::RingCoeffs{ℛ} where {ℛ, T<:CRTEncoded, S<:AbstractArray{T}}
    error("NNTT methods currently only implemented for CRTEncoded in StructArray format")
end

function sample_ring_array(rng::Random.AbstractRNG, ℛ::NegacyclicRing{<:CRTEncoded{N, moduli}}, coeff_distribution) where {N,moduli}
    StructArray(coefftype(ℛ)(rand(rng, coeff_distribution)) for _ in 1:degree(modulus(ℛ)))
end
