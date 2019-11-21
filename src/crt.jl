using Hecke
using AbstractAlgebra
using GaloisFields
using GaloisFields: PrimeField
using Random
using StructArrays

struct CRTEncoded{N, M <: NTuple{N, PrimeField}} <: Number
    c::M
end

"""
    CRTExpand{T<:PrimeField}

A singleton representing expansion of the CRT basis by adding `T` as an
additional CRT component. In the integer domain this represents multiplication
by `modulus(T)`, so this type is used by multiplication with an existing CRT
encoded value.

# Examples:
```
julia> x = CRTEncoded{2, Tuple{𝔽₅, 𝔽₇}}(3)
CRTEncoded{2,Tuple{𝔽₅,𝔽₇}}((3, 3))

julia> y = CRTExpand{𝔽₁₁}();

julia> x*y
CRTEncoded{3,Tuple{𝔽₅,𝔽₇,𝔽₁₁}}((3, 5, 0))

julia> convert(Integer, x*y)
333
```
"""
struct CRTExpand{T<:PrimeField} <: Number; end
Base.convert(::Type{Integer}, c::CRTExpand{T}) where {T} = modulus(T)

function Base.:*(a::CRTEncoded{N, M}, b::CRTExpand{T}) where {T,N,M}
    CRTEncoded{N+1, Tuple{M.parameters..., T}}(tuple((convert(Integer, b) * a).c..., zero(T)))
end

"""
    CRTResidual{T<:PrimeField}

Represents the value c*[(q/q_i)^{-1}]_{q_i} q/q_i where `q_i=modulus(c)` and `q`
is the modulus of some larger CRT representation. That CRT representation is
not represented in this type and may thus be.

# Example:
```julia
julia> r = CRTResidual(𝔽₃(3))

julia> convert(Integer, CRTEncoded{3, Tuple{𝔽₃, 𝔽₅, 𝔽₇}}(r))
308

julia> mod(3*invmod(77, 5), 5)*77
308
```
"""
struct CRTResidual{T<:PrimeField}
    c::T
end

function CRTEncoded{N, M}(c::CRTResidual{T}) where {N,M<:NTuple{N, PrimeField},T}
    found = false
    ret = CRTEncoded{N, M}(tuple(map(M.parameters) do F
        if F == T
            @assert !found
            found = true
            c.c
        else
            F(0)
        end
    end...))
    @assert found
    ret
end

Base.getproperty(c::CRTEncoded, i::Int64) = c.c[i]
NTT.modulus(::Type{CRTEncoded{N,mods}}) where {N,mods} = prod(x->BigInt(NTT.modulus(x)), fieldtypes(mods))
NTT.modulus(c::CRTEncoded) = NTT.modulus(typeof(c))
moduli(::Type{CRTEncoded{N,mods}}) where {N,mods} = mods
moduli(::Type{T}) where {T<:NTT.RingElement} where {N,mods} = moduli(coefftype(T))
moduli(::T) where {T<:NTT.RingElement} where {N,mods} = moduli(T)
moduli(::NegacyclicRing{BaseRing}) where {BaseRing} = moduli(BaseRing)

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
    if length(x.c) == 1
        convert(Integer, x.c[1])
    else
        # TODO: Do this ourselves?
        AbstractAlgebra.crt(collect(map(y->BigInt(convert(Integer, y)), x.c)), collect(map(x->BigInt(modulus(x)), fieldtypes(moduli))))
    end
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
    CRTEncoded{N, moduli}(map(*, a.c, b.c))
end

function Base.:-(a::CRTEncoded{N, moduli}) where {N,moduli}
    CRTEncoded{N, moduli}(.-(a.c))
end

function Base.:-(a::CRTEncoded{N, moduli}, b::CRTEncoded{N, moduli}) where {T,N,moduli}
    CRTEncoded{N, moduli}(a.c .- b.c)
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


# Modswitch to next modulus in CRT chain

struct DropLastParams{P<:SHEShemeParams} <: SHEShemeParams
    params::P
end
modswitch_drop(params::SHEShemeParams) = DropLastParams(params)
scheme_name(::Type{DropLastParams{P}}) where P = scheme_name(P)
ℛ_plain(p::DropLastParams) = ℛ_plain(p.params)
ℛ_cipher(p::DropLastParams) = drop_last(ℛ_cipher(p.params))
relin_window(p::DropLastParams) = relin_window(p.params)

π⁻¹(p::DropLastParams, plaintext) = modswitch_drop(π⁻¹(p.params, plaintext))
π(p::DropLastParams, b) = b

struct DropLastSampler{Ring} <: Random.Sampler{Ring}
    s::Union{RingSampler{Ring}, DropLastSampler{Ring}}
end

function Random.rand(rng::Random.AbstractRNG, r::DropLastSampler)
    # TODO: We could just directly generate the appropriate ring element here.
    modswitch_drop(rand(rng, r.s))
end

𝒩(p::DropLastParams) = DropLastSampler(𝒩(p.params))
𝒢(p::DropLastParams) = DropLastSampler(𝒢(p.params))

function crtselect(::Type{CRTEncoded{N,M}}, which) where {N,M}
    CRTEncoded{length(which),Tuple{M.parameters[which]...}}
end

function crtselect(ℛ::NegacyclicRing{T, N}, which) where {T<:CRTEncoded, N}
    T′ = crtselect(T, which)
    NegacyclicRing{T′, N}(T′(ℛ.ψ.c[which]))
end

function crtselect(::Type{NTT.RingElement{ℛ, Field, Storage}}, which) where {ℛ, Field, Storage}
    ℛnew = crtselect(ℛ, which)
    NTT.RingElement{ℛnew, coefftype(ℛnew), similar(Storage, coefftype(ℛnew))}
end

function crtselect(x::NTT.RingElement{ℛ, Field, Storage}, which) where {ℛ, Field, Storage}
    ℛnew = crtselect(ℛ, which)
    newA = similar(Storage, coefftype(ℛnew))
    newAConstruct = StructArray{newA.parameters[1:3]...}
    primaly = dualy = nothing
    if x.primal !== nothing
        primaly = OffsetArray(newAConstruct(fieldarrays(parent(x.primal))[which]), axes(x.primal))
    end
    if x.dual !== nothing
        dualy = OffsetArray(newAConstruct(fieldarrays(parent(x.dual))[which]), axes(x.dual))
    end
    NTT.RingElement{ℛnew, coefftype(ℛnew), newA}(primaly,dualy)
end

drop_last(x) = crtselect(x, 1:(length(moduli(x).parameters)-1))

function modswitch(crt::CRTEncoded)
    ct_qk = crt.c[end]
    CRTEncoded(map(crt.c[1:end-1]) do cc
        inv(typeof(cc)(modulus(ct_qk))) * (cc - typeof(cc)(convert(Integer, ct_qk)))
    end)
end

function modswitch_drop(crt::CRTEncoded)
    CRTEncoded(crt.c[1:end-1])
end

function modswitch(re::NTT.RingElement{ℛ, Field}) where {ℛ, Field <: CRTEncoded}
    NTT.RingElement{drop_last(ℛ)}(map(modswitch, NTT.coeffs_primal(re)), nothing)
end

function modswitch_drop(re::NTT.RingElement{ℛ, Field}) where {ℛ, Field <: CRTEncoded}
    NTT.RingElement{drop_last(ℛ)}(map(modswitch_drop, NTT.coeffs_primal(re)), nothing)
end

function modswitch_drop(c::CipherText{Enc}) where {Enc}
    CipherText{drop_last(Enc)}(modswitch_drop(c.params), map(modswitch_drop, c.cs))
end

function downswitch_keyelement(params, key::KeyComponent, elt::NTT.RingElement{<:Any, <:CRTEncoded})
    length(moduli(key.mask).parameters) == length(moduli(elt).parameters) && return key
    KeyComponent(
        crtselect(key.mask, 1:length(moduli(elt).parameters)),
        crtselect(key.masked, 1:length(moduli(elt).parameters))
    )
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
