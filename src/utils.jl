module Utils

using Base.Meta
using GaloisFields
using GaloisFields: PrimeField
using Distributions
using Nemo
using AbstractAlgebra
using ..NTT
using ..NTT: modulus, degree
using Primes
using Mods
using StructArrays
using OffsetArrays
using BitIntegers
using Hecke

macro fields_as_locals(a)
    @assert isexpr(a, :(::))
    sym, T = a.args
    T = getfield(__module__, T)
    quote
        $([:(local $(esc(fname))) for fname in fieldnames(T)]...)
        let s = $(esc(sym))
            $([:($(esc(fname)) = s.$fname) for fname in fieldnames(T)]...)
        end
    end
end

# HACK - DiscreteUniform is the default for types, but it'd be nice to
# allow this.
Distributions.DiscreteUniform(T::Type) = T

# HACK
Base.BroadcastStyle(::Type{<:OffsetArray{<:Any, <:Any, T}}) where T<:StructArray = Base.BroadcastStyle(T)

# HACK
Base.convert(::Type{Integer}, x::PrimeField) = x.n

# HACK
Base.convert(::Type{<:StructArray{T}}, A::AbstractArray{T}) where {T} =
    StructArray(A)

# Hack
Base.convert(::Type{Integer}, x::AbstractAlgebra.Generic.Res{fmpz}) = Nemo.lift(x)

# Hack
Base.div(a::Nemo.fmpz, b::BigInt, r::RoundingMode{:NearestTiesAway}) =
    div(BigInt(a), b, r)

# Hack
function Base.convert(T::Type{<:PrimeField}, x::PrimeField)
    x = convert(Integer, x)
    x > modulus(T) && throw(InexactError(:convert, T, x))
    T(x)
end

# Hack
Primes.isprime(x::Int256) = Primes.isprime(big(x))

# Hack
(::NmodPolyRing)(x::AbstractAlgebra.Generic.Res{fmpz_mod_poly}) = x

# Hack
Base.gcdx(x::nmod_poly, y::nmod_poly) = Hecke.rresx(x, y)

# Hack
Nemo.nmod_poly(a::Integer, b::Integer) = Nemo.nmod_poly(UInt64(a), UInt64(b))

# Hack (https://github.com/Nemocas/Nemo.jl/issues/668)
function (a::FqNmodFiniteField)(b::Nemo.nmod)
    a(data(b))
end

# Hack
function Base.oftype(R::AbstractAlgebra.Generic.Res{fmpz_mod_poly}, e::AbstractAlgebra.Generic.Res{nmod_poly})
    ℤx = Nemo.ZZ["x"][1]
    @assert Nemo.lift(ℤx, modulus(parent(R))) == Nemo.lift(ℤx, modulus(parent(e)))
    parent(R)(Nemo.lift(ℤx, data(e)))
end

struct UsageError
    msg::String
end

# For now
Base.widen(::Type{Int128}) = Int256

end
