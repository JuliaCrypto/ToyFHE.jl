import Base: +, -, *

"""
Represents `mod m` modular arithmetic when identified as (-m/2,m/2]. The
identification does not change the way arithmetic is performed, but does
change conversions to other types.
"""
struct SignedMod{T} <: Number
    x::T
end

function Base.convert(::Type{Integer}, x::SignedMod)
    n = convert(Integer, x.x)
    if n > div(modulus(x.x), 2)
        return n - modulus(x.x)
    else
        return n
    end
end

Base.oftype(a::SignedMod, b) = SignedMod(oftype(a.x, b))
Base.promote_rule(::Type{SignedMod{T}}, ::Type{S}) where {T, S} = SignedMod{promote_type(T, S)}

for f in (:+, :-, :*)
    @eval $f(a::SignedMod{T}, b::SignedMod{T}) where {T} = SignedMod{T}($f(a.x, b.x))
end

Base.:*(a::SignedMod, b::Integer) = a*oftype(a, b)

function Base.div(e::SignedMod, x::Integer, r::RoundingMode)
    oftype(e, div(convert(Integer, e), x, r))
end

function Base.mod(e::SignedMod, x::Integer)
    oftype(e, mod(convert(Integer, e), x))
end

function Base.show(io::IO, s::SignedMod)
    show(io, convert(Integer, s))
end
