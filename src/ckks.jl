export CKKSParams, FixedRational

################################################################################
#                        CKKS Scheme definition
################################################################################

struct CKKSParams <: SHEShemeParams
    # The Cypertext ring over which operations are performed
    ℛ
    # The big ring used during multiplication
    ℛbig
    relin_window
    σ
end
scheme_name(p::Type{CKKSParams}) = "CKKS"

# From the RLWE perspective, ℛ is both the plain and ciphertext. The encoder
# takes care of the conversion to/from complex numbers.
ℛ_plain(p::CKKSParams) = p.ℛ
ℛ_cipher(p::CKKSParams) = p.ℛ

π⁻¹(params::CKKSParams, plaintext) = params.ℛ(plaintext)
π(params::CKKSParams, b) = b

𝒩(params::CKKSParams) = RingSampler(params.ℛ, DiscreteNormal(0, params.σ))
𝒢(params::CKKSParams) = RingSampler(params.ℛ, DiscreteNormal(0, params.σ))

#=
mul_expand(params::CKKSParams, c::CipherText) = map(c->switch(params.ℛbig, c), c.cs)
function mul_contract(params::CKKSParams, c)
    @fields_as_locals params::CKKSParams
    map(c) do e
        switch(ℛ, e)
    end
end
=#

################################################################################
#                        CKKS Scheme definition
################################################################################

"""
    FixedRational{den, T<:Integer} <: Real
Rational number type, with numerator of type T and fixed denominator `den`.
"""
struct FixedRational{denom, T}
    x::T
    Base.reinterpret(::Type{FixedRational{denom, T}}, x::T) where {T,denom} =
        new{denom, T}(x)
    Base.reinterpret(::Type{FixedRational{denom}}, x::T) where {T,denom} =
        new{denom, T}(x)
    function FixedRational{denom, T}(x::Real) where {T, denom}
        n = round(BigInt, big(x)*denom)
        if n < 0
            n = modulus(T) + n
        end
        new{denom, T}(convert(T, n))
    end
end
drop_last(::Type{FixedRational{denom, T}}) where {T, denom} = FixedRational{denom, drop_last(T)}
drop_last(::Type{FixedRational{denom}}) where {denom} = FixedRational{denom}


function Base.convert(::Type{Float64}, fr::FixedRational{denom, T}) where {T, denom}
    n = convert(Integer, fr.x)
    if n > div(modulus(T), 2)
        n = n - modulus(T)
    end
    Float64(n/denom)
end

function Base.show(io::IO, fr::FixedRational{denom}) where {denom}
    print(io, convert(Float64, fr))
end

function maybe_wide_mul(a::Integer, b::Integer)
    T = promote_type(typeof(a), typeof(b))
    c = widemul(a, b)
    isa(T, BigInt) && return Float64(c)
    c < typemax(T) && return T(c)
    isa(c, BigInt) && return Float64(c)
    return c
end
maybe_wide_mul(a, b) = a*b

function maybe_wide_sq(a::Integer, b::Integer)
    T = promote_type(typeof(a), typeof(b))
    c = big(a)^b
    isa(T, BigInt) && return Float64(c)
    c < typemax(T) && return T(c)
    isa(c, BigInt) && return Float64(c)
    return c
end
maybe_wide_sq(a, b) = a^b

Base.:^(a::Type{FixedRational{denom, T}}, b::Int64) where {T,denom} = FixedRational{maybe_wide_sq(denom, b), T}
Base.:^(a::Type{FixedRational{denom}}, b::Int64) where {T,denom} = FixedRational{maybe_wide_sq(denom, b)}
Base.:*(a::Type{FixedRational{denom, T}}, b::Number) where {T,denom} = FixedRational{b/denom, T}
Base.:*(a::Type{FixedRational{denom}}, b::Number) where {T,denom} = FixedRational{b/denom, T}
Base.:*(a::Type{FixedRational{denom1, T}}, b::Type{FixedRational{denom2, T}}) where {T,denom1,denom2} = FixedRational{maybe_wide_mul(denom1, denom2), T}
Base.:*(a::Type{FixedRational{denom1}}, b::Type{FixedRational{denom2}}) where {denom1,denom2} = FixedRational{maybe_wide_mul(denom1, denom2)}
StructArrays.createinstance(FR::Type{FixedRational{denom, T}}, x::T) where {T,denom} =
    reinterpret(FR, x)
