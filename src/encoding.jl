using Primes

export SlotEncoding, ScalarEncoding

struct ScalarEncoding{T, R} <: AbstractArray{T, 0}
    plain::R
end

ScalarEncoding{T}(plain::R) where {T,R} =
    ScalarEncoding{T,R}(plain)

Base.axes(a::ScalarEncoding) = ()
Base.size(a::ScalarEncoding) = 1
Base.getindex(a::ScalarEncoding) = NTT.coeffs_primal(a.plain)[0]
Base.setindex!(a::ScalarEncoding, v) =  NTT.coeffs_primal(a.plain)[0] = v

function Base.convert(T::Type{<:NTT.RingElement}, s::ScalarEncoding)
    convert(T, s.plain)
end

function Base.showarg(io::IO, v::ScalarEncoding{T, R}, toplevel) where {T,R}
    # In the default case where we're just encoding into a ring
    # element, no need to make special note of it.
    if R <: NTT.RingElement
        print(io, "ScalarEncoding{", T, "}")
    else
        print(io, typeof(v))
    end
end

struct SlotEncoding{T, S, R} <: AbstractVector{T}
    plain::R
end

function SlotEncoding(r::R) where R<:NTT.RingElement{ℛ} where ℛ
    N = degree(ℛ)
    p = modulus(eltype(R))
    if Primes.isprime(p) && gcd(p-1, 2N) == 2N
        # TODO: Should allow CRTEncoded here also
        return SlotEncoding{eltype(R), eltype(R), R}(R(NTT.coeffs_dual(r), nothing))
    end
    error("Not supported yet")
end

Base.axes(s::SlotEncoding) = Base.axes(s.plain)
Base.size(s::SlotEncoding) = Base.size(s.plain)
Base.getindex(s::SlotEncoding, idxs...) = getindex(s.plain, idxs...)
Base.setindex!(s::SlotEncoding, args...) = setindex!(s.plain, args...)

function Base.convert(::Type{NTT.RingElement}, s::SlotEncoding)
    typeof(s.plain)(nothing, NTT.coeffs_primal(s.plain))
end

function Base.convert(T::Type{<:NTT.RingElement}, s::SlotEncoding)
    T(nothing, NTT.coeffs_primal(s.plain))
end

function SlotEncoding(r::R) where R<:NmodPolyRing
    f = modulus(r)
    p = modulus(base_ring(r))
end
