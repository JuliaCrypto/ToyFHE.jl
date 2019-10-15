using Primes

export SlotEncoding

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

function Base.convert(T::Type{<:NTT.RingElement}, s::SlotEncoding)
    T(nothing, NTT.coeffs_primal(s.plain))
end
