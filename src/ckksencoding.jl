using FFTW

struct CKKSEncoding{ScaleT} <: AbstractVector{Complex{Float64}}
    # TODO: Should this just be looked up by type?
    PlainT::Type
    # We permute the elements of the encoding in order to make the galois group
    # act cyclically on the plaintext elements.
    data::OffsetVector{Complex{Float64}}
end
Base.size(e::CKKSEncoding) = Base.size(e.data)
Base.axes(e::CKKSEncoding) = Base.axes(e.data)
Base.getindex(e::CKKSEncoding, i::Integer) = Base.getindex(e.data, i)
Base.setindex!(e::CKKSEncoding, v, i::Integer) = Base.setindex!(e.data, v, i)
scale(::Type{CKKSEncoding{ScaleT}}) where {ScaleT} = ScaleT
drop_last(::Type{CKKSEncoding{ScaleT}}) where {ScaleT} = CKKSEncoding{drop_last(ScaleT)}

"""
    ℤmstarPermutation

Recall that in a packed encoding the individual plaintext slots `c_i` correspond
to evaluation of the polynomial `p` at the i-th power of some m-th root of unity
``ζ_m``` for ``i ∈ ℤ_m^\\star`` where ``ℤ_m^\\star`` is the multiplicative group
of integers modulo ``m`` [1]. Alternatively, for m a power of 2, recall that
the ``c_i`` are the length ϕ(n) negacyclic fft of the coefficents of `p`
(the correspondence between these two is explained in the manual). This
identification is useful, because `ℤ_m^\\star` generally decomposes as a product
of cyclic groups. For example for `m=2^k` a power of two, we have:

```math
    \\mathbb{Z}_m^\\star \\simeq C_2 \\times C_{2^{k-2}}
```

Another useful representation of this cyclic group is a `2x(k-2)` matrix of
plaintext slots where the group action acts as as [`circshift`](@ref) on the
rows and columns respectively (corresponding to the action in the ``C_2``
and ``C_{2^{k-2}}`` factors respectively).

This type implements that representation, arranging the elements of
``ℤ_m^\\star`` (as integers) in the mentioned matrix representation.

[1] https://en.wikipedia.org/wiki/Multiplicative_group_of_integers_modulo_n
"""
struct ℤmstarPermutation <: AbstractArray{Int, 2}
    M::Int
    ℤmstarPermutation(m::Integer) = (@assert ispow2(m); new(convert(Int,m)))
end

Base.size(p::ℤmstarPermutation) = (2, div(p.M, 4))
function Base.getindex(p::ℤmstarPermutation, row::Integer, col::Integer)
    ℤmstar_generator = 3
    rs, cls = size(p)
    (row in 1:rs && col in 1:cls) || Base.throw_boundserror(p, (row, col))
    mod(row == 2 ? p.M-ℤmstar_generator^col : ℤmstar_generator^col, p.M)
end

function CKKSEncoding{ScaleT}(plain::PlainT) where {ScaleT, PlainT}
    # Decode
    scaled = map(x->reinterpret(ScaleT, x), NTT.coeffs_primal(plain))
    scaled = Float64[convert(Float64, x) for x in scaled]
    # Undo the root of unity premul
    multed = map(x->convert(Float64, x), scaled) .* [Complex{Float64}(exp(big(-2*k/(2*length(scaled))*pi*im))) for k in eachindex(scaled)]
    # FFT it and take only the non-conjugated coefficients
    plain = OffsetArray(fft(collect(multed)), 0:length(multed)-1)[ℤmstarPermutation(2length(multed))[1,:] .>> 1]
    plain = OffsetArray(plain, 0:(length(multed)÷2-1))
    CKKSEncoding{ScaleT}(PlainT, plain)
end

function Base.convert(::Type{NTT.RingElement}, s::CKKSEncoding{<:Any})
    convert(s.PlainT, s)
end

function Base.convert(T::Type{<:NTT.RingElement}, s::CKKSEncoding{ScaleT}) where {ScaleT}
    data = collect(s.data)

    cmplx = OffsetArray(zeros(Complex{Float64}, 2*length(s.data)), 0:2*length(s.data)-1)
    idxs = ℤmstarPermutation(4length(s.data)) .>> 1
    for i in axes(s.data, 1)
        cmplx[idxs[1, i+1]] = s.data[i]
        cmplx[idxs[2, i+1]] = conj(s.data[i])
    end

    ipoints = OffsetArray(ifft(cmplx.parent), 0:2*length(s.data)-1)
    # Make it negacyclic
    nipoints = ipoints .* [Complex{Float64}(exp(big(2*k/(2*length(ipoints))*pi*im))) for k in eachindex(ipoints)]

    # Project to real coefficents
    nipointsr = map(nipoints) do p
        @assert isapprox(imag(p), 0, atol=10^-10)
        real(p)
    end

    # Encode into the fixed fraction representation
    encoded = map(ScaleT{eltype(s.PlainT)}, nipointsr)

    plaintext = map(x->x.x, encoded)
    s.PlainT(plaintext, nothing)
end


# Multiplication with plaintext scalars
function Base.:*(a::CipherText{CKKSEncoding{Tscale}, P, T, N}, b::AbstractFloat) where {Tscale, P, T, N}
    scaled = convert(Integer, Tscale(b).x)
    CipherText{CKKSEncoding{Tscale^2}, P, T, N}(a.params, map(c->c*scaled, a.cs))
end

function Base.broadcasted(::typeof(*), a::Array{<:AbstractFloat, 1}, c::CipherText{CKKSEncoding{Tscale}, P, T, N}) where {Tscale, P, T, N}
    plain = CKKSEncoding{Tscale}(zero(c.cs[1]))
    plain .= OffsetArray(a, 0:length(a)-1)
    re = convert(NTT.RingElement, plain)
    CipherText{CKKSEncoding{Tscale^2}, P, T, N}(c.params, map(c->c*re, c.cs))
end

function Base.broadcasted(::typeof(+), a::CipherText{CKKSEncoding{Tscale}, P, T, N}, b::AbstractFloat) where {Tscale, P, T, N}
    plain = CKKSEncoding{Tscale}(zero(T))
    plain .= b
    CipherText{CKKSEncoding{Tscale}, P, T, N}(a.params, (a.cs[1] + convert(T, plain), a.cs[2:end]...))
end

# Modswitching a ciphertext
function modswitch(c::CipherText{CKKSEncoding{Tscale}, P, T, N}) where {scale, Tscale<:FixedRational{scale}, CC<:CRTEncoded, P, T<:NTT.RingElement{<:Any, CC}, N}
    CipherText{CKKSEncoding{drop_last(FixedRational{scale/modulus(moduli(CC).parameters[end])})}}(modswitch_drop(c.params),
        map(modswitch, c.cs))
end

# Multiplication of ciphertexts
function *(c1::CipherText{C1, Enc, T}, c2::CipherText{C2, Enc, T}) where {C1 <: CKKSEncoding, C2 <: CKKSEncoding, Enc, T}
    CipherText{CKKSEncoding{scale(C1)*scale(C2)}}(c1.params, enc_mul(c1, c2))
end
