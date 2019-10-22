using FFTW

struct CKKSEncoding{ScaleT, PlainT} <: AbstractVector{Complex{Float64}}
    # TODO: Should this just be looked up by type?
    data::OffsetVector{Complex{Float64}}
end
Base.size(e::CKKSEncoding) = Base.size(e.data)
Base.axes(e::CKKSEncoding) = Base.axes(e.data)
Base.getindex(e::CKKSEncoding, args...) = Base.getindex(e.data, args...)
Base.setindex!(e::CKKSEncoding, args...) = Base.setindex!(e.data, args...)

function CKKSEncoding{ScaleT}(plain::PlainT) where {ScaleT, PlainT}
    # Decode
    scaled = map(x->reinterpret(ScaleT, x), NTT.coeffs_primal(plain))
    # Undo the root of unity premul
    multed = map(x->convert(Float64, x), scaled) .* [Complex{Float64}(exp(big(-2*k/(2*length(scaled))*pi*im))) for k in eachindex(scaled)]
    # FFT it and take only the non-conjugated coefficients
    CKKSEncoding{ScaleT, PlainT}(OffsetArray(fft(collect(multed))[1:2048], 0:2047))
end

function Base.convert(::Type{NTT.RingElement}, s::CKKSEncoding{<:Any, PlainT}) where {PlainT}
    convert(PlainT, s)
end

function Base.convert(T::Type{<:NTT.RingElement}, s::CKKSEncoding{ScaleT, PlainT}) where {ScaleT, PlainT}
    data = collect(s.data)
    ipoints = OffsetArray(ifft([data; map(conj, reverse(data))]), 0:2*length(s.data)-1)
    # Make it negacyclic
    nipoints = ipoints .* [Complex{Float64}(exp(big(2*k/(2*length(ipoints))*pi*im))) for k in eachindex(ipoints)]

    # Project to real coefficents
    nipointsr = map(nipoints) do p
        @assert isapprox(imag(p), 0, atol=10^-10)
        real(p)
    end

    # Encode into the fixed fraction representation
    encoded = map(ScaleT, nipointsr)

    plaintext = map(x->x.x, encoded)
    PlainT(plaintext, nothing)
end
