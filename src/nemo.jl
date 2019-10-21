using Nemo
using AbstractAlgebra
using Distributions
using OffsetArrays

coefftype(ring::AbstractAlgebra.ResRing{T}) where T<:Union{nmod_poly,fmpz_mod_poly} = coefftype(ring.base_ring)
coefftype(ring::PolyRing) = ring.base_ring
Distributions.DiscreteUniform(R::Union{Nemo.NmodRing, AbstractAlgebra.Generic.ResRing}) = R
function (R::AbstractAlgebra.Generic.ResRing{T})(o::OffsetVector) where T<:Union{nmod_poly,fmpz_mod_poly}
    @assert first(axes(o)[1]) == 0
    R(base_ring(R)(collect(o)))
end
Base.oftype(x::AbstractAlgebra.Generic.Res, y::AbstractArray) = parent(x)(y)
function (R::NmodPolyRing)(o::OffsetVector)
    @assert first(axes(o)[1]) == 0
    R(collect(o))
end


struct CoeffView <: AbstractVector{Any}
    poly
    degree::Int
end
Base.getindex(p::CoeffView, i::Int) = coeff(p.poly, i)
Base.axes(p::CoeffView) = (0:p.degree,)
NTT.coeffs_primal(p::AbstractAlgebra.Generic.Res) = CoeffView(lift(p), degree(modulus(p)))
