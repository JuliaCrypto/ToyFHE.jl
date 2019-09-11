module Utils

using Base.Meta
using GaloisFields
using GaloisFields: PrimeField
using Distributions

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

# TODO: This isn't fully generic


# HACK - DiscreteUniform is the default for types, but it'd be nice to
# allow this.
Distributions.DiscreteUniform(T::Type) = T

end
