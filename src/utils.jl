module Utils

using Base.Meta
using GaloisFields
using GaloisFields: PrimeField

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
function fqmod(e::PrimeField, modulus::Integer)
    halfq = GaloisFields.char(e)/2
    if e.n > halfq
        mod(e.n - GaloisFields.char(e), modulus)
    else
        mod(e.n, modulus)
    end
end

end
