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


plaintext_space(r::ResRing, p) = PolynomialRing(ResidueRing(Nemo.ZZ, p), "x")[1]
function plaintext_space(r::NegacyclicRing, p)
    coefft = Primes.isprime(p) ? GaloisField(p) :
        p == 256 ? UInt8 :
        Mod(p)
    if Primes.isprime(p) && p > 2degree(modulus(r))
        # TODO: Also needs to check here if the prime admits 2n-th roots of
        # unities.
        NegacyclicRing{coefft, degree(modulus(r))}(
            GaloisFields.minimal_primitive_root(coefft, 2degree(modulus(r))))
    else
        NegacyclicRing{coefft, degree(modulus(r))}()
    end
end

# HACK - DiscreteUniform is the default for types, but it'd be nice to
# allow this.
Distributions.DiscreteUniform(T::Type) = T

function fqmod(e::Union{PrimeField, Nemo.nmod, AbstractAlgebra.Generic.Res{fmpz}}, nq::Integer)
    q = modulus(e)
    halfq = q >> 1
    e = lift(e)
    if e > halfq
        mod(e - q, nq)
    else
        mod(e, nq)
    end
end


end
