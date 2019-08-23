"""
    R-LWE is performed over rings of the form ℛ = 𝔽q[x]/Φₘ(x)
    where Φₘ (Φ_m) is the m-th cyclotomic polynomial (the security of RLWE
    relies on cyclotomicity). In practical implementations, for each of
    computation, we choose m such that Φₘ(x) = xⁿ+1. This map contains
    a pre-computed selection of such polynomials, mapping the cylcomatc order
    to a suitable RLWE ring.

    In addition, after fixing `n` (or equivalently the cyclotomic order `m`),
    we need to choose `q`, such that.

    This module consists of various standard parameter choices to match
    other pupular libraries.
"""
module CryptParameters

    using GaloisFields
    import ..LWERing

    """
        These parameters match PALISADE's defaults. See PALISADE's
        elementfactory.cpp.
    """
    const palisade_parameters = Dict(
        16 => LWERing{GaloisField(1099511627873), 8}(108163207722),
        1024 => LWERing{GaloisField(525313), 512}(513496),
        2048 => LWERing{GaloisField(34359724033), 1024}(7225104974),
        4096 => LWERing{GaloisField(1152921504606830593), 2048}(811032584449645127)
    )

end
