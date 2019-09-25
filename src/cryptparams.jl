"""
    R-LWE is performed over rings of the form ℛ = 𝔽q[x]/Φₘ(x)
    where Φₘ (Φ_m) is the m-th cyclotomic polynomial (the security of RLWE
    relies on cyclotomicity). In practical implementations, for each of
    computation, we choose m such that Φₘ(x) = xⁿ+1.

    This module consists of various standard parameter choices to match
    other popular libraries and standards.
"""
module CryptParameters

    using GaloisFields
    import ..NegacyclicRing
    export StdDistribution, HEStd_uniform, HEStd_error, HEStd_ternary,
        StdSecurity, HEStd_128_classic, HEStd_192_classic, HEStd_256_classic

    """
        These parameters match PALISADE's BGV defaults. See PALISADE's
        elementfactory.cpp.
    """
    const palisade_parameters = Dict(
        # m => ℛ
        16 => NegacyclicRing{GaloisField(1099511627873), 8}(108163207722),
        1024 => NegacyclicRing{GaloisField(525313), 512}(513496),
        2048 => NegacyclicRing{GaloisField(34359724033), 1024}(7225104974),
        4096 => NegacyclicRing{GaloisField(1152921504606830593), 2048}(811032584449645127)
    )

    @enum(StdDistribution,
        HEStd_uniform,
        HEStd_error,
        HeStd_ternary)

    @enum(StdSecurity,
        HEStd_128_classic = 1,
        HEStd_192_classic = 2,
        HEStd_256_classic = 3
    )

    struct StdParam
        dist::StdDistribution
        n::Int
        security::StdSecurity
        ringDim::Int
    end
    StdParam((a,b,c,d),) = StdParam(a,b,c,d)

    # Table from http://homomorphicencryption.org/wp-content/uploads/2018/11/HomomorphicEncryptionStandardv1.1.pdf
    const std_n = [1024, 2048, 4096, 8192, 16384, 32768]
    const std_params = Dict(
        HEStd_uniform => Dict(
             1024 => ( 29,  21,  16),
             2048 => ( 56,  39,  31),
             4096 => (111,  77,  60),
             8192 => (220, 154, 120),
            16384 => (440, 307, 239),
            32768 => (880, 612, 478)
        ),
        HEStd_error => Dict(
             1024 => ( 29,  21,  16),
             2048 => ( 56,  39,  31),
             4096 => (111,  77,  60),
             8192 => (220, 154, 120),
            16384 => (440, 307, 239),
            32768 => (883, 613, 478)
        ),
        HeStd_ternary => Dict(
            1024 => ( 27,  19,  14),
            2048 => ( 54,  37,  29),
            4096 => (109,  75,  58),
            8192 => (218, 152, 118),
           16384 => (438, 305, 237),
           32768 => (881, 611, 476)
        )
    )

    function std_ring_dim(dist_type, minSecLevel, curLogQ)
        std_n[findfirst(n->std_params[dist_type][n][Int(minSecLevel)] >= curLogQ, std_n)]
    end
end
