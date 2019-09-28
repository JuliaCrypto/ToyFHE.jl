
# Trivial BFV, scalar encoding small plaintext space, pow2 cyclotomic ring
include("bfv_triv.jl")

# Trivial BGV
include("bgv_triv.jl")

# BFV with p=65537 SIMD over cyclotomic ring
include("bfv_simd.jl")

# BFV with p=256 SIMD over non-cyclotomic ring
include("bfv_uint8.jl")

# Keyswitching for BFV
include("bfv_keyswitch.jl")

# Noise measurement for bfv
include("bfv_noise.jl")
