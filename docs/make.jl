using Documenter
using ToyFHE

makedocs(sitename="ToyFHE.jl", modules=[ToyFHE, ToyFHE.BFV],
pages = [
    "index.md",
    "Manual" => [
        "Background" => [
            "man/background.md",
            "man/background/rlwe.md",
            "man/background/fhe.md",
        ],
        "man/encoding.md"
    ]
])
