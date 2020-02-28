using Test

@testset "examples" begin
    src = joinpath(dirname(dirname(@__FILE__)), "examples")
    dst = joinpath(mktempdir(), "examples")
    cp(src, dst; force = true)
    cd(dst)

    @testset "encrypted_mnist" begin
        dir = joinpath(dst, "encrypted_mnist")
        cd(dir)
        rm(joinpath(dir, "mnist_conv.bson"); force = true)
        include(joinpath(dir, "train.jl"))
        include(joinpath(dir, "infer.jl"))
    end
end
