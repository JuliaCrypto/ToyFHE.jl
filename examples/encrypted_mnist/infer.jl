# Classifies MNIST digits with a convolutional network.
# Writes out saved model to the file "mnist_conv.bson".
# Demonstrates basic model construction, training, saving,
# conditional early-exit, and learning rate scheduling.
#
# This model, while simple, should hit around 99% test
# accuracy after training for approximately 20 epochs.

using Flux, Flux.Data.MNIST, Statistics
using Flux: onehotbatch, onecold, crossentropy, throttle
using Base.Iterators: repeated, partition
using Printf, BSON
using Base.Iterators
using LinearAlgebra
using ToyFHE
using Primes
using GaloisFields
using Test

using OffsetArrays

# Prepare test set in batches of 64
function make_minibatch(X, Y, idxs)
    X_batch = Array{Float32}(undef, size(X[1])..., 1, length(idxs))
    for i in 1:length(idxs)
        X_batch[:, :, :, i] = Float32.(X[idxs[i]])
    end
Y_batch = onehotbatch(Y[idxs], 0:9)
    return (X_batch, Y_batch)
end
batch_size = 64
test_imgs = MNIST.images(:test)
test_labels = MNIST.labels(:test)
mb_idxs = partition(1:length(test_labels), batch_size)
test_set =  [make_minibatch(test_imgs, test_labels, i) for i in mb_idxs]

# Load up the model we trained in train.jl
BSON.@load joinpath(dirname(@__FILE__), "mnist_conv.bson") model

function reshape_and_vcat(x)
    let y=reshape(x, 64, 4, size(x, 4))
        vcat((y[:,i,:] for i=axes(y,2))...)
    end
end

model = Chain(model.layers[1], reshape_and_vcat, model.layers[3:end]...)

function encrypted_matmul(weights, x)
    sum(diag(circshift(weights, (0,(k-1)))).*circshift(x, (k-1,0)) for k = 1:size(x, 1))
end

function naive_rectangular_matmul(weights, x)
    @assert size(weights, 1) < size(weights, 2)
    weights = vcat(weights, zeros(eltype(weights), size(weights, 2)-size(weights, 1), size(weights, 2)))
    encrypted_matmul(weights, x)
end

function public_preprocess(batch)
    ka = OffsetArray(0:7, 0:7)
    # Create feature extracted matrix
    I = [[batch[i′*3 .+ (1:7), j′*3 .+ (1:7), 1, k] for i′=ka, j′=ka] for k = 1:64]

    # Reshape into the ciphertext
    Iᵢⱼ = [[I[k][l...][i,j] for k=1:64, l=product(ka, ka)] for i=1:7, j=1:7]
end

function do_encrypted_inference(model, batch)
    Iᵢⱼ = public_preprocess(batch)

    # Evaluate the convolution
    weights = model.layers[1].weight
    conv_weights = reverse(reverse(weights, dims=1), dims=2)
    conved = [sum(Iᵢⱼ[i,j]*conv_weights[i,j,1,channel] for i=1:7, j=1:7) for channel = 1:4]
    conved = map(((x,b),)->x .+ b, zip(conved, model.layers[1].bias))
    sqed1 = map(x->x.^2, conved)
    sqed1 = map(x->reshape(x, 64, 64)', sqed1)
    fq1_weights = model.layers[3].W
    fq1 = sum(enumerate(partition(1:256, 64))) do (i,range)
        encrypted_matmul(fq1_weights[:, range], sqed1[i])
    end
    fq1 .+= model.layers[3].b
    sqed2 = fq1.^2
    fq2_weights = model.layers[4].W
    result = naive_rectangular_matmul(fq2_weights, sqed2)[1:10, :]
    result .+= model.layers[4].b
    result
end

batch = test_set[1][1]

# First test that the plaintext implementation above is correct
@test model(batch) ≈ do_encrypted_inference(model, batch)

# Set up the crypto parameters
# For now, we match the parameters of [JKLS19] and see how far we get

N = 2^13
q₀ = nextprime(Int128(2)^60 + 1, 1; interval=2N)
ps = nextprime(q₀ + 2N, 1; interval=2N)

q₁ = nextprime(2^40 + 1, 1; interval=2N)
q₂ = nextprime(q₁ + 2N, 1; interval=2N)
q₃ = nextprime(q₂ + 2N, 1; interval=2N)
q₄ = nextprime(q₃ + 2N, 1; interval=2N)
q₅ = nextprime(q₄ + 2N, 1; interval=2N)

ℛ = let CT = ToyFHE.CRTEncoded{7, Tuple{GaloisField.((q₀,q₁,q₂,q₃,q₄,q₅,ps))...}}
    ζ₂n = GaloisFields.minimal_primitive_root(CT, 2N)
    ToyFHE.NegacyclicRing{CT, N}(ζ₂n)
end

ckks_params = ModulusRaised(CKKSParams(ℛ, 0, 3.2))
kp = keygen(ckks_params)

Iᵢⱼ = public_preprocess(batch)

scale = Int128(2^40)
Tscale = FixedRational{scale}

C_Iij = map(Iᵢⱼ) do Iij
    plain = CKKSEncoding{Tscale}(zero(plaintext_space(ckks_params)))
    plain .= OffsetArray(vec(Iij), 0:(N÷2-1))
    encrypt(kp, plain)
end

weights = model.layers[1].weight
conv_weights = reverse(reverse(weights, dims=1), dims=2)
conved3 = [sum(C_Iij[i,j]*conv_weights[i,j,1,channel] for i=1:7, j=1:7) for channel = 1:4]
conved2 = map(((x,b),)->x .+ b, zip(conved3, model.layers[1].bias))
conved1 = map(ToyFHE.modswitch, conved2)

ek = keygen(EvalMultKey, kp.priv)
gk = keygen(GaloisKey, kp.priv; steps=64)

Csqed1 = map(x->x*x, conved1)
Csqed1 = map(x->keyswitch(ek, x), Csqed1)
Csqed1 = map(ToyFHE.modswitch, Csqed1)

decrypt(kp, Csqed1[1])

function encrypted_matmul(gk, weights, x::ToyFHE.CipherText)
    result = repeat(diag(weights), inner=64).*x
    rotated = x
    for k = 2:64
        @show k
        rotated = ToyFHE.rotate(gk, rotated)
        result += repeat(diag(circshift(weights, (0,(k-1)))), inner=64) .* rotated
    end
    result
end

decrypt_matrix(kp, xx) = reshape(collect(decrypt(kp, xx)), (64, 64))'

fq1_weights = model.layers[3].W
reshape(collect(decrypt(kp, encrypted_matmul(gk, fq1_weights[:, 1:64], Csqed1[1]))), (64, 64))'
fq1_weights[:, 1:64] * reshape(collect(decrypt(kp, Csqed1[1])), (64, 64))'

Cfq1 = sum(enumerate(partition(1:256, 64))) do (i,range)
    encrypted_matmul(gk, fq1_weights[:, range], Csqed1[i])
end

Cfq1 = Cfq1 .+ OffsetArray(repeat(model.layers[3].b, inner=64), 0:4095)
Cfq1 = modswitch(Cfq1)

Csqed2 = Cfq1*Cfq1
Csqed2 = keyswitch(ek, Csqed2)
Csqed2 = modswitch(Csqed2)

function naive_rectangular_matmul(gk, weights, x)
    @assert size(weights, 1) < size(weights, 2)
    weights = vcat(weights, zeros(eltype(weights), size(weights, 2)-size(weights, 1), size(weights, 2)))
    encrypted_matmul(gk, weights, x)
end

fq2_weights = model.layers[4].W
Cresult = naive_rectangular_matmul(gk, fq2_weights, Csqed2)
Cresult = Cresult .+ OffsetArray(repeat(vcat(model.layers[4].b, zeros(54)), inner=64), 0:4095)

plain_result = model(batch)
enc_result = real.(decrypt_matrix(kp, Cresult))

# N.B.: To see the images in your terminal, make sure to have `Images` and `TerminalExtensions` loaded,
# and make sure your terminal supports inline images (currently iTerm2 on OS X).
function compare_models(test_imgs, test_labels, plain_result, enc_result, plain_label, enc_label)
    print(" "^14)
    display(reduce(hcat, test_imgs[1:64]))
    println()
    print("Ground truth:  ")
    print(join(test_labels[1:64], " "))
    println()
    print(rpad(plain_label, 15))
    for idx = 1:64
        label = argmax(plain_result[:,idx])-1
        printstyled(label, color = label == test_labels[idx] ? :green : red)
        print(" ")
    end
    println()
    print(rpad(enc_label, 15))
    for idx = 1:64
        label = argmax(enc_result[:,idx])-1
        printstyled(label, color = label == test_labels[idx] ? :green : red)
        print(" ")
    end
    println()
end

compare_models(test_imgs, test_labels, model(batch), real.(decrypt_matrix(kp, Cresult)), "Model (Plain):", "Model (Enc):")
