# Machine Learning on encrypted data

This example shows the evaluation of a simple convolutional neural network on
encrypted input data. It roughly follows the presentation in [JKLS19], but uses
a slightly different model (it adds back the weight biases) to be more consistent
with the standard MNIST example. Additionally, we only consider the case where
the data is encrypted, but the model is public. The intended application here
is that a cloud provider has a machine learning model they would like to offer
as an API, with the users sending their queries encrypted and receiving encrypted
answers (with the cloud provider not learning any information about the query).
In this scenario, there is little reason to also encrypt the model weights, so
we don't in this scenario. That more general scenario is however discussed in
[JKLS19] and the interested reader is encouraged to read that paper.

## Running the example

```julia
julia> include("train.jl")

julia> include("infer.jl")
```

[JKLS19] Xiaoqian Jiang, Miran Kim, Kristin Lauter, Yongsoo Song
         "Secure Outsourced Matrix Computation and Application to Neural Networks"
         https://eprint.iacr.org/2018/1041.pdf
