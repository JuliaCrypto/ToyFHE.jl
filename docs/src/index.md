# ToyFHE.jl

*A pure-julia implementation of homomorphic encryption.*

## Overview

This package provides clean implementations of a variety of homomorphic
encryption protocols. Special attention is paid to clarity of presentation.
In particular, a goal of this package is that the implementations can be used
as references for the particular cryptographic protocol in question. This gives
rise to a couple of design decisions:

    1. Where possible code is shared between the various encryption protocols.
       Multiple dispatch makes it relatively easy to select an appropriate
       implementation for each protocol where they differ. However, as much
       as possible common code should be preferred to highlight similarities
       between the protocols.

    2. Performance optimizations are hidden behind type abstractions. As much
       as possible the main implementations for each protocol should match the
       textbook pseudocode for each scheme. Where optimizations are done for
       performance, these should be separated out into methods on a particular
       implementation type corresponding to that optimizations (or set of
       optimizations).

    3. The code should be kept small and modular.

In addition to modular arithmetic, the primary computational primitive needed
by the schemes implemented in this package is polynomial multiplication modulo
a certain cyclotomic polynomial. Two backends are available that provide this
primitive:

- A native julia one implemented in this package, based on
(FourierTransforms.jl)[https://github.com/JuliaComputing/FourierTransforms.jl].
- `libflint`, via the (Nemo.jl)[https://github.com/Nemocas/Nemo.jl] stack of packages.

The former is written in pure julia all the way down and
expected to be more performant and more easily portable to new hardware
architectures, but is restricted to power of two cyclotomics, while the latter
is ultimately implemented in `C`, but has support for arbitrary cyclotomics.
The implications of this will be discussed later and functionality that requires
arbitrary cyclotomics should be appropriatly marked in the documentation.

# Usage

Users simply wishing to try out the package are encouraged to check the
[Quickstart Guide](@ref quickstart) or the [Examples](@ref quickstart).
However, those intending to seriously use or extend the package are encouraged
to read the full manual:
```@contents
Pages = [
    "man/prerequisites.md",
    "man/background.md",
    "man/encoding.md",
    "man/ckks.md"
]
Depth = 1
```
