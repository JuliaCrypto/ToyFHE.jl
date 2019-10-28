# ToyFHE.jl - A toy implementation of FHE algorithms

> I frequently hear music in the heart of noise. - George Gershwin


| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![][docs-dev-img]][docs-dev-url] | [![][travis-img]][travis-url] [![][codecov-img]][codecov-url] |

## Installation

This package currently depends on a number of modifications
to upstream packages. The included `Manifest.toml` lists
known working versions for these packages. To use those
versions, first clone this package to a location of your
choice:

```
$ git clone https://github.com/JuliaComputing/ToyFHE.jl ToyFHE
```

Then load up the project within Julia:
```
$ julia --project=ToyFHE
```

If you do not have the correct versions of the dependencies installed, you may be asked to install them via `instantiate`.

## Documentation

- [**DEVEL**][docs-dev-url] &mdash; *documentation of the in-development version.*

## Features

**WARNING**: The code in this repository is in an extremely alpha quality. You
may want to check back in a little bit once things have been cleaned up.

This repository is a playground for exploring homomorphic encryption protocols.
The design goal is ease of use and ease of readability over absolute performance
or suitability for production HE applications. The goal of this code is to aid
in research and the quick exploration of HE applications.

This package currently contains (partial) implementations of the following HE
schemes:

- Brakerski/Fan-Vercauteren (BFV)
- Brakerski-Gentry-Vaikuntanathan (BGV)
- Cheon-Kim-Kim-Song (CKKS)

Both power-of-two and general cyclotomic rings are supported for homomorphic
operations. The former is based on a pure Julia FFT implementation and thus
likely suitable for multi-threading and GPU applications with little additional
effort (those this has not been done so far). The latter is based on the
[Nemo](http://nemocas.org/) stack of Julia packages, which are ultimately using
[FLINT](http://www.flintlib.org/) as the execution engine.

# Disclaimer
## Performance notice

This package has not been optimized for performance. The only implementated
performance optimizations are those that were absolutely required to perform
the desired algorithmic exploration. PRs are welcome to improve performance
(as long as readability is preserved), but such work is not currently on the
roadmap.

## Security Notice

This package currently has known issues (weak RNG, known timing side channels)
that make it unsuitable for use other than for algorithmic research. DO NOT USE
FOR PRODUCTION APPLICATIONS (I mean it). As with the performance consideration,
these issues are addressable, but not currently on the roadmap.

## Omnibus

I am not a cryptographer. I am not your cryptographer. This code is a toy.
This code is not intended for production use. This code has not been audited
or validated. Consult a professional cryptographer, before using cryptography.
No warranty (see LICENSE). Don't sue me.


[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://juliacomputing.github.io/ToyFHE.jl/dev

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://juliacomputing.github.io/ToyFHE.jl/stable

[travis-img]: https://travis-ci.org/JuliaComputing/ToyFHE.jl.svg?branch=master
[travis-url]: https://travis-ci.org/JuliaComputing/ToyFHE.jl

[codecov-img]: https://codecov.io/gh/JuliaComputing/ToyFHE.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/JuliaComputing/ToyFHE.jl

[issues-url]: https://github.com/JuliaComputing/ToyFHE.jl/issues
