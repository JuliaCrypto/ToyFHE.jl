> I frequently hear music in the heart of noise. - George Gershwin

# ToyFHE.jl - A toy implementation of FHE algorithms

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
