# Ring LWE and ℛq

## The ring of m-th cyclotomic integers mod q

All schemes implemented by this package make use of the RLWE setting.
In particular, the security assumption underlying the protocols is the general
[Ring Learning With Errors](https://en.wikipedia.org/wiki/Ring_learning_with_errors)
problem and the primary mathematical object of interest is the ring

```math
    \mathscr{R}_q = \frac{(\mathbb{Z}/q \mathbb{Z})[x]}{\Phi_m(x)}
```

where $\Phi_m$ is the m-th cyclotomic polynomial. In particular, we may represent
this ring as polynomials with coefficients in ``\mathbb{Z}/q \mathbb{Z}``, of
degree less than ``n = deg(\Phi_m(x)) = \phi(m)`` where ``\phi`` is the Euler [`totient`](@ref Primes.totient) function. Alternatively we may think of ``\mathscr{R}_q`` as the residue ring mod q of the
m-th cyclotomic ring of integers ``\mathscr{R} = \mathbb{Z} / \Phi_m(x)``.
We shall call ``q`` the ciphertext modulus, and ``n`` the ring dimension. In
general we will take the ciphertext modulus to be a large prime or a product of
somewhat large primes (the correctness of the crytopgraphic scheme depends on
choosing `q` such that the sum of the number of bits of the prime factors is
sufficiently large).

This definition may seem a bit abstract at first, but since it is the primary
mathematical object of interest, it is import that we get acquianted with it a
bit. Luckily for us, it is quite easy to play with these definitions using
the [Nemo.jl](https://github.com/Nemocas/Nemo.jl) stack of packages:

```jldoctest
using Nemo

# ℤ, the ring of integers
ℤ = Nemo.ZZ

# ℤ[x], the ring of polynomials with integer coefficients
ℤx, x = PolynomialRing(ℤ, "x")

# Nemo provides a function to compute m-th cyclotomics for us
Φ(m, x) = Nemo.cyclotomic(m, x)

# Let's take a look at first 10 cyclotomics
[Φ(m, x) for m = 1:10]

# output
10-element Array{fmpz_poly,1}:
 x-1
 x+1
 x^2+x+1
 x^2+1
 x^4+x^3+x^2+x+1
 x^2-x+1
 x^6+x^5+x^4+x^3+x^2+x+1
 x^4+1
 x^6+x^3+1
 x^4-x^3+x^2-x+1
```

There are two practical observations worth taking away from this exercise.

1. The degree of the m-th cyclotomic polynomial is indeed given by the Euler totient function:

```jldoctest
using Primes

ϕ(n) = Primes.totient(n)

[ϕ(m) for m = 1:10]

# output
10-element Array{Int64,1}:
 1
 1
 2
 2
 4
 2
 6
 4
 6
 4
```

2. When ``m`` is a power of two, the polynomial will have the form ``x^{m/2} + 1``. This holds generally, but let's just do the first couple of powers of two as examples:

```jldoctest
julia> [Φ(2^i, x) for i = 1:13]
13-element Array{fmpz_poly,1}:
 x+1
 x^2+1
 x^4+1
 x^8+1
 x^16+1
 x^32+1
 x^64+1
 x^128+1
 x^256+1
 x^512+1
 x^1024+1
 x^2048+1
 x^4096+1
```

For reasons having to do with the security properties of the resulting scheme,
`m` will generally be either prime or a power of two.

## Power of two cyclotomic rings

### Using Nemo.jl / Negacyclic convolutions
Because of their importance to practical implementations, this package contains
a custom implementation of power-of-two cyclotomic rings. However, before we get
there, let's build this ring manually using Nemo. To be concrete, we will chose
`m=2^3=8` (i.e. `n=2^2=4`) and `q=97`[1]:

```jldoctest
m = 2^3
q = 97

# ℤq = ℤ/qℤ, the ring of integers mod `q`
# Note that in this case `q` is a prime, so ℤq = 𝔽q is a finite field
ℤq = ResidueRing(ℤ, q)

# (ℤ/qℤ)[x], the ring of polynomials with integer coefficients mod `q`
ℤqx, xq = PolynomialRing(ℤq, "x")

# Finally our ring of interest (ℤ/qℤ)[x]/Φ_m(x)
ℛ = ResidueRing(ℤqx, Φ(m, x))

# output
Residue ring of Univariate Polynomial Ring in x over Integers modulo 31 modulo x^4+1
```

Now that we have the ring, let's construct a few polynomials in it.
```jldoctest
p1 = ℛ(x+1)
p2 = ℛ(x^3)
p3 = ℛ(4)
p4 = ℛ(5)
```

Note in particular that we have an embedding of ``\mathbb{F}_q`` into ``\mathscr{R}``
by simply embedding into the `x^0 = 1` term of the polynomial. Multiplying these
polynomials works as we'd expect:
```jldoctest
julia> p3*p4
20
```

and similarly polynomials whose product is of degree smaller than the degree of
the cyclotomic polynomial will behave as expected:
```jldoctest
julia> p1^2
x^2+2*x+1
```

but once we reach the degree of the cyclotomic, we will wrap around
```jldoctest
julia> p2*p1
x^3+96
```

which is expected since ``x^4 \equiv -1 (mod x^4 + 1)`` and
```math
(x+1)(x^3) = x^4 + x^3 \equiv x^3 - 1 \equiv x^3 + 96 mod (x^4 + 1, 97)
```

As we can see the behavior (at least for power of two cyclotomics) is fairly
simple: Once our multiplication wraps around, we simply subtract the resulting
coefficients starting at `x^0` again. The resulting operation is thus essentially
a cyclic convolution, except that wrapping around introduces an extra `-` sign.
We thus refer to this operation as a **negacyclic** convolution.

### The `NegacyclicRing` type

Because of the importance of power of two cyclotomics rings (i.e. negacyclic rings),
this package comes with a specialized implementation of these rings that provides
improved performance over the generic version using Nemo (as well as being written
in pure Julia). The [`NegacyclicRing`](@ref) type provides the entrypoint for
this functionality:

```
using GaloisFields
𝔽₉₇ = GaloisField(q)
ℛ = NegacyclicRing{𝔽₉₇, m ÷ 2}()

# output
NegacyclicRing{𝔽₉₇,4}(33)
```

We can repeat our experiment from above
```jldoctest
using OffsetArrays
O(a) = OffsetArray(𝔽₉₇.(a), 0:3)
p1 = ℛ(O([1, 1, 0, 0])) # x + 1
p2 = ℛ(O([0, 0, 0, 1])) # x^3
p3 = ℛ(O([4, 0, 0, 0])) # 4
p4 = ℛ(O([5, 0, 0, 0])) # 5
```

We can already see from the construction that `NegacyclicRing` ring uses a
different representation of the ring. Instead of polynomials, we think of ring
elements as arrays with a funny multiplication law (of course the two
representations are completely isomorphic). This is purely for convenience of
presentation when `m` is large (large arrays of numbers print more conveniently
than large sums of coefficients times ``x^i``). The results are analogous to
the what we obtained above using Nemo:
```jldoctest
julia> [p3*p4, p1^2, p1*p2]
3-element Array{ToyFHE.NTT.RingElement{NegacyclicRing{𝔽₉₇,4}(33),𝔽₉₇,Array{𝔽₉₇,1}},1}:
 [20, 0, 0, 0]
 [1, 2, 1, 0]
 [96, 0, 0, 1]
```

Internally, the `NegacyclicRing` type uses an analogue of the fast fourier transform
algorithm over finite fields to turn the (negacyclic) convolutions into pointwise
multiplications, thus performing these computations in `\mathcal{O}(n log n)` time
rather than `\mathcal{O}(n^2)` time. This transformation is known as the
negacyclic number theoretic transform or fermat number transform. See the corresponding [`section in the mathematical prerequisites chapter`](@ref nntt) for details.

[1] Note that `m` was chosen arbitrarily, but
`q` was chosen to be the first prime that such that `m` divies `q-1` (we will
see later why this criterion is important to allow a performant specialized
implementation).
