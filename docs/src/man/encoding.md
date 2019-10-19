# Ciphertext Encodings

In the [background section](@ref background), we learned how to homomorphically perform operations over the ring ``\mathscr{R}_q``.
However, in general we may want to perform computations over mathematical spaces other than ``\mathscr{R}_q``. To do, we must devise a scheme and selection of parameters `p,q,m` that allow us
to efficiently **encode** our object of interested into an element
of ``\mathscr{R}_q`` such that arithmetic in ``\mathscr{R}_q``
corresponds to whatever computation we're interested in for
the encoded object. In this section we shall see several strategies for encoding numbers into single ciphertexts, while in the next section we shall consider the encoding of more complex higher level objects into multiple ciphertexts. To make the above mentioned parameterization explicit, we will use the notation ``\mathscr{R}_{q,p,m}`` in this section rather than leaving the parameters ``p,m`` implicit.

## Scalar Encodings
As noted in the [background section](@ref rlwe), there is a natural embedding of ``\mathbb{Z}/p \mathbb{Z}`` into ``\mathscr{R}_{q,p,m}`` (for arbitrary `q`,`m`, though of course we require `q>p` and security considerations will in general force larger values).

```jldoctest
ℛ = NegacyclicRing{GaloisField(7), 2}(nothing)
p1, p2 = zero(ℛ), zero(ℛ)
p1[0] = 3
p2[0] = 4
p1*p2

# output
2-element ToyFHE.NTT.RingElement{NegacyclicRing{𝔽₇,2}(0),𝔽₇,Array{𝔽₇,1}} with indices 0:1:
 5
 0
```

To make the use of this encoding explicit, there is an encoding wrapper type that simply perfoms this encoding (in this example with UInt8, which in Julia is defined as ``\mathbb{Z}/256\mathbb{Z}``)

```jldoctest
ℛ = NegacyclicRing{UInt8, 2}(nothing)
a, b = ScalarEncoding{UInt8}(zero(ℛ)), ScalarEncoding{UInt8}(zero(ℛ))
a[] = 3
b[] = 15
convert(RingElement, a) * convert(RingElement, b)

# output
2-element RingElement{NegacyclicRing{UInt8,2}(0x00),UInt8,Array{UInt8,1}} with indices 0:1:
 0x2d
 0x00
```

## SIMD Encoding

In general, the ring dimension ``n`` will be quite large (e.g. 2048, 4096 for small examples, or even larger for programs of practical interest). As such, it seems wasteful to only use the lowest coefficient for computation, while the others just sit around warming the planet. And of we can always use the other coefficients if the operations if we're ok with ``\mathscr{R}_q`` arithmetic. However,
in practical operations we are much more interested in pointwise operations on vectors rater than negacyclic convolutions of vectors.
There are several techniques to remidy this disconnect and let us treat a single ``\mathscr{R}_q`` as multiple independent plaintexts
with both addition and multiplication being applied elementwise,
similar to SIMD lanes in a modern CPU or GPU. In the literature,
these independent plaintexts are generally referred to as plaintext "slots". In practice, how to go between the vector of slots and an element in ``\mathscr{R}_q`` depends on the which parameter domain we're in, though they are all just applications of the chinese remainder theorem for polynomials.

### SIMD Encoding using NTT
You may recall from the [background section](@ref background) that we are using a mathematical trick known as the Number Theory Transform (NTT) to efficiently perform negacyclic convolutions as pointwise multiplications. Naturally,
we can also apply this trick in reverse and use the fact that we
can perform negacyclic convolutions to obtain a method to perform
element-wise multiplications. The catch here is that for this to work,
we need to be able to do negacyclic convolutions in *both* the plaintext space and the ciphertext space (we always need to be able to do the latter, but in the general scheme, we can chose the plaintext modulus `p`). Since the ring dimension is often quite large, this
method may not be applicable to values of `p` of interest, or may force us to use a value of `p` larger than we intended. Recall thus the condition for being able to apply the NTT:

    If the prime factors of p are `p = p₁ ⋯ pⱼ`, the NTT exists iff the `pᵢ` are all distrinct and 2n divides `pᵢ` for all `i` in `1:j`.

Since in practice `N` is generally in the thousands to tens of thousands and is a power of two, the available values of `p` is fairly
restricted. Nevertheless if the plaintext space of interest is e.g.
machine integers (of 16, 32 or 64 bit width) and lack of overflow is
guaranteed, this encoding may be an attractive option. In particular note that $65537 = 2^16 + 1$ is prime and can thus be encoded with minimal overhead. Note further that this methods gives us exactly ``n`` plaintext slots, so is a very efficient method of encoding such integers into a ciphertext.

This package provides a convenience wrapper that acts as a view
into a plaintext, but implicity performs the NTT:

```jldoctest
ℛ = NegacyclicRing{GaloisField(65537), 2048}()

a = SlotEncoding(zero(ℛ))
a[0:9] = 1:10

b = SlotEncoding(zero(ℛ))
b[:] .= 10

SlotEncoding(convert(RingElement, a) * convert(RingElement, b))[0:10]

# output
11-element Array{𝔽₆₅₅₃₇,1}:
  10
  20
  30
  40
  50
  60
  70
  80
  90
 100
   0
```

### SIMD Encoding over general cyclotomics

The example in the previous section is a special case of a more general technique: The polynomial chinese remainder theorem.
Just like its analogue for integers (which is discussed in the [mathematical background](@ref crt) and used for the [CRT encoding of integer coefficients](@ref crt-encoding)), if a polynomial f has coprime factors `f = f₁ ⋯ fⱼ`, then there is an ismorphism between
polynomials mod `f` and polynomials mod each of the factors. As a
result, if ``\Phi_m(x) \in \mathbb{Z}_p[x]`` splits as ``\Phi_m(x) = f_1(x) \cdots f_j(x)``, this induces an ismorphism

```math
\mathbb{Z}_p[x]/\Phi_m(x) \cong \mathbb{Z}_p[x]/f_1(x) \times \cdots \times \mathbb{Z}_p[x]/f_j(x)
```

We can now see the correspondence between this definition and
the method of the previous subection. Recall that one of the definitions of the m-th cyclotomic polynomial is

```math
\Phi_m(x) = \Pi_{k=1}^{\phi(m)} (x - \zeta_{m,k})
```

where the ``\zeta_{m,k}`` are the m-th [primitive roots of unity](@ref primitive-roots). The definition of the previous section simply requires that ``\zeta_{m,k} \in \mathbb{Z}_p`` (i.e. that ``\mathbb{Z}_p`` contains the m-th primitive roots of unity), and so by the chinese remainder theorem our ring factors into ``\phi(m) = n`` copies of ``\mathbb{Z}_p``. So now let us consider what happens when the setting of the previous section does not apply.

For concreteness, let us suppose we choose `m = 7`,`p=2`. We have (for a more general introduction to cyclotomics using Nemo, see the [RLWE section of the background](@ref background-rlwe)):

```jldoctest
using Nemo
ℤ = Nemo.ZZ
ℤx, x = PolynomialRing(ℤ, "x")
ℤpx, xp = PolynomialRing(ResidueRing(ℤ, 2), "x")
Φ(m, x) = Nemo.cyclotomic(m, x)

Nemo.factor(ℤpx(Φ(7, x)))

# output

1 * (x^3+x+1) * (x^3+x^2+1)
```

Thus, we see that ``\mathbb{F}_2/\Phi_7(x) \cong \mathbb{F}_2/(x^3 + x + 1) \times \mathbb{F}_2/(x^3 + x^2 + 1) \cong \mathbb{F}_{2^3}^2`` (i.e. two copies of the finite field ``\mathbb{F}_{2^3}``). It is these two copies that we can independently address
in a SIMD manner. Of course we may also use any subfields of this factorization,
and in particular we may use the subfield ``\mathbb{F}_2`` if we want to
perform binary arithmetic.

In summary: We get one copy of ``\mathbb{F}_{p^r}`` for every *distinct* irreducible
factor of ``\Phi_m(x)`` over ``\mathbb{F}_p``.

The `PolyCRTEncoding` type knows how to perform these isomorphisms and will
automatically apply them when possible:

```jldoctest
const ℤ = ZZ
const ℤx, x = ZZ["x"]
const ℤpx, xp = ResidueRing(ZZ, 2)["xp"]
ℛ = ResidueRing(ℤpx, cyclotomic(7, x))

plain = PolyCRTEncoding(zero(ℛ))
β = GaloisFields.gen(eltype(plain))

plain[1] = β+1
plain[2] = β^2+1

p = convert(ℛ, plain)

# output
xp^4+xp^3+xp^2+xp
```

An of course as before it works the inverse way also:
```jldoctest
julia> PolyCRTEncoding(ℛ(xp^4+xp^3+xp^2+xp))
2-element PolyCRTEncoding{𝔽₈}:
   β + 1
 β^2 + 1
```

### SIMD Encoding from Hensel lifting

### SIMD Summary
