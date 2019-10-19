# Cryptographic Schemes for FHE

## FHE generally

A general (public-key) FHE scheme provides the following functions:

 - `(pub, priv) = keygen(...)` - Generates public and private keys
 - `cipher = encrypt(pub, plain)` - Encrypts `x` using the public key `pub`
 - `plain = decrypt(priv, cipher)` - Decrypts the ciphertext using private key `priv`
 - `cipher_new = eval(f, ciphers...)` - Compute the function `f` on the encrypted values

One particular important distinction between different FHE schemes is what
functions `f` are allowed to be homomorphically. In particular, in order to be
considered a **somewhat** or **fully** homomorphic encryption scheme (abbreviated
SHE and FHE respectively), `f` needs to range over at least two different values
that together allow some sort of universality (e.g. `+` and `*` on booleans
is sufficient). The difference between SHE and FHE schemes is that SHE schemes
may have a limit on the number of times that `eval` may be called, while FHE
schemes are defined not to have such a limit. It is possible to turn many SHE
schemes into FHE schemes using [`bootstrapping`](@ref)

## FHE using RLWE

Since we are in the RLWE setting, the basic primitives we have available are
polynomial multiplication (mod `q` and some cyclotomic) and polynomial addition
or equivalently negacyclic convolution and addition for length `n` arrays. All
three supported schemes `BGV`, `BFV` and `CKKS` allow evaluating negacyclic
convolution and addition on encrypted length `n` arrays, though they differ in
precisely how they embed their plaintext space into each polynomial.
Broadly speaking BGV encodes the plaintext into the (masked) lower bits of the
ciphertext, BFV encodes bits into the upper bits of the ciphertext and CKKS is,
well complicated, but somewhat similar to BFV. Nevertheless, they do share a lot
in common, so let's go over those parts and then revisit.

### Key generation for RLWE FHE schemes

The key structure is the same among all three schemes:
```
const ℛE = RingElement{ℛ}

struct PrivKey
    s::ℛE
end

struct PubKey
    mask::ℛE
    masked::ℛE
end
```

We now define three probability distributions ``\mathscr{U}`` the uniform
distribution over ``\mathscr{R}E``, ``\mathscr{G}`` whose probability has
discrete support near `0` and "small" standard deviation (common choices are
discrete gaussian distributions with standard deviation ``\frac{8}{\sqrt{2\pi}}``$
or the uniform distribution on ``\{-1, 0, 1\}``), as well as a noise distribution
``\mathscr{N}`` that is scheme specific and may or may not be the same as
``\mathscr{G}``.

We then compute
```julia
s = rand(𝒢)
mask = rand(𝒰)
masked = -(mask * s + rand(𝒩))
```

Note that that the problem of computing the private key `s` from the public key
`(mask, masked)` is precisely the ring learning with error problem (i.e.
the added noise makes computing `s` from `masked` hard, even if we do know
the `mask`).

!!! note

    We use the opposite sign for `masked` from what is found in the literature
    on the BGV scheme. The reason for this deviation is to unify the presentation
    for all three schemes (The papers on BFV and CKKS generally follow the same
    sign convention as were are using here).


### Encryption in RLWE FHE

A fresh cipher text consits of two ring elements `c₀` and `c₁` computed as
follows:

```julia
u = rand(𝒢)
c₀ = u * masked + rand(𝒩) + π⁻¹(plain)
c₁ = u * mask + rand(𝒩)
```

where `γ` is a scheme specific embedding function that maps the plaintext space
into the cipher text space.

### Homomorphic addition

Homomorphic addition is simple. For two ciphertext `(c₀¹, c₁¹)` and `(c₀², c₁²)`,
we simply add the ciphertexts elementwise and compute:

```julia
c₀′ = c₀¹ + c₀²
c₁′ = c₁¹ + c₁²
```

### Homomorphic multiplication

Homomorphic multiplication is generally similar between the schemes but differs
in the details. This presentation most closely matches BGV. For two ciphertext
`(c₀, c₁)` and `(d₀, d₁)`, we compute

```julia
c₀′ = c₀ * d₀
c₁′ = c₀ * d₁ + d₀ + c₁
c₂′ = c₁ * d₁
```

Note that the number of elements in the ciphertext increases. As a result,
we generally need to perform *relinerization* which is described below to
get the number of elements of the ciphertext back down to 2. However, it is
possible to decrypt ciphertexts with arbitrarily many elements. In general the
number of elements in the ciphertexts is the sum of the number of elements in
the input ciphertexts minus 1. Whether or not to do relinerization depends on
the tradeoff between the cost of relinerization and the cost of performing
multiplications with large numbers of ciphertext elements. The general formula
for multiplying two ciphertexts is:

```math
c_k\prime = \sum_{i+j == k} c_i * d_k
```

### Decryption

For a ciphertext `c = (c₀, ..., cⱼ)` the decryption formula is:

```julia
plain = π(sum(cᵢ * s^i for (i,cᵢ) in enumerate(c)))
```

where `π` is once again the scheme specific embedding we are familiar with from
the encryption step. Why does this work? Well, lets take it one step at a time.
Suppose we have a fresh encryption `(c₀, c₁)`. Then we have:

```
c₀ = u * masked + rand(𝒩) + π⁻¹(plain) = - u * mask * s - s * rand(𝒩) + π⁻¹(plain)
s * c₁ = u * mask * s + s * rand(𝒩)
```
and thus
```
c₀ + c₁ = - s * rand(𝒩) + s * rand(𝒩) + π⁻¹(plain)
```
We note that the two noise terms do not cancel, because they are independent
draws (they do however cancel in expectation of course). Nevertheless, it is
evident that one could choose `π` and `𝒩` appropriately such that `π` projects
out the noise terms (recall that `s` was chosen to have small coefficients). In
particular note that we have no term proportional to `mask`, which was sampled
from the uniform distribution.

Repeating this exercise for the ciphertext after multiplication we obtain:
```
c₀′ = (u * masked + e₁ + π⁻¹(m₁)) * (v * masked + e₂ + π⁻¹(m₂))
    = u * v * masked ^ 2 +
      masked * (u * e₂ + v * e₁ + u * π⁻¹(m₂) + v * π⁻¹(m₁)) +
      e₁ * e₂ +
      π⁻¹(m₂)π⁻¹(m₂)
    =  u * v * mask^2 * s^2 + u * v * e^2
     - mask * s * (u * e₂ + v * e₁ + u * π⁻¹(m₂) + v * π⁻¹(m₁))
     + e * (u * e₂ + v * e₁ + u * π⁻¹(m₂) + v * π⁻¹(m₁))
     + e₁ * e₂
     + π⁻¹(m₂)π⁻¹(m₂)
s * c₁′ = s * (u * masked + e₁ + π⁻¹(m₁)) * v * mask +
          s * (u * masked + e₁ + π⁻¹(m₁)) * f₂ +
          s * (v * masked + e₂ + π⁻¹(m₂)) * u * mask +
          s * (v * masked + e₂ + π⁻¹(m₂)) * f₁
        = 2 s * u * v * masked * mask
        + mask * s * (u * e₂ + v * e₁ + u * π⁻¹(m₂) + v * π⁻¹(m₁))
        + s * (u * masked + e₁ + π⁻¹(m₁)) * f₂
        + s * (v * masked + e₂ + π⁻¹(m₂)) * f₁
        = - 2 * u * v * s^2 * mask^2
          + mask * s * (u * e₂ + v * e₁ + u * π⁻¹(m₂) + v * π⁻¹(m₁))
          - s^2 * u * mask * f₂
          - s^2 * u * mask * f₁
          + s * (u * e + e₁ + π⁻¹(m₁)) * f₂
          + s * (v * e + e₂ + π⁻¹(m₂)) * f₁
          + 2 * u * v * s * e
s^2 * c₂′ = (u * mask + f₁) * (v * mask + f₂)
          = u * mask * mask^2 * s^2 +
            s^2 * u * mask * f₂ +
            s^2 * v * mask * f₁ +
            s^2 * f₁ * f₂
```

Summing these up, we get:
```
c₀′ + s * c₁′ + s^2 * c₂′ =
    π⁻¹(m₁)π⁻¹(m₂) +
    + e * (u * e₂ + v * e₁ + u * π⁻¹(m₂) + v * π⁻¹(m₁))
    + e₁ * e₂
    + s * (u * e + e₁ + π⁻¹(m₁)) * f₂
    + s * (v * e + e₂ + π⁻¹(m₂)) * f₁
    + 2 * u * v * s * e
    + s^2 * f₁ * f₂
```

where `e₁, e₂, e, f₁, f₂` are the various values samples from `𝒩`. Clearly
this decryption is more complicated. All the various noise terms are "small",
but we're multiplying a lot of them (and some by the embedded message), so the
noise is certainly significant. However, the important thing to notice is that
we have canceled all terms proportional to `mask`.

### Keyswitching / Modulus Switching

*TODO*

### BGV

We now turn to the three schemes themselves in order to specify the various
operations we had left open above. For, BGV, we have:

```
𝒩 = p * 𝒢
π⁻¹(plain) = plain
π(c) = mod(c, p)
```

It is easy to see that this works. In our expression above, ever value sampled
from `𝒩` is proportional to `p` and thus our decryption will be equal to
``m₁ * m₂ + p * noise \equiv m₁ * m₂ (mod p)`` (as long as the total magnitude
of ``p*noise`` is less than the ciphertext modulus ``q``).
We have thus constructed a scheme that can homomorphically evaluate
operations in `ℤ/pℤ`. Recall from the introduction that we said BGV embeds the
plaintext in the "low order bits" of the ciphertext. Perhaps a more precise
statement would be that in BGV the plaintext is the (masked) last digit of
the ciphertext in base `p`.

### BFV

BFV is essentially the opposite of BGV. We set:
```
𝒩 = 𝒢
π⁻¹(plain) = q/p * plain
π(c) = div(c, q/p, RoundToNearest)
```

Because of this encoding we make one additional modification to the scheme.
After every multiplication we divide by `Δ = q/p` (multiply by `p/q`) to rescale
the plaintext to be once again be proportional to `Δ`. In order to do this, we
amend the multiplication procedure as follows:
Before performing the multiplication, we switch the coefficients to some larger
modulus `q*t` where `t` is of same magnitude as (in number of bits) but slightly
larger `q` than q. As a result, we know that `(Δ * m₁) * (Δ * m₂) < Δ^2 m₁ m₂ < q^2 < q * t`
and thus no modular wraparound occurs (which would lose the top bits of our
plaintext).

Of course we could instead go all the way to ℤ. However, using `q*s` lets us
continue using the fast NTT-based arithmetic to the multiplication.

### CKKS

CKKS is special among the three schemes we support in that it is not exact, i.e.
we only ensure that the decrypted value approximately correct (for values of
approximately that depend on the parameters of the scheme). As such, we use
the perhaps paradoxically simple setup:

```
𝒩 = 𝒢
π⁻¹(plain) = plain
π(c) = c
```

We thus see that even for fresh encryptions, the low order digits lose precision.
Note that we do often still perform a BFV-like rescaling operation during
multiplication in CKKS as we tend to think of the plaintext as some fixed point
integer. Because the details are tricky and require additional functionality
we have not yet discussed, working with CKKS will be explained in its own
manual section later.
