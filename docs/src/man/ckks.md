# The CKKS encryption scheme

The CKKS scheme may seem a bit strange at first, but it is actually a quite
straightforward extension of the ideas we have already seen. Rather than
describing what it is, let us explore how one might come to making the same
choices.

Let's say you've just finished writing for wonderful implementation of BGV,
or BFV. Basking in the glory of your accomplishments you proudly show off
encrypted integer multiplies to your boss. "Look", you say excitedly, "it
multiplied those numbers even though they were encrypted. It even did
[a few thousand of those at once](encoding.md)." "That's nice", your boss
replies "now make it do ℝ; ℤ is so last century". Discouraged you return to
your desk - there goes the weekend.

## Fixed point encoding

What now? A faint memory dawns of the dark ages before IEEE754. Back when the
world was black and white and the real numbers were fixed point. You dust off
the trusty fixed point number library and try it:

```
using ToyFHE
using FixedPointNumbers

params = BFVParams(
    65537, # plaintext modulus
    ; eval_mult_count = 1
)
kp = keygen(params)

a = ScalarEncoding{Int16}(zero(plaintext_space(params)))
a[] = reinterpret(Int16, 1.5Q11f4)

b = ScalarEncoding{Int16}(zero(plaintext_space(params)))
b[] = reinterpret(Int16, 2.0Q11f4)

c1, c2 = encrypt.(kp, (a, b))
y = c1*c2
reinterpret(Fixed{Int16,8}, convert(Integer, decrypt(kp, y)[0]) % Int16)

# output

3.0Q7f8
```

Hey that worked! One particular thing to note is that multiplication increased
the size of the fraction. We started off with a `Q11f4` (i.e. 11 quotient bits,
4 fraction bits) and ended up with a `Q7f8`. Managing the size of the fraction
will become important. We'll see that below.

Alright, let's keep going. Can we do SIMD?

```
a = SlotEncoding(zero(plaintext_space(params)))
b = SlotEncoding(zero(plaintext_space(params)))
plain_range = LinRange(0.0, 2.0, 4096)
a[:] .= reinterpret.(Int16, Fixed{Int16,6}.(OffsetArray((plain_range, 0:4095)))
b[:] .= reinterpret(Int16, Fixed{Int16,6}(2.0Q11f4))


c1, c2 = encrypt.(kp, (a, b))
y = c1*c2
# Let's just look at every fourth output
reinterpret.(Fixed{Int16,12}, convert.(Integer, SlotEncoding(decrypt(kp, y))) .% Int16)[1:4:end]

# output
512-element Array{Q3f12,1} with eltype Fixed{Int16,12}:
 0.0Q3f12
 0.0Q3f12
 0.0312Q3f12
 0.0312Q3f12
 0.0312Q3f12
 0.0312Q3f12
 ⋮
 3.9688Q3f12
 3.9688Q3f12
 4.0Q3f12
 4.0Q3f12
```

Ok, that kinda works, but even with 6 bits in the fraction, we don't really have
enough precision to take advantage of the whole 4096 plaintext slots we have.
Can we come up with something better?

## Noisy encoding

So here comes one of the main ideas of CKKS: What if instead of letting all
those lower bits just sit around and accumulate noise, we used them in the
computation? Sure, they'll get gradually overwritten by noise, but at least
we can take advantage of the higher precision in the first couple intermediate
multiplications. Let's see what we can do. We define our new scheme (partially
extracted from the ckks.jl source file):

```
struct CKKSParams <: SHEShemeParams
    # The Cypertext ring over which operations are performed
    ℛ
    # The big ring used during multiplication
    ℛbig
    relin_window
    σ
end
scheme_name(p::Type{CKKSParams}) = "CKKS"

# Just, like, put the plaintext right there
π⁻¹(params::CKKSParams, plaintext) = params.ℛ(plaintext)
π(params::CKKSParams, b) = b
```

so let's try it:
```
ckksp = CKKSParams(params.ℛ, params.ℛbig, params.relin_window, params.σ)
kp = keygen(ckksp)
a = ScalarEncoding{eltype(params.ℛ)}(zero(plaintext_space(ckksp)))
a[] = reinterpret(Int128, Fixed{Int128, 35}(plain_range[2]))

b = ScalarEncoding{eltype(params.ℛ)}(zero(plaintext_space(ckksp)))
b[] = reinterpret(Int128, Fixed{Int128, 35}(2.0))

c1, c2 = encrypt.(kp, (a, b))
y = c1*c2
reinterpret(Fixed{Int128, 70}, convert(Integer, ScalarEncoding{eltype(params.ℛ)}(decrypt(kp, y))[]))

# output
0.0009768583554942668Q57f70
```

Alright, we're getting there. Recall from before that when we tried this in
exact arithmetic, we just got it rounded to `0.0Q3f12`. Here we're much closer
to something usable. The correctly rounded answer would be:

```
julia> plain_range[2]*2.0
0.0009768009768009768
```

so we've got about four significant digits of precision here while still using
the exact same ring as before (so the performance is the same) - not bad. Let's
just quickly try the same thing in SIMD and then we're off to the beach:

```
a = SlotEncoding(zero(plaintext_space(ckksp)))
b = SlotEncoding(zero(plaintext_space(ckksp)))

a[:] .= reinterpret.(Int128, Fixed{Int128, 35}.(OffsetArray(plain_range, 0:4095)))
b[:] .= reinterpret(Int128, Fixed{Int128, 35}(2.0Q11f4))

c1, c2 = encrypt.(kp, (a, b))
y = c1*c2
reinterpret.(Fixed{Int128, 70}, convert.(Integer, SlotEncoding(decrypt(kp, y))))

# output
4096-element OffsetArray(::Array{Q57f70,1}, 0:4095) with eltype Fixed{Int128,70} with indices 0:4095:
   403.04509382624775Q57f70
  9155.89125414592Q57f70
  6508.182113138049Q57f70
  4818.847233626708Q57f70
  3260.3910792678466Q57f70
 10573.82362550682Q57f70
      ⋮
```

Wait what happend? We did all the same things as before. Everything should work...
Noooo!!!! *drops picnic basket*. Well, recall from the [encoding](encoding.md)
section how our slot encoding works: We're doing a finite field analogue of the
FFT. And it turns out that that procedure is particularly sensitive to errors,
so while we were ok with the round off errors in the simple scalar example above,
once we start doing ffts everything falls apart. What can we do?

## Complex FFT

We get to the second major idea in CKKS: We permute the order of the fixed point
encoding and the FFT. In particular, recall the operation that is native to the
scheme: A negacyclic convolution of integers. With our fixed point encoding,
we get a negacyclic convolution of reals. What can we do with that? A lot it
turns out. In the same say way, we can turn a convolution in a finite field
negacyclic (by multiplying every element by a power of the square root of the
appropriate root of unity), we can do the same in the complex domain and obtain
the *negacyclic convolution theorem* (for complex numbers):

```math
f \cdot g = \psi^{-1}(\mathscr{F}^{-1}(\mathscr{F}(\psi(f)) *_{neg} \mathscr{F}(\psi(f))))
```

The insight here is that the complex fft is more robuts to roundoff errors than
doing the fft in the finite field. There's one small catch. We said we can do
a *real* negacyclic convolution, but in general, our convolution needs to
operate on complex coefficients. Can we work around that? As it turns out, we can.
Suppose we take image of ``\mathbb{R}^n`` under ``\psi^{-1} \circ \mathscr{F}^{-1}``.
What do we get? As it turns out get get the subspace ``\mathbb{H} \subset \mathbb{C}^n``,
where ``\mathbb{H} = { c_i \in \mathbb{C} | c_{i} = \overline{c_{n-i}} }`` (see
the [mathematical prerequisites](prerequisites) sections for a proof).
As such, as long as the value we try to encod is in ``\mathbb{H}`` our coefficients
will be real and we can use our primitive to compute them.

So let's try that:

```julia
using FFTW

ψ(a) = a .* [Complex{Float64}(exp(big(-2*k/(2*length(a))*pi*im))) for k in eachindex(a)]
ψ⁻¹(a) = a .* [Complex{Float64}(exp(big(2*k/(2*length(a))*pi*im))) for k in eachindex(a)]

plain1 = LinRange(0.0, 2.0, 2048)
plain2 = [2.0 for i=0:2047]

# Make it an element of ℍ
plain1 = [plain1; map(conj, reverse(plain1))]
plain2 = [plain2; map(conj, reverse(plain2))]

to_fixed(x) = reinterpret.(Int128, Fixed{Int128, 35}.(x))
a = zero(plaintext_space(ckksp))
b = zero(plaintext_space(ckksp))
a .= to_fixed(real.(ψ⁻¹(OffsetArray(ifft(collect(plain1)), 0:4095))))
b .= to_fixed(real.(ψ⁻¹(OffsetArray(ifft(collect(plain2)), 0:4095))))

c1, c2 = encrypt.(kp, (a, b))
y = c1*c2
dec_raw = Float64.(reinterpret.(Fixed{Int128, 70}, convert.(Integer, ToyFHE.SignedMod.(ToyFHE.NTT.coeffs_primal(decrypt(kp, y))))))
fft(collect(ψ(dec_raw)))[1:2048]

# output
2048-element Array{Complex{Float64},1}:
 -1.5783516373879536e-6 + 8.360154309303663e-7im
   0.001953265729577769 - 6.287625729356563e-7im
     0.0039068079904514 - 2.7258837087105554e-6im
   0.005860982470267118 + 1.540880824177161e-7im
                        ⋮
```

which compares quite favourably to the true answer:

```julia
julia> LinRange(0.0, 2.0, 2048) * 2.0
2048-element LinRange{Float64}:
 0.0,0.00195408,0.00390816,0.00586224,0.00781632,0.00977040,0.0117245,0.0136786,0.0156326,…,3.98046,3.98241,3.98437,3.98632,3.98828,3.99023,3.99218,3.99414,3.99609,3.99805,4.0
```

(to a precision of about 10^-6).

Of course, this is all a bit complicated, so the library provides utilities
that do all of this encoding and decoding for you:

```julia
scale = 2^35
Tscale = FixedRational{scale}
Tscalesq = FixedRational{Int128(scale)^2}

plain1 = CKKSEncoding{Tscale}(zero(plaintext_space(ckksp)))
plain1 .= OffsetArray(LinRange(0.0, 2.0, 2048), 0:2047)

plain2 = CKKSEncoding{Tscale}(zero(plaintext_space(ckksp)))
plain2 .= 2.0

c1, c2 = encrypt.(kp, (plain1, plain2))
y = c1*c2
CKKSEncoding{Tscalesq}(decrypt(kp, y))
```

## Rescaling by modswitch

TODO

## Summary

So that's it. It is often stated the the plaintext space of CKKS is ``\mathbb{C}^{n/2}``,
which is not wrong of course, but by itself also doesn't grant much insight. Having
gone through this construction, we see that arriving at the CKKS construction is
quite natural when one sits down and attempts to do arithmetic on encrypted reals.
Let us Nevertheless recall the two core ideas:
    - Using the otherwise unused noise bits in exchange for a noisy output (but
      recovering lots of intermediate precision).
    - Doing the FFT in the complex domain instead of the finite field to get
      improved noise properties.
