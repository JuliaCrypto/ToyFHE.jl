export ModulusRaised

"""
    ModulusRaised{P<:SHEShemeParams}

Modifies the underlying scheme to treat the last prime in the CRT basis as
a "special prime" to be used for modulus raising. This technique is used in
Microsoft SEAL (as of version 3.3) and described in https://eprint.iacr.org/2019/524.pdf.
Using the special prime reduces the noise introduced by keyswitching, without
requiring an (expensive) more fine grained relinerization radix.
"""
struct ModulusRaised{P<:SHEShemeParams} <: PassthroughParams{P}
    params::P
end

scheme_name(::Type{ModulusRaised{P}}) where P = string(scheme_name(P), " (with special prime)")

# The special prime is reserved for keys, so the ciphertext ring is the subring
# with the last crt component droppped.
ℛ_cipher(params::ModulusRaised) = drop_last(ℛ_cipher(parent_params(params)))
ℛ_plain(params::ModulusRaised{<:CKKSParams}) = drop_last(ℛ_plain(parent_params(params)))

function encrypt(rng::AbstractRNG, key::PubKey{P}, ::Zero) where {P<:ModulusRaised}
    c = encrypt(rng, PubKey(parent_params(key.params), key.key), Zero())
    CipherText{Zero}(key.params, map(modswitch_drop, c.cs))
end

function make_eval_key(rng::AbstractRNG, (old, new)::Pair{<:Any, P}) where {P<:PrivKey{<:ModulusRaised}}
    ps = modulus(moduli(eltype(eltype(new.secret))).parameters[end])
    parent_new = PrivKey(parent_params(new.params), new.secret)
    KeySwitchKey(new.params, make_eval_key(rng, ps*old=>parent_new).key)
end


function keyswitch_expand(ek::KeySwitchKey{<:ModulusRaised}, c)
    ℛkey = typeof(ek.key[1].mask)
    key_moduli = moduli(coefftype(ℛkey)).parameters
    𝔽ps = key_moduli[end]
    ℛexpanded = crtselect(typeof(ek.key[1].mask), [1:length(moduli(c).parameters); length(key_moduli)])
    ℛexpanded(c .* CRTExpand{𝔽ps}(), nothing)
end
keyswitch_contract(ek::KeySwitchKey{<:ModulusRaised}, c) = modswitch(c)
function downswitch_keyelement(params::ModulusRaised, key::KeyComponent, elt)
    which = [1:length(moduli(elt).parameters)-1; length(moduli(key.mask).parameters)]
    KeyComponent(
        crtselect(key.mask, which),
        crtselect(key.masked, which)
    )
end
