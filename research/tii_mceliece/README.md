# TII McEliece challenge keys — vendored from `mjosaarinen/tii-solved`

This directory is a verbatim import of the public data and verification code
from <https://github.com/mjosaarinen/tii-solved> (Markku-Juhani O. Saarinen),
pinned at commit `fe91133007ffcb09ac0112f35055667fde6a8faf`.

The upstream project publishes recovered secret keys for eight TII Track-2
McEliece key-recovery challenges — six new (TII-129, 213, 240, 246, 252, 254)
and two reproduction controls (TII-83, 248) whose keys Hemmert had already
published. See `UPSTREAM_README.md` for the original text and the timing
tables, and `../../RESEARCH_TII_MCELIECE.md` for our notes on the attack ideas
and how they relate to this repository.

Upstream carries no license file, so treat everything here as the authors'
work, kept for study and reproduction with attribution. The public keys are
themselves copies of the TII challenge repository (`ElenaKirshanova/tii_decoding_challenge`)
and the Hemmert archive; `SOURCES.md` records the pinned upstream commits and
the SHA-256 of every public-key file.

## Contents

```
tii_public_keys/            original public parity-check matrices (verbatim)
tii_secret_keys/
  secret_key_tii_<N>.json   canonical recovered key: GF(2^m) integers + params + provenance
  sk_McEliece_<N>.txt       TII Track-2 two-line text view
verify_recovered_key.sage   standalone SageMath verifier
export_secret_keys.sage     regenerates the .txt / .pckl views from canonical JSON
KEY_FORMATS.md              exact field-element mapping and interoperability spec
SOURCES.md                  upstream URLs, pinned commits, public-key digests
timings/                    upstream recovery and verification run records
```

### Deliberate omission: no pickles

Upstream also ships `secret_key_tii_<N>.pckl` (Sage `[support, g]` pickles for
Hemmert's verifier). Those are **not** vendored here: unpickling executes
arbitrary code, and they are a redundant view of the canonical JSON. If a
consumer needs them, regenerate locally with

```bash
sage export_secret_keys.sage all pickle
```

## Verifying

```bash
cd research/tii_mceliece
sage verify_recovered_key.sage 252     # one challenge
sage verify_recovered_key.sage all     # every key present
```

For each key the verifier rebuilds `GF(2^m)` from the stored modulus, re-checks
the public key's SHA-256 against the digest the recovery was bound to,
reconstructs the Goppa parity check

```
H_rec[j*r + k, l] = (y_l · x_l^k)^(2^j),   y_l = 1/g(x_l),
0 ≤ k < r = deg g,   0 ≤ j < m,   0 ≤ l < n
```

and accepts iff `RowSpace(H_rec) == RowSpace(H)` over `GF(2^m)`. Comparing row
spaces rather than a fixed echelon form makes the check independent of how `H`
was stored. `(m, r, n)` are checked against a verifier-owned table of published
challenge parameters, not re-derived from the candidate; all 17 checks must
pass. SageMath is the only dependency.
