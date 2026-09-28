# Exact n131 rotated-arity screen: balanced m10 alone passes the frozen counting gate

The [pre-outcome protocol](PROTOCOL.md) and source/input hashes were
frozen in draft [PR #778](https://github.com/aburan28/crypto/pull/778) at
`2b565b6`, before any of these new-arity outcomes were read. The frozen-file
SHA-256 is `3e930522ef1102770b7a748477959a108f8656fa2aaaa0b22a17598b2856eedb`.
All four complete mask enumerations and independent natural-mask replays passed.
The producer and verifier row SHA-256 digests match for every mask in each arm;
the 64 SHA-selected full-group-law samples per arm also passed.

| Arm | Masks | Liftable nonzero x, L | Physical points F=1+2L | Nonzero signed ±[4] columns C | Exact F^m/q | Frozen N≥q and C≤16,384 screen |
|:--|--:|--:|--:|--:|--:|:--|
| m7,d18 | 262,144 | 130,854 | 261,709 | 130,854 | 0.1235552391 | FAIL both |
| m8,d16 | 65,536 | 32,873 | 65,747 | 32,873 | 0.5130244782 | FAIL both |
| m9,d14 | 16,384 | 8,371 | 16,743 | 8,371 | 0.1519252832 | FAIL tuple count |
| m10,d13 | 8,192 | 3,988 | 7,977 | 3,988 | 1.5329446704 | **PASS** |

Every liftable nonzero x yielded a distinct signed projected column within
its own F0 in these four arms; x=0 supplied the identity column. The ratios
are **necessary ordered-physical-tuple ceilings**, `F^m/q`, before support
collisions, valid-point equations, the PDP search, rank, or cost. For m10 the
ceiling clips at 100%; `1.533` is tuple capacity per subgroup point, not a
hit probability. The original m6,d21 arm in [#773](https://github.com/aburan28/crypto/pull/773)
had F=2,096,269 and C=1,048,134, giving exact `F^6/q=0.1246845468`.
Even four such equal-size bases have a union support upper bound of
`4F^6/q=0.4987381872`; the n19 four-base toy gain in #775 cannot overcome
that n131 counting bound at the same m6,d21 size.

The balanced m10 arm is selected only for the next *implicit PDP feasibility*
gate. Its raw affine-chain layout has 130 factor bits and eight 131-bit
intermediate x states, or 1,178 variables before exceptional branches and
solver auxiliaries. The m7/m8/m9 raw layouts have 781/914/1,043 variables.
The small m10 factor base reduces the projected column count, but the longer
chain may dominate SAT or Gröbner cost. No branch-complete m10 exporter,
measured relation yield, independent row rank, target decomposition, end-to-end
ECDLP work, or comparison with rho exists yet.

The immediate alternative worth testing is an **unequal-dimension slot
allocation** that uses the unused normal-basis coordinates without adding as
many S3 links: for m9, take 15-bit slots i=0..4 and 14-bit slots i=5..8.
Their conjugate indices are exactly 0..130, so their combined field-coordinate
rank is 131. Their physical tuple count is `F15^5 * F14^4`, not one balanced
arm's `F^9`. After inverse-Frobenius normalization the 14-bit base is a
subset of the 15-bit base, so the global signed-column set is the 15-bit
set. Exact F14, F15 and C15 require a separate preregistered census. This
field-coordinate rank is not relation-matrix rank. The m9 raw affine-chain
floor would be 1,048 variables, 130 fewer than balanced m10. This is a
parameter lead, not an admitted solver or an ECC2K-130 attack claim.

The cold measured producer wall times were 10.506/2.835/0.921/0.532 s for
m7–m10; independent full replay took 124.310/31.522/8.336/4.107 s. Each
stage passed its 600/1,200 s and 256 MiB acceptance cap. The maximum recorded
producer/independent-verifier peak RSS was 42,270,720/39,501,824 bytes. These are census
costs, not PDP or full-attack timings. [Raw evidence](evidence/receipt.json),
including all per-arm count JSON, chunk digests, stdout/stderr, code/input
hashes, and the independent results, is committed. The receipt SHA-256 is
`3e902102ccf3e97dddf0adcaa86140f4b19d79ea2f568f2e6679be69ab155670`.
Run `python3 ci_replay.py --evidence evidence` to repeat the full verifier
against the archived producer data.
