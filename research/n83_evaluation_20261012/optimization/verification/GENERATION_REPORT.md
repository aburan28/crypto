# Deterministic native degree-83 training and holdout freeze

2026-10-12. Eight public synthetic fixtures were generated and independently
replayed using standalone C++ field/curve arithmetic. This is a target/input
freeze and correctness result; no benchmark was executed by these tools.

The question was to prepare three training and five holdout points for the
exact retained degree-83, a=1, K600 comparison while avoiding target orbits
already used by the prior observation or the retained base. The exact model
is `icv1-f2m83-t6151469093347-9eba6cd5`, field
`z^83+z^7+z^4+z^2+1`, curve `y²+xy=x³+x²+1`, subgroup
`r=8569786107849059`. The retained G was checked explicitly.

![The scalar sidecar is for the independent validator only; measured children receive public point files.](generation.pdf)

## Generation law and checks

The generator uses SplitMix64 seed `202610128301`. For bound `r-1`, it
rejects 64-bit draws below `2^64 mod(r-1)` and maps accepted draws to
`d=draw mod(r-1)+1`. This produces uniform nonzero subgroup scalar
proposals without modulo bias. It computes `Q=[d]G`, rejects a candidate
whose signed Frobenius orbit meets the prior fixed Q, K600 base or earlier
accepted fixtures, and freezes the first three accepted points as training
followed by five holdouts.

The canonical orbit key is the minimum polynomial-coordinate x among the
83 Frobenius shifts; negation leaves x unchanged. The native checker
verified the field, subgroup and scalar eigenvalue and all retained base
points/labels before using this exclusion key. The resulting accepted law
is conditioned on the declared exclusions and avoids replacement of an
orbit; this report does not call it an unconditioned IID sample.

Recorded counts:

| Quantity | Recorded count |
|---|---:|
| Training fixtures accepted | 3 |
| Holdout fixtures accepted | 5 |
| SplitMix64 draws | 8 |
| Uniform-range rejections | 0 |
| Excluded-orbit rejections | 0 |
| On-curve, subgroup and scalar point replays | 8/8 |
| Distinct signed orbits versus prior/base/one another | 8/8 |

`check_optimization --fixtures-only` reproduced the full seed/rejection
stream, all eight point/scalar identities, signed-orbit exclusions,
individual point-file digests and exact five-row combined holdout file.
The new receipt status is `PASS_FIXTURES_ONLY`; zero invalid rows.
The independent executable uses the same standalone arithmetic helpers as
the generator and imports no research solver. It is independent of the
subsequent Rust solver runs, not an independent implementation of every
arithmetic primitive or an independent host.

## Frozen files and remaining work

The eight individual JSONL files contain one `[x,y]` array each, with decimal
coordinate strings and no scalar. The five-row holdout file uses the same
point-only encoding. Scalars are retained separately in
`validation_sidecar.json` and must not enter measured child inputs.
The manifest records source/executable and retained-input hashes; external
`sha256sum -c` confirmed every listed artifact. Generation and receipt
creation fail closed on existing paths rather than overwriting evidence.

The generator receipt SHA-256 is:

```text
c93edc32a77411443229191e2db4cc3eeb7e99976b7f79a8a753ba103bce9909
```

The independent fixture receipt SHA-256 is:

```text
fb263efde55e8177466ced20df6f09bd3b74c60b7c7505b62f91d6a650860315
```

Native compile/replay commands and file layout are in `README.md`. The
parent task will freeze the protocol before timings and later invoke the
prepared output verifier after all timed children finish. No future scalar
recovery, witness verification, timing ratio, statistical winner or full
rank/factor-log transcript is claimed by this fixture-generation receipt.
The existing historical model-alias correction remains explicit; no old
receipt or target file was changed.
