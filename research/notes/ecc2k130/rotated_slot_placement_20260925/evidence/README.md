# Six-position raw evidence and replay

The preregistered source/input/workflow freeze is `FROZEN.json` SHA-256
`ff4bd870a0c21e823fc7efa117630876118e0685da92079232c48cb2f07c5fcf`.
Draft [PR #793](https://github.com/aburan28/crypto/pull/793) and its exact-head
hash-only `parse` and `replay` checks passed before the first six-sum run.
The earlier draft head `b61d8f1` had no triggered workflow; it was corrected
pre-outcome at `5b2b577`, and its initial freeze SHA is retained inside the
final `FROZEN.json`. No outcome was computed under the initial head.

The single frozen run began `2026-09-25T16:01:21.297832+00:00` and ended
`2026-09-25T16:02:34.244616+00:00`. All six producer and six independent
verifier children, plus the exact pair analysis, exited zero. Every
producer/verifier stayed below its 300/600-second and 512-MiB caps.
The archive `raw.tar.gz` is 745,073 bytes with SHA-256
`495e522064f440f0317640d28d812c5ef012f3dd12580a07e3b941c806b92d68`.
`raw.tar.gz.sha256` checks the archive; `receipt.json` records every command,
exit status, external cap, start/end time, stdout/stderr hash, source freeze,
and each raw member hash. The archive contains six `counts.u32le` files of
exactly 130,873 uint32 entries, six result files with all 64 shared-target
witness rows, six independent verification receipts, exact pair analysis,
and all child logs. `analysis.json` is a duplicate of the archived file for
quick review. A local temporary run directory is not required to replay it.

From the repository root, use:

```sh
python3 research/notes/ecc2k130/rotated_slot_placement_20260925/ci_replay.py
python3 research/notes/ecc2k130/rotated_slot_placement_20260925/ci_replay.py --evidence research/notes/ecc2k130/rotated_slot_placement_20260925/evidence
```

The first command reads only frozen source/input hashes; the second extracts
the committed archive, independently rebuilds every factor and q-point table
under bit-serial/Fermat arithmetic, recomputes all six complete multiplicity
arrays, verifies the 64 same-target witnesses and modular ranks per arm, and
reapplies all 15 pair comparisons. It passed locally after the run and is run again by PR CI when this
evidence is committed. `accounting.py` and `accounting.json` are explicitly
*post-outcome cost reconstruction*: they rederive exact sparse-convolution
integer update counts and prefix supports from the fixed factors and compare
the final support to the archived result. They do not change the frozen
decision, inputs, or census program.
