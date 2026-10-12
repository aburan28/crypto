# Degree-263 all-line dual transport

The earlier [normalized-dual certificate](../../../ecc2k130_dual_transport_20260925/README.md)
checked eight saved quotient maps. This package derives the generator and
complement for any of the 264 canonical rational degree-263 kernel directions
from a saved twist-torsion basis, then checks the normalized forward and dual
maps. It is variable-time research arithmetic, not a PDP extractor.

The [protocol](PROTOCOL.md) was first pushed at `c3d1172a` and amended after a
developmental timing pilot, before any all-line result. Final source was
pushed at `f41cd690`; [LOCK.json](LOCK.json) was pushed at `0292ceba` before
the campaign. The lock covers all production, reference, workflow, protocol,
and input bytes. The frozen run uses all 264 directions of the first archived
basis plus six fixed directions of the second. Each tested direction checks
both public full-point dual identities, infinity, all 524 nonzero forward and
reverse twist-kernel inputs, and two invalid kernel-abscissa inputs. Twelve
fixed directions also receive independent bit-polynomial direct Vélu replay.

To reproduce from a checkout with the locked bytes:

```sh
python3 research/notes/ecc2k130/dual_all_lines_20260930/certificate.py \
  --out /tmp/dual_all_lines_result.json
python3 research/notes/ecc2k130/dual_all_lines_20260930/verify.py \
  --result /tmp/dual_all_lines_result.json \
  --out /tmp/dual_all_lines_verify.json
```

Both commands refuse to overwrite an existing receipt. The frozen
[producer](RESULT.json) returned `PASS_ALL_LINES` and the [independent
replay](VERIFY.json) returned `PASS`; the [decision](DECISION.md) gives the
exact counts, certificate timing, and claim boundary. Supplied torsion bases
exclude cold discovery. This establishes map correctness on the frozen
rational degree-263 lines, not natural PDP yield, relation rank, or an
ECC2K-130 index-calculus crossover.
