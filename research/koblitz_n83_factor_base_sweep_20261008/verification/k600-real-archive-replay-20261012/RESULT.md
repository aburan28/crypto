# K=600 exact-archive constructed replay

The K=600 natural-search auditor passed one constructed known-answer partial
relation against the exact retained compressed K=600 and K=2,048 archives.
The [machine receipt](replay.json) binds both object SHA-256 digests, the
fixture and replay-source digests, the selected point indices and orbit
columns, the constructed target, and the results of three negative controls.
The K=600 object is already stored at
`s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/a0/objects/b4dda45e6c3edfca673db6f3f82c16cce8c860f5e62c092155258df399a0363e.jsonl.gz`;
the envelope remains the pinned K=2,048 S3 object named in its manifest.

The existing public known-answer fixture contributes three small points at
global indices 0, 166 and 332, all in K=600, and two residual points from
envelope columns 1,182 and 1,183, both in the K=600 suffix. The fixture's
five-point sum is `S`. The replay constructs a separate subgroup point
`Q = b^{-1}(S - aG)` using the worker's frozen first-trial coefficients `a`
and `b`, so the auditor checks the exact relation `aG + bQ = S`. This `Q` is
a constructed control input; it is not the natural public-corpus target. The
relation is a graph pivot with two residuals, not a completed factor-base row.

The auditor accepted the correct archived relation and rejected (1) a
different, still archived small point because its group equation changed,
(2) a small-domain point substituted for a residual because it was absent
from the exact envelope suffix, and (3) a 415,685-clause event because the
K=600 three-small/two-residual domain must contain 415,686 clauses. These
checks exercised actual point rows, orbit labels, subgroup arithmetic and
the new supervisor's event replay, beyond its synthetic unit controls.
`replay.py` provides the source and can regenerate `replay.json` from local
copies of the same S3 objects:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/k600-real-archive-replay-20261012/replay.py \
  --small /path/to/k600.jsonl.gz \
  --envelope /path/to/k2048.jsonl.gz \
  --output /path/to/new-replay.json
```

The event and its counters were constructed for an auditor test. No SAT
model, natural-target relation search, graph cycle, rank growth or complete
index-calculus runtime was measured. The conservative retained-N83 local
charge remains 15,000 seconds; no new S3 write or billed resource was used.

Validation: the replay passed twice against the exact local copies, and
`cmp` found the two JSON receipts byte-for-byte identical. The full release
library suite passed 4,000 tests with 113 ignored; the N83 study Python suite
passed 51 tests and the boundary suite passed 43. Python compilation, JSON
parsing and staged whitespace checks passed.
