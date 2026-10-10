# Cache and batch inversion beyond eleven bits

This campaign compares the original hybrid, a cached hybrid and the matched
elementary-symmetric S4 circuit on fresh targets at (n,d)=(11,5),(23,6),(29,6).
All 144 trials use five-second all-phase budgets and separately measure first
relation and complete projected enumeration. Deterministically shuffled
variant order reduces fixed-order bias. Each instance pays its own setup.

The cached solver computes L_V once, moves the a=0 residual norm outside the
inner root loop, and uses one batch inversion for all eligible z and r+z.
For each nonzero b it computes one inverse, then uses

```
alpha_z = z/(r+z)
beta_z  = (r+z)/z²
t       = b*alpha_z
AS input = K*beta_z/b².
```

This is the same input K/((r+z)*t²) used by the original solver. Prefix
products and a backward pass compute every inverse with one inversion and
linear multiplication work. This changes neither candidate order nor the
full cubic support predicate. Setup and elimination work remain charged.

Validation checks the full candidate sequence and complete projected set
against the original on all 43 affine five-bit targets. The new independent
group-law oracle indexes sums of pairs of signed factor points, then matches
R-P3 against that table. Its answers equal the old signed-triple oracle on
all 43 targets. Larger-field correctness uses this pair oracle, avoiding
exhaustive enumeration of all field elements or all target points.

The cached solver resolves all four supported-target first-relation slots
at both 23 and 29 bits. It resolves seven of eight full enumerations at
23 bits and none at 29 under five seconds. The symmetric S4 control resolves
none of the larger-field slots under that budget. All 144 trials have zero
validation errors. This clears some larger cases but fails the campaign's
all-slots-complete scaling gate.

The four uniform targets at each larger size have no restricted decomposition.
Their completion is evidence of a correct negative decision, not relation
yield. Supported targets are sampled from random distinct signed factor
triples and are biased toward targets with more decompositions. All three
variants see the same targets. The first-d-coordinate subspaces have only
2^d x-values; at d=6 their density falls with field degree. No ECDLP exponent,
natural-yield gain, or full attack speedup is inferred.

Field API operation vectors are comparable between the two hybrids only for
matched completed work. SAT operations have no calibrated conversion. Actual
wall time and within-budget flags are retained because watchdog deadlines
are soft. Raw timeouts remain unknown even when a partial run found all
oracle tuples. See ../scaling_23_29.md for the combined follow-up and the
image-space successor, which was frozen after this profiling evidence.
