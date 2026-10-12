# Frozen replay protocol: EXP-TOY-BINARY-FORMULA-001

Status: stage diagnostic; accounting classification; no performance claim.
The original exploratory run predates this protocol. This protocol freezes
the subsequent repository replay; it is not a retrospective preregistration.

Hypothesis: for fixed binary coordinate encoding and subspace V, varying
nonzero a6 changes only the constant terms of the Boolean ANF of f3.
Reference: exact identity f_c = f_0 + c. Conjugacy control:
f_(c^2)(x^2,y^2,z^2) = f_c(x,y,z)^2, transporting V by squaring.

Frozen inputs: fields of degrees 3, 4, 5; polynomial moduli 0b1011, 0b10011,
0b100101; all coordinate subspaces of dimensions 2 and 3; every nonzero
coefficient; every ordered triple in each subspace. No seeds or sampling.
The complete inputs are encoded in audit.py and expanded in results.json.

Success: zero assertion failures, exactly 798 systems, 34 subspaces,
219968 pointwise conjugacy checks and 219968 symmetry checks, and byte-identical
replay of results.json. Stop on any mismatch; preserve the output for diagnosis.
Do not extend the field sizes as part of this experiment.

Cost accounting: assertion and assignment counts are validation workload,
not calibrated cryptanalytic operation counts. Setup includes complete tiny
field multiplication tables; encoding uses exhaustive truth tables and ANF
transforms; verification includes conjugacy and permutation checks. No solver,
relation collection, or scalar recovery is performed. End-to-end operations,
S, ratios to rho, speedup, and degree of regularity are all unset (null).
The WDSat performance regression is inapplicable because no solver or
production algorithm is changed. No runtime comparison is claimed.

Replay from this directory:

```sh
python3 audit.py --output /tmp/binary-formula-replay.json
cmp results.json /tmp/binary-formula-replay.json
```

Decision: retain the coefficient/conjugacy control as an algebraic baseline.
Do not infer an isogeny advantage or verified rational-point relation yield
from these polynomial-zero counts. See README.md for results and limits.
