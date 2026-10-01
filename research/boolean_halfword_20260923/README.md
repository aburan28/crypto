# Complete Boolean solves with 16-bit equation syndromes

This standalone experiment tests the representation proposed after the negative
implicit-minor diagnostic. Its inputs are generated Boolean quadratic systems
with 12, 16, 20 and 24 variables. Here n counts Boolean variables, not a curve's
subgroup order. There are no curve inputs, imported targets, scalar recovery or
production solver calls. Full-IC, operation-calibrated and rho costs remain null.

## Exact representation and recovery

Each bit of the retained 32-bit syndrome names an equation. The new representation
projects onto the first 16 equation coordinates and evaluates eight such syndromes
per native 128-bit vector. XOR updates commute exactly with this projection.
A zero projected syndrome is only a candidate: every hit is checked on the complete
original polynomial system, and a failed check continues the scan.

The fused encoder validates the full input domain while directly constructing the
first sixteen equation coordinates. It retains coefficient cancellations and checks
omitted equations for unsupported monomials too. Unsupported input returns UNKNOWN.
The experiment does not cache coefficient-dependent output across different systems.

Both 16-point and 64-point blocks reuse the established quadratic Gray identities
and four-bit compiled transition schedule. The 64-point order is four original
16-point Gray blocks; its reflected low-group order changes with high-block parity.
The candidate preserves the model, completed point count and batch count of the
corresponding retained full-word policy. Each scalar-source and native arm must also
agree on every projected hit, original check and rejection. Scalar Rust source may
be auto-vectorized; no scalar-instruction or native-only causal claim is assumed.

Caps are charged in complete blocks. An insufficient point cap stops before the
block and returns UNKNOWN, even if the unchecked block contains a projected zero.
Exhaustion of the whole domain is required for UNSAT. Tiny domains retain the old
fallback order and block sizes.

## Frozen comparison

The same binary includes all 52 prior complete solvers, plus scalar-source and
native 16/64-point treatments, a fixed native dispatcher, and matched 32-bit/16-bit
EOR3 treatments. Runtime feature detection selects the target-feature path only
when supported; portable fallbacks preserve the two-XOR policy. Its initial policy
uses 16-point blocks through n20 and 64-point blocks above, matching the retained
dispatcher. New variants are compared with the strongest pre-existing reference;
native/scalar ratios are reported separately to distinguish representation gains
from any native-intrinsic effect.

Every arm pays for a fresh complete cold solve: validation and encoding of its
input domain, narrowing, schedule construction, scan, original-system hit checks,
temporary-workspace release and external result validation. Common fixture and
reference preparation, diagnostic serialization and returned-record destruction
are outside the arm clocks and inside worker receipts. No phase is subtracted.
"Cold" means fresh solver state, not an empty processor or operating-system cache.

The current discovery grid has 24 fixtures and 61 A/B arms, with eight random/reverse
repetitions and a preceding A/A calibration. The full grid retains all 216 predecessor
fixtures and adds 24 fresh holdouts: 240 systems and 117,120 A/B observations, plus
3,840 A/A observations. The 57-arm historical probes remain unchanged and unqualified
under the current isolation policy. Its source must match a sealed
discovery run before it can start. The previous inputs are hash-bound in
`REFERENCE_FIXTURES.json` and independently regenerated during verification.

The dramatic criterion remains a paired 95% median-ratio lower bound above 2.0
against the pointwise fastest retained complete reference in every n16/n20/n24 ×
family × regression/holdout group, with all cells complete and verified. Incremental
improvement above 1.0 is separate. A positive primary result still requires a new
unchanged-source confirmation with unused holdouts. No source tuning follows
holdout timing. Discovery never promotes a performance claim.

This is a complete generic mathematical solver comparison. There is no measured
conversion into cryptanalytic operations or a rho boundary; wall times are finite
solver diagnostics. The WDSat/full-IC suite is inapplicable to this input-only
Boolean program. The retained equivalent Boolean corpus and new holdouts provide
the matched regression here, without supplying any missing cryptanalytic costs.

## Run and verification

New timing runs require Linux CPU affinity and pressure interfaces and run through
the repository isolation controller. The local macOS host supplies correctness
checks only. Use new output directories from the repository root with Python 3.10+
and stable Rust:

```sh
python3 research/boolean_halfword_20260923/isolated_run.py --phase discovery --out research/boolean_halfword_20260923/qualified_probe_01
python3 research/boolean_halfword_20260923/isolated_run.py --phase full --discovery research/boolean_halfword_20260923/qualified_probe_01 --out research/boolean_halfword_20260923/qualified_full_01
python3 tools/isolated_bench.py busy -- python3 -m unittest discover -s research/boolean_halfword_20260923 -p 'test_*.py' -v
```

The runner snapshots source and protocol, records compiler/executable hashes,
retains raw output and fresh-worker receipts, and seals success or failure with a
manifest. The verifier reads sealed evidence; do not rerun the writing analyzer
inside an archived directory. A full run binds its discovery manifest and unchanged
timed Rust sources.

The correctness suite covers all initial blocks and Gray transitions on bounded
systems, the highest retained bit, the 16/17-equation boundary, false and multiple
hits, every small cap, all late models on a complete eight-variable fixture family,
all pairs of three-variable quadratics, unsupported omitted rows, and the retained
solver suite. Original-equation verification is never inferred from a projected hit.

See `RESOURCE_PLAN.md` for resource refusal, A/A noise, retained binaries and the Linux ARM64 CI route. The legacy `run.py` refuses current protocols. `QUALIFIED_RUNS.json` records accepted qualified evidence; null entries remain unmeasured.

Schema 3 executes A/A and A/B for one fixture in a single reserved worker and retains the combined stdout plus exact phase slices. n12 uses sixteen repetitions; other sizes use eight. The resource record covers the paired fixture. This adds real calibration/comparison work without padding or changing the 10% resource threshold; see `ISOLATION_ATTEMPTS.md`. The current discovery therefore has 14,640 comparison and 480 A/A observations; the full grid has 146,400 comparison and 4,800 A/A observations. Earlier counts describe their frozen predecessor protocols.
