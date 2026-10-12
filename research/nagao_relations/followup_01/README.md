# Follow-up 01: denser bases and larger SAT instances

Executed the frozen [contract](contract.json) using the existing encodings and
independent curve arithmetic. Sources, input targets, exact solutions, failures,
and counts are in [raw.jsonl](raw.jsonl); aggregates are in [summary.json](summary.json).
No mathematical mismatch or infrastructure exception was recorded.

## SAT completion under the same five-second limit

Eight uniformly sampled affine targets per field, seed 2026091301, three
summands, weight at most two. Both encodings use identical nonzero/distinct-x
restrictions and full projected solution enumeration. Completion requires
solver UNSAT after blocking all solutions, followed by equality with the
exhaustive group-law oracle.

| Field degree | Variant | Complete / attempted | Verified x-tuples in complete cases | Verified x-tuples in incomplete cases | Full-DLP S / rho ratio |
|---:|---|---:|---:|---:|---|
| 9 | Chained S3 | 0 / 8 | — | 128 | Unmeasured |
| 9 | RR norm | 8 / 8 | 363 | — | Unmeasured |
| 11 | Chained S3 | 0 / 8 | — | 4 | Unmeasured |
| 11 | RR norm | 0 / 8 | — | 9 | Unmeasured |

All eight completed norm solution sets match exactly. Incomplete instances are
watchdog outcomes, never UNSAT or mathematical rejection. The practical
completion advantage at n=9 is observed on this panel. It supplies no SAT
operation-count speedup, asymptotic conclusion, or prediction for ECC2K-130.
Circuit size and diagnostic times remain in the raw records. Class:
**engineering observation**, with no measured full-attack improvement.

## Function-first support search

Over GF(32), V is the span of the first d normal-basis vectors. The d=4 panel
uses every affine target; d=5 uses eight uniformly sampled targets with the
same fixed seed. Each target exhausts all 32×31 coefficient pairs a,b≠0.
Every signed point triple is independently verified and the entire solution
set agrees with direct point enumeration.

| dim(V) | Signed base points | Target cases | Supported cases | Verified signed triples |
|---:|---:|---:|---:|---:|
| 4 | 20 | 43 | 43 | 864 |
| 5 | 42 | 8 | 8 | 1,780 |

The following is a deliberately explicit **unit-cost field-API model**, adding
one unit per field addition, multiplication or squaring. Inversions and curve
calls are expanded rather than added again. It includes each implementation's
setup, unsuccessful candidates, extraction and verification. It excludes
Python overhead and is not a measured conversion to curve additions or machine
instructions. Separate primitive categories and phases are retained in the raw
data. Ratios compare identical inputs within each dimension; do not compare
the two dimensions as though the factor base or target count were fixed.

| dim(V) | Variant | Field-API units | Ratio to point-enumeration reference | Correct |
|---:|---|---:|---:|---|
| 4 | Point enumeration | 1,379,967 | 1.000 | Yes |
| 4 | Function-coefficient enumeration | 4,399,218 | 3.188 | Yes |
| 5 | Point enumeration | 2,835,805 | 1.000 | Yes |
| 5 | Function-coefficient enumeration | 1,531,283 | 0.540 | Yes |

The dense d=5 case is the entire field as an abscissa domain. Its favorable
ratio is against a brute-force reference on a tiny, almost-full-group factor
base, not the best Semaev solver or a scalable index-calculus attack. Counting
different primitives equally is a model assumption, not hardware calibration.
There is no full-DLP or lower-bound ratio for this microbenchmark.

Membership tests consume 67.0% of RR field-API units in d=4 and 45.5% in d=5;
root extraction rises from 3.0% to 32.5%. These identify where to experiment
next without hiding the increasing cost of extracting dense solution sets.

## Next experiments justified by these observations

1. Compile H-divides-L_V into the shared IR, replacing exhaustive coefficient
   search, and compare against RR norm on exactly the same subspace base.
2. At n=11, distinguish first verified relation from complete enumeration and
   test coefficient conditioning on a new frozen target panel. Charge all
   branches and unsuccessful targets; retain the current timed-out baselines.
3. Before any full-attack claim, price rank growth, relation collection,
   linear algebra and descent against rho on the same subgroup.

Reproduce in a new evidence directory: copy `run.py` and `contract.json` to a
sibling of `followup_01`, then run the copied script using Python with
`pycryptosat`. Existing evidence paths are refused rather than overwritten.
The [publication mapping](publication.json) identifies the exact source
snapshot corresponding to the local commit recorded by the run.
