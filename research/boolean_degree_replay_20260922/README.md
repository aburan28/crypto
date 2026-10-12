# Frozen Boolean degree evidence replay

This is a bounded, offline correctness corpus extracted from the point
decomposition investigation. It contains 1,152 Boolean systems with four,
six, or eight variables: 768 initial S3 cases, 320 dimension-four follow-ups,
and 64 S4 cases. Rebased and membership variants are paired observations,
not independent statistical trials. Classification: **accounting**.

## Run

From the repository root, with Python 3.10 or later:

```sh
python -m unittest discover -s research/boolean_degree_replay_20260922 -p test_check.py -v
python research/boolean_degree_replay_20260922/check.py
python -m pip install -r research/boolean_degree_replay_20260922/requirements.txt
python research/boolean_degree_replay_20260922/check.py --sympy
```

All commands are read-only. They do not overwrite the frozen corpus. The CLI
reads only its sibling corpus; it has no curve or external-target interface.

## What is checked

1. The corpus digest, exact record count, and unique identifiers.
2. Every input solution, by exhaustively evaluating at most 256 assignments.
3. The reduced Boolean basis: reduced form, vanishing at every input root,
   and the standard-monomial count equaling the number of input roots.
4. The recorded degree-limited closure ranks and the first tested bound at
   which that closure contains the complete basis.
5. With `--sympy`, equality with an independently computed ordinary reduced
   Gröbner basis over GF(2), with the Boolean field equations included.

The completeness check in step 3 uses the finite Boolean ring: its ideals
are determined by their zero sets. The proposed leading ideal bounds the
quotient dimension above; evaluation at the known roots bounds it below.
Equality establishes the claimed basis. Step 5 provides a separate oracle.

The bound in step 4 has a precise restricted meaning. Start with the input
equations, reduce rows over GF(2), and multiply each derived row of degree
strictly less than D by each variable, reducing modulo x_i^2=x_i. Repeat to
closure. Test D from max(1, input degree) through the reported bound. This
is the least completing bound for this rule and range, not the least degree
of every possible Boolean proof, an F4/F5 trace, or homogeneous regularity.
The full basis is used only to check completion, never as an input row.

## Replayed result

| Check | Records | Outcome |
|---|---:|---|
| Complete input solution sets | 1,152 | 338 SAT; 814 UNSAT |
| Reduced Boolean basis and dimension | 1,152 | All match |
| Degree-limited closure and first completion | 1,152 | All match |
| Independent SymPy reduced basis | 1,152 | All match |
| Negative and boundary tests | 13 | All pass |

| Truncated-closure bound | Records |
|---:|---:|
| 2 | 559 |
| 3 | 433 |
| 4 | 96 |
| 6 | 64 |

These counts are validation coverage, not performance comparisons. In
particular, the 64 S4 cases already have input degree six. A bound of six
there does not demonstrate growth beyond the input degree.

## Provenance and limits

`manifest.json` records the hashes and counts of the three original raw and
certificate files. `corpus.jsonl` is a normalized subset: equations, claimed
bases, roots, and profiles up to the claimed bound. Hexadecimal coefficients
use the manifest's exact monomial ordering. The source ideal-membership
certificates are not needed for this exhaustive replay and are not included.

The original curve metadata, rational-point membership construction,
decomposition yield, field-product ranks, and timings are outside this
checker's claims. It validates the supplied Boolean equations, not their
derivation from curve arithmetic. It does not reproduce the entire original
experiment. No ratio to a cost boundary, end-to-end speedup, or conclusion
at extension degree 131 is inferred from this corpus.

Tamper tests reject removed roots, a fabricated root, an omitted or duplicate
basis row, false completion, altered rank or input degree, missing lower-degree
profiles, invalid sizes/types, changed corpus bytes, missing records, and
duplicate identifiers. CI repeats both replay paths and the tests.
