# Source-only capacity screen for paired N83 factor bases

The frozen Koblitz mixed-orbit worker has a 150,000-variable cap and a
12-million orbit-domain-clause cap. [This screen](screen.json) applies the
exact degree-83 formulas used by its chained-S3 and Frobenius-orbit-domain
builders to all **55 same-stream size pairs** in the v3 paired design. Each
size pair occurs in 18 signed-Frobenius arm/policy/seed specifications, for
990 pair specifications and four worker-shaped 5/6-summand, 1/2-residual
cases each. The 1,980 other-closure pair specifications are outside this
worker's signed-orbit domain contract.

| Summands | Residuals | Size pairs within both caps | Size pairs over variable cap | Over domain-clause cap |
| ---: | ---: | ---: | ---: | ---: |
| 5 | 1 | 54 of 55 | 1 | 0 |
| 5 | 2 | 45 of 55 | 10 | 0 |
| 6 | 1 | 48 of 55 | 7 | 0 |
| 6 | 2 | 42 of 55 | 13 | 0 |

The saved receipt gives the exact variable, domain-clause and XOR counts for
each of the 220 size-pair/shape combinations. A size pair above the variable
cap cannot start the current worker's first, all-affine identity pattern
under that frozen configuration; increasing the cap would require separate
memory and runtime validation. None of these 220 combinations exceeds the
domain-clause cap. These counts apply even where the v3 design has no
constructed archive or executable source adapter; passing a source cap is
not evidence that a base exists or a solver will find a relation.

The formulas reproduce three counters from **both** independently audited
retained 3-small/2-residual trials. K=600/K=2,048 required 100,472 variables,
415,686 domain clauses and 4,980 XOR rows; K=1,182/K=2,048 required 101,636,
465,738 and 4,980. The [generator](screen.py) pins the v3 design, both raw
trial events and audits, and the three frozen worker/chain source files by
SHA-256. The eight v3 source-available recipes, across K=600 or K=1,182,
five/six summands and one/two residuals, all pass these two source caps.

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/paired-capacity-screen-20261013/screen.py --check
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008/verification/paired-capacity-screen-20261013 -p 'test_*.py' -v
```

The screen does not build a SAT instance, measure memory peak, run a solver,
construct or upload a base, or estimate complete index-calculus runtime. It
launched no retained N83 computation, so the conservative local ledger
remains 17,400 of 17,400 authorized seconds. A best-runtime factor base still
requires natural useful-yield, rank, graph closure, linear algebra,
individual-log and matched-rho measurements.
