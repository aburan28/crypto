# Recount of the instance and prime counts (sheriff case SHF-20261016-3caa6f)

Raw outputs behind three numbers in `../structural_completeness.tex`, re-run on
2026-10-16 at commit `f25d5e520` (branch `claude/sheriff-major-findings`).
~~the four scripts are unchanged since `8b6b8c8ba`~~ Three of the four cited
scripts (`exp_zprime_l_residue.py`, `exp_zprime2_generalized_lemma.py`,
`thread12_cm73_sweep.gp`) are unchanged since `8b6b8c8ba`. The fourth,
`exp_zdoubleprime_lemma_extended.py`, was edited in `11aa00591` (the commit
that added this directory) to stop printing the stale `16562` total, after its
2026-10-16 output had been captured; that output therefore came from the
`8b6b8c8ba` version. It was regenerated on 2026-10-17 with the current script
(Python 3.13.16, `python3 -I -B exp_zdoubleprime_lemma_extended.py`, exit 0,
3.4 s, at `b59187e55`); the only change is line 63, which now reads
`(Exp Z''): 1,802 instances.` instead of `Grand total: 16562 + 1,802 = 18,364
instances.` (Correction 2026-10-17, sheriff case SHF-20261017-e4d7f6, finding
P-N5.)

| File | Command (from `secp256k1_cm_audit/` unless noted) | Exit | Wall time |
|---|---|---|---|
| `exp_zprime_l_residue.txt` | `python3 -I -B exp_zprime_l_residue.py` | 0 | 62 s |
| `exp_zprime2_generalized_lemma.txt` | `python3 -I -B exp_zprime2_generalized_lemma.py` | 0 | 11 s |
| `exp_zdoubleprime_lemma_extended.txt` | `python3 -I -B exp_zdoubleprime_lemma_extended.py` | 0 | 4 s |
| `thread12_cm73_sweep.txt` | `gp --stacksize 32000000 -q thread12_cm73_sweep.gp` | 0 | 255 s |
| `count.txt` | `python3 -I -B papers/structural-completeness/recount/count.py` (from the repository root) | 0 | 5 s |

Python 3.13.16; PARI/GP 2.15.4 (x86-64, GMP 6.3.0). Each `.txt` output is the
script's stdout and stderr, followed by an `exit=` line.

`count.py` imports the cited scripts' own `count_points_j0` and `is_prime`,
repeats their enumeration, and splits each instance `(l, p, b)` with
`l | #E` for `y^2 = x^3 + b` by curve type (ordinary iff `p = 1 mod 3`) and
by residue class (`p mod l` inside the lemma's hypothesis or equal to `l-1`).
Its totals agree with the per-class tables the cited scripts print.

## What changed in the paper

- **Enhanced Lemma check** (remark after Corollary "B5 universality"). The
  paper said 16 562 ordinary instances for `l in {3,5,7}` and 18 364 in
  total. 16 562 was a hard-coded literal in
  `exp_zdoubleprime_lemma_extended.py`; the cited
  `exp_zprime2_generalized_lemma.py` enumerates 24 476 instances (8 734
  ordinary, 15 742 supersingular). All 8 156 instances in the lemma's
  hypothesis classes are ordinary and none violates the lemma. With
  `l in {11,13}` (1 802 instances, ordinary only) the total is 26 278, of
  which 10 536 are ordinary and 9 864 ordinary instances satisfy the
  hypothesis: 0 violations.
- **CM-73 sweep** (Remark "Toy example: split Jacobian"). The paper said
  621 primes `p = 1 mod 6` with `p <= 5000`. `thread12_cm73_sweep.gp` loops
  `forprime(p=7,5000)` and returns early unless `p % 6 == 1`, so it processed
  330 primes (7 to 4999), out of 669 primes `p <= 5000`. Its output (4 CM-73
  primes, 11 norm-form primes, the seven `a2 - 2p` values) matches the
  paper.
- **Abstract** now states the theorem's hypotheses (cofactor
  `h in {1,2,4,8}`, bound `sqrt(r / |Aut(E)|)`, the attack family only), the
  standing assumptions of Section 8, and the bounded ranges of every
  empirical check.

~~The `exp_zdoubleprime_lemma_extended.py` output still prints the stale
literal (`16562 + 1,802 = 18,364`) in its coverage summary; the paper no
longer uses it.~~ That was true of the output captured before the script was
edited; the regenerated output (2026-10-17, above) no longer prints the
literal. All counts are unaffected.

## Report

`python3 make_report.py && pdflatex report.tex` renders `recount.svg` and
`report.pdf` (counts table and stacked bar chart) from `count.txt`.
