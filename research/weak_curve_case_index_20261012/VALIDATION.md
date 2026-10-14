# Validation receipt — 12 October 2026

Worktree parent commit at import start:
`c9a6767a81e4c060339ceb2b4a569e1fecc53029`. The only new repository
paths in this follow-on are under `research/weak_curve_case_index_20261012/`.
No source curve/model entry in the main database was relabeled.

| Check | Result |
| --- | --- |
| `sh build.sh` | Exit 0, `case-index-ok`, `PRAGMA integrity_check=ok`; 36 consecutive numbered cases and 5 minimum-genus pattern rows |
| Rebuild the SQLite file twice and compare SHA-256 | Identical digest `b74448a12c09ad6f1f733593d12e0795da83f494b59cad5ec91eca5402ae6f88` |
| `sqlite3 -header -csv published_cases.sqlite 'SELECT d,COUNT(*) cases,SUM(e=0) minimal,SUM(e>0) higher,MAX(cover_genus) max_genus FROM published_case GROUP BY d ORDER BY d;'` | Rows `(2,4,1,3,5)`, `(3,8,2,6,17)`, `(5,11,1,10,129)`, `(7,13,1,12,769)` |
| `rustc --edition 2021 --test ramification_orbits.rs -o /private/tmp/weak_ramification_orbits_test && /private/tmp/weak_ramification_orbits_test` | Exit 0; 3/3 native tests passed: factor divisibility, cubic mask, factor-orbit counts |
| `rustc --edition 2021 ramification_orbits.rs -o /private/tmp/weak_ramification_orbits && /private/tmp/weak_ramification_orbits > /private/tmp/weak_ramification_orbits.csv && cmp ramification_orbits.csv /private/tmp/weak_ramification_orbits.csv` | Exit 0; committed 27 orbit rows match native replay |
| `lake build PublishedCases` in `lean/` | Exit 0; Lean 4.34.1, mathlib `d13f23b723b8a846827a245b89c10fc7d3f11612` |
| `lake env lean AxiomAudit.lean` | Exit 0; both declarations use only `propext`, `Classical.choice`, and `Quot.sound` |
| `lake env leanchecker --fresh --verbose PublishedCases` | Exit 0; `replaying PublishedCases with --fresh` |
| `xmllint --noout case_profile.svg construction_pipeline.svg` | Exit 0; both vector diagrams parse |
| `tectonic --only-cached --outdir /private/tmp REPORT.tex` | Exit 0 outside the filesystem sandbox using existing cached resources; final `REPORT.pdf` is one US-letter page, searchable, and visually inspected after the last source edit |
| `git diff --check` and Lean source search for `sorry`, `admit`, `sorryAx`, `axiom` | Exit 0; no matches in the two Lean modules |

The additional primary-source screen used [Costa et al., Theorem 2.8](https://arxiv.org/pdf/2002.02067).
The exact read-only query on the existing model/class database was:

```sql
SELECT mechanism_scope, trace % 2 AS odd, COUNT(*) AS rows,
       SUM(CASE WHEN combined_support=1 OR representative_count>0
                THEN 1 ELSE 0 END) AS positive
FROM class_observation
WHERE trace_status='ordinary'
GROUP BY mechanism_scope, odd ORDER BY mechanism_scope, odd;
```

It returned `two_branch_union,0,87486,85412` and
`cubic_norm_one,0,246512,122058`. The second row is a probabilistic
diagnostic population; the first is a complete exact two-branch population.
No odd trace occurs in either, so the source's odd-trace exclusion does not
split the observed high-depth zero rows.

The repository-wide required check `cargo test --release --lib --offline`
ran against the unchanged root library source. It exited 101 after 192.23 s:
3,932 passed, 11 failed, 112 ignored. Nine failures were sandbox socket
`Operation not permitted` errors. The N83 relation-query scalar assertion
failed again in a focused replay; the Z3 `blocked model returned again`
failure passed in a focused replay. The complete output is frozen in
[`root_cargo_test.log.gz`](root_cargo_test.log.gz). These failures are outside
the added case-index code and remain reported rather than counted as passing
root verification.
