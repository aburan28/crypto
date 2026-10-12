# Native replay of the completed degree-83 panel

2026-10-12. Correct exact model:
`icv1-f2m83-t6151469093347-9eba6cd5`. This descriptive report records
independent arithmetic replay of completed outputs, not new measurements of
algorithm performance.

![The native verifier checked the common instance and every recorded answer after all timed children completed. Counts include warmups and A/A controls and do not estimate success on independent targets.](fresh_verification.pdf)

## Observed verification result

The parent evaluation task ran `check_fresh` after all timed children
finished. The native receipt reports `PASS`, **22 independently verified
answers**, zero invalid rows and 22 retained process dispositions with status
`COMPLETED_VERIFIED`:

| Producer arm | Retained rows | Scalar replay to exact public Q | Four-point witnesses |
|---|---:|---:|---:|
| Strong signed-Frobenius rho | 16 | 16/16 | Not applicable |
| IC | 6 | 6/6 | 6/6 re-add to Q |

The counts include one warmup for each arm, ten strong A/A runs, and five
formal matched repeats for each arm. All rows concern the same public point
and are preserved individually. These are archive-replay counts, not an
independent-target completion estimate or a speedup claim. Failed outcomes,
caps and missing answers would remain explicit dispositions; none occurred
in the checked 22 rows. The parent reports checker execution of approximately
0.62 seconds, outside all producer timing intervals. The checker was not
rerun during this source/receipt review.

## Exact arithmetic scope

The retained field is defined by `z^83+z^7+z^4+z^2+1`; the curve is
`y²+xy=x³+x²+1`. The native checker revalidates field irreducibility,
prime subgroup `r=8569786107849059`, curve order and cofactor,
generator/target membership, Frobenius eigenvalue and the complete retained
factor base. The validated base has 600 representatives and 99,600
distinct signed-orbit points with correct scalar coefficient labels.

Every fresh recovered scalar passes a standalone polynomial-field affine
replay `[d]G=Q` against these exact coordinates:

```
G_x = 6614244651398184860292015
G_y = 8868260410029249391539529
Q_x = 2655898345537166180532960
Q_y = 3881526307991673265292615
```

For every IC row the checker also validates four point indices against the
independently checked base, compares x codes and intermediate abscissae,
and re-adds the actual signed points to Q. The validation-only fixture
scalar is compared only after a recovered answer passes scalar replay.
This verifier did not execute either research solver and was not a solver
input. It uses an independent C++ implementation of the standard field and
curve arithmetic, importing no research solver. The fresh check ran on the
same host after the fresh panel, not on an independent performance host.

## Identity correction and evidence boundaries

The correct full explicit-field model identity is:

```
ICV1:f2m-83-a5fd0d39:6151469093347:
9671406556910881928556062:0x1:-7:unk:r:9eba6cd5d359
```

The line break above is for rendering only. The older manually inserted
alias `icv1-f2m83-t6151469093347-cdcc5432` names a different modulus
`[0,1,2,45]`. The original archive receipt and archive report remain unchanged
and must be read with `CURVE_ID_ERRATUM.md`. Their explicit arithmetic inputs
were correct; the identity-label correction does not change their point checks.
The fresh checker checks the explicit modulus, subgroup, G and Q and depends
on neither alias.

The 600 coefficient coordinates are formal factor-log bookkeeping, not
600 independent group generators. The rational subgroup is cyclic, generated
by G. Each `R_j=[L_j]G`; the guided inhomogeneous system recovers 600 scalar
values L_j relative to G's known logarithm one. `INVARIANTS.md` explains the
coefficient map, orbit phase transport, root-index witness and strong rho's
two-coefficient collision equation with implementation references.

Correct point answers do not independently certify the complete rank
construction. The retained target and summary rows do not expose all 600
rank equations, their decomposition witnesses, the final matrix and solved
factor-log vector. Independent replay of that preparation transcript remains
unverified. The native replay also does not reproduce operation counts,
validate host-isolation promotion gates or certify timings on another host.

Both fresh producers retain positive timer residuals, separate from
canonicalization's runs of zero bits. In formal pair 1, rho's full interval
exceeds its phase sum by 401 ns; IC's full target interval exceeds its five
phase sum by 0.029956 ms. Full online intervals include those residuals. No
exact phase partition is asserted by this verification report.

## Durable evidence and reproduction

- Fresh native receipt: `fresh_replay_receipt.json`.
- Independent native sources: `check_fresh.cpp` and included `check.cpp`.
- Original archive verification: `receipt.json`, `REPORT.md`, `REPORT.pdf`.
- Identity correction: `curve_id_erratum.json`, `CURVE_ID_ERRATUM.md`.
- Existing algorithm explanation: `INVARIANTS.md`.
- Input/source/executable provenance: the three preserved SHA-256 files.

Fresh receipt SHA-256:

```
980a3d7d52db96a559dd36ad12e317c9b12be9f1c161a2d0d86803aafcc72441
```

Compile/replay instructions are in `README.md`. The editable verification
diagram is `fresh_verification.svg`, with its vector PDF rendering.
The parent evaluation task owns the main study report, timed records and
statistical analysis. This copy adds source/receipt verification documentation
without altering those artifacts or claiming repository publication.
