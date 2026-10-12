# P-192 and P-224 large prime-degree isogeny search

Frozen before execution, 2026-10-08. Requested scope: use the new native
isogeny algorithms on P-192 and P-224, including prime degrees beyond 1009.

## Question and inputs

Can the new Hecke/Newton modular-polynomial and BMSS fastElkies' construction
produce explicit rational prime-degree maps beyond the older walker's practical
degree limit of approximately 61? Extend the structural candidate search to
every prime in [67, 4093]. This is an isogeny construction experiment.

Use the standard models and published prime group orders in the existing
curve registry. Source revision: c70c32d486a3ac7531fe27f7193d9f09caa58344;
isogeny_algos revision: 7819ee8c0. Native source-order checks use seed 1;
construction uses seed 1; replay uses a different field/kernel implementation
and new public point samples. Preserve the protocol SHA-256 in the run manifest.

## Frozen search and stopping rule

1. For each prime degree in [67, 4093], compute the Frobenius discriminant
   modulo the degree. Record split, inert, or repeated-eigenvalue status,
   eigenvalues and their multiplicative orders. An inert degree is excluded
   from rational prime-degree construction. A repeated eigenvalue remains
   unresolved until its kernel is constructed; do not infer a kernel count.
2. Attempt every split degree in [67, 257]. Also attempt, separately for
   each source curve, the first split prime at or above each of 509, 1009,
   2003, and 4000, within the frozen upper bound. Do not replace failures
   with easier degrees. All other degrees are structural screening only.
3. Give each construction process 180 seconds and an 8 GiB resident-memory
   limit, including modular-polynomial
   setup, root finding, kernel construction, and the CLI's checks. Kill the
   child at the deadline, preserve its partial stdout/stderr, and record a
   timeout. Execute one construction at a time. No speed comparison is made;
   process elapsed times are operational receipts on a shared macOS host,
   not isolated benchmark measurements.

The memory limit was recorded before the second launch, after the preserved
first launch stopped at the P-192 dispatcher preflight. The upper-degree
Hecke/Newton tables otherwise grow beyond this host's available memory.
macOS rejected `ulimit -v`, so the native parent samples child RSS every
25 milliseconds with `proc_pidinfo` and kills it if it exceeds the limit.
An unavailable monitor fails closed. Sampled peak RSS and stop reasons are
retained. Brief growth between samples is possible; this is a supervised cap.
4. A construction is accepted only if the CLI exits successfully, all its
   checks pass, and the number of maps equals the two distinct Frobenius
   eigenlines for a split degree. Retain any discrepancy as a failure.
5. Independently replay every completed kernel using the existing walker's
   squarefreeness, division-polynomial torsion, subgroup-closure, and Vélu
   codomain checks. Replay the rational map on fresh known public points,
   including subgroup-order checks. Report only replayed maps as certified.

## Evidence and reporting

Keep exact commands, public parameters, source revisions, protocol/source
hashes, structural screening, raw CLI records, exit codes, failed attempts,
timeouts, independent replay receipts, and artifact SHA-256 digests. Record
kernel degree, map coefficient counts, and target model identities. Register
new target ICV1 identities before citing them. Keep full-map construction and
structural candidate counts distinct. Construction beyond 1009 is an explicit
acceptance obligation; timeouts leave that obligation partial.

Produce a dated report, editable vector diagram/coverage plot, and visually
checked PDF. Check the existing scoreboard, progress timeline, curve diagrams,
leaderboard and browser for affected views, documenting updates or why their
contracts are unaffected. No ECDLP-cost comparison or ledger promotion is
part of this experiment.

References: [BMSS](https://arxiv.org/abs/cs/0609020),
[the new CLI](../../isogeny_algos/docs/USAGE.md),
[the existing walker](../../docs/isogeny-walk/README.md).
