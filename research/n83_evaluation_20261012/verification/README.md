# Native independent verification artifacts

The fresh panel replay passed for **16 strong-rho producer rows and six IC
target rows**, including warmups and A/A controls: all 22 answers replay to
the same exact public point, and all six IC four-point witnesses re-add to it.
All 22 process dispositions remain `COMPLETED_VERIFIED`. This verification
does not reproduce the timings or the 600-relation preparation transcript.

Read [FRESH_VERIFICATION_REPORT.md](FRESH_VERIFICATION_REPORT.md) and its
[PDF](FRESH_VERIFICATION_REPORT.pdf) for the new verification scope.
[fresh_replay_receipt.json](fresh_replay_receipt.json) is the complete native
receipt. Its SHA-256 is
`980a3d7d52db96a559dd36ad12e317c9b12be9f1c161a2d0d86803aafcc72441`.

The exact model is `icv1-f2m83-t6151469093347-9eba6cd5`, using modulus low
terms `[0,2,4,7]`. The original [archive receipt](receipt.json) and
[archive report PDF](REPORT.pdf) incorrectly print historical alias
`icv1-f2m83-t6151469093347-cdcc5432`, which represents low terms
`[0,1,2,45]`. Those files are preserved unchanged; their explicit numerical
field, subgroup, generator, target and arithmetic checks were correct.
Read the [identity erratum](CURVE_ID_ERRATUM.md) with either original report.
The [invariants note](INVARIANTS.md) uses the corrected identity and explains
the existing root-index and rho algorithms without a novelty claim.

Native sources are [check.cpp](check.cpp) and [check_fresh.cpp](check_fresh.cpp).
The fresh checker includes the former for JSON/field/curve helpers, renames
its archived entry point, and never invokes that entry point. Sources and
historical executable hashes are preserved; compiled executables are excluded
from this directory. `checker_hashes.sha256`, `fresh_checker_hashes.sha256`
and `source_hashes.sha256` retain the original execution-path provenance.

To build and replay, from this directory:

```sh
g++ -std=c++17 -O2 -Wall -Wextra -Werror check_fresh.cpp -o /tmp/n83-fresh-check
/tmp/n83-fresh-check ARCHIVED_FIXTURE_BASE_DIR COMPLETED_PANEL_DIR NEW_RECEIPT_JSON
```

`ARCHIVED_FIXTURE_BASE_DIR` is the repository's
`experiments/koblitz-single-target-n83-20261006` directory.
`COMPLETED_PANEL_DIR` contains the completed fresh process directories and
their `producer.json`, `target.jsonl` and `execution.json` files.
Use a new receipt path: the checker refuses to overwrite a receipt. Run the
verifier after timed children finish. It recursively retains failures,
timeouts, cap outcomes and missing answers rather than counting them as wins.

`receipt.json` is the earlier archived-data verification, not the fresh panel
receipt. Its original SHA-256 remains
`7d4e290524ea9b83ff75d48d20b17449fb60a6991bae787fdfee53ef28d461d5`.
No frozen input or original receipt was rewritten when preparing this copy.

The files originated outside the checkout under
`/workspace/n83_evaluation/independent_check`; preserved receipt paths name
that execution workspace. The linked source locations in the invariants
note likewise identify the reviewed local checkout. Shared study source and
records are managed by the parent evaluation task; this verification copy
does not alter the runner, protocol or timed records.

No independent host run or repository publication is implied. Conductor
tooling and `.conductor/WORKFLOW.md` were unavailable in the provided
checkout; this copy was prepared only within the parent's explicitly assigned
verification path. No git commit or external write was performed here.
