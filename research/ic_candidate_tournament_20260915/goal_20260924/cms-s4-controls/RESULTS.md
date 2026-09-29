# Wide-S4/CryptoMiniSat control pilot: no solver verdict

The [registered four-point panel](panel.json) has ended. All four source-bound
exporter runs succeeded. The first dispatch then failed at the manifest-reader
gate (`KeyError: 'mitm'`) before starting any CryptoMiniSat job. The
[versioned gate repair](GATE-REPAIR.md) checked the sealed source exports and
started each of the four originally registered SAT jobs exactly once. Every
child exited with code `-6` before reading its instance: macOS `dyld` could
not load `@rpath/libcryptominisat5.5.14.dylib` from the relocated executable.
All four SAT stdout files are empty. This is a binary packaging failure; the
pilot has **no SAT/UNSAT/timeout verdict, model, or point witness**. Neither
positive query is a discovered relation and neither negative query is proved
infeasible by SAT. The exact feasibility labels come only from the parent
independent group oracle.

| Trial | Parent exact label | Export | Original adapter | Continued SAT process | SAT finding |
| --- | --- | --- | --- | --- | --- |
| 0 | no three-point relation | valid | `INVALID_EXPORT` | `SOLVER_ERROR`, `dyld`, `-6` | unknown |
| 1 | no three-point relation | valid | `INVALID_EXPORT` | `SOLVER_ERROR`, `dyld`, `-6` | unknown |
| 3 | relation exists | valid | `INVALID_EXPORT` | `SOLVER_ERROR`, `dyld`, `-6` | unknown |
| 10 | relation exists | valid | `INVALID_EXPORT` | `SOLVER_ERROR`, `dyld`, `-6` | unknown |

The original [stage-A raw archive](stage-a-evidence.tar.gz) is SHA-256
`2b1fa7abb093c94e741cebecdce35ea8fd59dbc0fb301e9aac9f1b95e0b195ef`.
The [stage-B raw archive](stage-b-evidence.tar.gz), including all four
commands, process metrics, stdout/stderr, original exports, source/build
receipts, and continuation code, is SHA-256
`c9246c617ef5d60ecc5281631decb6721e985ca3604514b6546c6177f05b064b`.
Its summary is SHA-256
`ffb5b334d62c6eeb8c807f7083e049534e1f5e0c70ce748f4f8b01a8e7402b07`.
The archive is independently replayed by `test_cms_s4_evidence.py`. The
original results remain unchanged and both failed stages remain visible.

The panel's no-retry rule ends this pilot. A future SAT test needs a **new,
pre-registered panel with different disclosed controls** and a loader preflight
of the executable at the exact path from which it will be run. The executable
and its transitive dynamic libraries need source/dependency build receipts
before any external-SAT result can be called source-bound. These controls
cannot estimate natural relation yield. No complete SAT IC solve or online
speedup is established.
