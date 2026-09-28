# First frozen attempt: relation PASS, archive replay path failure

Frozen PR head: `c7c8a1049c38561694d0f794baa4b42d37c3395e`.
Freeze SHA-256: `db93bd6185255a65ffd5b0d838a09db9a89f0b25a742c1f6705745851a60e383`.

The single cold producer and independent verifier both exited 0 and found zero
truth-table discrepancies. The subsequent archive-only `ci_replay.py --evidence`
call with the relative workflow path failed: it launched `verify.py` with a
relative producer directory but changed `cwd` to the source directory. The
verifier therefore could not locate `producer/result.json`. This is an archive
replay path bug, not a failed full-point relation control; the failed process
and original source/freeze are retained here without changing its receipt.

The correction on the next PR head resolves the supplied evidence path before
spawning the independent verifier. Its new freeze is SHA-linked back to this
one, and a new cold outcome is required; this archived attempt is not relabelled
as a passing final CI replay. See `failed_replay.json`, captured stdout/stderr,
`receipt.json`, and the original `FROZEN.json` and `ci_replay.py`.
