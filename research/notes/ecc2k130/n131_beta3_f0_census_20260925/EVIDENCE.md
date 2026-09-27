# Exact F0 census evidence

The frozen runner began at 2026-09-25 13:36:07.773774 UTC and finished at 13:52:40.891906 UTC on the local macOS arm64 host. Its measured source commit is `a07e5e932a653a509fcd89b46d4d9b25e5acdac2`; the pre-outcome source/input freeze is `FROZEN.json` (SHA-256 `bfb69239607652ae71ccfc52cabe939d49ae998596440270396cc08091641dd0`). All four children exited zero, the two independent replays returned `pass`, and the runner status is `complete`. There were no measured failed attempts or caps. The prior `f0b973d` preregistration and `a07e5e9` review amendment are preserved in Git history.

`evidence/receipt.json` (SHA-256 `335d68c319ca9a76fe6ba62cd11107f26ecc27d925b95a3ba44314061d70e631`) records commands, UTC boundaries, exit codes, local Python/psutil versions, source commit and hashes, sampled child RSS, and the SHA-256 of every raw file. `evidence/pilot/` and `evidence/full/` each contain a producer summary, independent replay report and 8,192-row chunk hashes/counts. The top-level stdout/stderr files preserve exact child output, including empty stderr. The per-row data stream is intentionally not committed: the deterministic Gray traversal and independent natural-mask replay reconstruct it byte-for-byte, and both full digests are `b8b03e28e783ebec12be04c82542d74448896e3a128c0a5784a77a0a688da30d`. This compact raw evidence occupies about 104 KiB rather than a large two-million-row serialization.

From repository root, the archive-only check is:

```sh
python3.12 research/notes/ecc2k130/n131_beta3_f0_census_20260925/ci_replay.py --evidence research/notes/ecc2k130/n131_beta3_f0_census_20260925/evidence
```

That check verifies the frozen source and independent n13/n19 and 64-mask controls, every committed raw file hash, stage/count invariants and the archived full independent replay result. It does not rerun the 902-second local full verifier in CI. For a fresh full re-enumeration of all mask rows, run `replay.py --input .../INPUT.json --summary .../evidence/full/summary.json --output /private/tmp/fresh-f0-replay.json` with Python 3.12 and the frozen sources; the preregistered verifier wall/RSS caps still apply. The original child command lines and resource costs are in `evidence/receipt.json`.
