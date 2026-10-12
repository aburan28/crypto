# Frozen rotated PDP corpus evidence

This directory is the durable record of the single preregistered run of
[`PROTOCOL.md`](../PROTOCOL.md). The source/input freeze was committed at
`39d7327` in [PR #767](https://github.com/aburan28/crypto/pull/767) before
the run. Its protocol-anchored `FROZEN.json` SHA-256 is
`7e0ed24783547ac8a80a83cd62a27976f72798730483121db0c2e55341e45e6b`.
The frozen input manifest SHA-256 is
`6c457c4f216ef942ee6a9e6a3c7832de990d20cd8bbefc518c206788c9d28222`.
The parent #762 merge is `2de218f583d2dddeafdd0180dd02329dc57d53e4`;
its `raw.tar.gz` SHA-256 is
`fad4b4cc48f980349d1ef472d18327368dd0f93f3d1284fc8bbdc20f0d3086d7`.
The prerequisite #763 merge is `ec990d85b61eacccefa9dee97f2a835e7a20b785`.

The exact frozen run command, issued from the repository root with Python
3.13.1 on macOS arm64, was:

```sh
python3 research/notes/ecc2k130/rotated_pdp_corpus_20260925/run.py --out /private/tmp/rotated-pdp-corpus-run-20260925
```

`receipt.json` records the four process commands, UTC start/finish, successful
exit codes, platform, SHA-256 of every raw file, source and input hashes, and
parent hashes. It has SHA-256
`5addca04a9b880caedf5b140086fec35b9f58cecf826b62f600277eb28ddd4fa`.
The archive `raw.tar.gz` is 2,058,524 bytes with SHA-256
`39f16990213e2c7525c6dae92c131127b12f914bbfd6182264249da72bc9546c`;
the adjacent `.sha256` file pins this hash. The eight per-child stdout/stderr
files are retained; all are empty because each child wrote structured JSON
directly and exited zero. There was one run and no failed or substituted arm.

The archive contains, per arm, factor points, full-point and projected
histograms with integer multiplicities, one witness per supported full sum,
every SHA selection attempt, sixteen total targets, producer summary and the
independent verifier receipt. Both complete group-law histograms, not just
the selected target labels, are checked by CI. Read the archive with
`tar -tzf` or extract it with:

```sh
tar -xzf research/notes/ecc2k130/rotated_pdp_corpus_20260925/evidence/raw.tar.gz -C /tmp
```

From the repository root, independently replay the committed archive using
the separately implemented bit-serial verifier:

```sh
python3 research/notes/ecc2k130/rotated_pdp_corpus_20260925/ci_replay.py --evidence research/notes/ecc2k130/rotated_pdp_corpus_20260925/evidence
```

That command first checks the protocol's literal freeze SHA, every frozen
source and input hash, the merged parent source and archive hashes, receipt,
tar member set and individual raw-file hashes. It then rebuilds both full
histograms by direct tuple enumeration with independent field/group law and
checks every target/selection stream. The committed CI workflow uses this
archive-only path; it does not create replacement outcomes or retime a solver.
