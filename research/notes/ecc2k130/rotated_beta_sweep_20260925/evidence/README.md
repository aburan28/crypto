# Immutable n19 normal-beta sweep evidence

This directory preserves the one frozen four-beta support run behind
[PR #769](https://github.com/aburan28/crypto/pull/769). The pre-outcome
protocol was committed at `a6519f6`; frozen source/selection at `8553b62`.
The protocol-anchored `FROZEN.json` SHA-256 is
`44ac5784e3544c3bc32ed5b63f0c93ff27029c6a12b4598a26160cf492516543`.
The hash-selected `selection.json` SHA-256 is
`b88d2ca8c33b90650ac466937a01e82d428aa6bd1ef53f6da7f5958c843a5d6d`.
The merged #767 reference archive SHA-256 is
`39f16990213e2c7525c6dae92c131127b12f914bbfd6182264249da72bc9546c`.

The exact measured command from the repository root was:

```sh
python3 research/notes/ecc2k130/rotated_beta_sweep_20260925/run.py --out /private/tmp/rotated-beta-sweep-run-20260925
```

`receipt.json` has SHA-256
`8ce12a489a605745115a7018dfadd10d89059d75b2d0fb1417d38a0914a06f5c`.
It records source/input/selection/reference hashes, platform, UTC times,
process commands and exit status, and SHA-256 for each raw file. All nine
children exited zero: independent selection replay, then four producer and
four independent verifier processes. There was no substitute beta, failed
measured child or cap exceedance. Every child stdout/stderr file is retained;
only selection replay emits a short JSON result.

The final `raw.tar.gz` is 8,931,901 bytes with SHA-256
`fe84aef6a2cf7f6f4c950245c9c8e870354fb750666f5482b81f9a997d107140`.
`raw.tar.gz.sha256` pins it and `SHA256SUMS` covers every evidence file except
this README and the ledger itself. `SHA256SUMS` SHA-256 is
`a36de809220fc3386622589468a5fb412941d55d4dae5e6735a9bd4545e64ba4`.
For each beta, the archive includes all factor points, complete full-point
and projected histograms with multiplicity and witnesses, an exact q-entry
little-endian `target_counts.u32le` array in `[4]H` scalar order, the eight
fixed #767 Q target decisions with all four torsion-coset counts/witnesses,
producer summary and independent verifier receipt.

The first *packaging* attempt (not a measured support run) failed the
fail-closed archive member-set check because macOS `tar` inserted AppleDouble
`._` members. `archive_preflight_failure.json` retains the initial archive
hash and error. The same raw run directory was repacked with
`COPYFILE_DISABLE=1`; no source, selection, target or measured child was
changed. The corrected archive contains exactly the expected raw files.

From the repository root, list or extract the archive with:

```sh
tar -tzf research/notes/ecc2k130/rotated_beta_sweep_20260925/evidence/raw.tar.gz
tar -xzf research/notes/ecc2k130/rotated_beta_sweep_20260925/evidence/raw.tar.gz -C /tmp
```

Run the committed **archive-only** replay with Python 3.13:

```sh
python3 research/notes/ecc2k130/rotated_beta_sweep_20260925/ci_replay.py --evidence research/notes/ecc2k130/rotated_beta_sweep_20260925/evidence
```

It checks the independent protocol freeze anchor, every frozen source and
input hash, the #767 archive hash, selection stream, raw member set and
per-file hashes. Its independent bit-serial/Fermat verifier directly
re-enumerates all `4 * 7^6 = 470596` labelled tuples and all
`4 * 130873 = 523492` subgroup-target counts, witnesses, energies,
overlap and eight fixed Q decisions per beta. The CI workflow runs this
replay on committed evidence only. It does not produce new sweep outcomes.
