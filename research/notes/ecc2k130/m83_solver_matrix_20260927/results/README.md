# Raw evidence and reproduction

The full receipts are committed as
`raw_receipts.tar.gz.base64` (base64 text containing a gzip-compressed tar).
The compressed bytes are **47,248 bytes**, SHA-256
`6792ab6016437fe2f04a18ebed0d244d5ffb503526db1f2c23f25c1ea98e12a7`.
`raw_manifest.json` records the size and SHA-256 of each of the **131**
embedded JSON files. It is a durable in-repository archive; no local absolute
path is needed to retrieve it. The archive keeps original subprocess stdout,
stderr, return code, elapsed time, source/configuration hashes, the four
quotient-size preflight failures, full-orbit group-operation counts and every
solver timeout, alongside complete successful cases.

In a clean checkout, from this directory:

```bash
base64 -d raw_receipts.tar.gz.base64 | tar -xz -C .
python ../verify_archive.py
python ../summarize.py
```

`verify_archive.py` checks the compressed hash, every file against the
manifest, and all pinned source hashes; `summarize.py` verifies equal input
and equation digests, algebraic model sets and independent group-check
results for all completed matched solver rows. It also matches the
full-orbit target hashes to the original fixed-phase inputs. The derived
`validated_summary.json` is committed separately for easy review.

To repeat the bounded measurements, install Python 3.12, Z3 4.15.3 and
pycryptosat 5.16.0. From the parent directory, on a machine that enforces
`RLIMIT_AS`:

```bash
PYTHONHASHSEED=0 python run_v1.py --out results/fresh_original
PYTHONHASHSEED=0 python run_extended.py --out results/fresh_extension
PYTHONHASHSEED=0 python run_orbit.py --out results/fresh_orbit
```

The original archive hashes correspond exactly to `run_v1.py`. For the
second command, `run_extended.py` imports `run.py`. The orbit runner's
target comparison reads `results/run_20260927` after archive extraction.
Fresh runs must use new output directories to avoid overwriting evidence.
