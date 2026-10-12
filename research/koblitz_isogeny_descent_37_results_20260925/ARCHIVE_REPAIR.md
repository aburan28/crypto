# n37 degree-73 raw archive repair

The original result document named a 40,311-byte `raw.json.gz` with SHA-256
`eb4773d556886672e8b80735fc486c16b73765fd5069ff81133100bfc3eafa90`.
The bytes previously committed at this path were 30,060 bytes, had SHA-256
`e16d66c0d3d529329f87c40991899159d7f6e0aa0656d5e6400b9b08868857be`,
and did not decompress as gzip. That Git blob was
`c820c107751592571889f0d023e7ff867ddaa4f9`. The prior archive could not
support the published summary even though the summary itself was present.

The restored file was extracted from the original
[GitHub Actions run 36087300313](https://github.com/aburan28/crypto/actions/runs/36087300313),
artifact `10843942663` (`koblitz-isogeny-descent-37`), using:

```sh
gh run download 36087300313 -n koblitz-isogeny-descent-37 -D /tmp/n37-archived-run
cp /tmp/n37-archived-run/__w/_temp/koblitz-37-run/raw.json.gz \
  research/koblitz_isogeny_descent_37_results_20260925/raw.json.gz
```

The restored gzip file is byte-for-byte consistent with the original
published size and SHA-256. Its embedded producer-source hashes match the
current frozen `experiment.py`, `contract.json`, and `matrix.py`. The existing
`summary.json` was retained unchanged. The manifest [EVIDENCE.json](EVIDENCE.json)
records the old and restored identities.

Run the standard-library archive check with:

```sh
python3 research/koblitz_isogeny_descent_37_20260925/verify.py \
  research/koblitz_isogeny_descent_37_results_20260925
python3 research/koblitz_isogeny_descent_37_results_20260925/verify_archive.py \
  --out /tmp/n37-archive-replay.json
cmp /tmp/n37-archive-replay.json \
  research/koblitz_isogeny_descent_37_results_20260925/REPLAY.json
```

The second check verifies hashes, the complete paired case grid, matching
transported workloads, matrix replay digests, stage-only accounting, and a
fresh reconstruction of every published aggregate from the raw cases. CI
performs these checks on changes to the archive or its frozen inputs. This is
an archive custody and internal-consistency repair; it is not a new Sage run
or an independent proof of the map. The original registered criterion still
fails, and no end-to-end discrete-logarithm speedup is established.
