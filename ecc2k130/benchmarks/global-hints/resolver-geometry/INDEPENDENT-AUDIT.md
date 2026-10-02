# Resolver geometry independent audit

The post-run audit is separate from the GPU producer and does not allocate a
GPU.  It reopens the raw archive extracted by `modal_job.py` and fails closed
on source, binary, artifact, schedule, work, resource, corpus, checkpoint or
decision drift.

The audited inputs were:

- measured source commit: `5c91ec05166342bea08251a4ce6351f926f9fc6c`;
- raw archive:
  `/private/tmp/ecc2k-resolver-geometry-gpu-20261002/results.tgz`;
- raw archive SHA-256:
  `02edcd5e654ee16a2077bc4f3dd69518644fde2f9752fc1f970afe7d3580060d`;
- launch receipt SHA-256:
  `e5b7690f299138610dbe6aeb0c5e380dcfc11ebeace66331d8ac4f842c3c7b11`;
- producer result SHA-256:
  `b9e8d2aa388b15c63ee04b2677d91a384f5707c7bf68b2d1f89854ae4e9ac9b9`;
- audit source SHA-256:
  `59b938279294fbc278ae9d481c92aac90599ce84f45810e3704b367e94d8c75a`;
- audit output SHA-256:
  `35f5d10d6333fdbaf23b75cae02b34564f22430be81e14563c15bf97bf1b7b87`.

The exact native invocation was:

```
g++ -O2 -std=c++17 -Wall -Wextra -Werror \
  /private/tmp/resolver_geometry_independent_audit.cpp \
  -o /private/tmp/resolver_geometry_independent_audit
/private/tmp/resolver_geometry_independent_audit --self-test
/private/tmp/resolver_geometry_independent_audit \
  /private/tmp/ecc2k-resolver-geometry-gpu-20261002/results \
  /private/tmp/ecc2k130-resolver-threads-20261002/ecc2k130 \
  /private/tmp/ecc2k-resolver-geometry-gpu-20261002/launch.json \
  /private/tmp/ecc2k-resolver-geometry-gpu-20261002/results.tgz \
  /private/tmp/resolver-geometry-independent-audit.json
```

The same source also passed a C++17 UBSAN self-test build with
`-fsanitize=undefined -fno-sanitize-recover=undefined`.

Audit result: **PASS / DO_NOT_PROMOTE**.  The audit independently checked the
exact 27-file source inventory, 11 binary hashes, 178 in-job artifact hashes,
34 sample-log hashes, all 66 runtime logs (32 correctness and 34 timing), 32
v3 corpus files and 20 checkpoints.  Every runtime log has exactly one
finished row and its registered final iteration count; numeric ledger fields
must parse in full without trailing data.
It recomputed the A/A drift, all five ratio families, four medians, paired
geometry decisions, fixed tie rule, parent-map gate and strict 26,000 M it/s
goal from `samples.tsv`, then required exact agreement with `result.json`.

The archive listing contains 182 safe entries.  A sorted scan rejected any
absolute path or `..` traversal member before the audit used the extraction.
The committed `archive-list.txt` and `archive-receipt.json` bind that check.

`artifact-files.sha256` intentionally excludes `job.log` and `exit-code`:
the Modal wrapper finishes those only after the producer creates its in-job
manifest.  The launch receipt, exit code, raw archive hash and independent
audit bind those wrapper-owned files after collection.
