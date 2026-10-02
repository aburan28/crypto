# Hosted replay import fix before measurement

[PR validation run 36762769428](https://github.com/aburan28/crypto/actions/runs/36762769428)
passed the frozen-input replay but failed before the exact-source build.
Its broader sparse checkout contained a different experiment's `prepare.py`.
`verify_cold.py` had put that experiment directory ahead of its own script
directory before importing `run_cold`, so Python selected the wrong
`prepare` module and raised `ImportError: cannot import name 'CELLS'`.

The verifier now imports `run_cold` before putting the older experiment
directory on `sys.path`. The runner, solver source, frozen Q, source hashes,
decision rules and local smoke raw bytes are unchanged. The measure job was
skipped by the one-shot dispatch guard; no timed child started. The failed
validation remains in GitHub's log and this note rather than being
overwritten or treated as an experiment result.
