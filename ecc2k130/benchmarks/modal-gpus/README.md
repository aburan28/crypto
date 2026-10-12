# Modal GPU survey receipts

See [SURVEY.md](SURVEY.md) for the boundary, the target, and the table.
This directory holds the catalog, the freezer, and the frozen Modal logs
the note cites.

```sh
bash benchmarks/modal-gpus/run-all.sh
bash benchmarks/modal-gpus/run-one.sh L40S
make bench-modal-gpu SURVEY_GPU=H100!
make bench-modal-gpu SURVEY_GPU=T4 SURVEY_CLMAD=0
```

Workers stay automatic. T4 is CLMAD=0; every other catalog GPU is
CLMAD=1. Logs go to `/tmp/ecc2k130-survey` until Modal returns, then
copy into this directory.

`python3 benchmarks/modal-gpus/summarize.py` rebuilds `summary.json`
from the frozen receipts. It does not invent a rate.

Every pinned SKU returned three valid repeats. No SKU beat the 6000
reference 15.115792 B/s. See [SURVEY.md](SURVEY.md) §3.

