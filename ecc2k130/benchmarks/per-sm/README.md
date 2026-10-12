# Per-SM occupancy receipts

See [WAVES.md](WAVES.md) for the boundary, the target, and the table.

```sh
bash benchmarks/per-sm/run.sh                  # RTX-PRO-6000, waves 1,4,6,8
WAVES_GPU=L40S bash benchmarks/per-sm/run.sh
make bench-waves-modal
```

Logs go to `/tmp/ecc2k130-per-sm` until Modal returns. `freeze.py` refuses
a log that is missing a wave or whose GPU name does not match `WAVES_GPU`.
