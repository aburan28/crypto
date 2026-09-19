# B200 Modal receipts

See [../../B200.md](../../B200.md) for the boundary, the target, and the
table. This directory holds the frozen Modal logs and the slim JSON the
note cites.

```sh
bash benchmarks/b200/run-modal.sh            # validate, then CLMAD=1 bench
B200_CLMAD=0 bash benchmarks/b200/run-modal.sh   # software arm only (skip validate)
```

`make bench-b200-modal` is the bench half; `make validate-b200-modal` is
the planted-log half. Workers stay automatic.
