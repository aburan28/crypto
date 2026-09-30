# Frobenius quotient follow-up: independent SAT reproduction (in progress)

Independent re-implementation of the target-only four-point decomposition
formulations recorded in
[`../frobenius_quotient_followup_20260926`](../frobenius_quotient_followup_20260926/README.md),
whose original solver code is not in this repository. Stage diagnostic only;
`S`, end-to-end cost and speedup are unset.

- `encode.py` builds an explicit field circuit (AND gates, native XOR or
  clause-expanded XOR) from the target and factor-space basis only.
- `run.py` runs one fresh single-thread CryptoMiniSat 5.16.0 solver per
  attempt. It verifies each proposed tuple on the curve, blocks rejected
  tuples and continues within the same budget.
- `test_encode.py` checks the S3 formula, the XOR expansion, the field-product
  circuit and planted-witness satisfiability.
- `summarize.py` builds the results tables.

```sh
pip install pycryptosat==5.16.0
python3 -m unittest test_encode -v
python3 run.py --n 13 --s 4 --budget 2 --out /tmp/n13.jsonl
python3 summarize.py /tmp/n13.jsonl
```

Status: 2 s grids at n = 13 and 19 are committed in `results/`. The 300 s
grids are running, and the write-up will follow.
