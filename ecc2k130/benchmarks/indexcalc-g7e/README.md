# G7e pair-enumeration index calculus

Engineering measurement of the cheapest three-summand decomposition
oracle on one RTX PRO 6000 Blackwell (`g7e`). It does not change the
ECC2K-130 index-calculus exponent.

Boundaries, table, and class:
[`RESEARCH_ECC2K130_G7E_INDEX_CALCULUS.md`](../../../RESEARCH_ECC2K130_G7E_INDEX_CALCULUS.md).
Frozen numbers: [`summary.json`](summary.json).

```sh
make -C ../.. indexcalc-cuda ARCH='-gencode arch=compute_120,code=sm_120'
python3 run.py
```

`--skip-gpu` runs only the degree-5/9 pair-enumeration discrete logs.
