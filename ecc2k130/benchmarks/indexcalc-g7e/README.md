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

SAT (CryptoMiniSat on the host CPU, not the GPU) is a separate receipt.
System Python on this host is PEP 668, so use a venv:

```sh
python3 -m venv ~/ic-venv
~/ic-venv/bin/pip install -r requirements-sat.txt
cd ../../codegen
~/ic-venv/bin/python testdecomp.py
~/ic-venv/bin/python -m unittest test_indexcalc_pairs.py test_indexcalc_e2e.py
cd ../benchmarks/indexcalc-g7e
~/ic-venv/bin/python run_sat.py
```

`run_sat.py` writes `sat.json` and patches `summary.json['sat']`. It does
not change the product-law floor or the GPU pair-scan row.
