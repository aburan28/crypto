"""Nonallocating m83 phase-model Gröbner/FES backend eligibility checks."""
import importlib.util
import json
from pathlib import Path
import shutil

HERE = Path(__file__).resolve().parent
phase = HERE / "results" / "run_20260927"
rows = []
for case in sorted(phase.glob("phase-*.json")):
    wrapper = json.loads(case.read_text())
    equation = next(json.loads(line) for line in wrapper["stdout"].splitlines()
                    if json.loads(line)["stage"] == "equations")
    nvars = equation["nvars"]
    for engine in ("f4-boolean", "f5b-boolean", "fes"):
        rows.append({"case": case.stem, "engine": engine,
                     "status": "unsupported_memory", "boolean_variables": nvars,
                     "squarefree_columns": str(1 << nvars),
                     "reason": "full squarefree monomial/assignment layout exceeds 1 GiB"})
for name, binary in (("msolve-F4", "msolve"), ("WDSat", "WDSat"),
                     ("Sage", "sage"), ("CryptoMiniSat-CLI", "cryptominisat5")):
    rows.append({"engine": name, "status": "unavailable" if not shutil.which(binary) else "available",
                 "path": shutil.which(binary)})
rows.append({"engine": "pycryptosat-native-XOR", "status": "available" if
             importlib.util.find_spec("pycryptosat") else "unavailable"})
(HERE / "results" / "engine_preflight.json").write_text(json.dumps(rows, sort_keys=True, indent=2) + "\n")
