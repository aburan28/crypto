#!/usr/bin/env python3
"""Sparse sequential-one-hot, allowed-output CNF from #781 frozen inputs."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import resource
import signal
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
DENSE = HERE.parent / "rotated_s3_o_branch_20260925"
PANELS = (("n2-m4", 2, 4, 0x7), ("n3-m5", 3, 5, 0xb),
          ("n4-m4", 4, 4, 0x13), ("n13-m5", 13, 5, 0x201b))
CAP_SECONDS = 180
CAP_RSS = 512 * 1024 * 1024


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value):
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def rss() -> int:
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw if sys.platform == "darwin" else raw * 1024


def dense_module():
    spec = importlib.util.spec_from_file_location("sparse_pinned_dense_export", DENSE / "export.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def frozen_input(frozen: dict):
    assert sha(Path(__file__)) == frozen["export_sha256"]
    assert sha(HERE / "PROTOCOL.md") == frozen["protocol_sha256"]
    assert sha(DENSE / "export.py") == frozen["dense_export_sha256"]
    assert sha(DENSE / "verify.py") == frozen["dense_verify_sha256"]
    assert sha(DENSE / "FROZEN.json") == frozen["dense_frozen_sha256"]
    assert sha(DENSE / "evidence/receipt.json") == frozen["dense_receipt_sha256"]
    assert sha(DENSE / "evidence/producer/result.json") == frozen["dense_result_sha256"]
    for relative, digest in frozen["dense_dependency_sha256"].items():
        assert sha(HERE.parent / relative) == digest, relative
    for name, *_ in PANELS:
        for file in ("base.cnf", "schema.json", "paths.jsonl.gz"):
            assert sha(DENSE / "evidence/producer" / name / file) == frozen["panels"][name][file]["sha256"]


def sequential_group(primary: list[int], auxiliary: list[int]):
    k = len(primary)
    assert k >= 1 and len(auxiliary) == k - 1
    rows = [tuple(primary)]
    if k == 1:
        return rows
    rows.append((-primary[0], auxiliary[0]))
    for index in range(1, k - 1):
        rows.extend(((-primary[index], auxiliary[index]),
                     (-auxiliary[index - 1], auxiliary[index]),
                     (-primary[index], -auxiliary[index - 1])))
    rows.append((-primary[-1], -auxiliary[-1]))
    assert len(rows) == 3 * k - 3
    return rows


def export_panel(dense, frozen: dict, panel: tuple, out: Path):
    name, n, m, poly = panel
    started = time.perf_counter()
    out.mkdir(parents=True, exist_ok=False)
    original_path = DENSE / "evidence/producer" / name / "schema.json"
    baseline = json.loads(original_path.read_text())
    assert (baseline["panel"], baseline["field_degree"], baseline["m"], baseline["field_poly"]) == panel
    assert [g["name"] for g in baseline["groups"]] == ([f"F{i}" for i in range(m)] +
                                                      [f"S{i}" for i in range(2, m + 1)])
    field = dense.GF(n, poly)
    assert dense.validate_factors(field, [[tuple(p) for p in slot]
                                          for slot in baseline["factor_points"]]) == [
                                              g["values"] for g in baseline["groups"][:m]]
    primary_count = baseline["variables"]
    next_var = primary_count + 1
    aux_groups = []
    clauses = []
    for group in baseline["groups"]:
        vars = group["vars"]
        auxiliaries = list(range(next_var, next_var + len(vars) - 1))
        aux_groups.append({"name": group["name"], "vars": auxiliaries})
        next_var += len(auxiliaries)
        clauses.extend(sequential_group(vars, auxiliaries))
    onehot_count = len(clauses)
    by_name = {group["name"]: group for group in baseline["groups"]}
    implication_count = 0
    allowed_size_counts = {}
    for stage in baseline["stages"]:
        left, right, end = (by_name[stage[key]] for key in ("left", "right", "out"))
        for u, vu in zip(left["values"], left["vars"]):
            for a, va in zip(right["values"], right["vars"]):
                permitted = set(dense.local(field, u, a))
                outputs = [vv for v, vv in zip(end["values"], end["vars"])
                           if v in permitted]
                assert len(outputs) == len(permitted & set(end["values"]))
                clauses.append(tuple([-vu, -va, *outputs]))
                implication_count += 1
                allowed_size_counts[len(outputs)] = allowed_size_counts.get(len(outputs), 0) + 1
    assert len(clauses) == onehot_count + implication_count
    assert implication_count == sum(len(by_name[s["left"]]["vars"]) *
                                    len(by_name[s["right"]]["vars"])
                                    for s in baseline["stages"])
    total_vars = next_var - 1
    schema = {**baseline,
              "dense_schema_sha256": frozen["panels"][name]["schema.json"]["sha256"],
              "dense_base_sha256": frozen["panels"][name]["base.cnf"]["sha256"],
              "dense_paths_sha256": frozen["panels"][name]["paths.jsonl.gz"]["sha256"],
              "encoding": "sinz-sequential-amo-and-allowed-output-v1",
              "primary_variables": primary_count,
              "auxiliary_groups": aux_groups,
              "variables": total_vars,
              "clauses": len(clauses),
              "onehot_clauses": onehot_count,
              "implication_clauses": implication_count}
    save(out / "schema.json", schema)
    with (out / "base.cnf").open("w") as stream:
        stream.write("c frozen sparse O-aware rational recursive S3 base\n")
        stream.write(f"p cnf {total_vars} {len(clauses)}\n")
        for clause in clauses:
            stream.write(" ".join(map(str, clause)) + " 0\n")
    domains = [group["values"] for group in baseline["groups"][:m]]
    paths = dense.enumerate_paths(field, domains)
    dense.write_rows(out / "paths.jsonl.gz", paths)
    assert sha(out / "paths.jsonl.gz") == frozen["panels"][name]["paths.jsonl.gz"]["sha256"]
    dense_count = baseline["clauses"]
    dense_bytes = frozen["panels"][name]["base.cnf"]["bytes"]
    row = {"panel": name, "variables": total_vars, "primary_variables": primary_count,
           "auxiliary_variables": total_vars - primary_count,
           "clauses": len(clauses), "onehot_clauses": onehot_count,
           "implication_clauses": implication_count,
           "allowed_output_size_counts": dict(sorted(allowed_size_counts.items())),
           "bytes": (out / "base.cnf").stat().st_size,
           "dense_variables": primary_count, "dense_clauses": dense_count,
           "dense_bytes": dense_bytes,
           "clause_ratio": [len(clauses), dense_count],
           "byte_ratio": [(out / "base.cnf").stat().st_size, dense_bytes],
           "model_paths": len(paths), "target_labels": len(baseline["targets"]),
           "factor_x_tuples": math.prod(len(group["values"]) for group in baseline["groups"][:m]),
           "schema_sha256": sha(out / "schema.json"),
           "base_sha256": sha(out / "base.cnf"),
           "paths_sha256": sha(out / "paths.jsonl.gz"),
           "panel_wall_seconds": time.perf_counter() - started}
    assert allowed_size_counts.get(0, 0) == (3 if name == "n13-m5" else 0)
    if name == "n13-m5":
        assert (row["variables"], row["clauses"]) == (2517, 4731)
        assert row["bytes"] <= 250000 and row["clauses"] * 100 <= dense_count
    return row


def run(out: Path):
    started, cpu = time.perf_counter(), time.process_time()
    def expired(_signal, _frame):
        raise TimeoutError(f"{CAP_SECONDS}s sparse exporter cap")
    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, CAP_SECONDS)
    out.mkdir(parents=True, exist_ok=False)
    try:
        frozen = json.loads((HERE / "FROZEN.json").read_text())
        frozen_input(frozen)
        dense = dense_module()
        panels = [export_panel(dense, frozen, panel, out / panel[0]) for panel in PANELS]
        result = {"domain": frozen["domain"], "panels": panels,
                  "wall_seconds": time.perf_counter() - started,
                  "cpu_seconds": time.process_time() - cpu,
                  "peak_rss_bytes": rss()}
        assert result["wall_seconds"] <= CAP_SECONDS and result["peak_rss_bytes"] <= CAP_RSS
        save(out / "result.json", result)
    except Exception as error:
        save(out / "failure.json", {"error": repr(error),
                                    "wall_seconds": time.perf_counter() - started,
                                    "cpu_seconds": time.process_time() - cpu,
                                    "peak_rss_bytes": rss()})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args.out)


if __name__ == "__main__":
    main()
