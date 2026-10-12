#!/usr/bin/env python3
"""Independent exact-point checks for the frozen n13 O-aware CNF solver panel."""
from __future__ import annotations

import hashlib
import importlib.util
import itertools
import json
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
CNF_DIR = NOTES / "rotated_s3_o_branch_20260925/evidence/producer/n13-m5"
CORPUS = NOTES / "rotated_pdp_corpus_20260925/evidence/raw.tar.gz"
CORPUS_VERIFY = NOTES / "rotated_pdp_corpus_20260925/verify.py"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def curve_module():
    spec = importlib.util.spec_from_file_location("independent_point_corpus", CORPUS_VERIFY)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def oracle_rows():
    with tarfile.open(CORPUS, "r:gz") as archive:
        stream = archive.extractfile("raw/n13-m5/targets.json")
        assert stream is not None
        rows = json.load(stream)
    assert len(rows) == 8
    return rows


def schema_and_truth():
    schema = json.loads((CNF_DIR / "schema.json").read_text())
    assert (schema["field_degree"], schema["field_poly"], schema["m"]) == (13, 0x201b, 5)
    assert (schema["variables"], schema["clauses"]) == (1263, 1195344)
    targets = schema["targets"]
    assert len(targets) == 33 and targets[-1]["id"] == "O"
    parent = curve_module()
    curve = parent.parent_verify.E(parent.parent_verify.GF(13, 0x201b))
    torsion = (None, (0, 1), (1, 0), (1, 1))
    rows = oracle_rows()
    truth = {}
    for qi, row in enumerate(rows):
        assert row["class"] == ("planted" if qi < 4 else "negative")
        q = tuple(row["Q"])
        assert curve.on(q)
        for ti, shift in enumerate(torsion):
            target = targets[4 * qi + ti]
            point = curve.add(q, shift)
            assert target["id"] == f"Q{qi}T{ti}"
            assert target["point"] == (None if point is None else list(point))
            assert target["class"] == row["class"]
            assert 1 <= target["assumption_literal"] <= schema["variables"]
            truth[target["id"]] = row["coset_multiplicities"][ti] > 0
    assert [group["name"] for group in schema["groups"][:5]] == [f"F{i}" for i in range(5)]
    assert sum(len(group["values"]) for group in schema["groups"][:5]) == 15
    return schema, truth, curve, parent


def query_bytes(base: bytes, literal: int, variables: int, clauses: int) -> bytes:
    first, header, rest = base.split(b"\n", 2)
    assert first.startswith(b"c ")
    assert header == f"p cnf {variables} {clauses}".encode()
    assert rest.endswith(b"\n")
    assert 1 <= literal <= variables
    return first + b"\n" + f"p cnf {variables} {clauses + 1}\n".encode() + rest + f"{literal} 0\n".encode()


def parse_solver_output(raw: bytes, exit_code: int, variables: int):
    status = None
    assignments = {}
    for line in raw.decode("utf-8", errors="replace").splitlines():
        line = line.strip()
        if line.startswith("s "):
            verdict = line[2:].strip()
            if verdict not in ("SATISFIABLE", "UNSATISFIABLE", "UNKNOWN"):
                continue
            if status is not None and status != verdict:
                raise ValueError("conflicting status lines")
            status = verdict
        elif line.startswith("v "):
            for token in line[2:].split():
                lit = int(token)
                if lit == 0:
                    continue
                if not 1 <= abs(lit) <= variables:
                    raise ValueError("out-of-range model literal")
                var, value = abs(lit), lit > 0
                if var in assignments and assignments[var] != value:
                    raise ValueError("contradictory model literal")
                assignments[var] = value
    if status == "SATISFIABLE":
        if exit_code != 10 or not assignments:
            raise ValueError("SAT lacks code 10 or model")
        return "SAT", assignments
    if status == "UNSATISFIABLE":
        if exit_code != 20 or assignments:
            raise ValueError("UNSAT lacks code 20 or has model")
        return "UNSAT", {}
    if status == "UNKNOWN":
        return "UNKNOWN", {}
    raise ValueError("no recognized SAT status")


def model_certificate(schema, target, assignments, curve, parent):
    xs = []
    for group in schema["groups"][:5]:
        selected = [i for i, var in enumerate(group["vars"]) if assignments.get(var) is True]
        if len(selected) != 1:
            raise ValueError("factor group lacks exactly one positive literal")
        xs.append(group["values"][selected[0]])
    fibres = [parent.lifts(curve, x) for x in xs]
    if any(not fibre or len(fibre) > 2 for fibre in fibres):
        raise ValueError("factor x is not a complete rational fibre")
    target_point = None if target["point"] is None else tuple(target["point"])
    assert curve.on(target_point)
    combinations = 1
    for fibre in fibres:
        combinations *= len(fibre)
    assert combinations <= 32
    for choice in itertools.product(*fibres):
        total = None
        for point in choice:
            total = curve.add(total, point)
        if total == target_point:
            return {"x": xs, "witness": [list(point) for point in choice],
                    "sign_assignments_checked_upper": combinations}
    raise ValueError("model x tuple has no exact-point signed lift")
