#!/usr/bin/env python3
"""Frozen O-aware recursive-S3 one-hot DIMACS exporter; see PROTOCOL.md."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import itertools
import json
import math
import resource
import signal
import sys
import tarfile
import time
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
CORPUS = NOTES / "rotated_pdp_corpus_20260925/evidence/raw.tar.gz"
PRIOR = NOTES / "rotated_s3_candidate_20260925/evidence/n13-m5/producer/result.json"
PARENT = NOTES / "rotated_subspace_support_20260925/gate.py"
GATE_INPUTS = NOTES / "rotated_m56_export_gate_20260925/INPUTS.json"
GATE_FROZEN = NOTES / "rotated_m56_export_gate_20260925/FROZEN.json"
GATE_N13 = NOTES / "rotated_m56_export_gate_20260925/evidence/n13-m5.json"
PANELS = (("n2-m4", 2, 4, 0x7), ("n3-m5", 3, 5, 0xb),
          ("n4-m4", 4, 4, 0x13), ("n13-m5", 13, 5, 0x201b))
CAP_SECONDS = 180
CAP_RSS = 512 * 1024 * 1024


class InputDomainError(ValueError):
    pass


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value):
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def rss() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def remainder(a: int, modulus: int) -> int:
    while a.bit_length() >= modulus.bit_length():
        a ^= modulus << (a.bit_length() - modulus.bit_length())
    return a


def irreducible(poly: int, n: int) -> bool:
    if poly.bit_length() != n + 1 or not poly & 1:
        return False
    for degree in range(1, n // 2 + 1):
        for low in range(1, 1 << degree, 2):
            if remainder(poly, (1 << degree) | low) == 0:
                return False
    return True


class GF:
    def __init__(self, n: int, poly: int):
        self.n, self.poly, self.size = n, poly, 1 << n
        self.ops = Counter()
        self.inverses = {}
        self.as_roots = None
        assert irreducible(poly, n)

    def mul(self, a: int, b: int) -> int:
        self.ops["field_mul"] += 1
        out = 0
        while b:
            if b & 1:
                out ^= a
            b >>= 1
            a <<= 1
            if a & self.size:
                a ^= self.poly
        return out

    def square(self, a: int) -> int:
        self.ops["field_square"] += 1
        return self.mul(a, a)

    def inv(self, a: int) -> int:
        assert a
        if a not in self.inverses:
            self.ops["field_inverse"] += 1
            exponent, base, answer = self.size - 2, a, 1
            while exponent:
                if exponent & 1:
                    answer = self.mul(answer, base)
                base = self.square(base)
                exponent >>= 1
            assert self.mul(a, answer) == 1
            self.inverses[a] = answer
        return self.inverses[a]

    def as_map(self):
        if self.as_roots is None:
            roots = {}
            for z in range(self.size):
                roots.setdefault(self.square(z) ^ z, []).append(z)
            assert all(len(row) == 2 for row in roots.values())
            self.as_roots = roots
        return self.as_roots

    def sqrt(self, a: int) -> int:
        result = a
        for _ in range(self.n - 1):
            result = self.square(result)
        assert self.square(result) == a
        return result

    def fibre(self, x: int):
        if x == 0:
            return ((0, 1),)
        h = x ^ self.square(self.inv(x))
        return tuple(sorted((x, self.mul(x, z)) for z in self.as_map().get(h, ())))

    def roots(self, a: int, b: int):
        product = self.mul(a, b)
        A, B, C = self.square(a ^ b), product, self.square(product) ^ 1
        if A == 0:
            return () if B == 0 else (self.mul(C, self.inv(B)),)
        if B == 0:
            return (self.sqrt(self.mul(C, self.inv(A))),)
        h = self.mul(self.mul(A, C), self.inv(self.square(B)))
        scale = self.mul(B, self.inv(A))
        return tuple(sorted(self.mul(scale, z) for z in self.as_map().get(h, ())))


def parent_module():
    spec = importlib.util.spec_from_file_location("s3_o_basis_parent", PARENT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def state_sort(value):
    return (-1 if value is None else value)



def target_state(point):
    return None if point is None else point[0]


def validate_factors(f: GF, factors):
    if not factors or any(not slot for slot in factors):
        raise InputDomainError("empty factor slot")
    domains = []
    for slot in factors:
        actual = {tuple(p) for p in slot}
        for x, y in actual:
            if not (0 <= x < f.size and 0 <= y < f.size):
                raise InputDomainError("point outside field")
        xs = sorted({p[0] for p in actual})
        for x in xs:
            lift = set(f.fibre(x))
            if not lift:
                raise InputDomainError(f"nonlift factor x={x}")
            if {p for p in actual if p[0] == x} != lift:
                raise InputDomainError(f"sign-incomplete or invalid factor x={x}")
        domains.append(xs)
    return domains


def toy_input(f: GF, m: int):
    points = sorted(p for x in range(f.size) for p in f.fibre(x))
    factors = [points[:] for _ in range(m)]
    targets = [{"id": "O", "point": None, "class": "identity"}]
    targets += [{"id": f"p-{x}-{y}", "point": [x, y], "class": "all_points"}
                for x, y in points]
    return factors, targets, None


def n13_input(f: GF):
    with tarfile.open(CORPUS, "r:gz") as tar:
        content = tar.extractfile("raw/n13-m5/factors.json").read()
    factors = [[tuple(p) for p in slot] for slot in json.loads(content)]
    prior = json.loads(PRIOR.read_text())
    targets = [{"id": f"Q{row['Q_index']}T{row['T_index']}",
                "point": row["target"], "class": row["target_class"]}
               for row in prior["targets"]]
    assert len(targets) == 32 and targets[12]["point"] == [7256, 3272]
    targets.append({"id": "O", "point": None, "class": "identity"})
    parent = parent_module()
    field = parent.Field(13, [0, 1, 3, 4])
    bases = parent.subspace_basis(parent.normal_conjugates(field, 3), 5, 2)
    mask_by_x = []
    for basis in bases:
        mapping = {}
        for mask in range(4):
            x = 0
            for bit, element in enumerate(basis):
                if mask >> bit & 1:
                    x ^= element
            assert x not in mapping
            mapping[x] = mask
        mask_by_x.append(mapping)
    assert [mask_by_x[i][x] for i, x in enumerate((0, 0, 0, 6433, 217))] == [0, 0, 0, 2, 1]
    assert all(all(p in f.fibre(p[0]) for p in slot) for slot in factors)
    return factors, targets, mask_by_x


def local(f: GF, left, right: int):
    if left is None:
        return (right,)
    roots = list(f.roots(left, right))
    if left == right:
        roots.append(None)
    return tuple(sorted(roots, key=state_sort))


def make_groups(f: GF, domains, targets):
    m = len(domains)
    reachable = None
    states = []
    for length in range(2, m + 1):
        if length == 2:
            values = {v for a in domains[0] for b in domains[1] for v in local(f, a, b)}
        else:
            values = {v for u in reachable for a in domains[length - 1]
                      for v in local(f, u, a)}
        reachable = values
        admitted = set(values) | {None}
        if length == m:
            admitted.update(target_state(row["point"]) for row in targets)
        states.append(sorted(admitted, key=state_sort))
    groups = []
    next_var = 1
    for i, domain in enumerate(domains):
        variables = list(range(next_var, next_var + len(domain)))
        groups.append({"name": f"F{i}", "values": domain, "vars": variables})
        next_var += len(domain)
    for length, domain in enumerate(states, start=2):
        variables = list(range(next_var, next_var + len(domain)))
        groups.append({"name": f"S{length}", "values": domain, "vars": variables})
        next_var += len(domain)
    return groups, next_var - 1


def stages(m: int):
    result = [{"left": "F0", "right": "F1", "out": "S2"}]
    result += [{"left": f"S{i}", "right": f"F{i}", "out": f"S{i+1}"}
               for i in range(2, m)]
    return result


def clauses(f: GF, groups, stage_rows):
    by_name = {g["name"]: g for g in groups}
    out = []
    case_count = Counter()
    for group in groups:
        vars = group["vars"]
        out.append(tuple(vars))
        out.extend((-a, -b) for a, b in itertools.combinations(vars, 2))
    for stage in stage_rows:
        left, right, end = (by_name[stage[key]] for key in ("left", "right", "out"))
        for u, vu in zip(left["values"], left["vars"]):
            for a, va in zip(right["values"], right["vars"]):
                allowed = set(local(f, u, a))
                case_count["O_plus_factor" if u is None else
                           "equal_x" if u == a else "distinct_x"] += 1
                for v, vv in zip(end["values"], end["vars"]):
                    if v not in allowed:
                        out.append((-vu, -va, -vv))
    return out, case_count


def write_cnf(path: Path, variables: int, clauses_list):
    with path.open("w") as stream:
        stream.write("c frozen O-aware sign-complete rational recursive S3 base\n")
        stream.write(f"p cnf {variables} {len(clauses_list)}\n")
        for row in clauses_list:
            stream.write(" ".join(map(str, row)) + " 0\n")


def enumerate_paths(f: GF, domains):
    rows = []
    for x_tuple in itertools.product(*domains):
        paths = [(s,) for s in local(f, x_tuple[0], x_tuple[1])]
        for factor_x in x_tuple[2:]:
            paths = [path + (next_state,) for path in paths
                     for next_state in local(f, path[-1], factor_x)]
        rows.extend({"factor_x": list(x_tuple), "states": list(path)} for path in paths)
    rows.sort(key=lambda row: (row["factor_x"], [state_sort(v) for v in row["states"]]))
    return rows


def write_rows(path: Path, rows):
    with path.open("wb") as file:
        with gzip.GzipFile(filename="", mode="wb", fileobj=file, mtime=0) as zipped:
            for row in rows:
                zipped.write((json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n").encode())


def target_rows(targets, groups, paths, mask_by_x):
    final = groups[-1]
    literal = {value: var for value, var in zip(final["values"], final["vars"])}
    output = []
    for row in targets:
        state = target_state(row["point"])
        selected = [path for path in paths if path["states"][-1] == state]
        masks = {tuple(path["factor_x"]) for path in selected}
        affine = {tuple(path["factor_x"]) for path in selected
                  if None not in path["states"][:-1]}
        exceptions = sorted(masks - affine)
        result = {**row, "assumption_literal": literal[state],
                  "model_path_count": len(selected), "model_x_tuple_count": len(masks),
                  "o_prefix_path_count": sum(None in path["states"][:-1] for path in selected),
                  "exceptional_only_x_tuples": [list(x) for x in exceptions]}
        if mask_by_x is not None:
            result["exceptional_only_masks"] = [[mask_by_x[i][x]
                                                  for i, x in enumerate(xs)]
                                                 for xs in exceptions]
        output.append(result)
    return output


def run_panel(name: str, n: int, m: int, poly: int, out: Path):
    out.mkdir(parents=True, exist_ok=False)
    f = GF(n, poly)
    factors, targets, mask_by_x = n13_input(f) if n == 13 else toy_input(f, m)
    domains = validate_factors(f, factors)
    assert len(domains) == m
    assert all(row["point"] is None or tuple(row["point"]) in f.fibre(row["point"][0])
               for row in targets)
    groups, variables = make_groups(f, domains, targets)
    stage_rows = stages(m)
    clause_rows, cases = clauses(f, groups, stage_rows)
    paths = enumerate_paths(f, domains)
    targets_out = target_rows(targets, groups, paths, mask_by_x)
    if n == 13:
        assert [0, 0, 0, 2, 1] in targets_out[12]["exceptional_only_masks"]
        assert any(path["factor_x"] == [0, 0, 0, 6433, 217] and
                   path["states"][0] is None and
                   path["states"][-1] == 7256 for path in paths)
    schema = {"panel": name, "field_degree": n, "field_poly": poly,
              "m": m, "groups": groups, "stages": stage_rows,
              "factor_points": [[list(p) for p in slot] for slot in factors],
              "targets": targets_out,
              "mask_by_x": [{str(x): mask for x, mask in row.items()}
                            for row in mask_by_x] if mask_by_x else None,
              "variables": variables, "clauses": len(clause_rows)}
    save(out / "schema.json", schema)
    write_cnf(out / "base.cnf", variables, clause_rows)
    write_rows(out / "paths.jsonl.gz", paths)
    return {"panel": name, "field_degree": n, "m": m,
            "factor_x_tuples": math.prod(map(len, domains)),
            "signed_point_tuples": math.prod(map(len, factors)),
            "model_paths": len(paths), "variables": variables,
            "clauses": len(clause_rows), "targets": len(targets),
            "target_rows": targets_out, "transition_cases": dict(sorted(cases.items())),
            "field_operations": dict(sorted(f.ops.items())),
            "base_sha256": sha(out / "base.cnf"),
            "schema_sha256": sha(out / "schema.json"),
            "paths_sha256": sha(out / "paths.jsonl.gz")}


def input_negative_controls():
    f3 = GF(3, 0xb)
    try:
        validate_factors(f3, [[(2, 0)]])
    except InputDomainError as error:
        nonlift = str(error)
    else:
        raise AssertionError("nonlift factor input accepted")
    f2 = GF(2, 0x7)
    try:
        validate_factors(f2, [[f2.fibre(1)[0]]])
    except InputDomainError as error:
        sign = str(error)
    else:
        raise AssertionError("sign-incomplete factor input accepted")
    assert nonlift == "nonlift factor x=2"
    assert sign == "sign-incomplete or invalid factor x=1"
    return {"nonlift": nonlift, "sign_incomplete": sign}


def run(out: Path):
    started, cpu = time.perf_counter(), time.process_time()
    def expired(_signal, _frame):
        raise TimeoutError(f"{CAP_SECONDS}s exporter cap")
    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, CAP_SECONDS)
    out.mkdir(parents=True, exist_ok=False)
    try:
        frozen = json.loads((HERE / "FROZEN.json").read_text())
        assert sha(Path(__file__)) == frozen["export_sha256"]
        assert sha(HERE / "PROTOCOL.md") == frozen["protocol_sha256"]
        assert sha(HERE / "PROOF.md") == frozen["proof_sha256"]
        assert sha(CORPUS) == frozen["corpus_sha256"]
        assert sha(PRIOR) == frozen["prior_sha256"]
        assert sha(PARENT) == frozen["parent_sha256"]
        assert sha(GATE_INPUTS) == frozen["gate_inputs_sha256"]
        assert sha(GATE_FROZEN) == frozen["gate_frozen_sha256"]
        assert sha(GATE_N13) == frozen["gate_n13_sha256"]
        controls = input_negative_controls()
        panel_rows = [run_panel(name, n, m, poly, out / name)
                      for name, n, m, poly in PANELS]
        result = {"domain": frozen["domain"], "panels": panel_rows,
                  "negative_input_controls": controls,
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
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    run(args.out)


if __name__ == "__main__":
    main()
