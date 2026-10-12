"""Preregistered native-XOR SAT and s=3,4 comparison; see PROTOCOL_EXTENSION.md."""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import sys
import time

import run as base
from test_cubic_quotient import NormalCoordinates, make_curve

HERE = Path(__file__).resolve().parent
ENGINES = base.ENGINES + ("sat-cms-xor",)
old_solve = base.solve
archived_setup = base.Setup


class PlainSpace:
    @staticmethod
    def encode(payload):
        return payload


class PlainSetup:
    """Nested normal-coordinate linear factor base, distinct from the quotient."""

    def __init__(self, n, s, seed, arm, modulus):
        self.n, self.s, self.seed, self.arm = n, s, seed, arm
        self.curve = make_curve(n, modulus)
        self.normal = NormalCoordinates(self.curve.f, seed)
        self.space = PlainSpace()
        self.phases = (0, 1, 2)
        self.counts = collections.Counter()
        self.phase_maps = [[{1 << (i * s + bit): self.normal.basis[(bit + phase) % n]
                             for bit in range(s)} for phase in self.phases]
                           for i in range(3)]


def cms_solve(equations):
    import pycryptosat
    nvars = 3 * base.S
    deadline = time.perf_counter() + base.LIMIT
    solver = pycryptosat.Solver(threads=1)
    monomials = sorted(set().union(*equations))
    auxiliary = {mask: nvars + i + 1 for i, mask in enumerate(
        mask for mask in monomials if mask.bit_count() > 1)}
    clauses = xors = 0
    for mask, out in auxiliary.items():
        inputs = [i + 1 for i in range(nvars) if mask >> i & 1]
        solver.add_clause([out] + [-variable for variable in inputs])
        clauses += 1
        for variable in inputs:
            solver.add_clause([-out, variable])
            clauses += 1
    for eq in equations:
        literals = [auxiliary[mask] if mask in auxiliary else mask.bit_length()
                    for mask in sorted(eq) if mask]
        rhs = (0 in eq)
        if not literals:
            if rhs:
                solver.add_clause([])
                clauses += 1
            continue
        solver.add_xor_clause(literals, rhs)
        xors += 1
    if time.perf_counter() > deadline:
        raise TimeoutError("CryptoMiniSat setup deadline")
    roots = []
    while True:
        remaining = deadline - time.perf_counter()
        if remaining <= 0:
            raise TimeoutError(f"CryptoMiniSat enumeration deadline after {len(roots)} models")
        sat, model = solver.solve(time_limit=remaining)
        if sat is None:
            raise TimeoutError(f"CryptoMiniSat timed out after {len(roots)} models")
        if not sat:
            break
        bits = sum(int(bool(model[i + 1])) << i for i in range(nvars))
        assert all(base.evaleq(eq, bits) == 0 for eq in equations)
        roots.append(bits)
        solver.add_clause([-(i + 1) if bits >> i & 1 else i + 1
                           for i in range(nvars)])
        clauses += 1
    return sorted(roots), {"cms_version": pycryptosat.__version__,
                           "seed_configurable": False, "threads": 1,
                           "auxiliary_variables": len(auxiliary),
                           "cnf_clauses": clauses, "xor_clauses": xors,
                           "models_enumerated": len(roots),
                           "decisions": None}


def extension_solve(equations, engine):
    if engine == "sat-cms-xor":
        return cms_solve(equations)
    return old_solve(equations, engine)


base.solve = extension_solve


def child(payload_bits, seed, index, engine, policy):
    base.S = payload_bits
    base.Setup = archived_setup if policy == "quotient" else PlainSetup
    result = base.fixed_child(seed, index, engine)
    result["payload_bits"] = payload_bits
    result["factor_base_policy"] = policy
    result["protocol_extension"] = True
    return result


def parent(out):
    out.mkdir(parents=True, exist_ok=False)
    files = sorted((HERE / "source").rglob("*.py")) + [HERE / "run.py", HERE / "run_extended.py",
            HERE / "PROTOCOL.md", HERE / "PROTOCOL_EXTENSION.md", HERE / "AMENDMENT.md",
            HERE / "PROTOCOL_PLAIN_SUBSPACE.md"]
    versions = {}
    for pkg in ("z3", "pycryptosat"):
        try:
            m = __import__(pkg)
            versions[pkg] = m.get_version_string() if pkg == "z3" else m.__version__
        except ImportError:
            versions[pkg] = None
    manifest = {"payload_sizes": (2, 3, 4), "seeds": base.SEEDS,
                "fixed_phase": (0, 1, 2), "engines": ENGINES,
                "limits": {"address_space_bytes": 1 << 30, "child_wall_seconds": 70,
                           "equation_seconds": 25, "solver_seconds": base.LIMIT},
                "versions": versions, "python": sys.version, "platform": platform.platform(),
                "source_sha256": {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in files}}
    (out / "manifest.json").write_text(json.dumps(manifest, sort_keys=True, indent=2) + "\n")
    preregistered_preflight(out)
    for policy, sizes in (("quotient", (2,)), ("normal_subspace", (2, 3, 4))):
        destination = out / ("quotient_s2" if policy == "quotient" else "plain_subspace")
        destination.mkdir()
        for s in sizes:
            for seed in base.SEEDS:
                indices = range(6) if policy == "quotient" else (0, 2)
                engines = ("sat-cms-xor",) if policy == "quotient" else ENGINES
                for index in indices:
                    for engine in engines:
                        key = f"s{s}-seed{seed}-i{index}-{engine}"
                        cmd = [sys.executable, str(HERE / "run_extended.py"), "--child",
                               "--s", str(s), "--seed", str(seed), "--policy", policy,
                               "--index", str(index), "--engine", engine]
                        base.run_child(destination, key, cmd)


def preregistered_preflight(out):
    """Reproduce the archived quotient construction failure on both seeds."""
    records = []
    for s in (3, 4):
        for seed in base.SEEDS:
            try:
                archived_setup(base.N, s, seed, "ternary_inline", base.MODULUS)
            except Exception as error:
                records.append({"payload_bits": s, "seed": seed, "status": "setup_failed",
                                "type": type(error).__name__, "message": str(error)})
            else:
                records.append({"payload_bits": s, "seed": seed, "status": "unexpected_success"})
    (out / "quotient_preflight_replay.json").write_text(
        json.dumps(records, sort_keys=True, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path)
    p.add_argument("--child", action="store_true")
    p.add_argument("--s", type=int, choices=(2, 3, 4))
    p.add_argument("--policy", choices=("quotient", "normal_subspace"))
    p.add_argument("--seed", type=int, choices=base.SEEDS)
    p.add_argument("--index", type=int)
    p.add_argument("--engine", choices=ENGINES)
    a = p.parse_args()
    if a.child:
        start = base.clock()
        try:
            print(json.dumps(child(a.s, a.seed, a.index, a.engine, a.policy), sort_keys=True))
        except (TimeoutError, MemoryError) as error:
            print(json.dumps({"status": "timeout" if isinstance(error, TimeoutError) else "oom",
                              "type": type(error).__name__, "message": str(error),
                              "payload_bits": a.s, "seed": a.seed, "index": a.index,
                              "engine": a.engine, "factor_base_policy": a.policy,
                              "timing": {"total_until_failure": base.delta(start)},
                              "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024}),
                  flush=True)
        except Exception as error:
            print(json.dumps({"status": "exception", "type": type(error).__name__,
                              "message": str(error), "payload_bits": a.s,
                              "seed": a.seed, "index": a.index, "engine": a.engine,
                              "factor_base_policy": a.policy}), flush=True)
            raise
    else:
        parent(a.out)
