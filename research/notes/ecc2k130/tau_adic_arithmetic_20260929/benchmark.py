#!/usr/bin/env python3
"""Standalone, variable-time known-scalar arithmetic experiment; not a solver.

Run with sage -python. Uses only deterministic, synthetic inputs.
"""
import argparse
import ctypes
import ctypes.util
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import shutil
import statistics
import subprocess
import sys
import time
import traceback

METHODS = ("binary", "binary_naf", "tau_naf", "reduced_tau_naf", "sage_native")
TAPS = {5: (2, 0), 31: (3, 0), 83: (45, 2, 1, 0), 131: (13, 2, 1, 0)}
SEEDS = (20260929, 20260930)


def sha(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()


def command(args):
    try:
        p = subprocess.run(args, capture_output=True, text=True, timeout=10)
        return {"returncode": p.returncode, "stdout": p.stdout, "stderr": p.stderr}
    except (OSError, subprocess.TimeoutExpired) as e:
        return {"error": str(e)}


def read(path):
    try:
        return Path(path).read_text().strip()
    except OSError:
        return None


def host():
    allowed = sorted(os.sched_getaffinity(0))
    cpu = allowed[0]
    os.sched_setaffinity(0, {cpu})
    nodes = sorted(Path(f"/sys/devices/system/cpu/cpu{cpu}").glob("node*"))
    node = int(nodes[0].name[4:]) if nodes else None
    numa = {"selected_node": node, "bind_succeeded": False}
    lib = ctypes.util.find_library("numa")
    if node is not None and lib:
        try:
            n = ctypes.CDLL(lib, use_errno=True)
            n.set_mempolicy.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_ulong), ctypes.c_ulong]
            mask = ctypes.c_ulong(1 << node)
            rc = n.set_mempolicy(2, ctypes.byref(mask), ctypes.sizeof(mask) * 8)
            numa.update(returncode=rc, errno=ctypes.get_errno(), bind_succeeded=(rc == 0))
            mode = ctypes.c_int()
            actual = ctypes.c_ulong()
            n.get_mempolicy.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_ulong), ctypes.c_ulong, ctypes.c_void_p, ctypes.c_ulong]
            rc2 = n.get_mempolicy(ctypes.byref(mode), ctypes.byref(actual), 64, None, 0)
            numa.update(query_returncode=rc2, observed_policy=mode.value if rc2 == 0 else None,
                        observed_mask=actual.value if rc2 == 0 else None)
        except (OSError, AttributeError) as e:
            numa["error"] = str(e)
    return {
        "platform": platform.platform(), "python": sys.version,
        "lscpu": command(["lscpu", "-J"]),
        "allowed_cpus_before": allowed, "affinity_after": sorted(os.sched_getaffinity(0)),
        "core_exclusive": False, "numa": numa,
        "visible_nodes": {p.name: read(p / "cpulist") for p in Path("/sys/devices/system/node").glob("node[0-9]*")},
        "memory_type": "unknown: DIMM type not exposed by this VM",
        "memory_limit": read("/sys/fs/cgroup/memory.max"),
        "cpu_limit": read("/sys/fs/cgroup/cpu.max"),
        "meminfo": read("/proc/meminfo"),
        "load_before": os.getloadavg(),
        "processes_before": command(["ps", "-eo", "comm,pcpu,pmem"]),
        "nvcc": shutil.which("nvcc"), "nvidia_smi": shutil.which("nvidia-smi"),
    }


def ring_mul(x, y, mu):
    a, b = x
    c, d = y
    return a*c - 2*b*d, a*d + b*c + mu*b*d


def norm(z, mu):
    a, b = z
    return a*a + mu*a*b + 2*b*b


def annihilator(m, mu):
    z = (1, 0)
    for _ in range(m):
        z = ring_mul(z, (0, 1), mu)
    return z[0] - 1, z[1]


def reduced(k, delta, mu):
    a, b = delta
    den = norm(delta, mu)
    # Round the exact coefficients of k/delta, then search a fixed neighborhood.
    q0 = (2*k*(a + mu*b) + den) // (2*den)
    q1 = (-2*k*b + den) // (2*den)
    candidates = []
    for u in range(q0 - 1, q0 + 2):
        for v in range(q1 - 1, q1 + 2):
            c, d = ring_mul((u, v), delta, mu)
            z = k-c, -d
            candidates.append(z)
    return min(candidates, key=lambda z: (norm(z, mu), z))


def tau_digits(z, mu):
    a, b = z
    digits = []
    while a or b:
        u = 2 - ((a - 2*b) % 4) if a & 1 else 0
        digits.append(u)
        even = a-u
        a, b = b + mu*(even//2), -(even//2)
        if len(digits) > 4096:
            raise ArithmeticError("tau recoding did not terminate")
    return digits


def binary_digits(k, naf=False):
    sign = -1 if k < 0 else 1
    k = abs(k)
    out = []
    while k:
        u = (2-k%4 if naf else 1) if k & 1 else 0
        out.append(sign*u)
        k = (k-u)//2
    return out


def point_json(p):
    # Both Sage's Givaro (small fields) and NTL GF2E expose polynomial().
    return None if p is None else [
        hex(sum(int(bit) << i for i, bit in enumerate(c.polynomial().list())))
        for c in p
    ]


class Arithmetic:
    def __init__(self, field, a, counted=False):
        self.f, self.a = field, field(a)
        self.counts = {k: 0 for k in ("additions", "doublings", "frobenius", "field_M", "field_S", "field_I")} if counted else None

    def count(self, **kwargs):
        if self.counts is not None:
            for k, v in kwargs.items():
                self.counts[k] += v

    @staticmethod
    def neg(p):
        return None if p is None else (p[0], p[0]+p[1])

    def double(self, p):
        if p is None or not p[0]:
            return None
        self.count(doublings=1, field_M=2, field_S=2, field_I=1)
        x, y = p
        t = x + y*(~x)
        xx = t*t + t + self.a
        return xx, x*x + (t+1)*xx

    def add(self, p, q):
        if p is None:
            return q
        if q is None:
            return p
        x, y = p
        u, v = q
        if x == u:
            return self.double(p) if y == v else None
        self.count(additions=1, field_M=2, field_S=1, field_I=1)
        t = (y+v)*(~(x+u))
        xx = t*t+t+x+u+self.a
        return xx, t*(x+xx)+xx+y

    def frob(self, p):
        if p is None:
            return None
        self.count(frobenius=1, field_S=2)
        return p[0]*p[0], p[1]*p[1]

    def evaluate(self, p, ds, use_tau):
        r = None
        neg = self.neg(p)
        advance = self.frob if use_tau else self.double
        for u in reversed(ds):
            r = advance(r)
            if u:
                r = self.add(r, p if u == 1 else neg)
        return r


def make_curve(m, a):
    ring = PolynomialRing(GF(2), "z")
    z = ring.gen()
    polynomial = z**m + sum(z**j for j in TAPS[m])
    assert polynomial.is_irreducible()
    f = GF(2**m, name="z", modulus=polynomial)
    e = EllipticCurve(f, [1, a, 0, 0, 1])
    mu = 2*a-1
    delta = annihilator(m, mu)
    # #E(F_2^m) = Norm(tau^m-1), checked against Sage on the exhaustive curves.
    order = norm(delta, mu)
    return {"m": m, "a": a, "mu": mu, "f": f, "e": e, "delta": delta,
            "order": order, "field_modulus": str(polynomial)}


def sage_result(c, p, k):
    s = c["e"](0) if p is None else c["e"](p)
    r = Integer(k)*s
    return None if r.is_zero() else (r[0], r[1])


def multiply(c, p, k, method, counted=False):
    if method == "sage_native":
        return sage_result(c, p, k), None, None
    if method.startswith("binary"):
        ds = binary_digits(k, method == "binary_naf")
        use_tau = False
    else:
        z = reduced(k, c["delta"], c["mu"]) if method == "reduced_tau_naf" else (k, 0)
        ds = tau_digits(z, c["mu"])
        use_tau = True
    ops = Arithmetic(c["f"], c["a"], counted)
    answer = ops.evaluate(p, ds, use_tau)
    return answer, ops.counts, {"digits": len(ds), "nonzero": sum(u != 0 for u in ds)}


def identities(c, p):
    ops = Arithmetic(c["f"], c["a"])
    t = ops.frob(p)
    lhs = ops.add(ops.frob(t), ops.double(p))
    assert lhs == (t if c["mu"] == 1 else ops.neg(t))
    r = p
    for _ in range(c["m"]):
        r = ops.frob(r)
    assert r == p


def check_case(c, p, k):
    expected = sage_result(c, p, k)
    for method in METHODS[:-1]:
        got, _, _ = multiply(c, p, k, method)
        assert got == expected, (c["m"], c["a"], point_json(p), k, method)
    for z in ((k, 0), reduced(k, c["delta"], c["mu"])):
        ds = tau_digits(z, c["mu"])
        acc = (0, 0)
        for d in reversed(ds):
            acc = ring_mul(acc, (0, 1), c["mu"])
            acc = acc[0]+d, acc[1]
        assert acc == z
        assert all(not(ds[i] and ds[i+1]) for i in range(len(ds)-1))


def fixture(c, rng):
    f, m = c["f"], c["m"]
    ops = Arithmetic(f, c["a"])
    while True:
        x = f.from_integer(rng.getrandbits(m))
        if not x:
            continue
        v = x + f(c["a"]) + ~(x*x)
        if v.trace():
            continue
        h, t = v, v
        for _ in range(1, (m+1)//2):
            t = t*t
            t = t*t
            h += t
        assert h*h+h == v
        p = ops.double(ops.double((x, x*h)))
        if p is not None:
            assert sage_result(c, p, c["order"]//4) is None
            return p


def timing(c, cases, method, expected):
    start = time.perf_counter_ns()
    outputs = [multiply(c, p, k, method)[0] for p, k in cases]
    elapsed = time.perf_counter_ns()-start
    assert outputs == expected
    return {"ns": elapsed, "ns_per_scalar": elapsed/len(cases),
            "output_sha256": sha([point_json(p) for p in outputs])}


def run(out):
    t0 = time.perf_counter()
    total_cases = 0
    for a in (0, 1):
        c = make_curve(5, a)
        assert int(c["e"].cardinality()) == c["order"]
        for q in c["e"].points():
            p = None if q.is_zero() else (q[0], q[1])
            identities(c, p)
            for k in range(-17, 18):
                check_case(c, p, k)
                total_cases += 1
    out["exhaustive"] = {"cases": total_cases, "candidate_equalities": total_cases*4,
                          "seconds": time.perf_counter()-t0, "status": "passed"}
    out["panels"] = []
    for m in (31, 83, 131):
        setup = time.perf_counter()
        c = make_curve(m, 0)
        edge_p = fixture(c, random.Random(SEEDS[0]+m))
        edges = [0, 1, 2, 3, 2**(m-1), 2**m-1, 2**m, 2**m+1]
        edge_cases = 0
        for p in (None, (c["f"](0), c["f"](1)), edge_p):
            identities(c, p)
            for k in sorted(set(edges + [-v for v in edges])):
                check_case(c, p, k)
                edge_cases += 1
        for seed in SEEDS:
            rng = random.Random(seed+m)
            cases = [(fixture(c, rng), rng.randrange(1, c["order"]//4)) for _ in range(24)]
            assert len({tuple(point_json(p)) for p, _ in cases}) == len(cases)
            expected = []
            for p, k in cases:
                check_case(c, p, k)
                identities(c, p)
                expected.append(sage_result(c, p, k))
            frozen = [{"point": point_json(p), "scalar": str(k)} for p, k in cases]
            panel = {"m": m, "a": 0, "seed": seed, "holdout": seed == SEEDS[1],
                     "curve_id": f"binary-koblitz:a=0;m={m};basis=polynomial;taps={TAPS[m]};cofactor_clear=4",
                     "field_modulus": c["field_modulus"], "group_order": str(c["order"]),
                     "quarter_order_is_prime": bool(Integer(c["order"]//4).is_prime()),
                     "edge_cases_verified": edge_cases, "cases": frozen,
                     "input_sha256": sha(frozen), "output_sha256": sha([point_json(p) for p in expected]),
                     "shared_fixture_setup_and_validation_seconds": time.perf_counter()-setup,
                     "aa": [], "rounds": [], "summary": {}}
            out["panels"].append(panel)
            for method in METHODS:
                timing(c, cases, method, expected)
            for _ in range(5):
                a = timing(c, cases, "binary_naf", expected)
                b = timing(c, cases, "binary_naf", expected)
                panel["aa"].append({"a": a, "b": b, "ratio_a_over_b": a["ns"]/b["ns"]})
            for rep in range(7):
                order = METHODS if rep%2 == 0 else tuple(reversed(METHODS))
                panel["rounds"].append({"index": rep, "order": order,
                    "measurements": {method: timing(c, cases, method, expected) for method in order}})
            for method in METHODS:
                samples = [r["measurements"][method]["ns_per_scalar"] for r in panel["rounds"]]
                ratios = [r["measurements"][method]["ns"]/r["measurements"]["binary_naf"]["ns"] for r in panel["rounds"]]
                counts, digits = [], []
                for (p, k), want in zip(cases, expected):
                    got, count, ds = multiply(c, p, k, method, counted=True)
                    assert got == want
                    if count is not None:
                        counts.append(count)
                        digits.append(ds)
                panel["summary"][method] = {
                    "median_ns": statistics.median(samples), "min_ns": min(samples),
                    "median_paired_cost_ratio_to_binary_naf": statistics.median(ratios),
                    "counted_per_case": counts, "digits_per_case": digits,
                    "mean_counts": {key: statistics.mean(r[key] for r in counts) for key in counts[0]} if counts else None,
                    "mean_digits": {key: statistics.mean(r[key] for r in digits) for key in digits[0]} if digits else None,
                }
            noise = max(abs(r["ratio_a_over_b"]-1) for r in panel["aa"])
            ratio = panel["summary"]["reduced_tau_naf"]["median_paired_cost_ratio_to_binary_naf"]
            panel["aa_max_deviation"] = noise
            panel["arithmetic_screen_passed"] = ratio <= 0.9 and 1-ratio > noise
            print(json.dumps({"m": m, "seed": seed, "summary": {k: {a:b for a,b in v.items() if a not in ("counted_per_case", "digits_per_case")} for k,v in panel["summary"].items()}, "aa": noise}), flush=True)
            setup = time.perf_counter()
    out["end_to_end_speedup"] = None
    out["gpu_throughput"] = None
    out["whole_walk_cost"] = None
    out["ecdlp_operations"] = None
    out["eligible_runtime_fraction"] = None
    out["status"] = "passed"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    # Reserve the path before doing any work; never overwrite an earlier run.
    with args.output.open("x") as f:
        f.write("{}\n")
    out = {"schema": 1, "status": "started", "utc_started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "classification": "known-scalar arithmetic stage diagnostic",
           "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           "protocol_sha256": hashlib.sha256(Path(__file__).with_name("PROTOCOL.md").read_bytes()).hexdigest(),
           "threads": 1, "seeds": SEEDS, "cases_per_panel": 24, "repetitions": 7}
    code = 0
    try:
        out["host"] = host()
        out["sage_version"] = sage_version
        run(out)
    except Exception:
        out["status"] = "failed"
        out["error"] = traceback.format_exc()
        print(out["error"], file=sys.stderr)
        code = 1
    finally:
        out["load_after"] = os.getloadavg()
        out["processes_after"] = command(["ps", "-eo", "comm,pcpu,pmem"])
        out["utc_finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        args.output.write_text(json.dumps(out, indent=2)+"\n")
    return code


if __name__ == "__main__":
    from sage.all import GF, EllipticCurve, Integer, PolynomialRing
    from sage.env import SAGE_VERSION as sage_version
    sys.exit(main())
