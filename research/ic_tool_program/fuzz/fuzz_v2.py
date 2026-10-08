#!/usr/bin/env python3
"""B6's generator (rounds/B6-fuzzing/PROTOCOL.md): seeded schema v2
documents, valid and corrupted, across every field kind, curve form,
target form and method, run through `ic check` and `ic price` with a
wall budget of 1-3 s.

A quarter of the documents are random, with corruptions. The rest are
the conformance cases' parameter files mutated: lightly (the target or
the method only, so that most still run) or anywhere.

Each run must:
- exit with a documented status: 0-4, never a panic's 101 or a signal;
- write a JSON report when it exits 0, 2, 3 or 4 (exit 1 is a budget
  stop or a failed run);
- never print "panicked" on stderr;
- for a `complete` report, replay every certificate as [d]G = Q in this
  file's own Python arithmetic, which shares nothing with the Rust, and
  match a known logarithm;
- for a quarter of the complete reports, give the same scalar when the
  report's `resolved` document is run again (design §6).

    python3 fuzz_v2.py <ic binary> <seed> <count> <out dir>
"""
from __future__ import annotations

import json
import random
import subprocess
import sys
import tempfile
from pathlib import Path

IC, SEED, COUNT, OUT = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), Path(sys.argv[4])
rng = random.Random(SEED)

PRIMES = [5, 7, 11, 13, 1009, 10007, 65537, 1000003, 16777259, 2147483647, (1 << 61) - 1,
          (1 << 89) - 1, (1 << 127) - 1, 2**255 - 19, 2**256 - 2**32 - 977]
COMPOSITES = [4, 9, 15, 561, 16777217, 2**64 + 1, 3215031751]
BINARY = [(5, 0b100101), (7, 0b10000011), (13, (1 << 13) | 0b11011), (17, (1 << 17) | 0b1001),
          (23, (1 << 23) | 0b100001), (31, (1 << 31) | 0b1001), (32, 0x10000008D), (37, (1 << 37) | 0b111111),
          (61, (1 << 61) | 0b100111), (83, 0x800000000200000000007), (100, 0x10000000000000000000000065),
          (127, (1 << 127) | 0b11), (163, (1 << 163) | 0b11001001), (2, 0b111), (3, 0b1011)]


def num(v: int) -> str:
    return rng.choice([str(v), hex(v)]) if v >= 0 else str(v)


def element(bound: int) -> int:
    return rng.choice([0, 1, 2, bound - 1, bound, bound + 1, rng.randrange(max(bound, 1)),
                       rng.randrange(max(bound, 1)), rng.randrange(max(bound, 1))])


def junk():
    return rng.choice([None, True, -1, 1.5, "", "0x", "+5", "five", [], {}, [1, 2], {"a": 1}, "0" * 500,
                       str(2**1100), "0x" + "f" * 300])


def field() -> tuple[dict, str, int]:
    kind = rng.choice(["binary", "binary", "prime", "prime", "prime_extension"])
    if kind == "binary":
        n, f = rng.choice(BINARY)
        doc = {"kind": "binary", "degree": n}
        if rng.random() < 0.9:
            doc["modulus"] = num(f if rng.random() < 0.9 else f ^ (1 << rng.randrange(n)))
        return doc, kind, 1 << n
    if kind == "prime":
        p = rng.choice(PRIMES + COMPOSITES[:2])
        return {"kind": "prime", "p": num(p)}, kind, p
    p = rng.choice([3, 5, 7, 1009, 10007, 2147483647])
    k = rng.choice([2, 3, 4, 1, 65])
    mod = [num(element(p)) for _ in range(max(k, 0))]
    return {"kind": "prime_extension", "p": num(p), "degree": k, "modulus": mod}, kind, p


def value(kind: str, bound: int, k: int):
    if kind == "prime_extension":
        return [num(element(bound)) for _ in range(k)]
    return num(element(bound))


def curve(kind: str, bound: int, k: int) -> dict:
    forms = {"binary": ["koblitz", "binary_weierstrass", "general_weierstrass"],
             "prime": ["short_weierstrass", "montgomery", "twisted_edwards", "general_weierstrass"],
             "prime_extension": ["short_weierstrass", "general_weierstrass"]}[kind]
    if rng.random() < 0.05:
        forms = ["koblitz", "binary_weierstrass", "general_weierstrass", "short_weierstrass",
                 "montgomery", "twisted_edwards", "hessian"]
    form = rng.choice(forms)
    keys = {"koblitz": ["a"], "binary_weierstrass": ["a", "b"], "short_weierstrass": ["a", "b"],
            "general_weierstrass": ["a1", "a2", "a3", "a4", "a6"], "montgomery": ["A", "B"],
            "twisted_edwards": ["a", "d"], "hessian": ["a"]}[form]
    doc = {"form": form}
    for key in keys:
        doc[key] = rng.choice([0, 1]) if form == "koblitz" else value(kind, bound, k)
    return doc


def point(kind: str, bound: int, k: int):
    return {"x": value(kind, bound, k), "y": value(kind, bound, k)}


def document() -> dict:
    if rng.random() < 0.08:
        doc = {"schema_version": 2, "name": "fuzz", "named": rng.choice(["secp256k1", "sect163k1", "ecc2k-95",
                                                                         "ecc2k-130", "nope"])}
        kind, bound, k = "prime", 1, 1
    else:
        fdoc, kind, bound = field()
        k = fdoc.get("degree", 1) if kind == "prime_extension" else 1
        k = k if isinstance(k, int) and 0 < k < 70 else 1
        doc = {"schema_version": 2, "name": "fuzz", "field": fdoc, "curve": curve(kind, bound, k)}
        r = rng.choice([2, 3, 7, 1439393, 5591479, bound, bound // 3 + 1, 2**81 + 1])
        sub = {"generator": point(kind, bound, k) if rng.random() < 0.9 else {"rule": "koblitz_search_v1"}}
        if rng.random() < 0.8:
            sub["order"] = num(r)
            sub["cofactor"] = num(rng.choice([1, 2, 3, 4, 1492, bound // max(r, 1) + 1]))
        doc["subgroup"] = sub
    tform = rng.choice(["point", "identity", "known_log", "random_seed", "public_hash_seed"])
    doc["target"] = {"point": point(kind, bound, k)} if tform == "point" else \
        {"point": "identity"} if tform == "identity" else \
        {tform: num(rng.randrange(1, 10**6)) if tform == "known_log" else rng.randrange(10**6)}
    method = {"solve": rng.choice(["paired", "rho", "index_calculus", "check"])}
    if rng.random() < 0.5:
        method["fidelity"] = rng.choice(["auto", "F0", "F2", "F1"])
    if rng.random() < 0.5:
        method["rho"] = {"pipeline": rng.choice(["auto", "rho-koblitz", "rho-negation", "rho-bignum", "trivial"]),
                         "seed": rng.randrange(2**32)}
    if rng.random() < 0.4:
        method["index_calculus"] = {"pipeline": rng.choice(["auto", "kic", "ic-binary-s4", "ic-prime-s3"]),
                                    "recipe": "auto"}
    doc["method"] = method
    doc["budget"] = {"wall_seconds": rng.choice([1, 2, 3])}
    # Corruptions.
    for _ in range(rng.choice([0, 0, 0, 1, 2])):
        path = rng.choice(["field", "curve", "subgroup", "target", "method", "budget", ""])
        node = doc.get(path, doc) if path else doc
        if isinstance(node, dict) and node:
            key = rng.choice(list(node) + ["extra"])
            node[key] = junk()
    return doc


# The conformance cases' parameter files, the mutations' seeds.
PROGRAMME = Path(__file__).resolve().parents[1]
SEEDS = []
for d in ("v2/params", "v2-b2/params", "v2-b3/params"):
    for f in sorted((PROGRAMME / "conformance" / d).glob("*.json")):
        SEEDS.append(json.loads(f.read_text()))


def bump(v):
    """A neighbouring value of a written integer or array."""
    if isinstance(v, list):
        return [bump(x) if rng.random() < 0.5 else x for x in v]
    if isinstance(v, str):
        try:
            n = int(v, 16) if v.startswith("0x") else int(v)
        except ValueError:
            return v
        return num(max(0, n + rng.choice([-1, 1, 2, -2, n, -n // 2])))
    if isinstance(v, int) and not isinstance(v, bool):
        return v + rng.choice([-1, 1, v])
    return v


def mutate(doc: dict, light: bool) -> dict:
    """A conformance document changed: lightly (its target or method,
    so that most still run) or anywhere."""
    doc = json.loads(json.dumps(doc))
    doc["budget"] = {"wall_seconds": rng.choice([1, 2, 3])}
    for _ in range(1 if light else rng.choice([1, 1, 2, 3])):
        what = rng.choice(["method", "target"] if light else
                          ["bump", "bump", "method", "target", "drop", "junk", "form"])
        if what == "bump":
            part = rng.choice([k for k in ("field", "curve", "subgroup", "target") if isinstance(doc.get(k), dict)] or ["name"])
            node = doc.get(part)
            if isinstance(node, dict) and node:
                key = rng.choice(list(node))
                if isinstance(node[key], dict) and node[key]:
                    sub = rng.choice(list(node[key]))
                    node[key][sub] = bump(node[key][sub])
                else:
                    node[key] = bump(node[key])
        elif what == "method":
            m = doc.setdefault("method", {})
            m["solve"] = rng.choice(["paired", "rho", "index_calculus", "check"])
            if rng.random() < 0.5:
                m["fidelity"] = rng.choice(["auto", "F0", "F2", "F1"])
            if rng.random() < 0.3:
                m["rho"] = {"pipeline": rng.choice(["auto", "rho-koblitz", "rho-negation", "rho-bignum"]),
                            "seed": rng.randrange(2**20)}
            if rng.random() < 0.3:
                # An earlier "junk" pass may have left a non-object here
                # (B6's amendment 1): replace it, drawing nothing more.
                if not isinstance(m.get("index_calculus"), dict):
                    m["index_calculus"] = {}
                m["index_calculus"]["pipeline"] = rng.choice(["auto", "kic", "ic-binary-s4", "ic-prime-s3"])
        elif what == "target":
            doc["target"] = rng.choice([{"known_log": num(rng.randrange(1, 10**7))}, {"random_seed": rng.randrange(1000)},
                                        {"public_hash_seed": rng.randrange(1000)}, {"point": "identity"}])
        elif what == "drop":
            part = rng.choice(list(doc))
            if part not in ("schema_version",):
                doc.pop(part)
        elif what == "junk":
            part = rng.choice(list(doc))
            node = doc[part]
            if isinstance(node, dict) and node:
                node[rng.choice(list(node))] = junk()
        elif what == "form" and isinstance(doc.get("curve"), dict):
            doc["curve"]["form"] = rng.choice(["koblitz", "binary_weierstrass", "general_weierstrass",
                                               "short_weierstrass", "montgomery", "twisted_edwards"])
    return doc


# ── An independent replay: [d]G = Q in Python arithmetic ─────────────


def bmul(a: int, b: int, f: int, n: int) -> int:
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a >> n & 1:
            a ^= f
    return out


def binv(a: int, f: int, n: int) -> int:
    # a^(2^n - 2) by square and multiply.
    out, base, e = 1, a, (1 << n) - 2
    while e:
        if e & 1:
            out = bmul(out, base, f, n)
        base = bmul(base, base, f, n)
        e >>= 1
    return out


def badd(P, Q, a, f, n):
    if P is None:
        return Q
    if Q is None:
        return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2:
        if y1 ^ y2 == x1:
            return None
        if x1 == 0:
            return None
        lam = x1 ^ bmul(y1, binv(x1, f, n), f, n)
        x3 = bmul(lam, lam, f, n) ^ lam ^ a
        return x3, bmul(x1, x1, f, n) ^ bmul(lam ^ 1, x3, f, n)
    lam = bmul(y1 ^ y2, binv(x1 ^ x2, f, n), f, n)
    x3 = bmul(lam, lam, f, n) ^ lam ^ x1 ^ x2 ^ a
    return x3, bmul(lam, x1 ^ x3, f, n) ^ x3 ^ y1


def padd(P, Q, a, p):
    if P is None:
        return Q
    if Q is None:
        return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None
        lam = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return x3, (lam * (x1 - x3) - y1) % p


def times(k: int, P, add):
    out = None
    while k:
        if k & 1:
            out = add(out, P)
        P = add(P, P)
        k >>= 1
    return out


def as_int(v) -> int:
    return int(v, 16) if isinstance(v, str) and v.startswith("0x") else int(v)


def replay(rep: dict) -> str | None:
    """None when every certificate replays, else what failed."""
    certs = rep.get("certificates") or {}
    for arm in ("ic", "rho"):
        c = certs.get(arm)
        if not c:
            continue
        d = as_int(c["scalar"])
        if c.get("schema") == "ic-v2-replay-v1":
            field, curve, sub = c["curve"]["field"], c["curve"]["curve"], c["curve"]["subgroup"]
            G = (as_int(sub["generator"]["x"]), as_int(sub["generator"]["y"]))
            Q = (as_int(c["target"]["x"]), as_int(c["target"]["y"]))
            if field["kind"] == "binary":
                n, f = field["degree"], as_int(field["modulus"])
                a = as_int(curve["a"]) if curve["form"] != "koblitz" else curve["a"]
                got = times(d, G, lambda P, R: badd(P, R, a, f, n))
            else:
                p = as_int(field["p"])
                got = times(d, G, lambda P, R: padd(P, R, as_int(curve["a"]), p))
        else:
            cv = c["curve"]
            n = cv["degree"]
            f = (1 << n) | sum(1 << t for t in cv["irreducible_low_terms"])
            G = tuple(int(v) for v in cv["generator"])
            Q = tuple(int(v) for v in c["target"])
            got = times(d, G, lambda P, R: badd(P, R, cv["curve_a"], f, n))
        if got != Q:
            return f"{arm}: [d]G != Q"
    known = (rep.get("resolved") or {}).get("target", {}).get("known_log")
    scalar = (rep.get("result") or {}).get("scalar")
    if known is not None and scalar is not None and as_int(known) % as_int(rep["resolved"]["subgroup"]["order"]) != as_int(scalar):
        return "the scalar is not the known logarithm"
    return None


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    bad = []
    outcomes: dict = {}
    with tempfile.TemporaryDirectory() as tmp:
        for i in range(COUNT):
            roll = rng.random()
            doc = (mutate(rng.choice(SEEDS), light=roll < 0.35) if SEEDS and roll < 0.75 else document())
            path = Path(tmp) / f"d{i}.json"
            path.write_text(json.dumps(doc))
            report = Path(tmp) / f"r{i}.json"
            command = rng.choice(["price", "price", "check"])
            argv = [IC, command, "--params", str(path), "--json", "--out", str(report)]
            if command == "price":
                argv += ["--repeats", "1", "--repeats-fast", "1"]
            try:
                p = subprocess.run(argv, capture_output=True, text=True, timeout=20,
                                   env={"RAYON_NUM_THREADS": "1", "PATH": "/usr/bin:/bin"})
                code, err = p.returncode, p.stderr
            except subprocess.TimeoutExpired:
                code, err = "timeout", ""
            ok = code in (0, 1, 2, 3, 4) and "panicked" not in err
            try:
                rep = json.loads(report.read_text()) if report.exists() else {}
            except json.JSONDecodeError:
                rep = {}
            key = f"{code}:{rep.get('status')}:{(rep.get('refusal') or {}).get('code', '')}"
            outcomes[key] = outcomes.get(key, 0) + 1
            if ok and code != 1 and not report.exists():
                ok = False
            why = None
            if ok and rep.get("status") == "complete":
                why = replay(rep)
                ok = why is None
                outcomes["replayed"] = outcomes.get("replayed", 0) + 1
                if ok and rep.get("resolved") and rng.random() < 0.25:
                    again = Path(tmp) / f"again{i}.json"
                    again.write_text(json.dumps(rep["resolved"]))
                    rep2_path = Path(tmp) / f"r2{i}.json"
                    subprocess.run([IC, "price", "--params", str(again), "--json", "--out", str(rep2_path),
                                    "--repeats", "1", "--repeats-fast", "1"], capture_output=True, timeout=30,
                                   env={"RAYON_NUM_THREADS": "1", "PATH": "/usr/bin:/bin"})
                    try:
                        rep2 = json.loads(rep2_path.read_text())
                    except (OSError, json.JSONDecodeError):
                        rep2 = {}
                    outcomes["resolved_rerun"] = outcomes.get("resolved_rerun", 0) + 1
                    if rep2.get("status") == "complete" and rep2["result"]["scalar"] != rep["result"]["scalar"]:
                        ok, why = False, "resolved gave another scalar"
            if not ok:
                case = {"i": i, "command": command, "exit": code, "why": why, "stderr": err[-800:], "doc": doc}
                bad.append(case)
                (OUT / f"bad-{SEED}-{i}.json").write_text(json.dumps(case, indent=1))
    print(json.dumps(dict(sorted(outcomes.items(), key=lambda kv: -kv[1])), indent=0))
    print(json.dumps({"seed": SEED, "count": COUNT, "bad": len(bad),
                      "first": [(b["i"], b["exit"], b["stderr"][-200:]) for b in bad[:5]]}, indent=1))


if __name__ == "__main__":
    main()
