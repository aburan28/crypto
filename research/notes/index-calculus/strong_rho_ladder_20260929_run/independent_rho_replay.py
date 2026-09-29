#!/usr/bin/env python3
"""Independent check of the strong-rho records: pure-Python GF(2^n), no Rust.

For each rho fixture record it checks, with the arithmetic of
`autolab_orbit_extract_20260924/independent_replay.py` (which shares no code with
any Rust producer):

  * the generator lies on the curve and has order r (once per run),
  * the published target Q lies on the curve,
  * [published scalar] G == Q,
  * the recovered scalar equals the published scalar,
  * for a sample of fixtures, [r] Q == infinity (Q lies in the prime-order
    subgroup).

usage: independent_rho_replay.py <rho.jsonl> <ic_records.jsonl> <n> <a> <order>
                                 [low_terms comma list | 'search']
`ic_records.jsonl` supplies the generator (its first record's "generator").
With 'search', the reduction polynomial is found as the low-term set under
which the generator lies on the curve (used when no header is retained).
"""
import gzip
import importlib.util
import itertools
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPLAY = os.path.join(
    HERE,
    "..", "..", "..",
    "sat_factor_base_review_20260908", "autolab_orbit_extract_20260924",
    "independent_replay.py",
)
spec = importlib.util.spec_from_file_location("independent_replay", REPLAY)
ir = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ir)


def _open(path):
    if not os.path.exists(path) and os.path.exists(path + ".gz"):
        return gzip.open(path + ".gz", "rt")
    return open(path)


def find_low_terms(n, a, gen):
    """Low-term set (always containing 0) making `gen` a point of the curve."""
    for extra in (1, 2, 3):
        for terms in itertools.combinations(range(1, n), extra):
            field = ir.Field(n, (0,) + terms)
            if ir.Curve(field, a).on_curve(gen):
                return (0,) + terms
    raise SystemExit("no reduction polynomial with <=4 low terms puts the generator on the curve")


def main():
    rho_path, ic_path, n, a, order = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    terms_arg = sys.argv[6] if len(sys.argv) > 6 else "search"
    with _open(ic_path) as fh:
        gen = tuple(json.loads(fh.readline())["generator"])
    if terms_arg == "search":
        low_terms = find_low_terms(n, a, gen)
    else:
        low_terms = tuple(int(t) for t in terms_arg.split(","))
    curve = ir.Curve(ir.Field(n, low_terms), a)

    report = {
        "rho_records": rho_path,
        "n": n,
        "a": a,
        "field_modulus_low_terms": list(low_terms),
        "subgroup_order": order,
        "generator": list(gen),
        "generator_on_curve": curve.on_curve(gen),
        "generator_has_order_r": curve.mul(order, gen) is None,
        "fixtures": 0,
        "pass": 0,
        "fail": 0,
        "order_r_sampled": 0,
        "failures": [],
    }
    with _open(rho_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            if rec.get("kind") != "rho_ks_batch_fixture":
                continue
            report["fixtures"] += 1
            q = tuple(rec["published_q"])
            d = rec["published_fixture_scalar"]
            checks = {
                "q_on_curve": curve.on_curve(q),
                "q_equals_d_times_g": curve.mul(d, gen) == q,
                "recovered_equals_published": rec["recovered_fixture_scalar"] == d,
            }
            if rec["fixture_index"] % 16 == 0:
                checks["q_has_order_r"] = curve.mul(order, q) is None
                report["order_r_sampled"] += 1
            if all(checks.values()):
                report["pass"] += 1
            else:
                report["fail"] += 1
                report["failures"].append({"fixture_index": rec["fixture_index"], "checks": checks})
    report["all_pass"] = (
        report["fail"] == 0
        and report["generator_on_curve"]
        and report["generator_has_order_r"]
        and report["fixtures"] > 0
    )
    print(json.dumps(report))
    sys.exit(0 if report["all_pass"] else 1)


if __name__ == "__main__":
    main()
