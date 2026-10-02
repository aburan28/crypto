#!/usr/bin/env python3
"""Independent arithmetic replay of the KS v2 rho records of one n=61 block.

Pure-Python GF(2^61) and Koblitz arithmetic from ``independent_replay.py``;
shares no code with either Rust producer. For every rho fixture record it
checks, against the base header, the eval scalar corpus and the IC records of
the same block:

  * the recovered scalar equals the published scalar and lies in [1, r),
  * the published scalar equals line ``fixture_index`` of the corpus,
  * published Q is on the curve and equals [recovered scalar] G,
  * published Q equals the IC arm's target for the same fixture index.

Usage: rho_replay.py BASE.jsonl SCALARS.txt IC.jsonl KS.jsonl REPORT.json
"""

import hashlib
import json
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from independent_replay import Curve, Field  # noqa: E402

CURVE = None
GEN = None
ORDER = None


def init(low_terms, a, generator, order):
    global CURVE, GEN, ORDER
    CURVE = Curve(Field(61, low_terms), a)
    GEN = tuple(generator)
    ORDER = order


def check(job):
    rec, corpus_scalar, ic_target = job
    d = rec["recovered_fixture_scalar"]
    q = tuple(rec["published_q"])
    checks = {
        "recovered_equals_published": d == rec["published_fixture_scalar"],
        "recovered_in_range": 0 < d < ORDER,
        "published_equals_corpus": rec["published_fixture_scalar"] == corpus_scalar,
        "q_on_curve": CURVE.on_curve(q),
        "recovered_times_g_is_q": CURVE.mul(d, GEN) == q,
        "q_equals_ic_target": ic_target is not None and q == tuple(ic_target),
        "producer_verified": rec.get("verified") is True,
    }
    return rec["fixture_index"], [k for k, v in checks.items() if not v]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    base, scalars, ic, ks, report = sys.argv[1:6]
    header = json.loads(Path(base).open("rb").readline())
    corpus = [int(line) for line in Path(scalars).read_text().split()]
    ic_targets, generator = {}, None
    for line in Path(ic).read_text().splitlines():
        rec = json.loads(line)
        if rec.get("kind") == "compact_orbit_dlp_target":
            ic_targets[rec["fixture_index"]] = rec["target"]
            generator = generator or rec["generator"]
    jobs, summary_lines = [], 0
    for line in Path(ks).read_text().splitlines():
        rec = json.loads(line)
        if rec.get("kind") == "rho_ks_batch_summary":
            summary_lines += 1
            continue
        assert rec["kind"] == "rho_ks_batch_fixture" and rec["n"] == 61
        i = rec["fixture_index"]
        jobs.append((rec, corpus[i], ic_targets.get(i)))
    init_args = (header["field_modulus_low_terms"], header["a"], generator,
                 header["subgroup_order"])
    init(*init_args)
    assert CURVE.on_curve(GEN) and CURVE.mul(ORDER, GEN) is None
    with Pool(initializer=init, initargs=init_args) as pool:
        results = pool.map(check, jobs, chunksize=256)
    failures = [{"fixture_index": i, "failed": f} for i, f in results if f]
    indices = sorted(i for i, _ in results)
    out = {
        "schema": "n61-L65536-ks-v2-rho-independent-replay-v1",
        "inputs": {"base_sha256": sha(base), "scalars_sha256": sha(scalars),
                   "ic_records_sha256": sha(ic), "rho_records_sha256": sha(ks)},
        "records": len(results),
        "pass": len(results) - len(failures),
        "fail": len(failures),
        "fixture_indices_complete": indices == list(range(len(corpus))),
        "summary_lines": summary_lines,
        "failures": failures[:100],
        "all_pass": not failures and indices == list(range(len(corpus))),
    }
    Path(report).write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: out[k] for k in ("records", "pass", "fail", "all_pass")}))
    return 0 if out["all_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
