"""ICMS command line.  Run from the repository root:

    python3 tools/icms validate docs/ic/measurement/specs/*.yaml
    python3 tools/icms capsule
    python3 tools/icms plan A.yaml B.yaml
    python3 tools/icms run A.yaml B.yaml A.yaml --cpus 3 --out docs/ic/measurement/sessions/<name>
    python3 tools/icms compare docs/ic/measurement/sessions/<name> --a 0 --b 1 --declare factor_base.params.dimension
    python3 tools/icms summarize docs/ic/measurement/sessions/<name>
    python3 tools/icms audit-sessions docs/ic/measurement/sessions

Listing a spec twice in ``run`` adds an A/A arm, which ``compare`` reports as
the session's noise floor beside every A/B ratio.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "icms"  # noqa: A001

from icms import adapters  # noqa: E402
from icms.environment import capsule as make_capsule, parse_cpu_list  # noqa: E402
from icms.registry import Registry  # noqa: E402
from icms.spec import SpecError, load as load_spec  # noqa: E402


def cmd_validate(args) -> int:
    reg = Registry.load()
    bad = 0
    for p in args.specs:
        try:
            sp = load_spec(p, reg)
        except SpecError as exc:
            bad += 1
            print(f"FAIL {p}")
            for prob in exc.problems:
                print(f"     {prob}")
            continue
        probs = adapters.get(sp["spec"]["execution"]["adapter"]).check(sp["spec"])
        if probs:
            bad += 1
            print(f"FAIL {p} ({sp['spec_id']}): adapter {sp['spec']['execution']['adapter']} refuses")
            for prob in probs:
                print(f"     {prob}")
            continue
        print(f"ok   {p}  {sp['spec_id']}  workloads {[w['workload_id'] for w in sp['workloads']]}")
    return 1 if bad else 0


def cmd_capsule(args) -> int:
    cap = make_capsule()
    if args.out:
        if os.path.exists(args.out):
            print(f"{args.out} exists; refusing to overwrite", file=sys.stderr)
            return 2
        with open(args.out, "w") as fh:
            json.dump(cap, fh, indent=1, sort_keys=True)
    st = cap["stable"]
    print(f"env_class_id  {cap['env_class_id']}")
    print(f"cpu           {st['cpu']['model_name']} (family {st['cpu']['family']} model {st['cpu']['model']} "
          f"stepping {st['cpu']['stepping']}, microcode {st['cpu']['microcode']}, {st['cpu']['logical_cpus']} logical)")
    print(f"kernel        {st['kernel']['release']}  clocksource {st['kernel']['clocksource']}")
    print(f"isolation     isolcpus={st['topology']['isolated']!r} nohz_full={st['topology']['nohz_full']!r} "
          f"cmdline flags {st['kernel']['cmdline_flags']}")
    govs = sorted({str(c['governor']) for c in st['topology']['cpus'].values()})
    print(f"frequency     governors {govs} no_turbo={st['topology']['intel_pstate_no_turbo']} "
          f"boost={st['topology']['cpufreq_boost']}")
    print(f"smt           control={st['topology']['smt_control']} active={st['topology']['smt_active']}")
    print(f"virt          {st['virtualization']['detect_virt']}  thp {st['kernel']['thp_enabled']}")
    print(f"psi           {'available' if st['kernel']['psi_available'] else 'absent'}")
    print(f"toolchain     python {st['toolchain']['python']} | {st['toolchain']['rustc']} | {st['toolchain']['gcc']}")
    return 0


def cmd_plan(args) -> int:
    from icms.session import plan
    reg = Registry.load()
    arms = []
    for p in args.specs:
        sp = load_spec(p, reg)
        for wl in sp["workloads"]:
            arms.append({"spec": sp["spec"], "path": p, "spec_id": sp["spec_id"], "wl": wl["workload_id"]})
    for i, a in enumerate(arms):
        print(f"arm {i}: {a['path']} {a['spec_id']} {a['wl']}")
    order = plan(arms, 0)
    print(" ".join(f"{i}{'w' if w else ''}" for i, r, w in order))
    return 0


def cmd_run(args) -> int:
    from icms.session import SessionError, auto_cpus, run_session
    try:
        cpus = auto_cpus() if args.cpus == "auto" else set(parse_cpu_list(args.cpus))
        run_session(args.specs, args.out, cpus, lock=args.lock, wait=args.wait,
                    settle=args.settle, max_other_cpu=args.max_other_cpu, max_psi=args.max_psi,
                    allow_busy=args.allow_busy, session_label=args.label)
    except SessionError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    return 0


def cmd_compare(args) -> int:
    from icms.compare import compare
    res = compare(args.session, args.a, args.b, args.declare or [])
    text = json.dumps(res, indent=1, sort_keys=True)
    if args.out:
        if os.path.exists(args.out):
            print(f"{args.out} exists; refusing to overwrite", file=sys.stderr)
            return 2
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w") as fh:
            fh.write(text + "\n")
    ops, wall = res["ops"], res["wall"]
    print(f"A arm {args.a} {res['a']['label']}\nB arm {args.b} {res['b']['label']}")
    for r in res["refusals"]:
        print(f"REFUSED: {r['reason']}: {r.get('fields') or r.get('record_ids') or ''}")
    if ops.get("admitted"):
        print(f"ops  {ops['unit']} over {ops['window']}: A {ops['a_median']:.6g}  B {ops['b_median']:.6g}  "
              f"B/A {ops['ratio_b_over_a']:.4g}{'' if ops['deterministic'] else '  [host-dependent]'}")
    elif "unit" in ops:
        print(f"ops  not admitted: {ops.get('reason')}")
    est = wall.get("estimate") or {}
    if wall.get("admitted"):
        print(f"wall B/A median {est['median_ratio']:.4g} 95% CI {est['ci95']} over {est['n_pairs']} pairs")
    else:
        for r in wall.get("refusals", []):
            print(f"wall not admitted: {r['reason']}")
        if est.get("median_ratio") is not None:
            print(f"     (descriptive only: median B/A {est['median_ratio']:.4g} over {est['n_pairs']} pairs)")
    for aa in res.get("aa_noise", []):
        e = aa["estimate"]
        print(f"A/A noise arms {aa['arms']}: median {e['median_ratio']} CI {e.get('ci95')}")
    # A refused comparison is a result (and is frozen when --out is given), but
    # it is not the comparison that was asked for, so the exit status says so.
    return 3 if res["refusals"] else 0


def cmd_audit_sessions(args) -> int:
    from icms.audit import audit
    n, problems = audit(args.paths)
    for p in problems:
        print(f"FAIL {p}")
    print(f"{n} session(s) audited, {len(problems)} problem(s)")
    return 1 if problems else 0


def cmd_summarize(args) -> int:
    from icms.compare import load_session
    session, records, _ = load_session(args.session)
    print(f"session {session['session_id']}  env {session['env_class_id']}  preflight quiet={session['preflight']['quiet']}")
    for a in session["arms"]:
        rs = [r for r in records if r["arm"] == a["index"] and not r["warmup"]]
        levels = sorted({r["isolation"]["earned_level"] for r in rs})
        walls = sorted(r["execution"]["wall_ns"] / 1e6 for r in rs)
        ok = sum(1 for r in rs if r["outcome"]["status"] == "complete")
        print(f"arm {a['index']}: {a['label']}  {a['spec_id']}  runs {len(rs)} complete {ok}  levels {levels}  "
              f"process wall ms min {walls[0]:.1f} median {walls[len(walls) // 2]:.1f}" if walls else f"arm {a['index']}: no runs")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="icms", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("validate"); p.add_argument("specs", nargs="+"); p.set_defaults(fn=cmd_validate)
    p = sub.add_parser("capsule"); p.add_argument("--out"); p.set_defaults(fn=cmd_capsule)
    p = sub.add_parser("plan"); p.add_argument("specs", nargs="+"); p.set_defaults(fn=cmd_plan)
    p = sub.add_parser("run")
    p.add_argument("specs", nargs="+")
    p.add_argument("--cpus", required=True, help="CPUs to reserve, e.g. 3 or 2-3, or 'auto' for the last whole core; CPU 0 never earns L1")
    p.add_argument("--out", required=True, help="new directory for the session")
    p.add_argument("--lock", default="/tmp/crypto-bench.lock")
    p.add_argument("--wait", action="store_true")
    p.add_argument("--settle", type=float, default=2.0)
    p.add_argument("--max-other-cpu", type=float, default=0.10)
    p.add_argument("--max-psi", type=float, default=5.0)
    p.add_argument("--allow-busy", action="store_true",
                   help="record a refused preflight and continue; every run then stays below L2")
    p.add_argument("--label", default="")
    p.set_defaults(fn=cmd_run)
    p = sub.add_parser("compare")
    p.add_argument("session"); p.add_argument("--a", type=int, required=True); p.add_argument("--b", type=int, required=True)
    p.add_argument("--declare", action="append", help="a spec field the comparison is about (repeatable)")
    p.add_argument("--out")
    p.set_defaults(fn=cmd_compare)
    p = sub.add_parser("summarize"); p.add_argument("session"); p.set_defaults(fn=cmd_summarize)
    p = sub.add_parser("audit-sessions", help="re-derive every hash, id, gate and comparison of committed sessions")
    p.add_argument("paths", nargs="+", help="session directories, or directories of sessions (a missing one holds none)")
    p.set_defaults(fn=cmd_audit_sessions)
    args = ap.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
