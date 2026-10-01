#!/usr/bin/env python3
"""Build the IC measurement-standard page from frozen files.

    python3 scripts/build_ic_measurement.py           # write docs/ic-measurement.html and its JSON
    python3 scripts/build_ic_measurement.py --check   # fail if either output is stale

The page cites and never computes (AGENTS.md section 7).  Every number on it
is read from a frozen file:

- docs/ic/measurement/registry.json            the vocabularies
- docs/ic/measurement/audit.json               the comparability findings
- docs/ic/measurement/sessions/*/session.json  each frozen session
- docs/ic/measurement/sessions/*/capsule.json  its host
- docs/ic/measurement/sessions/*/records.jsonl its runs
- docs/ic/measurement/sessions/*/comparisons/*.json  `icms compare --out` results

The only arithmetic here is counting records and taking medians and minima
of values the records already hold.  Ratios and intervals come from the
frozen comparison files.
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import statistics
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BASE = REPO / "docs" / "ic" / "measurement"
OUT_HTML = REPO / "docs" / "ic-measurement.html"
OUT_JSON = BASE / "page.json"
GITHUB = "https://github.com/aburan28/crypto/blob/main/"


def esc(s) -> str:
    return html.escape(str(s))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def g3(v) -> str:
    if v is None:
        return "—"
    if isinstance(v, (int,)) and not isinstance(v, bool):
        return f"{v:,}"
    if abs(v) >= 1000:
        return f"{v:,.0f}"
    return f"{v:.3g}"


def load_sessions() -> list[dict]:
    out = []
    for sdir in sorted((BASE / "sessions").glob("*/session.json")):
        d = sdir.parent
        session = json.loads(sdir.read_text())
        capsule = json.loads((d / "capsule.json").read_text())
        records = [json.loads(x) for x in (d / "records.jsonl").read_text().splitlines() if x.strip()]
        comparisons = []
        for c in sorted((d / "comparisons").glob("*.json")) if (d / "comparisons").exists() else []:
            comparisons.append({"name": c.stem, "path": c, "doc": json.loads(c.read_text())})
        out.append({"dir": d, "session": session, "capsule": capsule, "records": records, "comparisons": comparisons})
    return out


def arm_rows(s: dict) -> list[dict]:
    rows = []
    for arm in s["session"]["arms"]:
        rs = [r for r in s["records"] if r["arm"] == arm["index"] and not r["warmup"]]
        if not rs:
            continue
        first = rs[0]
        m = first.get("metrics") or {}
        fb, sysm, sol = m.get("factor_base") or {}, m.get("system") or {}, m.get("solver") or {}
        unit = first["unit"]
        totals = [((r.get("units") or {}).get(unit) or {}).get("total") for r in rs]
        totals = [t for t in totals if t is not None]
        walls = [r["execution"]["wall_ns"] / 1e6 for r in rs]
        win = [((r.get("windows") or {}).get(first["window"]) or {}).get("wall_ns") for r in rs]
        win = [w / 1e6 for w in win if w is not None]
        rows.append({
            "arm": arm["index"], "label": arm.get("label"), "spec_id": arm["spec_id"], "workload_id": arm["workload_id"],
            "runs": len(rs), "complete": sum(1 for r in rs if r["outcome"]["status"] == "complete" and r["outcome"].get("verified")),
            "levels": sorted({r["isolation"]["earned_level"] for r in rs}),
            "unit": unit, "window": first["window"],
            "deterministic": all(((r.get("units") or {}).get(unit) or {}).get("deterministic") for r in rs),
            "ops_median": statistics.median(totals) if totals else None,
            "ops_distinct": len(set(totals)),
            "S_median": statistics.median([((r.get("units") or {}).get(unit) or {}).get("S") for r in rs
                                           if ((r.get("units") or {}).get(unit) or {}).get("S") is not None] or [0]) or None,
            "window_wall_ms_median": statistics.median(win) if win else None,
            "window_wall_ms_min": min(win) if win else None,
            "process_wall_ms_median": statistics.median(walls), "process_wall_ms_min": min(walls),
            "usable_points": fb.get("usable_points"), "abscissae": fb.get("abscissae_with_points"),
            "columns": fb.get("columns"), "points_per_column": fb.get("points_per_column"),
            "n_vars": sysm.get("n_vars"), "n_equations": sysm.get("n_equations"),
            "cnf_clauses": sysm.get("cnf_clauses"), "xor_rows": sysm.get("xor_rows"),
            "solver": sol.get("name"), "solver_ops": sol.get("ops"), "solver_op_unit": sol.get("op_unit"),
            "priced_by": sol.get("priced_by"),
        })
    return rows


def gate_summary(s: dict) -> dict:
    checks: dict[str, dict[str, int]] = {}
    levels: dict[str, int] = {}
    for r in s["records"]:
        if r["warmup"]:
            continue
        levels[r["isolation"]["earned_level"]] = levels.get(r["isolation"]["earned_level"], 0) + 1
        for c in r["isolation"]["checks"]:
            d = checks.setdefault(c["id"], {"level": c["level"], "pass": 0, "fail": 0, "unknown": 0, "note": c["note"]})
            d[c["status"]] += 1
    return {"levels": levels, "checks": checks}


def capsule_facts(cap: dict) -> list[tuple[str, str]]:
    st = cap["stable"]
    t, k, c = st["topology"], st["kernel"], st["cpu"]
    govs = sorted({str(v.get("governor")) for v in t["cpus"].values()})
    feats = c.get("features") or {}
    return [
        ("env_class_id", cap["env_class_id"]),
        ("CPU", f"{c['model_name']} · family {c['family']} model {c['model']} stepping {c['stepping']} · microcode {c['microcode']} · {c['logical_cpus']} logical"),
        ("features", ", ".join(f for f in ("popcnt", "bmi2", "avx2", "avx512f", "pclmulqdq", "vpclmulqdq", "gfni", "aes", "hypervisor") if feats.get(f)) or "—"),
        ("kernel", f"{k['release']} · clocksource {k['clocksource']}"),
        ("isolation flags", f"isolcpus={t['isolated']!r} nohz_full={t['nohz_full']!r} · cmdline {json.dumps(k['cmdline_flags'], sort_keys=True)}"),
        ("frequency", f"governors {govs} · no_turbo {t['intel_pstate_no_turbo']} · boost {t['cpufreq_boost']}"),
        ("SMT", f"control {t['smt_control']} · active {t['smt_active']}"),
        ("virtualisation", str(st["virtualization"]["detect_virt"])),
        ("memory / THP", f"{st['memory'].get('MemTotal')} · THP {k['thp_enabled']}"),
        ("sysctls", ", ".join(f"{n.split('/')[-1]}={v}" for n, v in k["sysctl"].items() if v is not None)),
        ("toolchain", f"Python {st['toolchain']['python']} · {st['toolchain'].get('rustc')} · {st['toolchain'].get('gcc')}"),
        ("OS", f"{st['os'].get('PRETTY_NAME')} · {st['os'].get('libc')}"),
    ]


def build() -> dict:
    registry = json.loads((BASE / "registry.json").read_text())
    audit = json.loads((BASE / "audit.json").read_text())
    sessions = load_sessions()
    sources = {"registry": BASE / "registry.json", "audit": BASE / "audit.json"}
    doc_sessions = []
    for s in sessions:
        name = s["dir"].name
        for f in ("session.json", "capsule.json", "records.jsonl"):
            sources[f"{name}/{f}"] = s["dir"] / f
        for c in s["comparisons"]:
            sources[f"{name}/comparisons/{c['name']}.json"] = c["path"]
        doc_sessions.append({
            "name": name, "session_id": s["session"]["session_id"], "label": s["session"].get("label"),
            "env_class_id": s["session"]["env_class_id"],
            "commit": ((s["session"].get("implementation") or {}).get("repo") or {}).get("commit"),
            "dirty": ((s["session"].get("implementation") or {}).get("repo") or {}).get("dirty"),
            "preflight_quiet": s["session"]["preflight"]["quiet"],
            "reservation": s["session"].get("reservation"),
            "records": len(s["records"]), "measured": sum(1 for r in s["records"] if not r["warmup"]),
            "arms": arm_rows(s), "gates": gate_summary(s), "capsule": capsule_facts(s["capsule"]),
            "comparisons": [{"name": c["name"], **{k: c["doc"].get(k) for k in
                             ("a", "b", "aa", "declared_variables", "refusals", "ops", "wall", "aa_noise")}}
                            for c in s["comparisons"]],
        })
    return {
        "schema": "icms.page/v1",
        "generated_by": "scripts/build_ic_measurement.py",
        "what_this_is": "the ICMS v1 standard and its frozen demonstration sessions",
        "what_this_is_not": ["a speedup claim", "a measurement of any curve beyond the toy instances named"],
        "registry": {"windows": registry["windows"]["list"], "units": registry["units"]["list"],
                     "rho": registry["references"]["rho"], "size_fields": registry["factor_base_size_fields"]["fields"]},
        "audit": audit,
        "sessions": doc_sessions,
        "sources": {k: {"path": rel(p), "sha256": sha(p)} for k, p in sorted(sources.items())},
    }


STYLE = """
:root {
  --ground: #e7ebef; --surface: #fbfcfd; --sunk: #f0f3f6; --ink: #0e1922; --ink-2: #354855;
  --muted: #576976; --rule: #ccd5dc; --data: #1f6fd0; --data-soft: #d7e4f6; --bound: #b3372a;
  --bound-soft: #f3dcd8; --good: #227a4d; --good-soft: #d8efe2; --warn: #8a5a00; --warn-soft: #f6e7c8;
  --serif: "Newsreader", "Iowan Old Style", Georgia, serif;
  --sans: "IBM Plex Sans", system-ui, -apple-system, "Segoe UI", sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, "SF Mono", Menlo, monospace;
}
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --ground: #0d1318; --surface: #151d24; --sunk: #101820; --ink: #e3eaef; --ink-2: #b6c5cf;
  --muted: #8799a6; --rule: #25323a; --data: #4e92de; --data-soft: #1c2e42; --bound: #d9604e;
  --bound-soft: #38221f; --good: #4fc58a; --good-soft: #17301f; --warn: #e2b04a; --warn-soft: #33290f;
  color-scheme: dark; } }
:root[data-theme="dark"] {
  --ground: #0d1318; --surface: #151d24; --sunk: #101820; --ink: #e3eaef; --ink-2: #b6c5cf;
  --muted: #8799a6; --rule: #25323a; --data: #4e92de; --data-soft: #1c2e42; --bound: #d9604e;
  --bound-soft: #38221f; --good: #4fc58a; --good-soft: #17301f; --warn: #e2b04a; --warn-soft: #33290f;
  color-scheme: dark; }
* { box-sizing: border-box; }
body { margin: 0; padding-inline: 16px; padding-block: 36px 64px; background: var(--ground); color: var(--ink);
  font-family: var(--sans); font-size: 15px; line-height: 1.55; }
.wrap { max-width: 1080px; margin: 0 auto; display: flex; flex-direction: column; gap: 28px; }
header { display: flex; flex-direction: column; gap: 12px; }
.eyebrow { font-family: var(--mono); font-size: 11px; letter-spacing: .13em; text-transform: uppercase; color: var(--muted); }
h1 { font-family: var(--serif); font-weight: 500; font-size: clamp(30px, 5.2vw, 44px); line-height: 1.08; margin: 0; text-wrap: balance; }
h2 { font-family: var(--serif); font-weight: 500; font-size: 22px; margin: 0; text-wrap: balance; }
h3 { font-family: var(--sans); font-weight: 600; font-size: 14px; margin: 0; }
.verdict { font-family: var(--serif); font-size: clamp(17px, 2.2vw, 20px); line-height: 1.45; color: var(--ink-2); max-width: 66ch; margin: 0; }
.facts { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1px; background: var(--rule);
  border: 1px solid var(--rule); border-radius: 3px; overflow: hidden; margin: 0; }
.fact { background: var(--surface); padding: 14px 16px; display: flex; flex-direction: column; gap: 3px; min-width: 0; }
.fact dt { font-family: var(--mono); font-size: 10.5px; letter-spacing: .1em; text-transform: uppercase; color: var(--muted); }
.fact dd { margin: 0; font-family: var(--mono); font-size: 20px; font-weight: 500; font-variant-numeric: tabular-nums; }
.fact p { margin: 0; font-size: 12.5px; color: var(--muted); line-height: 1.4; overflow-wrap: anywhere; }
.panel { background: var(--surface); border: 1px solid var(--rule); border-radius: 4px; display: flex; flex-direction: column; gap: 12px; padding-block: 18px; }
.panel > * { padding-inline: 20px; }
.panel p, .panel li { margin: 0; max-width: 76ch; color: var(--ink-2); font-size: 14px; }
.panel ul, .panel ol { margin: 0; padding-inline-start: 40px; display: flex; flex-direction: column; gap: 4px; }
.scroll { overflow-x: auto; padding-inline: 0; }
table { border-collapse: collapse; width: 100%; font-size: 13px; font-variant-numeric: tabular-nums; }
th, td { padding: 7px 10px; text-align: left; border-bottom: 1px solid var(--rule); vertical-align: top; }
th { font-family: var(--mono); font-size: 10.5px; font-weight: 500; letter-spacing: .06em; text-transform: uppercase; color: var(--muted); white-space: nowrap; }
td.n, th.n { text-align: right; white-space: nowrap; }
code, .mono { font-family: var(--mono); font-size: 12px; overflow-wrap: anywhere; }
.chip { display: inline-block; font-family: var(--mono); font-size: 10px; letter-spacing: .08em; text-transform: uppercase;
  padding: 1px 6px; border: 1px solid var(--rule); border-radius: 2px; color: var(--ink-2); white-space: nowrap; }
.chip.pass { background: var(--good-soft); color: var(--good); border-color: transparent; }
.chip.fail { background: var(--bound-soft); color: var(--bound); border-color: transparent; }
.chip.unknown { background: var(--warn-soft); color: var(--warn); border-color: transparent; }
.chip.high { background: var(--bound-soft); color: var(--bound); border-color: transparent; }
.chip.medium { background: var(--warn-soft); color: var(--warn); border-color: transparent; }
.ratio { font-family: var(--mono); font-weight: 600; }
.muted { color: var(--muted); }
pre { margin: 0; background: var(--sunk); border: 1px solid var(--rule); border-radius: 3px; padding: 12px 14px;
  overflow-x: auto; font-family: var(--mono); font-size: 12px; line-height: 1.5; }
a { color: var(--data); }
a:focus-visible, summary:focus-visible { outline: 2px solid var(--data); outline-offset: 2px; }
details summary { cursor: pointer; font-family: var(--mono); font-size: 12px; color: var(--data); }
footer { font-size: 12.5px; color: var(--muted); display: flex; flex-direction: column; gap: 6px; }
footer p { margin: 0; max-width: 90ch; overflow-wrap: anywhere; }
"""


def chip(status: str) -> str:
    return f'<span class="chip {esc(status)}">{esc(status)}</span>'


def level_chip(level: str) -> str:
    cls = {"L3": "pass", "L2": "pass", "L1": "unknown"}.get(level, "fail")
    return f'<span class="chip {cls}">{esc(level)}</span>'


def link(path: str, text: str | None = None) -> str:
    return f'<a href="{esc(GITHUB + path)}">{esc(text or path)}</a>'


def page(doc: dict) -> str:
    sessions = doc["sessions"]
    audit = doc["audit"]["findings"]
    n_high = sum(1 for f in audit if f["severity"] == "high")
    measured = sum(s["measured"] for s in sessions)
    best = sorted({lvl for s in sessions for lvl in s["gates"]["levels"]})
    P = ['<div class="wrap"><header>',
         '<span class="eyebrow">ECDLP · index calculus · measurement standard v1 · no speedup claimed</span>',
         '<h1>Index Calculus Measurement Standard</h1>',
         '<p class="verdict">Two index-calculus figures are comparable only when they measured the same problem, in the same '
         'unit, over the same window, against the same reference, on the same class of host, under conditions that could '
         'not have moved the result. ICMS writes each of those down for every run and refuses a comparison when one '
         f'differs. The audit behind it found <strong>{len(audit)} ways the existing figures were not '
         f'apples to apples</strong>, {n_high} of them severe. The demonstration below runs on this repository\'s '
         'cloud container, which can be pinned but not isolated: wall time there is shown and refused, while operation '
         'counts are exact.</p>',
         '<dl class="facts">',
         f'<div class="fact"><dt>Audit findings</dt><dd>{len(audit)}</dd><p>{n_high} high severity; each with its evidence path</p></div>',
         f'<div class="fact"><dt>Runs recorded</dt><dd>{measured}</dd><p>in {len(sessions)} frozen session(s), warm-ups kept separately</p></div>',
         f'<div class="fact"><dt>Isolation earned here</dt><dd>{esc(", ".join(best) or "—")}</dd><p>L1 pinned, L2 quiet, L3 isolated; see the gate below</p></div>',
         f'<div class="fact"><dt>Units kept apart</dt><dd>{len(doc["registry"]["units"])}</dd><p>never compared with each other</p></div>',
         '</dl></header>']

    # What every figure states.
    P.append('<section class="panel" id="contract"><h2>What every figure must state</h2>'
             '<p>A record carries all of these, and a comparison checks each one before it computes anything. A field the '
             'host or producer did not expose is null and fails closed, never zero.</p><div class="scroll"><table>'
             '<thead><tr><th>element</th><th>what is recorded</th><th>a comparison refuses when</th></tr></thead><tbody>')
    rows = [
        ("Problem", "curve, subgroup, workload seed and target count, hashed into the workload id", "the workload ids differ"),
        ("Configuration", "the YAML spec: factor-base family and parameters, folds, large primes, partition; decomposition method, "
         "encoding, splitting, solver and every search-changing option; collector and stop rule; linear algebra", "any undeclared field differs"),
        ("Window", "online one target, cold end to end, whole process or stage only, with exact or derived conformance",
         "the windows differ, or a derived online window is offered as a wall-clock speedup"),
        ("Unit", "the unit id, whether the total was deterministic, and why not", "the units differ"),
        ("Reference", "the rho id, its method, automorphism count and runs", "the references differ"),
        ("Structure", "factor-base usable points, signed points, abscissae, columns; system variables, equations, degrees, "
         "clauses, XOR rows; solver conflicts and calls", "the producer's report contradicts the spec"),
        ("Host", "the capsule: CPU, microcode, features, kernel and command line, frequency policy, SMT, sysctls, THP, "
         "virtualisation, toolchain; hashed into env_class_id", "wall time on different env classes"),
        ("Conditions", "pinned CPUs read back, schedstat run-queue delay, steal, foreign runnable tasks, PSI, preemptions, "
         "the preflight; graded L0 to L3", "wall time below the declared level, or fewer than five paired rounds"),
        ("Versions", "commit, dirty flag and diff hash, binary sha256, compiled-in calibration hash, Cargo.lock, every ICMS "
         "source and the isolation tool", "—"),
    ]
    for a, b, c in rows:
        P.append(f'<tr><td><strong>{esc(a)}</strong></td><td>{esc(b)}</td><td>{esc(c)}</td></tr>')
    P.append('</tbody></table></div></section>')

    # Audit.
    P.append('<section class="panel" id="audit"><h2>Why the existing figures were not apples to apples</h2>'
             '<p>Found by reading the code and frozen evidence of this repository, cryptanalysis and crypto-autoresearcher. '
             'Each finding names where it was verified and how ICMS handles it.</p><div class="scroll"><table>'
             '<thead><tr><th>dimension</th><th>finding</th><th>severity</th><th>ICMS</th></tr></thead><tbody>')
    for f in audit:
        ev = "; ".join(link(e["path"], e["path"] + (f" ({e['where']})" if e.get("where") else "")) if e["repo"] == "crypto"
                       else esc(f"{e['repo']}: {e['path']}" + (f" ({e['where']})" if e.get('where') else "")) for e in f["evidence"])
        P.append(f'<tr><td><code>{esc(f["dimension"])}</code></td><td>{esc(f["finding"])}<br><span class="muted">{ev}</span></td>'
                 f'<td>{chip(f["severity"])}</td><td>{esc(f["icms"])}</td></tr>')
    P.append('</tbody></table></div></section>')

    # Levels.
    P.append('<section class="panel" id="levels"><h2>Isolation levels</h2>'
             '<p>A run earns the highest level whose checks, and all lower levels\' checks, passed. Unknown never passes. '
             'A spec may tighten a threshold, never loosen one. Wall time is admitted only at or above the spec\'s level, '
             'default L2. Operation counts need only L0, because they do not depend on contention.</p>'
             '<div class="scroll"><table><thead><tr><th>level</th><th>requires</th></tr></thead><tbody>'
             '<tr><td><strong>L0</strong> recorded</td><td>the run executed under a saved host capsule</td></tr>'
             '<tr><td><strong>L1</strong> pinned</td><td>child CPU mask read back equal to the reservation; CPU 0 excluded; other '
             'threads evicted; CPU time over wall at most the declared threads</td></tr>'
             '<tr><td><strong>L2</strong> quiet</td><td>child run-queue delay ≤ 0.5 %; zero hypervisor steal on the pinned CPUs; '
             'no foreign runnable task on them in any 50 ms sample; no memory stall; ≤ 50 preemptions/s; quiet preflight</td></tr>'
             '<tr><td><strong>L3</strong> isolated</td><td>pinned CPUs in isolcpus or nohz_full; performance governor, turbo off; '
             'SMT off or siblings reserved; bare metal</td></tr></tbody></table></div></section>')

    for s in sessions:
        P.append(f'<section class="panel" id="session-{esc(s["name"])}"><h2>Demonstration: {esc(s["label"] or s["name"])}</h2>'
                 f'<p>Session <code>{esc(s["session_id"])}</code> on host class <code>{esc(s["env_class_id"])}</code>, commit '
                 f'<code>{esc((s["commit"] or "")[:12])}</code>{" (dirty)" if s["dirty"] else ""}. Preflight '
                 f'{"quiet" if s["preflight_quiet"] else "refused: the machine was busy, so no run can reach L2"}; reserved CPUs '
                 f'{esc((s["reservation"] or {}).get("cpus"))}, {esc((s["reservation"] or {}).get("threads_moved"))} threads moved off. '
                 f'{s["measured"]} measured runs plus warm-ups, every one kept.</p>')
        P.append('<h3>Arms</h3><div class="scroll"><table><thead><tr><th class="n">arm</th><th>configuration</th>'
                 '<th class="n">runs ok</th><th>levels</th><th class="n">usable pts</th><th class="n">columns</th>'
                 '<th class="n">vars / eqs</th><th class="n">clauses / XOR</th><th class="n">ops (unit)</th><th class="n">S</th>'
                 '<th class="n">window ms med / min</th></tr></thead><tbody>')
        for a in s["arms"]:
            ops = f'{g3(a["ops_median"])}{"" if a["ops_distinct"] <= 1 else " (varies)"}'
            det = "" if a["deterministic"] else '<br><span class="muted">host-dependent</span>'
            P.append(f'<tr><td class="n">{a["arm"]}</td><td>{esc(a["label"])}<br><code class="muted">{esc(a["spec_id"])} · '
                     f'{esc(a["workload_id"])}</code>{("<br><span class=muted>solver " + esc(a["solver"]) + ", " + esc(a["priced_by"]) + "</span>") if a["solver"] else ""}</td>'
                     f'<td class="n">{a["complete"]}/{a["runs"]}</td><td>{" ".join(level_chip(lv) for lv in a["levels"])}</td>'
                     f'<td class="n">{g3(a["usable_points"])}</td><td class="n">{g3(a["columns"])}</td>'
                     f'<td class="n">{g3(a["n_vars"])} / {g3(a["n_equations"])}</td><td class="n">{g3(a["cnf_clauses"])} / {g3(a["xor_rows"])}</td>'
                     f'<td class="n">{ops}<br><span class="muted">{esc(a["unit"])}</span>{det}</td><td class="n">{g3(a["S_median"])}</td>'
                     f'<td class="n">{g3(a["window_wall_ms_median"])} / {g3(a["window_wall_ms_min"])}</td></tr>')
        P.append('</tbody></table></div>')
        if s["comparisons"]:
            P.append('<h3>Comparisons, as frozen by <code>icms compare</code></h3><div class="scroll"><table><thead><tr>'
                     '<th>pair</th><th>declared variable</th><th class="n">ops B/A</th><th>wall B/A</th><th>A/A noise</th></tr></thead><tbody>')
            for c in s["comparisons"]:
                ops, wall = c.get("ops") or {}, c.get("wall") or {}
                if c.get("refusals"):
                    ops_txt = "refused: " + "; ".join(r["reason"] for r in c["refusals"])
                elif ops.get("admitted"):
                    ops_txt = f'<span class="ratio">{ops["ratio_b_over_a"]:.4g}</span>' + ("" if ops["deterministic"] else ' <span class="muted">host-dependent</span>')
                else:
                    ops_txt = esc(ops.get("reason", "not admitted"))
                est = wall.get("estimate") or {}
                if wall.get("admitted"):
                    wall_txt = f'<span class="ratio">{est["median_ratio"]:.3g}</span> [{est["ci95"][0]:.3g}, {est["ci95"][1]:.3g}]'
                else:
                    reasons = "; ".join(r["reason"] for r in wall.get("refusals", [])) or "not admitted"
                    descr = (f'<br><span class="muted">descriptive only: median {est["median_ratio"]:.3g} over {est["n_pairs"]} pairs'
                             + (f', [{est["ci95"][0]:.3g}, {est["ci95"][1]:.3g}]' if est.get("ci95") else "") + '</span>') if est.get("median_ratio") else ""
                    wall_txt = f'refused: {esc(reasons)}{descr}'
                aa = "; ".join(f'arms {x["arms"]}: [{x["estimate"]["ci95"][0]:.3g}, {x["estimate"]["ci95"][1]:.3g}]'
                               for x in (c.get("aa_noise") or []) if (x.get("estimate") or {}).get("ci95")) or "—"
                decl = ", ".join(f'{d["field"]}: {d["a"]} → {d["b"]}' for d in (c.get("declared_variables") or [])) or ("none (A/A)" if c.get("aa") else "—")
                P.append(f'<tr><td>{esc((c.get("a") or {}).get("arm"))} vs {esc((c.get("b") or {}).get("arm"))}<br>'
                         f'<code class="muted">{esc(c["name"])}</code></td><td><code>{esc(decl)}</code></td>'
                         f'<td class="n">{ops_txt}</td><td>{wall_txt}</td><td class="mono">{esc(aa)}</td></tr>')
            P.append('</tbody></table></div>')
        P.append('<h3>What the gate saw</h3><div class="scroll"><table><thead><tr><th>check</th><th>level</th>'
                 '<th class="n">pass</th><th class="n">fail</th><th class="n">unknown</th><th>what it measures</th></tr></thead><tbody>')
        for cid, d in s["gates"]["checks"].items():
            P.append(f'<tr><td><code>{esc(cid)}</code></td><td>{esc(d["level"])}</td><td class="n">{d["pass"]}</td>'
                     f'<td class="n">{d["fail"]}</td><td class="n">{d["unknown"]}</td><td>{esc(d["note"])}</td></tr>')
        P.append('</tbody></table></div>')
        P.append('<details><summary>Host capsule</summary><div class="scroll"><table><tbody>')
        for k, v in s["capsule"]:
            P.append(f'<tr><th>{esc(k)}</th><td class="mono">{esc(v)}</td></tr>')
        P.append('</tbody></table></div></details></section>')

    # Registry.
    P.append('<section class="panel" id="registry"><h2>Windows and units</h2><div class="scroll"><table><thead><tr>'
             '<th>window</th><th>primary</th><th>definition</th></tr></thead><tbody>')
    for w in doc["registry"]["windows"]:
        P.append(f'<tr><td><code>{esc(w["name"])}</code></td><td>{"yes" if w["primary"] else "no"}</td><td>{esc(w["definition"])}</td></tr>')
    P.append('</tbody></table></div><div class="scroll"><table><thead><tr><th>unit</th><th>definition</th><th>deterministic</th></tr></thead><tbody>')
    for u in doc["registry"]["units"]:
        P.append(f'<tr><td><code>{esc(u["id"])}</code></td><td>{esc(u["definition"])}</td><td>{esc(u["deterministic"])}</td></tr>')
    P.append('</tbody></table></div><div class="scroll"><table><thead><tr><th>factor-base size field</th><th>means</th></tr></thead><tbody>')
    for k, v in doc["registry"]["size_fields"].items():
        P.append(f'<tr><td><code>{esc(k)}</code></td><td>{esc(v)}</td></tr>')
    P.append('</tbody></table></div></section>')

    P.append('<section class="panel" id="reproduce"><h2>Reproduce</h2><p>From the repository root, on the host whose '
             'capsule you want recorded. Building is heavy work: run it through <code>tools/isolated_bench.py busy</code> so '
             'it cannot overlap a timed run.</p><pre>python3 tools/isolated_bench.py busy -- cargo build --release --bin ic\n'
             'python3 tools/icms capsule\n'
             'python3 tools/icms validate docs/ic/measurement/specs/*.yaml\n'
             'python3 tools/icms run A.yaml B.yaml A.yaml --cpus 3 --out docs/ic/measurement/sessions/&lt;new&gt;\n'
             'python3 tools/icms compare docs/ic/measurement/sessions/&lt;new&gt; --a 0 --b 1 --declare &lt;field&gt; --out …/comparisons/&lt;name&gt;.json\n'
             'python3 scripts/build_ic_measurement.py</pre>'
             f'<p>The standard: {link("docs/ic/measurement/README.md")}. Implementation and tests: {link("tools/icms/__main__.py", "tools/icms/")}. '
             'Related pages: <a href="./">the cost ledger</a> and <a href="ic-leaderboard.html">the leaderboard</a>.</p></section>')

    P.append('<footer><p>Generated by <code>scripts/build_ic_measurement.py</code> from the frozen files below; this page '
             'computes nothing beyond counts and medians of recorded values.</p>')
    for k, v in doc["sources"].items():
        P.append(f'<p><code>{esc(v["path"])}</code> sha256 <code>{esc(v["sha256"][:12])}</code></p>')
    P.append('</footer></div>')
    head = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
            '<title>Index Calculus Measurement Standard</title>\n'
            '<meta name="description" content="ICMS v1: how index-calculus figures are made comparable: one spec, window, unit, '
            'reference and host record per run, an isolation gate, and comparisons that refuse a mismatch.">\n'
            '<link rel="preconnect" href="https://fonts.googleapis.com">\n'
            '<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>\n'
            '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,400;6..72,500'
            '&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500;600&display=swap">\n'
            f'<style>{STYLE}</style>\n</head>\n<body>\n<main id="main">\n')
    return head + "\n".join(P) + "\n</main>\n</body>\n</html>\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    doc = build()
    outs = {OUT_JSON: json.dumps(doc, indent=1, sort_keys=True) + "\n", OUT_HTML: page(doc)}
    if args.check:
        stale = [p for p, text in outs.items() if not p.exists() or p.read_text() != text]
        for p in stale:
            print(f"{rel(p)} is stale; run scripts/build_ic_measurement.py", file=sys.stderr)
        return 1 if stale else 0
    for p, text in outs.items():
        p.write_text(text)
        print(f"wrote {rel(p)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
