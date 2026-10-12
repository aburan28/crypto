#!/usr/bin/env python3
"""Regenerate the Hamming-ideal PDP panel of docs/index-calculus-scoreboard.html
between its begin/end markers from the frozen summary. The page cites and
never computes: every number here is read from results/main/summary.json.

    python3 scoreboard_panel.py
"""
import json, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
PAGE = os.path.join(REPO, "docs", "index-calculus-scoreboard.html")
BEGIN, END = "<!-- hamming-ideal-pdp:begin -->", "<!-- hamming-ideal-pdp:end -->"
ENC = ["SUB", "MONO", "C", "FC", "QFC"]

def main():
    s = json.load(open(os.path.join(HERE, "results", "main", "summary.json")))
    rows = {(r["n"], r["encoding"]): r for r in s["rows"]}
    sizes = sorted({r["n"] for r in s["rows"] if r["n"] >= 11})
    fits = {f["encoding"]: f for f in s["fits_calls"]}
    def cell(n, e, key, fmt):
        r = rows.get((n, e))
        return "pending" if r is None else fmt(r[key])
    head = "<th>encoding</th>" + "".join(f"<th>n = {n}: calls / target</th><th>n = {n}: tame depth</th>" for n in sizes) + "<th>exponent in |F|</th><th>class</th>"
    body = ""
    for e in ENC:
        body += f"<tr><td>{e}</td>"
        for n in sizes:
            body += f"<td>{cell(n, e, 'calls_mean', lambda v: f'{v:.1f}')}</td><td>{cell(n, e, 'tame_depth_mean', lambda v: f'{v:.1f}')}</td>"
        f = fits.get(e)
        body += f"<td>{f['exponent']:.2f} over n = {', '.join(map(str, f['sizes']))}</td>" if f else "<td>pending</td>"
        body += "<td>accounting</td></tr>"
    fsizes = "".join(f"<th>n = {n}</th>" for n in sizes)
    frow = "".join(f"<td>{rows[(n, 'FC')]['F_points'] if (n, 'FC') in rows else 'pending'}</td>" for n in sizes)
    panel = f'''{BEGIN}
    <div class="panel" id="hamming-ideal-pdp">
      <div class="panel-head">
        <h2>Hamming ideals as a membership ideal for a Frobenius-stable weight base <span class="chip">accounting</span></h2>
        <p>Frozen source: <code>research/hamming_ideal_pdp_20260930/results/main/summary.json</code> (derived from the per-target
          records beside it; protocol and result in the same directory). A normal-basis weight-bounded base
          <code>F_w = {{P : wt_NB(x(P)) &le; w}}</code> is Frobenius-stable at every prime <code>n</code>, needs no storage, and by
          La Scala&ndash;Marchesin&ndash;Tiwari has a bounded-degree membership ideal (C, FC, QFC). Their <code>MultiSolve</code>
          oracle, with a degree-truncated Boolean F4 inside, was run on <code>K_0/F_2^n</code> against the exhaustive oracle
          (<code>|F_w|</code> candidates per target) and a matched subspace base (SUB) through the same solver; MONO is the
          monomial-ideal control. Every target of every cell agrees with the exhaustive oracle. Unit: <code>GroebnerSafe</code>
          calls per target; the tame depth is the number of summand coordinates fixed before the basis became linear.
          No <code>S</code>, rho ratio or end-to-end figure: a stage diagnostic of one oracle, <code>m = 2</code>, <code>w = 2</code>.</p>
      </div>
      <div class="table-scroll"><table><caption>calls per target and mean tame depth (of n coordinates per summand); exponent = least-squares slope of log&#8322; calls against log&#8322; |F| over the prime sizes.</caption>
        <thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>
      <div class="table-scroll"><table><caption>the exhaustive oracle's candidate count |F_w| (points of the weight base) at each size.</caption>
        <thead><tr>{fsizes}</tr></thead><tbody><tr>{frow}</tr></tbody></table></div>
      <div class="panel-head"><p>Verdict: the paper&rsquo;s negative result on Classic McEliece transfers. The algebra resolves only after
        nearly a whole summand is fixed and the call count grows faster than the candidate count, so the oracle is the
        enumeration with a Gr&ouml;bner computation attached to each candidate. What the construction changes is the
        decomposition note&rsquo;s &sect;6 storage term for a Frobenius-stable set (removed); no operation count moves,
        and the verdict at the top of this page is unchanged. Budget sensitivity: 16&times; the per-call budget lowers the tame
        depth by one to two levels at <code>n = 11</code> (<code>results/tau/TABLE.md</code>).</p></div>
    </div>
    {END}'''
    page = open(PAGE).read()
    if BEGIN in page:
        page = re.sub(re.escape(BEGIN) + r".*?" + re.escape(END), lambda m: panel, page, flags=re.S)
    else:
        anchor = '''    <div class="panel">
      <div class="panel-head">
        <h2>The solver inside S: one base, four oracles, whole runs'''
        assert anchor in page
        page = page.replace(anchor, panel + "\n\n" + anchor, 1)
    open(PAGE, "w").write(page)
    print("panel written")

if __name__ == "__main__":
    main()
