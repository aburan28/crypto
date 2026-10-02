#!/usr/bin/env python3
"""Render a readable IC overview from frozen evidence; keep the historical ledger."""
import hashlib
import html
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT/'docs/index-calculus-scoreboard.html'
GOAL = ROOT/'research/ic_candidate_tournament_20260915/goal_20260924'
START = '<!-- BEGIN IC DASHBOARD OVERVIEW -->'
END = '<!-- END IC DASHBOARD OVERVIEW -->'
LIBRARY = '<!-- BEGIN IC EVIDENCE LIBRARY -->'
CSS_START = '/* BEGIN IC OVERVIEW STYLES */'
CSS_END = '/* END IC OVERVIEW STYLES */'


def source(relative):
    path = GOAL/relative
    raw = path.read_bytes()
    return json.loads(raw), dict(path=path.relative_to(ROOT).as_posix(),
                                sha256=hashlib.sha256(raw).hexdigest())


def link(relative, label):
    return f'<a href="../research/ic_candidate_tournament_20260915/goal_20260924/{relative}">{label}</a>'


def render():
    round3, round_pin = source('improvement/round3/RESULTS.json')
    f5, f5_pin = source('prepared-f5-v3-control-v1/TERMINAL.json')
    sat, sat_pin = source('prepared-one-target-controls-v1/outcome/result.json')
    confirmation = round3['decision']['confirmation']
    online = confirmation['online']
    assert round3['decision']['winner'] == 'incumbent'
    assert not round3['decision']['promotion_eligible']
    assert f5['scalar_verified'] and not f5['headline_online_admissible']
    assert not f5['fresh_paired_qualification'] and f5['native_executions'] == 1
    assert sat['arms']['sat']['status'] == 'ADMITTED_INCOMPLETE_PREPARED_SAT_CONTROL'
    assert sat['arms']['sat']['attempt_count'] == 8 and not sat['arms']['sat']['scalar_verified']
    rho_rows = {r['alias']: r for r in round3['stages']['confirmation']['table']}
    paired_rho = round3['decision']['winner_over_online_rho']['confirmation']['online']
    rho_ms = rho_rows['rho_online']['online_ms']
    rho_bars = []
    for i, (alias, name) in enumerate([('incumbent','Incumbent IC'), ('stop7_word','Challenger IC'), ('rho_online','Matched rho')]):
        value = rho_rows[alias]['online_ms']; y = 35+i*54
        assert rho_rows[alias]['complete'] and rho_rows[alias]['verified'] == rho_rows[alias]['scheduled']
        rho_bars.append(f'<text x="10" y="{y+17}">{name}</text><rect x="153" y="{y}" width="{value/rho_ms*340:.2f}" height="25" rx="3" fill="{"#738375" if alias == "rho_online" else "#438565"}"/><text x="635" y="{y+17}" text-anchor="end">{value:.6f} ms</text>')
    rho_graph = '<svg viewBox="0 0 665 205" role="img" aria-labelledby="rho-graph-title"><title id="rho-graph-title">Round three confirmation: paired online costs in milliseconds, after reusable preparation. Incumbent IC 0.033926 ms, challenger 0.033950 ms, matched rho 0.252643 ms.</title>'+''.join(rho_bars)+'<text x="153" y="200">0</text><text x="493" y="200" text-anchor="end">0.253 ms</text></svg>'
    data = dict(schema_version=1, scope='bounded IC autolab overview; not all repository research',
                sources=[round_pin, f5_pin, sat_pin],
                confirmation_online=online,
                confirmation_familywise_online_upper=confirmation['familywise']['metrics']['online_ns']['upper'],
                paired_targets=confirmation['paired_cases'], f5_control=f5,
                rho_online_comparison=paired_rho, rho_online_table=[{k:r[k] for k in ('alias','online_ms','rho_online_over_IC_online','verified','scheduled')} for r in rho_rows.values()])
    (ROOT/'docs/ic/dashboard-overview-data.json').write_text(json.dumps(data, sort_keys=True, indent=2)+'\n')
    rows = []
    # Coordinates transform recorded ratios only; the browser does no statistics.
    def x(value):
        return 150+(value-.85)/.25*430
    for index, (cell, value) in enumerate(online['per_cell'].items()):
        y = 54+index*34
        rows.append(f'<text x="15" y="{y+5}">{html.escape(cell)}</text>'
                    f'<line x1="{x(1):.2f}" x2="{x(value):.2f}" y1="{y}" y2="{y}" class="comparison-link"/>'
                    f'<circle cx="{x(value):.2f}" cy="{y}" r="6" class="{"slower" if value>1 else "faster"}"/>'
                    f'<text x="640" y="{y+5}" text-anchor="end">{value:.3f}×</text>')
    ratio, (lo, hi) = online['candidate_over_baseline'], online['ci95']
    forest = f'''<svg viewBox="0 0 680 310" role="img" aria-labelledby="ratio-title ratio-desc">
<title id="ratio-title">Round three challenger divided by the qualified online IC reference</title>
<desc id="ratio-desc">Six recorded curve-cell ratios. Four are below one and two above. Overall ratio {ratio:.6f}, descriptive 95 percent interval {lo:.6f} to {hi:.6f}. It crosses one. Incumbent retained.</desc>
<rect x="150" y="24" width="258" height="236" fill="var(--dash-green-soft)"/>
<text x="151" y="15">Less time</text><text x="579" y="15" text-anchor="end">More time</text>
<line x1="{x(1):.2f}" x2="{x(1):.2f}" y1="24" y2="265" class="reference-line"/>
{''.join(rows)}
<text x="15" y="265" class="chart-bold">Overall</text>
<line x1="{x(lo):.2f}" x2="{x(hi):.2f}" y1="260" y2="260" class="interval"/>
<circle cx="{x(ratio):.2f}" cy="260" r="6" class="overall-dot"/>
<text x="640" y="265" text-anchor="end">{ratio:.3f}×</text>
<text x="150" y="300">0.85×</text><text x="{x(1):.2f}" y="300" text-anchor="middle">1.00× · reference</text>
<text x="580" y="300" text-anchor="end">1.10×</text></svg>'''
    front = f'''{START}
<main class="dash" id="ic-overview">
  <nav class="dash-nav" aria-label="Dashboard"><a class="dash-brand" href="#ic-overview"><span class="brand-mark">IC</span> Research lab</a><div><a href="#lab-results">Results</a><a href="#lab-pipeline">Pipeline</a><a href="#lab-next">Next experiment</a><a href="#evidence-search">Evidence</a></div></nav>
  <header class="dash-hero"><p class="dash-kicker">Index calculus · bounded autolab · updated 2 October 2026</p><h1>What works.<br>What still needs proof.</h1><p class="dash-lead">The tournament has kept its incumbent. F5 now completes a source-bound prepared control; SAT still needs that gate. A fresh comparison of these solver families has not run.</p><a class="dash-primary" href="#lab-results">See the measured results <span aria-hidden="true">↓</span></a></header>
  <div class="dash-metrics" aria-label="Current status">
    <article><span class="status neutral">Tournament</span><h2>Incumbent retained</h2><p>Three closed improvement rounds. No challenger passed promotion.</p>{link('improvement/round3/README.md','Last round →')}</article>
    <article><span class="status good">F5 · disclosed control</span><h2>Complete &amp; verified</h2><p>Three target attempts. Two negatives and one valid decomposition; scalar independently checked.</p>{link('prepared-f5-v3-control-v1/RESULT.md','Original evidence →')}</article>
    <article><span class="status pending">SAT · prepared control</span><h2>Budget exhausted</h2><p>Eight inconclusive attempts. No recovered scalar. A source-valid feasible query is retained.</p>{link('prepared-one-target-controls-v1/RESULT.md','Failure &amp; diagnosis →')}</article>
  </div>
  <section class="dash-section" id="lab-results" aria-labelledby="results-title">
    <div class="dash-section-head"><div><p class="dash-kicker">01 / Measured results</p><h2 id="results-title">Did the last challenger win?</h2></div><span class="status neutral">No promotion</span></div>
    <div class="dash-result-grid"><figure class="dash-chart"><figcaption><strong>One target at a time, six curve cells</strong><span>Challenger time ÷ qualified IC online reference · smaller is better</span></figcaption>{forest}<p class="dash-chart-note">Dots are recorded cell estimates, not per-cell confidence intervals. The overall line is the descriptive 95% interval.</p></figure>
    <div class="dash-explanation"><h3>A small average gain is not enough.</h3><p>The <code>stop7_word</code> challenger used <strong>{ratio:.3f}×</strong> the qualified online reference time. The overall interval crosses 1, and two cells were slower.</p><p>The stricter familywise upper bound is <strong>{data['confirmation_familywise_online_upper']:.3f}×</strong>. The frozen acceptance gate failed, so the incumbent stays.</p><div class="dash-small-facts"><span><strong>{confirmation['paired_cases']}</strong> paired targets</span><span><strong>6</strong> fixed curve cells</span></div><p class="dash-footnote">Synthetic toy panel, round 3 confirmation only. The denominator is the qualified <code>ic_online</code> role. Cell labels identify field degree/model, not subgroup bits. {link('improvement/round3/README.md','Full table, reference and replay →')}</p></div></div>
    <details class="dash-details"><summary>Exact plotted values and what this graph does not establish</summary><div class="dash-detail-body"><p>Online interval: one supplied public point after reusable preparation through scalar replay. This graph compares IC implementations; it does not plot F5 or SAT, or establish a rho crossover. The three old confirmation sets are closed. New solver families require new targets and a new frozen protocol.</p><table><caption>Recorded round 3 online challenger/reference ratios</caption><thead><tr><th>Cell</th><th>Ratio</th></tr></thead><tbody>{''.join(f'<tr><td>{html.escape(k)}</td><td>{v:.6f}×</td></tr>' for k,v in online['per_cell'].items())}<tr><td>Overall descriptive 95% interval</td><td>{lo:.6f}–{hi:.6f}</td></tr></tbody></table>{link('improvement/round3/RESULTS.json','Frozen JSON')} · <a href="ic/dashboard-overview-data.json">Overview data and source hashes</a></div></details>
  </section>
  <section class="dash-section" aria-labelledby="rho-results-title"><div class="dash-section-head"><div><p class="dash-kicker">Same panel / Matched rho</p><h2 id="rho-results-title">Does the incumbent beat rho here?</h2></div><span class="status good">Online only · toy panel</span></div><div class="dash-result-grid"><figure class="dash-chart"><figcaption><strong>The same targets, the same online boundary</strong><span>Recorded online milliseconds · lower is better · linear scale</span></figcaption>{rho_graph}<p class="dash-chart-note">Equal-cell geometric means of three-process medians for each single target. No batch amortization. All plotted arms verified 216/216 repetitions.</p></figure><div class="dash-explanation"><h3>Yes, on this registered toy workload.</h3><p>The qualified rho online reference took <strong>{rho_rows['incumbent']['rho_online_over_IC_online']:.2f}×</strong> the incumbent’s online time. The paired IC/rho cost ratio was <strong>{paired_rho['candidate_over_baseline']:.3f}×</strong>, with descriptive 95% interval <strong>{paired_rho['ci95'][0]:.3f}–{paired_rho['ci95'][1]:.3f}</strong>.</p><p>This is the old pair-table incumbent, not F5 or SAT. Reusable preparation is excluded from both online intervals. It establishes neither a cold-start result nor an advantage at larger sizes.</p><p class="dash-footnote">One point per solve, same point and frozen resources across arms. Reference quality is bounded by its accepted qualification. {link('improvement/round3/README.md','Paired rho context, separate cold costs and limits →')}</p></div></div></section>
  <section class="dash-section" id="lab-pipeline" aria-labelledby="pipeline-title"><div class="dash-section-head"><div><p class="dash-kicker">02 / How it works</p><h2 id="pipeline-title">A complete pipeline, not a fast fragment.</h2></div></div><p class="dash-section-intro">Build reusable knowledge first. Then time everything needed to recover one previously unseen target. A fast solver call alone cannot win the tournament.</p>
    <div class="pipeline-zone"><div class="pipeline-zone-label">Reusable preparation <span>reported separately</span></div><ol class="pipeline"><li><b>1</b><strong>Choose a factor base</strong><span>Exact points, subgroup filtering, sign &amp; Frobenius orbits.</span></li><li><b>2</b><strong>Collect relations</strong><span>Generate ordinary queries. Solve PDP with F4/F5, SAT or other candidates; verify each relation.</span></li><li><b>3</b><strong>Solve the matrix</strong><span>Recover factor logs over the subgroup modulus. Check rank and every log.</span></li></ol></div>
    <div class="pipeline-zone online-zone"><div class="pipeline-zone-label">One-target online interval <span>headline comparison</span></div><ol class="pipeline"><li><b>4</b><strong>Decompose the target</strong><span>Query, encode and solve. Charge every failed and timed-out attempt.</span></li><li><b>5</b><strong>Recover its logarithm</strong><span>Verify the decomposition and perform complete target descent.</span></li><li><b>6</b><strong>Check the scalar</strong><span>Independently replay the recovered answer. Stop the clock here.</span></li></ol></div>
    <p class="dash-footnote">PDP means point decomposition. The solver’s internal Macaulay matrix is part of PDP; it is different from final relation linear algebra. The matched rho arm must solve the same one point under the same resources.</p>
  </section>
  <section class="dash-section" aria-labelledby="control-title"><div class="dash-section-head"><div><p class="dash-kicker">03 / Current solver gate</p><h2 id="control-title">F5 runs end to end. Its cost is in PDP.</h2></div><span class="status neutral">Correctness control</span></div>
    <div class="dash-control-grid"><div><p class="control-number">11.956 <span>seconds</span></p><p>Raw native online interval for one <strong>already disclosed</strong> n17 target. All three attempts and scalar replay are included.</p><div class="cost-bar" role="img" aria-label="Point decomposition accounts for 11.955937581 seconds out of the 11.955969917 second online interval"><span>PDP · virtually all the time</span></div><p class="dash-footnote">Uncalibrated physical macOS ARM64 control. Known successful fixture seed, one worker, no matched rho arm. This is not expected fresh-target performance.</p></div><div class="dash-control-ledger"><div><span>Target query generation</span><b>11,293 ns</b></div><div><span>PDP · both failures &amp; success</span><b>11,955,937,581 ns</b></div><div><span>Relation verification</span><b>251 ns</b></div><div><span>Target descent</span><b>11,875 ns</b></div><div><span>Scalar recovery check</span><b>8,917 ns</b></div><div class="unknown"><span>Matched rho speedup</span><b>Unknown</b></div></div></div><p class="dash-footnote">62 usable points before folding · 63 geometric points · 29 relation columns. No new ordinary-query yield estimate. {link('prepared-f5-v3-control-v1/RESULT.md','Timing boundary, failed attempts and archive replay →')}</p>
  </section>
  <section class="dash-section" id="lab-next" aria-labelledby="next-title"><div class="dash-section-head"><div><p class="dash-kicker">04 / Next experiment</p><h2 id="next-title">Earn a comparison before declaring a winner.</h2></div></div><ol class="next-list"><li><span class="next-state done">Done</span><div><strong>Close the F5 source gate</strong><p>Original claim, native output, two negatives, witness and scalar survive exact archive restoration.</p></div></li><li><span class="next-state active">Next</span><div><strong>Qualify a complete prepared SAT path</strong><p>Diagnose the retained feasible query, then freeze a separate version and budget. Keep the exhausted run closed.</p></div></li><li><span class="next-state">Pending</span><div><strong>Freeze a fair new tournament</strong><p>Finish historical target exclusions, source-bound reference, host calibration, resources and execution order.</p></div></li><li><span class="next-state">Pending</span><div><strong>Compare fresh paired targets</strong><p>Measure complete F5/SAT pipelines against the incumbent and rho. Preserve failures, uncertainty and alternative solver/base families.</p></div></li></ol><p class="dash-footnote">The goal is a verified improvement under a declared workload. Local wins cannot establish a globally fastest IC implementation.</p></section>
  <section class="dash-library-intro" id="evidence-search" aria-labelledby="library-title"><p class="dash-kicker">05 / Evidence library</p><h2 id="library-title">Details when you need them.</h2><p>Historical panels include different workloads, units and evidence levels. They remain intact; they are not one combined ranking.</p><label for="evidence-query">Find a report</label><input id="evidence-query" type="search" placeholder="Try Frobenius, F5, rho, factor base…" autocomplete="off" aria-controls="evidence-results"><p id="evidence-count" aria-live="polite"></p><ul id="evidence-results"></ul><noscript><p>Open the full ledger below and use your browser’s Find command.</p></noscript></section>
</main>
{END}
'''
    styles = '''
/* BEGIN IC OVERVIEW STYLES */
:root{--dash-green:#23664d;--dash-green-soft:#e8f1ec;--dash-red:#9e503a}html{scroll-behavior:smooth;scroll-padding-top:24px}body{padding:0 24px 60px;background:#f5f5f0;color:#232c28}.dash{max-width:1120px;margin:auto;font-family:system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}.dash-nav{display:flex;align-items:center;justify-content:space-between;gap:24px;padding:27px 0;border-bottom:1px solid #d8ded5}.dash-nav div{display:flex;gap:26px;flex-wrap:wrap}.dash a{color:#23664d;text-underline-offset:4px}.dash-nav a{text-decoration:none;font-size:14px;font-weight:550}.dash-brand{display:flex;gap:11px;align-items:center;color:#232c28!important}.brand-mark{display:grid;place-items:center;width:34px;height:34px;border-radius:8px;background:#23664d;color:white;font-size:13px;letter-spacing:.04em}.dash-hero{display:block;padding:64px 0 40px}.dash-kicker{font-size:11px;text-transform:uppercase;letter-spacing:.14em;color:#68766b;font-weight:650;margin:0 0 14px}.dash h1{font-family:inherit;font-size:clamp(38px,5.5vw,64px);font-weight:650;letter-spacing:-.055em;line-height:1.04;margin:0 0 22px}.dash-lead{max-width:740px;font-size:19px;line-height:1.65;color:#566257;margin:0 0 24px}.dash-primary{display:inline-flex;align-items:center;gap:26px;padding:12px 18px;background:#23664d;color:white!important;border-radius:7px;text-decoration:none;font-size:14px;font-weight:600}.dash-metrics{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px;margin-bottom:54px}.dash-metrics article{background:#fff;border:1px solid #d8ded5;border-radius:12px;padding:25px;display:flex;flex-direction:column;align-items:flex-start;gap:12px}.dash-metrics h2{font-size:21px;letter-spacing:-.025em;line-height:1.25;margin:2px 0;font-weight:650}.dash-metrics p{font-size:14px;line-height:1.6;margin:0;color:#59665c;flex:1}.dash-metrics a{font-size:13px}.status{display:inline-block;align-self:flex-start;white-space:nowrap;font-size:11px;font-weight:650;border-radius:5px;padding:5px 8px;background:#edf0eb;color:#586355}.status.good{background:#e3f0e8;color:#24664c}.status.pending{background:#f9eadc;color:#8d562a}.dash-section{border-top:1px solid #d8ded5;padding:35px 0 38px;display:block}.dash-section-head{display:flex;justify-content:space-between;gap:20px;align-items:center;margin-bottom:25px}.dash-section-head .dash-kicker{margin-bottom:10px}.dash h2{font-family:inherit;font-size:28px;letter-spacing:-.035em;line-height:1.23;font-weight:630;margin:0}.dash-metrics h2{font-size:21px}.dash h3{font-family:inherit;font-size:22px;letter-spacing:-.025em;font-weight:620;line-height:1.3;margin:0 0 14px}.dash-result-grid{display:grid;grid-template-columns:minmax(0,1.6fr) minmax(0,1fr);gap:32px}.dash-chart{margin:0;background:white;border:1px solid #d8ded5;padding:22px;border-radius:12px}.dash-chart figcaption{display:flex;flex-direction:column;gap:7px;font-size:14px;margin-bottom:22px}.dash-chart figcaption span{font-size:12px;color:#657365}.dash-chart svg{width:100%;height:auto;overflow:visible}.dash-chart svg text{font-family:inherit;fill:#5a685e;font-size:13px}.comparison-link{stroke:#c1cec3;stroke-width:2}.reference-line{stroke:#78837b;stroke-width:1.5;stroke-dasharray:4 5}.faster{fill:#23664d}.slower{fill:#9e503a}.interval{stroke:#232c28;stroke-width:3}.overall-dot{fill:#232c28}.dash-chart svg .chart-bold{font-weight:700;fill:#232c28}.dash-chart-note{font-size:11px;color:#667269;line-height:1.5;margin:13px 0 0}.dash-explanation{padding:12px 0;line-height:1.7;font-size:15px;color:#566257}.dash-explanation strong{color:#232c28}.dash-small-facts{display:flex;gap:32px;padding:12px 0}.dash-small-facts span{display:flex;flex-direction:column;font-size:12px}.dash-small-facts strong{font-size:25px;font-weight:600}.dash .dash-footnote{font-size:12px;color:#667269;line-height:1.65;margin:16px 0 0}.dash-details{margin-top:18px;border:1px solid #d8ded5;border-radius:8px;background:#f0f2ec}.dash-details summary{padding:14px 18px;cursor:pointer;font-size:13px;font-weight:550}.dash-detail-body{padding:0 18px 18px;font-size:13px;overflow:auto}.dash-detail-body table{width:100%;font-size:12px;border-collapse:collapse}.dash-detail-body th,.dash-detail-body td{text-align:left;padding:8px;border-bottom:1px solid #d8ded5}.dash-detail-body caption{text-align:left;padding:10px 0}.dash-section-intro{max-width:770px;font-size:16px;line-height:1.65;color:#566257;margin:0 0 25px}.pipeline-zone{border:1px solid #d8ded5;border-radius:10px;background:#ecefe7;padding:20px 22px;margin-bottom:16px}.online-zone{background:#e8f1ec;border-color:#c3d6c8}.pipeline-zone-label{font-size:12px;text-transform:uppercase;letter-spacing:.07em;font-weight:650;margin-bottom:20px}.pipeline-zone-label span{display:inline-block;margin-left:14px;text-transform:none;letter-spacing:0;font-weight:400;color:#68746b}.pipeline{list-style:none;padding:0;margin:0;display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:34px}.pipeline li{position:relative;display:grid;grid-template-columns:30px 1fr;column-gap:10px;align-content:start}.pipeline li:not(:last-child):after{content:"→";position:absolute;right:-25px;top:6px;color:#6e8472}.pipeline b{grid-row:1/3;display:grid;place-items:center;border:1px solid #bccbbe;width:28px;height:28px;border-radius:50%;font-size:12px;color:#23664d;background:#fff8}.pipeline strong{font-size:14px;font-weight:650;line-height:1.45}.pipeline span{font-size:12px;color:#607066;line-height:1.65;margin-top:7px}.dash-control-grid{display:grid;grid-template-columns:1fr 1.1fr;gap:45px}.control-number{font-size:48px;font-weight:600;letter-spacing:-.04em;margin:0 0 12px;line-height:1}.control-number span{font-size:18px;font-weight:400;letter-spacing:0;color:#6c776e}.dash-control-grid p:not(.control-number){font-size:14px;color:#59665c;line-height:1.65}.cost-bar{background:#23664d;color:white;height:40px;border-radius:5px;display:flex;align-items:center;padding:0 14px;font-size:12px;margin-top:22px}.dash-control-ledger div{display:flex;justify-content:space-between;gap:15px;padding:12px 0;border-bottom:1px solid #d8ded5;font-size:12px;line-height:1.5}.dash-control-ledger b{font-family:ui-monospace,monospace;white-space:nowrap;font-size:11px;font-weight:550}.dash-control-ledger .unknown{color:#7d6848}.next-list{list-style:none;padding:0;margin:0}.next-list li{display:grid;grid-template-columns:85px 1fr;gap:22px;padding:18px 0;border-bottom:1px solid #e2e6df}.next-list li:last-child{border-bottom:0}.next-list strong{font-size:15px}.next-list p{font-size:13px;line-height:1.6;margin:6px 0 0;color:#667269}.next-state{font-size:11px;border:1px solid #d8ded5;border-radius:5px;align-self:start;padding:5px 9px;text-align:center;color:#68746b}.next-state.done{background:#e3f0e8;border-color:#c5d9cd;color:#24664c}.next-state.active{background:#f9eadc;border-color:#e6cba8;color:#8d562a}.dash-library-intro{border-top:1px solid #d8ded5;padding:36px 0 16px;display:block}.dash-library-intro>p{color:#667269;font-size:14px;line-height:1.6}.dash-library-intro label{display:block;font-size:12px;font-weight:650;margin:23px 0 8px}.dash input{width:100%;max-width:590px;border:1px solid #bdcabe;background:white;border-radius:7px;padding:13px 15px;font:inherit;font-size:14px;color:#232c28}.dash #evidence-count{font-size:12px;min-height:20px}.dash #evidence-results{padding:0;margin:0;list-style:none;display:grid;grid-template-columns:1fr 1fr;gap:0 30px}.dash #evidence-results li{border-bottom:1px solid #e0e5dd;padding:11px 0}.dash #evidence-results a{font-size:13px;line-height:1.6;text-decoration:none}.dash #evidence-results a:hover{text-decoration:underline}.dash-library{max-width:1120px;margin:12px auto 0;border:1px solid #d8ded5;border-radius:10px;background:var(--surface)}.dash-library>summary{font-family:system-ui,sans-serif;padding:18px 22px;font-size:14px;font-weight:600;cursor:pointer;color:var(--ink)}#legacy-evidence{padding:22px;max-width:100%;overflow:hidden}#legacy-evidence:before{content:"Historical ledger · each panel keeps its original scope and date";display:block;padding:12px 0 24px;color:var(--muted);font-size:12px}#legacy-evidence section,#legacy-evidence .panel{scroll-margin-top:30px}#legacy-evidence h1{font-size:30px}.dash a:focus-visible,.dash input:focus-visible,summary:focus-visible{outline:3px solid #b58a42;outline-offset:4px}.dash code{font-size:.9em;overflow-wrap:anywhere}.dash section[id]{scroll-margin-top:25px}
@media(max-width:760px){body{padding:0 18px 32px}.dash-nav{align-items:flex-start;gap:15px}.dash-nav div{gap:12px;font-size:12px;justify-content:flex-end;max-width:215px}.dash-nav a{font-size:12px}.dash-hero{padding:42px 0 28px}.dash-lead{font-size:16px}.dash-metrics{grid-template-columns:1fr;gap:12px;margin-bottom:30px}.dash-metrics article{padding:20px}.dash-result-grid,.dash-control-grid{grid-template-columns:1fr;gap:22px}.dash h2{font-size:24px}.dash-section-head{align-items:flex-start}.dash-chart{padding:16px}.dash-chart svg text{font-size:26px}.dash-chart svg circle{r:8px}.pipeline{grid-template-columns:1fr;gap:28px}.pipeline li:not(:last-child):after{content:"↓";top:auto;bottom:-24px;left:8px;right:auto}.pipeline-zone-label span{display:block;margin:5px 0 0}.dash #evidence-results{grid-template-columns:1fr}.dash-control-ledger div{font-size:11px}.dash-control-ledger b{font-size:10px}.dash-section{padding:28px 0}.next-list li{grid-template-columns:68px 1fr;gap:14px}#legacy-evidence{padding:12px}.dash-section-head>.status{white-space:normal;text-align:center;max-width:100px}}
@media(prefers-reduced-motion:reduce){html{scroll-behavior:auto}}
@media(prefers-color-scheme:dark){body{background:#171e1a;color:#e4ebe6}:root{--dash-green-soft:#203d2c}.dash-nav,.dash-section,.dash-library-intro{border-color:#35473a}.dash-brand,.dash-explanation strong{color:#e4ebe6!important}.dash a{color:#9bcbae}.dash-lead,.dash-explanation,.dash-metrics p,.dash-section-intro,.dash-control-grid p:not(.control-number),.next-list p,.dash-library-intro>p{color:#a5b9ab}.dash-metrics article,.dash-chart,.dash input{background:#202a23;border-color:#35473a;color:#e4ebe6}.dash-kicker,.dash .dash-footnote,.dash-chart-note,.dash-chart figcaption span{color:#a1b6a7}.dash-chart svg text{fill:#b6c5bc}.dash-chart svg .chart-bold{fill:#e4ebe6}.overall-dot{fill:#e4ebe6}.interval{stroke:#e4ebe6}.faster{fill:#91c4a2}.slower{fill:#dfa187}.dash-details,.pipeline-zone{background:#232e25;border-color:#35473a}.online-zone{background:#20372a}.pipeline span,.pipeline-zone-label span{color:#b3c5b8}.pipeline b{background:#2e4134;color:#a6d3b5;border-color:#57735e}.dash-control-ledger div,.next-list li,.dash #evidence-results li{border-color:#35473a}.dash-library{border-color:#35473a}.dash-detail-body td,.dash-detail-body th{border-color:#35473a}}
/* END IC OVERVIEW STYLES */
'''
    script = '''<script id="ic-overview-interactions">
(function () {
  const library = document.getElementById('evidence-library');
  const ledger = document.getElementById('legacy-evidence');
  const query = document.getElementById('evidence-query');
  const list = document.getElementById('evidence-results');
  const count = document.getElementById('evidence-count');
  const reports = Array.from(ledger.querySelectorAll('h2')).map((heading, i) => {
    let target = heading.closest('[id]');
    if (!target || target === ledger) { heading.id = 'evidence-record-' + i; target = heading; }
    return { title: heading.textContent.trim(), id: target.id, text: (heading.closest('section') || heading.parentElement).textContent.toLowerCase() };
  });
  function search() {
    const words = query.value.toLowerCase().trim().split(/\\s+/).filter(Boolean);
    const matches = reports.filter(r => words.every(w => r.text.includes(w)));
    list.replaceChildren();
    matches.slice(0, 12).forEach(r => {
      const li = document.createElement('li'); const a = document.createElement('a');
      a.href = '#' + r.id; a.textContent = r.title;
      a.addEventListener('click', () => { library.open = true; });
      li.append(a); list.append(li);
    });
    count.textContent = matches.length ? matches.length + ' reports' + (matches.length > 12 ? ' · showing 12; narrow your search' : '') : 'No matching reports. Try a broader term, or open the full ledger.';
  }
  function revealHash() {
    let id; try { id = decodeURIComponent(location.hash.slice(1)); } catch (_) { return; }
    const target = document.getElementById(id);
    if (target && ledger.contains(target)) {
      library.open = true;
      requestAnimationFrame(() => target.scrollIntoView());
    }
  }
  query.addEventListener('input', search);
  window.addEventListener('hashchange', revealHash);
  search(); revealHash();
})();
</script>'''
    page = PAGE.read_text()
    if START in page:
        before, rest = page.split(START, 1)
        _, after = rest.split(END, 1)
        page = before + front.rstrip('\n') + after
        old_start = page.index('<script id="ic-overview-interactions">')
        old_end = page.index('</script>', old_start)+len('</script>')
        page = page[:old_start]+script+page[old_end:]
    else:
        assert page.count('<body>') == 1 and page.count('</body>') == 1
        page = page.replace('<body>', '<body>\n'+front+LIBRARY+'\n<details class="dash-library" id="evidence-library"><summary>Open the full historical ledger</summary><div id="legacy-evidence">', 1)
        page = page.replace('</body>', '</div></details>\n<!-- END IC EVIDENCE LIBRARY -->\n'+script+'\n</body>', 1)
    if CSS_START in page:
        start = page.index(CSS_START)
        end = page.index(CSS_END, start)+len(CSS_END)
        page = page[:start]+styles.strip()+page[end:]
    else:
        page = page.replace('</head>', '<style>'+styles+'</style>\n</head>', 1)
    page = page.replace('<title>Index Calculus Scoreboard</title>', '<title>IC Research Lab — Results, Pipeline &amp; Evidence</title>')
    PAGE.write_text(page)
    print('Rendered evidence-bound overview; original historical ledger retained.')


if __name__ == '__main__':
    render()
