#!/usr/bin/env python3
"""Render a sealed ecbench session into a descriptive graph and PDF source."""
import datetime as dt
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

session_dir = Path(sys.argv[1])
study_dir = Path(__file__).parent
session = json.loads((session_dir / 'session.json').read_text())
records = [json.loads(line) for line in (session_dir / 'records.jsonl').read_text().splitlines()]
measured = [row for row in records if not row['warmup']]
by_arm = defaultdict(list)
for row in measured:
    by_arm[row['arm']].append(row)
start = dt.datetime.fromtimestamp(session['started_unix_ms'] / 1000, dt.timezone.utc)
report_dir = study_dir / 'reports' / f"run-{start:%Y%m%dT%H%MZ}"
report_dir.mkdir(parents=True, exist_ok=False)

summary = {}
for arm, rows in sorted(by_arm.items()):
    costs = [row['cost']['total_gae'] for row in rows]
    invs = [row['field_ops']['invs'] for row in rows if row['field_ops']]
    summary[arm] = {
        'n': len(rows),
        'verified': sum(row['outcome']['status'] == 'verified' for row in rows),
        'gae_mean': statistics.mean(costs),
        'gae_min': min(costs),
        'gae_max': max(costs),
        'invs_mean': statistics.mean(invs) if invs else None,
        'isolation_levels': sorted({row['isolation']['level'] for row in rows}),
        'rss_kib_max': max(row['time']['max_rss_kib'] for row in rows),
    }
(report_dir / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')

arms = list(summary)
fig, axes = plt.subplots(1, 2, figsize=(11.5, max(3.2, 0.43 * len(arms) + 1.5)), layout='constrained')
y = list(range(len(arms)))
means = [summary[a]['gae_mean'] for a in arms]
low = [means[i] - summary[a]['gae_min'] for i, a in enumerate(arms)]
high = [summary[a]['gae_max'] - means[i] for i, a in enumerate(arms)]
axes[0].barh(y, means, color='#4c78a8')
axes[0].errorbar(means, y, xerr=[low, high], fmt='none', ecolor='#243b53', capsize=2)
axes[0].set_xlabel('Cold counted GAE; bars mean, whiskers min-max')
axes[0].set_yticks(y, arms)
axes[0].invert_yaxis()
axes[0].grid(axis='x', alpha=0.2)
invs = [summary[a]['invs_mean'] or 0 for a in arms]
axes[1].barh(y, invs, color='#f58518')
axes[1].set_xlabel('Native field inversions per execution, mean')
axes[1].set_yticks(y, arms)
axes[1].invert_yaxis()
axes[1].grid(axis='x', alpha=0.2)
fig.suptitle(f"{session['session_id']} | measured executions only | L0 counts")
fig.savefig(report_dir / 'counts.svg')
fig.savefig(report_dir / 'counts.png', dpi=150)
plt.close(fig)

lines = [
    '#set page(margin: 18mm)', '#set text(size: 10pt)',
    f"= {session['label']}",
    f"Session `{session['session_id']}`. Source commit `{session['git_commit']}`; binary SHA-256 `{session['binary_sha256']}`.",
    f"Status `{session['status']}`; {session['records_written']}/{session['executions_planned']} execution records. Host class `{session['env_class_id']}`.",
    'All local wall values are L0 diagnostics; the figure shows counted cold work and native inversion counts.',
    '#image("counts.svg", width: 100%)',
    'The bars are arithmetic means over measured executions; whiskers span their observed minimum and maximum. No confidence interval or wall advantage is inferred from the figure.',
]
for arm in arms:
    s = summary[arm]
    lines.append(f"- `{arm}`: {s['verified']}/{s['n']} verified; mean {s['gae_mean']:.3f} cold GAE; mean {s['invs_mean'] if s['invs_mean'] is not None else 'missing'} native inversions; peak {s['rss_kib_max']} KiB RSS.")
lines.append(f"Raw records: `{session_dir}/records.jsonl`. Plot data: `summary.json`. Audit receipt is linked from the study README.")
(report_dir / 'REPORT.typ').write_text('\n\n'.join(lines) + '\n')
(report_dir / 'REPORT.md').write_text('\n'.join([
    f"# {session['label']}", '',
    f"Session `{session['session_id']}`, binary `{session['binary_sha256']}`; {session['records_written']}/{session['executions_planned']} records.",
    '', '![Measured cold work and field inversions](counts.svg)', '',
    'L0 operation and inversion counts. All wall time is descriptive.', '',
    f"Raw records: [{session_dir.name}](../../sessions/{session_dir.name}/records.jsonl).",
]) + '\n')
print(report_dir)
