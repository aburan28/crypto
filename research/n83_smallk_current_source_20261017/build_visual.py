#!/usr/bin/env python3
"""Build the source-linked N83 exposed-point current-source control chart."""
import json
import os
import tempfile
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'matplotlib-n83-controls'))
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.hashsalt'] = 'n83-current-source-controls-20261017'
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / 'reports' / 'run-20261017-current-source-controls'
SOURCE_SHA256 = '8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12'
BINARY_SHA256 = 'd71a9662c9d4e2dde58443bfc646a05737c6b9ca6b3d5e831b632ce09384c906'
SOURCES = (
    ('K1248 root11', HERE / 'replay_K1248.json', HERE / 'receipt_K1248.json'),
    ('K1280 root11', HERE / 'replay_K1280.json', HERE / 'receipt_K1280.json'),
    ('K2543 root11', ROOT / 'research/n83_k2543_current_source_20261017/replay_K2543.json',
     ROOT / 'research/n83_k2543_current_source_20261017/receipt_K2543.json'),
    ('K2816 root11', ROOT / 'research/n83_k2816_current_source_20261017/replay_K2816.json',
     ROOT / 'research/n83_k2816_current_source_20261017/receipt_K2816.json'),
    ('K3072 root11', ROOT / 'research/n83_k3072_current_source_20261017/replay.json',
     ROOT / 'research/n83_k3072_current_source_20261017/receipt.json'),
    ('K3450 root11', ROOT / 'research/n83_k3072_current_source_20261017/replay_K3450.json',
     ROOT / 'research/n83_k3072_current_source_20261017/receipt_K3450.json'),
    ('K3450 wide3', ROOT / 'research/n83_k3072_current_source_20261017/replay_K3450_wide3.json',
     ROOT / 'research/n83_k3072_current_source_20261017/receipt_K3450_wide3.json'),
    ('K3500 root4', ROOT / 'research/n83_k3072_current_source_20261017/replay_K3500_root4.json',
     ROOT / 'research/n83_k3072_current_source_20261017/receipt_K3500_root4.json'),
)


def main():
    rows = []
    for label, replay_path, receipt_path in SOURCES:
        replay = json.loads(replay_path.read_text())
        receipt = json.loads(receipt_path.read_text())
        if (replay['status'] != 'PASS_INDEPENDENT_ARITHMETIC_REPLAY' or
                receipt['status'] != 'PASS_PRODUCER' or
                receipt['source_sha256'] != SOURCE_SHA256 or
                receipt['binary_sha256'] != BINARY_SHA256 or
                replay['candidate_id'] is None):
            raise ValueError(f'{label}: unverified input')
        rows.append({'label': label, 'online_ms': replay['one_target_online_ms'],
                     'cold_s': receipt['wall_seconds'],
                     'candidate_id': replay['candidate_id'],
                     'replay': str(replay_path.relative_to(ROOT)),
                     'receipt': str(receipt_path.relative_to(ROOT))})
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'plot_data.json').write_text(json.dumps({'schema': 'n83.current-source-control-visual/v1',
                                                    'population': 'one previously exposed public point per arm',
                                                    'rows': rows}, sort_keys=True, indent=2) + '\n')
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4), layout='constrained')
    labels = [row['label'] for row in rows]
    positions = list(range(len(rows)))
    colors = ['#0b6e75' if row['label'].startswith('K12') else '#345d9d' for row in rows]
    for axis, key, title, unit in ((axes[0], 'online_ms', 'One-target online interval', 'ms'),
                                   (axes[1], 'cold_s', 'Cold process duration', 's')):
        values = [row[key] for row in rows]
        axis.barh(positions, values, color=colors, height=0.65)
        axis.set_yticks(positions, labels)
        axis.invert_yaxis()
        axis.set_xlim(0, max(values) * 1.23)
        axis.set_title(title, fontsize=12)
        axis.set_xlabel(unit)
        axis.grid(axis='x', alpha=0.2)
        axis.set_axisbelow(True)
        for y, value in enumerate(values):
            axis.text(value + max(values) * 0.015, y, f'{value:.3f}', va='center', fontsize=8)
    fig.suptitle('N83 current-source factor bases · same exposed public point · shared macOS host', fontsize=13)
    fig.text(0.5, -0.01, 'One verified process per arm; uncertainty is not estimated from a single point. All online intervals exclude reusable setup.',
             ha='center', fontsize=8)
    for suffix in ('svg', 'png'):
        options = {'metadata': {'Date': None}} if suffix == 'svg' else {}
        fig.savefig(OUT / f'controls.{suffix}', dpi=150, bbox_inches='tight',
                    **options)
    svg_path = OUT / 'controls.svg'
    svg_path.write_text('\n'.join(line.rstrip() for line in
                                  svg_path.read_text().splitlines()) + '\n')
    plt.close(fig)
    print(json.dumps({'output': str(OUT), 'rows': len(rows)}, sort_keys=True))


if __name__ == '__main__':
    main()
