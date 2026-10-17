#!/usr/bin/env python3
"""Audit high-K fixture overlap against versioned N83 JSONL point rows."""

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EARLIER_FIXTURE = ROOT / 'research/n83_highk_tiebreak_20261017'


def main():
    targets = {tuple(json.loads(line)) for line in
               (HERE / 'targets.jsonl').read_text().splitlines()}
    if len(targets) != 10:
        raise ValueError('ten distinct high-K fixture points required')
    files, rows, overlaps = 0, 0, []
    for path in sorted((ROOT / 'research').glob('n83*/**/*.jsonl')):
        # generate.py separately binds the byte-identical earlier fixture.
        if path.parent in (HERE, EARLIER_FIXTURE):
            continue
        files += 1
        for line_number, line in enumerate(path.open(), start=1):
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            if (isinstance(row, list) and len(row) == 2 and
                    all(isinstance(value, str) for value in row)):
                rows += 1
                if tuple(row) in targets:
                    overlaps.append({'file': str(path.relative_to(ROOT)),
                                     'line': line_number})
    result = {'status': 'PASS_NO_OVERLAP' if not overlaps else 'OVERLAP',
              'fixture_points': len(targets), 'files_scanned': files,
              'two_string_point_rows': rows, 'overlaps': overlaps}
    print(json.dumps(result, sort_keys=True))
    return 0 if not overlaps else 1


if __name__ == '__main__':
    raise SystemExit(main())
