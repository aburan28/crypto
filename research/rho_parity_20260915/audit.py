#!/usr/bin/env python3
"""Fail closed if an instruction binary differs from its native revision."""
import gzip
import hashlib
import json
import lzma
from pathlib import Path

HERE = Path(__file__).resolve().parent


def rows(record):
    return {r['phase']: r for r in
            (json.loads(s) for s in record['stdout'].splitlines() if s.startswith('{'))}


def key(record):
    return (record.get('kind', record.get('mode')), record['case'],
            record['variant'], record['seed'], record['known_log'])


def main():
    mapping = {'reference': ('iteration-1', 'reference'),
               'confirmation-reference': ('confirmation', 'reference'),
               'confirmation-candidate': ('confirmation', 'candidate')}
    mapping.update({f'iteration-{i}': (f'iteration-{i}', 'candidate') for i in range(1, 6)})
    results = []
    for name, (native_name, revision) in mapping.items():
        directory = HERE/'counts'/name
        if not (directory/'complete.json').exists():
            continue
        assert not (directory/'REJECTED.json').exists()
        native_path = HERE/'results'/native_name/'processes.jsonl.gz'
        profiles_path = directory/'profiles.jsonl.xz'
        native = {}
        with gzip.open(native_path, 'rt') as f:
            for line in f:
                r = json.loads(line)
                if r['kind'] in ['dlp', 'rho'] and r['revision'] == revision and r['rep'] == 0:
                    native[key(r)] = rows(r)
        checks = fields = oracle_inputs = 0
        with lzma.open(profiles_path, 'rt') as f:
            for line in f:
                r = json.loads(line)
                assert r['status'] == 'finished', (name, r['status'])
                if r['mode'].startswith('calibrate'):
                    continue
                expected = native[key(r)]
                for phase, row in rows(r).items():
                    for field, value in row.items():
                        if field.endswith('_ns'):
                            continue
                        # The immutable cost harness predates the added residual
                        # assignment counter. All fields it does emit must match.
                        assert expected[phase][field] == value, (name, key(r), phase, field)
                        fields += 1
                    oracle_inputs += len(row.get('attempts', []))
                checks += 1
        assert checks == 200, (name, checks)
        results.append(dict(corpus=name, native_comparison=native_name, native_revision=revision,
                            matched_processes=checks, matched_fields=fields,
                            matched_oracle_inputs=oracle_inputs,
                            native_sha256=hashlib.sha256(native_path.read_bytes()).hexdigest(),
                            profiles_sha256=hashlib.sha256(profiles_path.read_bytes()).hexdigest()))
    output = dict(accepted=results, rejected_excluded=[p.parent.name for p in
                  sorted((HERE/'counts').glob('*/REJECTED.json'))])
    (HERE/'instruction-audit.json').write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
