#!/usr/bin/env python3
"""Merge native curve_standards output. Metadata only; no curve arithmetic."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / 'docs/curves/registry.json'
STANDARDS = ROOT / 'docs/curves/standards/registry.json'


def merge(catalog, imported):
    by_model = {row['model_json']: row for row in catalog['curves']}
    for new in imported['curves']:
        raw = new['model_json']
        digest = hashlib.sha256(raw.encode()).hexdigest()
        if not new['slug'].endswith(digest[:8]) or not new['icv1'].endswith(digest[:12]):
            raise ValueError('native model identity digest mismatch')
        if raw not in by_model:
            catalog['curves'].append(new)
            by_model[raw] = new
            continue
        old = by_model[raw]
        if old['slug'] != new['slug'] or old['order'] != new['order']:
            raise ValueError('conflicting model identity/order')
        for key in ('aliases', 'standard_names', 'standards_provenance', 'sources'):
            for value in new.get(key, []):
                if value not in old.setdefault(key, []):
                    old[key].append(value)
        reps = old.setdefault('representations', [])
        known = {r['curve_uid']: r for r in reps}
        for rep in new['representations']:
            if rep['curve_uid'] not in known:
                reps.append(rep)
                known[rep['curve_uid']] = rep
        if reps:
            old.pop('ec1_unresolved', None)
    # Ambiguous friendly labels do not establish a model identity.
    models = {}
    for row in catalog['curves']:
        for alias in row['aliases']:
            models.setdefault(alias.lower(), set()).add(row['slug'])
    for row in catalog['curves']:
        row['aliases'] = sorted([a for a in row['aliases'] if len(models[a.lower()]) == 1], key=str.lower)
    catalog['standards_source'] = 'docs/curves/standards/registry.json'
    return catalog


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--check', action='store_true')
    args = ap.parse_args()
    before = REGISTRY.read_text()
    result = merge(json.loads(before), json.loads(STANDARDS.read_text()))
    text = json.dumps(result, indent=1, ensure_ascii=False) + '\n'
    if args.check:
        if text != before:
            raise SystemExit('standard catalog join is stale')
    else:
        REGISTRY.write_text(text)
        # Alias index formatting only; no identity arithmetic is invoked.
        import curve_id
        curve_id.ALIAS_MAP.write_text(curve_id.alias_map_text(result))
    print(f"standard catalog join: {len(result['curves'])} models")


if __name__ == '__main__':
    main()
