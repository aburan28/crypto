#!/usr/bin/env python3
"""ICV1 curve identities and ICCAN1 candidate identities for a tournament round.

    python3 campaign_20260916/round24_identities.py runs/round-0024 [--json OUT]

The convention is research/isogeny_volcano_ic_20260924/identity_schema.json.
Each tournament cell is one curve, y^2 + xy = x^3 + a*x^2 + 1 over GF(2^n) with
the fixture's irreducible modulus; every case of a cell shares it, which this
script checks.

The schema leaves two encodings open. They are pinned here, and anyone
recomputing these identities must use the same choices:

* canonical JSON is `json.dumps(obj, sort_keys=True, separators=(',', ':'))`;
* `modhash8` is sha256 of the canonical JSON of the modulus
  `{"degree": m, "low_terms": [...]}`, first 8 hex characters;
* a field element is its polynomial-basis integer in lowercase hex, `0x`-prefixed;
* the curve model is the canonical JSON of MODEL below: field, modulus, basis,
  Weierstrass coefficients a2 and a6, #E, subgroup order and cofactor. The
  generator is excluded: it names a point, not the curve.

`end` and `level` are `unk` and `path` is `r`: this campaign certifies neither
an endomorphism order nor a volcano position. The display slug drops the
schema's `ell<l>-L<level>` segment for the same reason: there is no volcano.

Candidate identities read each arm's factor base off its verified receipts
(signed base size / 2n orbits, one value per cell or the round is refused).
"""
import argparse
import collections
import hashlib
import json
from pathlib import Path


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(',', ':'))


def sha(obj):
    return hashlib.sha256(canonical(obj).encode()).hexdigest()


def curve_identity(fixture):
    m, a = fixture['degree'], fixture['curve_a']
    modulus = {'degree': fixture['irreducible']['degree'],
               'low_terms': fixture['irreducible']['low_terms']}
    order = int(fixture['group_order'])
    trace = 2 ** m + 1 - order
    model = {'schema': 'ICV1-model-1', 'field': {'kind': 'f2m', 'degree': m, 'modulus': modulus,
                                                 'basis': 'polynomial'},
             'curve': {'form': 'y^2+xy=x^3+a2*x^2+a6', 'a2': hex(a), 'a6': '0x1'},
             'order': str(order), 'subgroup_order': fixture['subgroup_order'],
             'cofactor': fixture['cofactor']}
    field = f'f2m-{m}-{sha(modulus)[:8]}'
    j = '0x1'  # j = 1/a6 on an ordinary binary curve, and a6 = 1 here
    model_hash = sha(model)[:12]
    return {'curve_id': f'ICV1:{field}:{trace}:{order}:{j}:unk:unk:r:{model_hash}',
            'slug': f'icv1-f2m{m}-t{trace}-pr-{model_hash}', 'model': model}


RELATION = {('pair_table', 3): 'pair-s3', ('triple_counted', 4): 'triple-s4-counted'}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('round', type=Path)
    p.add_argument('--json', type=Path)
    args = p.parse_args()
    fixtures = json.loads((args.round / 'fixtures.json').read_text())
    arms = json.loads((args.round / 'candidates.json').read_text())
    curves = {}
    for cases in fixtures.values():
        for case in cases:
            ident = curve_identity(case['fixture'])
            assert curves.setdefault(case['cell'], ident)['curve_id'] == ident['curve_id'], case['cell']
    orbits = collections.defaultdict(set)
    for receipt in args.round.glob('runs/*/*/*/rep-*/receipt.json'):
        r = json.loads(receipt.read_text())
        if r['status'] == 'VERIFIED' and r['mode'] == 'ic':
            n = int(r['cell'][1:].split('a')[0])
            orbits[(r['arm'], r['cell'])].add(r['certificate']['signed_base_size'] // (2 * n))
    candidates = {}
    for arm in arms:
        cfg = arm['config']
        relation = RELATION[(cfg['solver'], cfg['summands'])]
        for cell, ident in sorted(curves.items()):
            seen = orbits.get((arm['id'], cell))
            if not seen:
                continue
            assert len(seen) == 1, (arm['id'], cell, seen)
            fb = f'orbit{seen.pop()}'
            candidates.setdefault(arm['id'], {})[cell] = (
                f"ICCAN1/{ident['slug']}/{fb}/{relation}/{cfg['linear_algebra']}/frob-neg/r01")
    print('| cell | ICV1 curve identity |\n|:--|:--|')
    for cell, ident in sorted(curves.items(), key=lambda kv: int(kv[0][1:].split('a')[0])):
        print(f"| `{cell}` | `{ident['curve_id']}` |")
    for arm, per_cell in candidates.items():
        print(f'\n{arm}:')
        for cell, cid in sorted(per_cell.items(), key=lambda kv: int(kv[0][1:].split('a')[0])):
            print(f'  {cell:6} {cid}')
    if args.json:
        args.json.write_text(json.dumps({'convention': 'research/isogeny_volcano_ic_20260924/identity_schema.json',
                                         'curves': curves, 'candidates': candidates}, indent=2) + '\n')


if __name__ == '__main__':
    main()
