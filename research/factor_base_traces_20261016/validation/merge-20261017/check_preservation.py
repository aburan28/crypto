#!/usr/bin/env python3
"""Read-only metadata/blob preservation audit; no cryptographic computation."""
import copy
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[4]
HEAD = '0678100c4e0660f69ef707724d07b22ed3a4b58e'
MAIN = 'db69b9bd58976dfdadab7623b9c75165b23ff21d'
PRIOR = 'd78eed09a48bec50b8fab3cafc9c0203fb7d5b99'
CREATION = '0b4a7451c7bb080642f03fccb641066dd11f795d'

def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)

def previous(ref, path):
    return json.loads(git('show', f'{ref}:{path}'))

def current(path):
    return json.loads((ROOT / path).read_text())

index = {}
for line in git('ls-files', '--stage', '-z').decode().split('\0'):
    if line:
        meta, path = line.split('\t', 1)
        mode, sha, stage = meta.split()
        assert stage == '0', f'unresolved merge: {path}'
        index[path] = (mode, sha)

payload = git('diff', '--name-only', f'{PRIOR}^2', PRIOR).decode().splitlines()
checked = 0
for path in payload:
    if path == 'store/index.json':
        continue
    assert index[path][1] == git('rev-parse', f'{PRIOR}:{path}').decode().strip(), path
    assert git('hash-object', '--', path).decode().strip() == index[path][1], path
    checked += 1
assert checked == 55
print('All 55 previous non-store PR payload blobs remain byte-identical.')

frozen = git('ls-tree', '-r', '-z', CREATION, '--',
             'research/factor_base_evaluation_20261015').decode().split('\0')
count = 0
for line in filter(None, frozen):
    meta, path = line.split('\t', 1)
    mode, kind, sha = meta.split()
    assert kind == 'blob' and index[path] == (mode, sha), path
    count += 1
print(f'All {count} historical evaluation blobs present at creation remain unchanged.')

reg_path = 'docs/curves/registry.json'
registry = current(reg_path)
models = {r['slug']: r for r in registry['curves']}
aliases = current('src/cryptanalysis/curve_aliases.json')
store = current('store/index.json')['entries']
for ref in (HEAD, MAIN):
    for row in previous(ref, reg_path)['curves']:
        assert models[row['slug']] == row, (ref, row['slug'])
    old = previous(ref, 'src/cryptanalysis/curve_aliases.json')
    assert all(aliases['aliases'][k] == v for k, v in old['aliases'].items()), ref
    assert all(aliases[k] == v for k, v in old.items() if k != 'aliases'), ref
    assert all(store[k] == v for k, v in previous(ref, 'store/index.json')['entries'].items()), ref
print(f"Both parents preserved: {len(models)} models, {len(aliases['aliases'])} aliases, {len(store)} store entries.")

def measured(doc):
    doc = copy.deepcopy(doc)
    doc.pop('roster')
    doc['sources']['registry'].pop('sha256')
    return doc

board = current('docs/ic/leaderboard.json')
assert measured(board) == measured(previous(MAIN, 'docs/ic/leaderboard.json'))
print('Every non-roster leaderboard field equals the pinned main canonical input.')

browser = current('docs/browser/data.json')
keys = {'curves': 'slug', 'methods': 'method_id', 'factor_bases': 'fb_id',
        'candidates': 'candidate_id', 'rounds': 'round', 'sessions': 'dir'}
for ref in (HEAD, MAIN):
    old = previous(ref, 'docs/browser/data.json')
    for field, key in keys.items():
        now = {r[key]: r for r in browser[field]}
        assert all(r[key] in now for r in old[field]), (ref, field)
        if field in ('rounds', 'sessions'):
            assert all(now[r[key]] == r for r in old[field]), (ref, field)
    assert all(r in browser['yields'] for r in old['yields']), ref
print('All browser identities, exact sessions/rounds and yields from both parents retained.')
print('Final browser counts: ' + json.dumps(browser['counts'], sort_keys=True))

card = current('backlog/XP-20261016-prime-base-oracle-tuning/card.json')
assert card['status'] == 'draft'
assert card['cairn']['funding']['units'] == 0
assert card['protocol'] == {'path': None, 'blob': None}
assert card['checkpoint'] == {'file': None, 'resume_flag': None}
assert card['run']['argv'] == ['TODO-freeze-native-ecbench-spec-and-command']
assert card == previous(MAIN, 'backlog/XP-20261016-prime-base-oracle-tuning/card.json')
print('Planning card unchanged: draft, zero-funded, no runnable protocol/checkpoint.')
