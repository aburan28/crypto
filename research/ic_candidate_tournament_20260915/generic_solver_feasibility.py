#!/usr/bin/env python3
"""Source-bound static layout preflight for the bounded generic IC worker.

This catches configurations the reviewed F4/F5 encoder cannot represent
*before* fresh target generation. Passing this check says nothing about
natural relation yield, runtime, or complete DLP success.
"""
import argparse
import json
from pathlib import Path
import re
import subprocess

from oracle import require

SOURCE_COMMIT = '765c3c5f19032bd852163805f257c56babef2040'
SOURCE_OBJECTS = {
    'examples/ic_tournament_worker.rs': '745247f3d89d46cdffb6d1f28778971caae7aa21',
    'src/cryptanalysis/koblitz_factor_base_search.rs': '8aacd0b9c0c4b4e8702c97d631319afdaa4e22e6',
    'src/cryptanalysis/koblitz_index_calculus.rs': '16a102a9be61630da2d1c49bcabcc55eaf5ab2f0',
    'src/cryptanalysis/koblitz_groebner.rs': '91ee5d6cfd99166697895aac62dbabb1f4d9778c',
    'src/cryptanalysis/polynomial_reuse.rs': 'f79411c0f2309090b4ba7ddfd9b80bb202859796',
}
MAX_VARS = 64
ALGEBRAIC = {'f4', 'f5', 'inherited_f4'}
SAT = {'sat_xor', 'sat_cnf'}
CELL = re.compile(r'n([1-9][0-9]*)a[01]\Z')


def check_source(repo):
    """Confirm the reviewed revision's source objects exist unchanged."""
    for path, expected in SOURCE_OBJECTS.items():
        actual = subprocess.check_output(
            ['git', 'rev-parse', f'{SOURCE_COMMIT}:{path}'], cwd=repo, text=True).strip()
        require(actual == expected, f'unreviewed generic source object: {path}')


def check_source_checkout(repo):
    """Preflight must inspect the exact checkout that will build the worker."""
    check_source(repo)
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                   cwd=repo, text=True).strip()
    require(head == SOURCE_COMMIT, 'generic source checkout is not the reviewed worker')
    for path, expected in SOURCE_OBJECTS.items():
        actual = subprocess.check_output(['git', 'hash-object', path],
                                         cwd=repo, text=True).strip()
        require(actual == expected, f'generic worker checkout has changed: {path}')
    status = subprocess.check_output(
        ['git', 'status', '--porcelain', '--untracked-files=all', '--',
         'src', 'examples/ic_tournament_worker.rs', 'Cargo.toml'],
        cwd=repo, text=True).strip()
    require(not status, 'generic worker source checkout has uncommitted changes')


def basis_length(recipe, degree):
    require(type(recipe) is dict, 'algebraic candidate has no explicit base policy')
    require(('recipe' in recipe) != ('kind' in recipe), 'ambiguous base policy')
    family = recipe.get('recipe', recipe.get('kind'))
    if family == 'subgroup_orbits':
        if 'recipe' in recipe:
            require(set(recipe) == {'recipe', 'seed', 'requested_points'},
                    'unknown sampled-orbit policy')
            requested = recipe['requested_points']
            require((type(requested) is int and requested > 0)
                    or (type(requested) is str
                        and re.fullmatch(r'[1-9][0-9]*n', requested)),
                    'invalid sampled-orbit point request')
        else:
            require(set(recipe) == {'kind', 'seed', 'points'}
                    and type(recipe['points']) is int and recipe['points'] > 0,
                    'unknown sampled-orbit policy')
        require(type(recipe['seed']) is int and 0 <= recipe['seed'] < 1 << 64,
                'invalid sampled-orbit seed')
        # build_explicit_frobenius_orbit_factor_base creates the full ambient
        # polynomial basis, regardless of its sampled point count.
        return degree
    if family == 'standard_subspace':
        require(set(recipe) == {('recipe' if 'recipe' in recipe else 'kind'), 'dimension'}
                and type(recipe['dimension']) is int
                and 1 <= recipe['dimension'] <= 20
                and recipe['dimension'] < degree,
                'invalid standard-subspace policy')
        return recipe['dimension']
    require(False, 'unreviewed algebraic factor-base policy')


def assess(panel):
    require(type(panel) is dict and type(panel.get('cells')) is list
            and type(panel.get('candidates')) is list, 'invalid candidate panel')
    cells = panel['cells']
    require(cells and len(set(cells)) == len(cells)
            and all(type(c) is str and CELL.fullmatch(c) for c in cells),
            'invalid or duplicate curve cells')
    require(panel['candidates'] and all(type(candidate) is dict
            and type(candidate.get('id')) is str and candidate['id']
            for candidate in panel['candidates']), 'invalid candidate IDs')
    ids = [candidate['id'] for candidate in panel['candidates']]
    require(len(set(ids)) == len(ids),
            'missing or duplicate candidate IDs')
    rows = []
    algebraic_arms = 0
    for candidate in panel['candidates']:
        require(type(candidate) is dict and type(candidate.get('id')) is str
                and type(candidate.get('config')) is dict,
                'invalid candidate declaration')
        solver = candidate['config'].get('solver')
        if solver not in ALGEBRAIC:
            rows.append(dict(candidate=candidate['id'], solver=solver,
                             status='SEPARATE_SAT_PATH_UNCHECKED' if solver in SAT
                                    else 'OUTSIDE_ALGEBRAIC_GATE'))
            continue
        m = candidate['config'].get('summands')
        require(type(m) is int and 2 <= m <= 4, 'invalid algebraic summand count')
        algebraic_arms += 1
        recipe = candidate['config'].get('factor_base')
        for cell in cells:
            n = int(CELL.fullmatch(cell).group(1))
            require(5 <= n <= 31 and n % 2 == 1, 'outside reviewed worker field range')
            ell = basis_length(recipe, n)
            variables = m * ell + (m - 2) * n
            rows.append(dict(candidate=candidate['id'], cell=cell, solver=solver,
                             summands=m, basis_length=ell, boolean_variables=variables,
                             max_boolean_variables=MAX_VARS,
                             status='ABOVE_LAYOUT_CAP' if variables > MAX_VARS
                                    else 'WITHIN_LAYOUT_CAP_ONLY'))
    require(algebraic_arms > 0, 'no algebraic arm to preflight')
    impossible = [row for row in rows if row['status'] == 'ABOVE_LAYOUT_CAP']
    return dict(schema_version=1, source_commit=SOURCE_COMMIT,
                source_objects=SOURCE_OBJECTS, layout='m*ell+(m-2)*n',
                max_boolean_variables=MAX_VARS,
                status='FAIL_STATIC_LAYOUT' if impossible else 'PASS_STATIC_LAYOUT_ONLY',
                impossible_algebraic_cells=len(impossible), rows=rows,
                interpretation='Static encoder check only; dispatch, natural yield, '
                               'complete recovery and cost remain unmeasured.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--panel', type=Path, required=True)
    parser.add_argument('--source-repo', type=Path,
                        default=Path(__file__).resolve().parents[2])
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--require-pass', action='store_true',
                        help='abort premeasurement if any algebraic cell exceeds the cap')
    args = parser.parse_args()
    check_source_checkout(args.source_repo)
    result = assess(json.loads(args.panel.read_text()))
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(f"{result['status']}: {result['impossible_algebraic_cells']} impossible cells")
    if args.require_pass:
        require(result['status'] == 'PASS_STATIC_LAYOUT_ONLY',
                'algebraic encoder cannot represent registered panel')


if __name__ == '__main__':
    main()
