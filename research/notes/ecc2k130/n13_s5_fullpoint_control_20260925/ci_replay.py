#!/usr/bin/env python3
"""Hash/syntax/input-only S5 preflight; never calls a SAT solver or census."""
from __future__ import annotations

import ast
import hashlib
import importlib.util
import json
import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CORPUS = ROOT / 'research/notes/ecc2k130/rotated_pdp_corpus_20260925/evidence/raw.tar.gz'
PRIOR = ROOT / 'research/notes/ecc2k130/rotated_s3_candidate_20260925/evidence/n13-m5/producer/result.json'
BASIS = ROOT / 'research/notes/ecc2k130/rotated_subspace_support_20260925/gate.py'
FERMAT = ROOT / 'research/notes/ecc2k130/rotated_subspace_support_20260925/verify.py'


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise AssertionError(f'missing module {path}')
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def main() -> None:
    if sys.version_info < (3, 10):
        raise SystemExit('Python >=3.10 required')
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    assert frozen['schema'] == 'ecc2k130_n13_s5_fullpoint_control_v1'
    assert frozen['required_merged_prs'] == [802, 804, 807, 811]
    assert frozen['semantic_wall_seconds'] == 3600
    assert frozen['semantic_rss_bytes'] == 4 * (1 << 30)
    assert frozen['cnf_byte_max'] == 64 * (1 << 20)
    for relative, expected in frozen['source_sha256'].items():
        path = HERE / relative
        assert path.is_file() and sha(path) == expected, relative
        if path.suffix == '.py':
            ast.parse(path.read_text(), filename=str(path))
    for category in ('merged_parent_sha256', 'input_source_sha256', 'workflow_sha256'):
        for relative, expected in frozen[category].items():
            path = ROOT / relative
            assert path.is_file() and sha(path) == expected, relative
    data = json.loads((HERE / 'INPUT.json').read_text())
    assert data['schema'] == 'ecc2k130_n13_s5_control_input_v1'
    assert (data['n'], data['field_polynomial'], data['m'], data['factor_points_per_slot']) == (13, 0x201B, 5, 5)
    assert data['source_corpus_sha256'] == sha(CORPUS)
    assert data['source_prior_result_sha256'] == sha(PRIOR)
    with tarfile.open(CORPUS, 'r:gz') as archive:
        factors = json.load(archive.extractfile('raw/n13-m5/factors.json'))
    assert data['factors'] == [[[0, *p] for p in slot] for slot in factors]
    prior = json.loads(PRIOR.read_text())['targets']
    assert len(prior) == len(data['targets']) == 32
    positive = 0
    for i, (row, old) in enumerate(zip(data['targets'], prior, strict=True)):
        assert (row['index'], row['Q_index'], row['T_index']) == (i, i // 4, i % 4)
        assert row['id'] == f'Q{i//4}T{i%4}'
        assert row['full_target'] == [0, *old['target']]
        assert row['class'] == old['target_class']
        assert row['archived_exact_tuple_count'] == old['true_point_tuple_count']
        positive += row['archived_exact_tuple_count'] > 0
    assert positive == 5
    exceptional = data['exceptional']
    assert exceptional['target_index'] == 12 and exceptional['target_id'] == 'Q3T0'
    assert exceptional['full_target'] == data['targets'][12]['full_target'] == [0, 7256, 3272]
    assert exceptional['factor_indices'] == [0, 0, 0, 4, 1]
    assert exceptional['factor_points'] == [data['factors'][i][j] for i, j in enumerate(exceptional['factor_indices'])]
    assert exceptional['prefix_after_first_two'] == [1, 0, 0]
    # Reproduce the published #770 x-mask mapping from its frozen basis,
    # independent of the #774 candidate path table.
    gate = module(BASIS, 's5_basis_preflight')
    field = gate.Field(13, [0, 1, 3, 4])
    bases = gate.subspace_basis(gate.normal_conjugates(field, 3), 5, 2)
    masks = []
    for i, basis in enumerate(bases):
        x = exceptional['factor_points'][i][1]
        found = [mask for mask in range(4)
                 if x == ((basis[0] if mask & 1 else 0) ^
                          (basis[1] if mask & 2 else 0))]
        assert len(found) == 1
        masks.append(found[0])
    assert masks == exceptional['x_mask'] == [0, 0, 0, 2, 1]
    # Confirm the already-published signed witness and O prefix with a
    # separate Fermat law; this does not run a new solver/query outcome.
    mod = module(FERMAT, 's5_fermat_preflight')
    curve = mod.E(mod.GF(13, 0x201B))
    total = None
    for i, p in enumerate(exceptional['factor_points']):
        total = curve.add(total, (p[1], p[2]))
        if i == 1:
            assert total is None
    assert [0, *total] == exceptional['full_target']
    toy = json.loads((HERE / 'evidence/toy_selftest.json').read_text())
    assert toy['decision'] == 'TOY_SELFTEST_PASS'
    assert toy['primary_bits'] == 78 and toy['cnf']['clauses'] == 6147
    print(json.dumps({'decision': 'HASH_INPUT_SYNTAX_PASS',
                      'source_files': len(frozen['source_sha256']),
                      'branch_labels': len(data['targets']),
                      'nonzero_exact_branches': positive,
                      'n13_outcome': 'UNRUN'}, sort_keys=True))


if __name__ == '__main__':
    main()
