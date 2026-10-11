#!/usr/bin/env python3
"""Reproduce a source-derived admission screen; never constructs SAT or reads S3."""
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
STUDY = ROOT / 'research/koblitz_n83_factor_base_sweep_20261008'
SOURCE = ROOT / 'src/cryptanalysis/binary_semaev_chain_sat.rs'
MANIFESTS = [
    STUDY / 'pilot-01/manifest.json',
    STUDY / 'verification/two-hour-extension-20261010/v2-object/manifest.json',
    STUDY / 'verification/remaining-budget-k2048-20261010/manifest.json',
    STUDY / 'verification/structured-d12-orbit-20261010/manifest.json',
]
ORDER = 2417851639230796216685689
N = 83
SIZES = [64, 256, 600, 1182, 2001, 2048, 4096, 8192, 16627]
VAR_CAP = 150000
DOMAIN_CLAUSE_CAP = 1000000


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_counts(k, m):
    # Transcribed from chain_s3_max_variables and frobenius_orbit_domain_size.
    bits = (N - 1).bit_length()
    core = N * (2 * m - 1) + 4 * N * (m - 1) + 3 * N * N * (m - 1)
    domain_vars = m * (2 * k - 1 + N + bits + 2 * N * bits)
    domain_clauses = m * (k * N + 3 * k - 3 + (1 << bits) - N + 4 * N * bits)
    domain_xors = m * N * (bits + 1)
    return dict(core_variables=core, domain_variables=domain_vars,
                total_variables=core + domain_vars, domain_clauses=domain_clauses,
                domain_xor_rows=domain_xors,
                within_150k_variable_and_1m_domain_clause_caps=(core + domain_vars <= VAR_CAP and domain_clauses <= DOMAIN_CLAUSE_CAP))


def main():
    source = SOURCE.read_text()
    assert 'pub fn chain_s3_max_variables' in source
    assert 'pub fn frobenius_orbit_domain_size' in source
    assert source_counts(2001, 5) == dict(core_variables=84743,
        domain_variables=26265, total_variables=111008,
        domain_clauses=872260, domain_xor_rows=3320,
        within_150k_variable_and_1m_domain_clause_caps=True)
    assert source_counts(2048, 5)['total_variables'] == 111478
    assert source_counts(2048, 5)['domain_clauses'] == 892470
    assert source_counts(1182, 6)['core_variables'] == 105908
    objects = []
    for path in MANIFESTS:
        data = json.loads(path.read_text())
        bases = data['bases'] if 'bases' in data else [data]
        for base in bases:
            if base.get('a', 0) != 0:
                continue
            k = base.get('columns', base.get('orbit_columns'))
            points = base.get('points', base.get('point_records'))
            assert k in SIZES and points == 166 * k
            assert base.get('s3_uri', '').startswith('s3://crypto-autoresearcher/factor-bases/icv1/etc/')
            objects.append(dict(columns=k, points=points,
                policy=base.get('policy', 'structured_d12'),
                compressed_blake3=base['compressed_blake3'],
                manifest=str(path.relative_to(ROOT))))
    rows = []
    for k in SIZES:
        stored = [obj for obj in objects if obj['columns'] == k]
        arities = {}
        for m in (5, 6):
            count = source_counts(k, m)
            numerator = math.comb(166 * k + m - 1, m)
            count['uniform_target_mean_multiset_count'] = f'{numerator / ORDER:.8g}'
            arities[str(m)] = count
        rows.append(dict(columns=k, point_count_if_constructed=166 * k,
                         stored_primary_objects=len(stored),
                         stored_object_hashes=[obj['compressed_blake3'] for obj in stored],
                         arities=arities))
    receipt = dict(schema='n83.orbit-domain-static-screen/v1',
                   measurement_scope='source and manifest arithmetic only; zero SAT constructions, searches or new S3 transfers',
                   inputs={str(p.relative_to(ROOT)): sha256(p) for p in [SOURCE, *MANIFESTS]},
                   subgroup_order=str(ORDER), field_bits=N,
                   capacity_limits=dict(variables=VAR_CAP, domain_clauses=DOMAIN_CLAUSE_CAP),
                   stored_primary_object_count=len(objects), rows=rows)
    path = Path(__file__).with_name('static-screen.json')
    path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(f'{path}: {len(objects)} stored primary objects; {len(rows)} sizes; 2 arities')


if __name__ == '__main__':
    main()
