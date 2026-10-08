#!/usr/bin/env python3
"""Read-only exact-group and source-model audit of the frozen SAT sample."""
import argparse
from collections import Counter
import hashlib
import io
from itertools import product
import json
import math
from pathlib import Path
import tarfile

from generic_bases import lifts
from oracle import require
from run_generic_exact_yield_audit import exact_three_sum, pair_index
from run_static_cms_s4_natural import PANEL, PANEL_SHA256, admit
from tournament import read, write


def digest(data):
    return hashlib.sha256(data).hexdigest()


def archive_files(path):
    data = path.read_bytes()
    with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as archive:
        members = archive.getmembers()
        names = [member.name for member in members]
        require(len(names) == len(set(names))
                and all(member.isfile() and not member.name.startswith('/')
                        and '..' not in Path(member.name).parts for member in members),
                'unsafe or duplicate natural evidence member')
        return digest(data), {member.name: archive.extractfile(member).read()
                              for member in members}


def complete_model(stdout, count):
    values = {}
    for line in stdout.decode().splitlines():
        if not line.startswith('v '):
            continue
        for token in line.split()[1:]:
            literal = int(token)
            if literal:
                index = abs(literal)
                require(1 <= index <= count
                        and (index not in values or values[index] == (literal > 0)),
                        'invalid, conflicting or extra SAT literal')
                values[index] = literal > 0
    require(set(values) == set(range(1, count+1)),
            'SAT assignment is not complete')
    return [values[index] for index in range(1, count+1)]


def check_formula(data, model):
    header = None
    clauses = []
    for line in data.decode().splitlines():
        line = line.strip()
        if not line or line.startswith('c'):
            continue
        if line.startswith('p '):
            parts = line.split()
            require(header is None and parts[:2] == ['p', 'cnf']
                    and len(parts) == 4, 'invalid XOR-DIMACS header')
            header = int(parts[2]), int(parts[3])
        else:
            clauses.append(line)
    require(header == (len(model), len(clauses)),
            'XOR-DIMACS variable or row count mismatch')
    for line in clauses:
        xor = line.startswith('x ')
        literals = [int(value) for value in line.split()[1 if xor else 0:]]
        require(len(literals) > 1 and literals[-1] == 0
                and all(1 <= abs(literal) <= len(model)
                        for literal in literals[:-1]),
                'invalid XOR-DIMACS clause')
        truth = [(model[abs(literal)-1] if literal > 0
                  else not model[abs(literal)-1]) for literal in literals[:-1]]
        require((sum(truth) % 2 == 1) if xor else any(truth),
                'SAT assignment does not satisfy the exported source system')


def group_lift(curve, base, model, point):
    xs = [sum(1 << bit for bit in range(6) if model[6*slot+bit])
          for slot in range(3)]
    allowed = set(base)
    for triple in product(*(tuple(p for p in lifts(curve, x) if p in allowed)
                            for x in xs)):
        if curve.add(curve.add(triple[0], triple[1]), triple[2]) == point:
            return xs, triple
    return xs, None


def wilson(success, count):
    z = 1.959963984540054
    p = success/count
    denominator = 1 + z*z/count
    center = (p + z*z/(2*count))/denominator
    radius = z*math.sqrt(p*(1-p)/count + z*z/(4*count*count))/denominator
    return [format(max(0, center-radius), '.6f'),
            format(min(1, center+radius), '.6f')]


def audit(archive):
    panel = read(PANEL)
    _, curve, base, built = admit(panel, require_local_binary=False)
    archive_sha256, files = archive_files(archive)
    require(files['registered-panel.json'] == PANEL.read_bytes()
            and digest(files['registered-panel.json']) == PANEL_SHA256
            and digest(files['cms-executable']) == panel['cms_executable_sha256']
            and digest(files['cms-build-receipt.json'])
                == panel['cms_build_receipt_sha256']
            and digest(files['cms-build-bundle-seal.json'])
                == panel['cms_build_bundle_seal_sha256']
            and '@rpath' not in files['cms-linkage.txt'].decode(),
            'registered source, binary, receipt or linkage changed')
    preflight = json.loads(files['cms-preflight.metrics.json'])
    require(preflight['returncode'] == 0 and not preflight['timed_out']
            and b'CryptoMiniSat version 5.14.7' in files['cms-preflight.stdout'],
            'copied static solver did not start in its final location')
    result = json.loads(files['summary.json'])
    require(result['status'] == 'NATURAL_STAGE_COMPLETE'
            and result['panel_sha256'] == PANEL_SHA256
            and result['cms_preflight'] == preflight
            and result['cms_binary_sha256'] == panel['cms_executable_sha256']
            and result['full_sat_ic_admission'] is False
            and result['natural_yield_estimate'] is None
            and result['online_speedup'] is None
            and len(result['rows']) == len(panel['schedule']) == 32,
            'natural stage summary is incomplete or overclaims')
    require([json.loads(line) for line in
             files['progress.jsonl'].decode().splitlines()] == result['rows'],
            'progress stream differs from final rows')
    pairs = pair_index(curve, base)
    statuses = Counter()
    rows = []
    total_export_wall = total_cms_wall = 0.0
    for item, measured in zip(panel['schedule'], result['rows']):
        trial = item['trial']
        prefix = f'trial-{trial:02d}/'
        require(measured['trial'] == trial
                and measured['probe_scalar'] == item['probe_scalar']
                and measured['public_point'] == item['point'],
                'ordinary public query changed')
        exact = exact_three_sum(curve, base, pairs, tuple(item['point']))
        export = json.loads(files[prefix+'export.metrics.json'])
        require(export == measured['exporter'], 'export process receipt changed')
        total_export_wall += export['metrics']['wall_seconds']
        status = measured['status']
        statuses[status] += 1
        cms = measured['cms']
        if cms is not None:
            require(cms == json.loads(files[prefix+'cms.metrics.json']),
                    'SAT process receipt changed')
            total_cms_wall += cms['metrics']['wall_seconds']
        if status not in {'EXPORT_FAILURE', 'INVALID_EXPORT'}:
            manifest = json.loads(files[prefix+'instance/manifest.json'])
            require(digest(files[prefix+'instance/manifest.json'])
                    == measured['manifest_sha256']
                    and manifest['representation'] == 'symmetrised_s4'
                    and [int(manifest['target'][key]) for key in ('x', 'y')]
                        == item['point']
                    and manifest['factor_base_geometry']['curve_points'] == len(base),
                    'source manifest or target differs')
            for name, descriptor in manifest['exports'].items():
                body = files[prefix+'instance/'+descriptor['path']]
                require(measured['exports'][name]
                        == {'bytes': len(body), 'sha256': digest(body)},
                        'source export bytes differ')
            stdout = files[prefix+'cms.stdout']
            if status == 'VALID_POINT_WITNESS':
                require(cms['returncode'] == 10 and not cms['timed_out']
                        and b's SATISFIABLE' in stdout,
                        'witness lacks a SAT verdict')
                count = manifest['exports']['cryptominisat_xor_dimacs']['variables']
                model = complete_model(stdout, count)
                check_formula(files[prefix+'instance/instance.xor.cnf'], model)
                xs, lifted = group_lift(curve, base, model, tuple(item['point']))
                require(lifted is not None
                        and exact is not None
                        and measured['source_model_valid'] is True
                        and measured['source_model_sha256'] == digest(bytes(model))
                        and measured['point_witness']['x_coordinates'] == xs
                        and measured['point_witness']['points']
                            == [list(point) for point in lifted]
                        and measured['point_witness']['point_indices']
                            == [base.index(point) for point in lifted],
                        'SAT witness fails formula or independent group replay')
            elif status == 'SOURCE_UNSAT':
                require(cms['returncode'] == 20 and not cms['timed_out']
                        and b's UNSATISFIABLE' in stdout and exact is None,
                        'reported source UNSAT contradicts exact group oracle')
            elif status == 'TIMEOUT':
                require(cms['timed_out'], 'timeout status lacks watchdog evidence')
            elif status == 'SOURCE_MODEL_NONLIFTING':
                require(cms['returncode'] == 10 and not cms['timed_out']
                        and b's SATISFIABLE' in stdout,
                        'nonlifting model lacks a SAT verdict')
                count = manifest['exports']['cryptominisat_xor_dimacs']['variables']
                model = complete_model(stdout, count)
                check_formula(files[prefix+'instance/instance.xor.cnf'], model)
                _, lifted = group_lift(curve, base, model, tuple(item['point']))
                require(lifted is None and measured['source_model_valid'] is True,
                        'nonlifting model actually contains a group witness')
            else:
                require(status in {'SOLVER_ERROR', 'UNKNOWN_INCONCLUSIVE',
                                   'INVALID_SOURCE_MODEL'},
                        'unknown solver status')
        rows.append(dict(trial=trial, probe_scalar=item['probe_scalar'],
                         public_point=item['point'], exact_feasible=exact is not None,
                         exact_witness_indices=exact, measured_status=status,
                         export_wall_seconds=format(export['metrics']['wall_seconds'], '.6f'),
                         cms_wall_seconds=None if cms is None else
                             format(cms['metrics']['wall_seconds'], '.6f')))
    found = statuses['VALID_POINT_WITNESS']
    feasible = sum(row['exact_feasible'] for row in rows)
    require(found <= feasible and sum(statuses.values()) == 32,
            'observed witnesses exceed exact feasible queries')
    return dict(schema_version=1, status='AUDITED_NATURAL_STAGE',
                panel_sha256=PANEL_SHA256,
                evidence_archive_sha256=archive_sha256,
                summary_sha256=digest(files['summary.json']),
                runner_sha256=digest(files['registered-runner.py']),
                source_build_bundle_status=built['status'],
                statuses=dict(sorted(statuses.items())),
                attempts=32, exact_feasible=feasible,
                verified_witnesses=found,
                feasible_missed=feasible-found,
                witness_rate_wilson95=wilson(found, 32),
                exact_feasible_rate_wilson95=wilson(feasible, 32),
                summed_export_process_wall_seconds=format(total_export_wall, '.6f'),
                summed_cms_process_wall_seconds=format(total_cms_wall, '.6f'),
                rows=rows, full_sat_ic_admission=False,
                online_speedup=None,
                scope='32 seeded ordinary queries only; no full-rank relation matrix or target DLP')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    require(not args.out.exists(), 'audit output already exists')
    write(args.out, audit(args.archive), exclusive=True)


if __name__ == '__main__':
    main()
