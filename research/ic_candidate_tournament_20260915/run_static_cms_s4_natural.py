#!/usr/bin/env python3
"""One-shot blinded 32-query natural-yield stage for source-receipted CMS."""
import argparse
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess

from generic_query_law import probe_scalar
from generic_solver_feasibility import check_source_checkout
from oracle import require
from run_cms_s4_controls import (
    EXPORTER, build_exporter, digest, lift_source_assignment, meter,
    parse_cms_model, sha_bytes, validate_xor_dimacs,
)
from run_static_cms_s4_controls import (
    BUNDLE, PANEL as STATIC_PANEL, SCRIPTS, admit as admit_static,
    validate_export,
)
from tournament import read, write

HERE = Path(__file__).resolve().parent
REGISTRATION = HERE/'goal_20260924/static-cms-s4-natural'
PANEL = REGISTRATION/'panel.json'
PANEL_SHA256 = 'f0655044e6d18053872dd33cd1255706861820975fbacd66f3b7331471158c33'
SHARED = (
    'parent_exact_result_sha256', 'parent_recovery_archive_file',
    'parent_recovery_archive_sha256', 'source_commit', 'exporter_source_sha256',
    'cms_path', 'cms_executable_sha256', 'cms_build_receipt_sha256',
    'cms_build_bundle_seal_sha256', 'cms_build_verifier_sha256',
    'process_meter_sha256', 'cms_parser_sha256', 'curve_cell',
    'field_degree', 'curve_a', 'irreducible_low_terms', 'factor_base',
    'summands', 'source_encoding', 'export_only', 'export_timeout_seconds',
    'cms_timeout_seconds', 'cms_conflict_budget', 'cms_threads',
    'cms_random_seed', 'cms_max_models_per_query',
)


def admit(panel, *, require_local_binary=True):
    require(digest(PANEL) == PANEL_SHA256, 'natural SAT panel changed')
    parent = read(STATIC_PANEL)
    files, curve, base, built = admit_static(
        parent, require_local_binary=require_local_binary)
    require(panel['schema_version'] == 1
            and panel['status'] == 'REGISTERED_BEFORE_EXECUTION'
            and panel['algorithm_seed'] == 2026092935
            and panel['export_nonce'] == 2026092934
            and all(panel[key] == parent[key] for key in SHARED)
            and len(panel['schedule']) == 32,
            'natural SAT source, schedule or limits changed')
    points = set()
    scalars = set()
    for trial, item in enumerate(panel['schedule']):
        require(set(item) == {'trial', 'probe_scalar', 'point'}
                and item['trial'] == trial
                and item['probe_scalar'] == probe_scalar(
                    panel['algorithm_seed'], trial, curve.r)
                and tuple(item['point']) == curve.mul(
                    curve.g, item['probe_scalar']),
                'query differs from frozen independently replayed law')
        points.add(tuple(item['point']))
        scalars.add(item['probe_scalar'])
    require(len(points) == len(scalars) == 32,
            'natural query schedule repeats a point')
    old = read(HERE/'goal_20260924/generic-exact-yield-audit/RESULT.json')
    prior = next(row for row in old['rows']
                 if row['cell'] == 'n17a1' and row['solver'] == 'f5')
    require(not scalars.intersection(item['a'] for item in prior['queries']),
            'new natural sample reuses a previous ordinary query')
    return files, curve, base, built


def one_query(panel, item, exporter, cms, curve, base, out):
    trial = item['trial']
    directory = out/f'trial-{trial:02d}'
    directory.mkdir()
    instance = directory/'instance'
    command = [exporter, '17', '6', 'standard',
               str(panel['export_nonce']), str(panel['cms_conflict_budget']),
               instance, '1', '0', '--target-x', str(item['point'][0]),
               '--target-y', str(item['point'][1]),
               '--blind-instance-id', f'control-{trial:02d}', '--export-only']
    exported = meter(command, directory, 'export',
                     panel['export_timeout_seconds'])
    row = dict(trial=trial, probe_scalar=item['probe_scalar'],
               public_point=item['point'], exporter=exported,
               cms=None, source_model_valid=None,
               point_witness=None, status='EXPORT_FAILURE')
    if exported['timed_out'] or exported['returncode'] != 0:
        return row
    try:
        manifest = read(instance/'manifest.json')
        row['exports'] = validate_export(manifest, item, panel, base, instance)
        row['manifest_sha256'] = digest(instance/'manifest.json')
    except (OSError, ValueError, KeyError, TypeError) as error:
        row.update(status='INVALID_EXPORT', reason=f'{type(error).__name__}: {error}')
        return row
    command = [cms, '--verb', '1', '--threads', '1', '--random', '1',
               '--maxsol', '1', '--maxconfl', str(panel['cms_conflict_budget']),
               str(instance/'instance.xor.cnf')]
    measured = meter(command, directory, 'cms', panel['cms_timeout_seconds'])
    row['cms'] = measured
    stdout = (directory/'cms.stdout').read_text()
    if measured['timed_out']:
        row['status'] = 'TIMEOUT'
    elif measured['returncode'] == 10 and 's SATISFIABLE' in stdout:
        maximum = manifest['exports']['cryptominisat_xor_dimacs']['variables']
        model = parse_cms_model(stdout, maximum)
        if model is None or not validate_xor_dimacs(instance/'instance.xor.cnf', model):
            row['status'] = 'INVALID_SOURCE_MODEL'
        else:
            row['source_model_valid'] = True
            row['source_model_sha256'] = sha_bytes(bytes(model))
            witness = lift_source_assignment(model, manifest, curve, base,
                                             tuple(item['point']))
            row['point_witness'] = witness
            row['status'] = ('VALID_POINT_WITNESS' if witness['group_replay']
                             else 'SOURCE_MODEL_NONLIFTING')
    elif measured['returncode'] == 20 and 's UNSATISFIABLE' in stdout:
        row['status'] = 'SOURCE_UNSAT'
    elif measured['returncode'] == 0:
        row['status'] = 'UNKNOWN_INCONCLUSIVE'
    else:
        row['status'] = 'SOLVER_ERROR'
    return row


def run(panel, source_root, out):
    files, curve, base, built = admit(panel)
    check_source_checkout(source_root)
    require(digest(source_root/EXPORTER) == panel['exporter_source_sha256'],
            'pinned exporter changed')
    require(not out.exists(), 'output exists; frozen query may not be retried')
    out.mkdir(parents=True)
    for source, name in ((PANEL, 'registered-panel.json'),
                         (REGISTRATION/'PROTOCOL.md', 'PROTOCOL.md'),
                         (Path(__file__), 'registered-runner.py'),
                         (SCRIPTS/'process_meter.py', 'process_meter.py'),
                         (SCRIPTS/'run_koblitz_pdp_matrix.py', 'cms_parser.py'),
                         (BUNDLE/'builds/cryptominisat-receipt.json',
                          'cms-build-receipt.json'),
                         (BUNDLE/'bundle-seal.json',
                          'cms-build-bundle-seal.json')):
        shutil.copy2(source, out/name)
    write(out/'cms-build-verification.json', built, exclusive=True)
    write(out/'host.json', dict(system=platform.system(),
                               machine=platform.machine(),
                               cpu_count=os.cpu_count(),
                               scope='one-host fresh ordinary-query SAT stage'),
          exclusive=True)
    from run_generic_exact_yield_audit import record
    parent_source = record(files, 'build/source-manifest.json')
    write(out/'parent-source-manifest.json', parent_source, exclusive=True)
    cms = out/'cms-executable'
    shutil.copy2(panel['cms_path'], cms)
    require(digest(cms) == panel['cms_executable_sha256'],
            'copied static solver changed')
    linkage = subprocess.check_output(['otool', '-L', str(cms)], text=True)
    (out/'cms-linkage.txt').write_text(linkage)
    require('@rpath' not in linkage and '/opt/homebrew' not in linkage,
            'copied SAT solver has a nonportable library path')
    version = meter([cms, '--version'], out, 'cms-preflight', 10)
    if (version['returncode'] != 0 or version['timed_out']
            or 'CryptoMiniSat version 5.14.7'
                not in (out/'cms-preflight.stdout').read_text()):
        write(out/'summary.json', dict(schema_version=1,
              status='PREFLIGHT_FAILURE', panel_sha256=PANEL_SHA256,
              rows=[], cms_preflight=version, natural_yield_estimate=None,
              full_sat_ic_admission=False, online_speedup=None), exclusive=True)
        return
    build_dir = out/'build'
    build_dir.mkdir()
    exporter, build_record = build_exporter(panel, source_root, build_dir,
                                             parent_source)
    rows = []
    for item in panel['schedule']:
        row = one_query(panel, item, exporter, cms, curve, base, out)
        rows.append(row)
        with (out/'progress.jsonl').open('a') as stream:
            stream.write(json.dumps(row, sort_keys=True)+'\n')
        print(json.dumps(dict(trial=row['trial'], status=row['status'])),
              flush=True)
    write(out/'summary.json', dict(schema_version=1,
          status='NATURAL_STAGE_COMPLETE', panel_sha256=PANEL_SHA256,
          source_commit=panel['source_commit'],
          exporter_build=build_record,
          cms_build_receipt_sha256=panel['cms_build_receipt_sha256'],
          cms_binary_sha256=panel['cms_executable_sha256'],
          cms_preflight=version, rows=rows,
          natural_yield_estimate=None, full_sat_ic_admission=False,
          online_speedup=None, scope=panel['claim_boundary']), exclusive=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run(read(PANEL), args.source_root.resolve(), args.out.resolve())


if __name__ == '__main__':
    main()
