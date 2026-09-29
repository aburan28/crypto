#!/usr/bin/env python3
"""Run the four frozen, disclosed wide-S4/CryptoMiniSat correctness controls."""
import argparse
import hashlib
from itertools import product
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

from generic_bases import lifts
from generic_build import command, controlled_environment, digest, source_manifest
from generic_solver_feasibility import check_source_checkout
from oracle import Curve, require
from run_generic_exact_yield_audit import (
    PANEL as EXACT_PANEL, audit as audit_exact, load_evidence, record,
)
from tournament import read, write

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCRIPTS = ROOT / 'scripts'
REGISTRATION = HERE / 'goal_20260924/cms-s4-controls'
PANEL = REGISTRATION / 'panel.json'
PANEL_SHA256 = '1735b401af3ea87a78a380c5cc2bae7eb12eb39117ae7b6038346fc810e44c23'
EXPORTER = 'examples/koblitz_pdp_export.rs'

sys.path.insert(0, str(SCRIPTS))
from run_koblitz_pdp_matrix import parse_cms_model, validate_xor_dimacs  # noqa: E402


def sha_bytes(data):
    return hashlib.sha256(data).hexdigest()


def preflight(panel, *, require_local_solver=True):
    require(digest(PANEL) == PANEL_SHA256, 'registered SAT control panel changed')
    require(panel['schema_version'] == 1
            and panel['status'] == 'REGISTERED_BEFORE_EXECUTION'
            and panel['source_commit'] == '765c3c5f19032bd852163805f257c56babef2040'
            and panel['curve_cell'] == 'n19a0'
            and panel['field_degree'] == 19 and panel['curve_a'] == 0
            and panel['irreducible_low_terms'] == [0, 1, 2, 5]
            and panel['factor_base'] == {'kind': 'standard_subspace', 'dimension': 6}
            and panel['summands'] == 3
            and panel['export_only'] is True
            and panel['export_nonce'] == 2026092931
            and panel['export_timeout_seconds'] == 60
            and panel['cms_timeout_seconds'] == 120
            and panel['cms_conflict_budget'] == 1000000
            and panel['cms_threads'] == panel['cms_random_seed'] == panel['cms_max_models_per_query'] == 1
            and [row['trial'] for row in panel['schedule']] == [0, 1, 3, 10],
            'registered SAT controls or resource policy changed')
    require(digest(SCRIPTS/'process_meter.py') == panel['process_meter_sha256']
            and digest(SCRIPTS/'run_koblitz_pdp_matrix.py') == panel['cms_parser_sha256'],
            'meter or SAT model parser changed')
    # Execution requires the registered local binary. Portable evidence replay
    # instead verifies its sealed archived bytes (and can run on Linux CI).
    if require_local_solver:
        cms = Path(panel['cms_path'])
        require(cms.is_file() and digest(cms) == panel['cms_executable_sha256'],
                'external SAT executable changed')
    exact_panel = read(EXACT_PANEL)
    require(digest(EXACT_PANEL.parent/'RESULT.json') == panel['parent_exact_result_sha256']
            and digest(HERE/panel['parent_recovery_archive_file'])
                == panel['parent_recovery_archive_sha256'],
            'parent exact or recovery evidence changed')
    exact = read(EXACT_PANEL.parent/'RESULT.json')
    replay = audit_exact(exact_panel)
    replay['runner_sha256'] = exact['runner_sha256']
    require(replay == exact, 'original exact-yield/source admission failed on replay')
    files = load_evidence(exact_panel)
    report = record(files, 'jobs/n19a0/f5/stdout.json')
    fixture = report['fixture']
    curve = Curve(fixture)
    require(fixture['degree'] == panel['field_degree']
            and fixture['curve_a'] == panel['curve_a']
            and fixture['irreducible']['low_terms'] == panel['irreducible_low_terms'],
            'registered curve representation changed')
    base = tuple(curve.decode(p) for p in report['factor_base'])
    require(len(base) == 65 and len(set(base)) == 65, 'geometric base changed')
    source_row = next(row for row in exact['rows']
                      if row['cell'] == 'n19a0' and row['solver'] == 'f5')
    for item in panel['schedule']:
        label = source_row['queries'][item['trial']]
        require(label['trial'] == item['trial']
                and label['a'] == item['probe_scalar']
                and label['exact_feasible'] is item['exact_relation_exists']
                and curve.mul(curve.g, item['probe_scalar']) == tuple(item['point']),
                'registered point, query or exact label changed')
    require([item['exact_relation_exists'] for item in panel['schedule']]
            == [False, False, True, True], 'registered control classes changed')
    return files, curve, base


def meter(command_line, directory, name, seconds):
    prefix = directory / name
    command = [sys.executable, str(SCRIPTS/'process_meter.py'),
               '--cwd', str(directory), '--timeout', str(seconds),
               '--stdout', str(prefix)+'.stdout', '--stderr', str(prefix)+'.stderr',
               '--metrics', str(prefix)+'.metrics.json', '--exclusive-create',
               '--', *map(str, command_line)]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    require(completed.returncode == 0,
            f'process meter failed for {name}: {completed.stderr}')
    return read(str(prefix)+'.metrics.json')


def build_exporter(panel, source_root, out, parent_source):
    check_source_checkout(source_root)
    require(digest(source_root/EXPORTER) == panel['exporter_source_sha256'],
            'exporter source changed')
    environment = controlled_environment(source_root)
    metadata = json.loads(command(['cargo', 'metadata', '--locked', '--offline',
        '--format-version', '1'], source_root, environment))
    require(source_manifest(source_root, metadata) == parent_source,
            'library/dependency source differs from admitted worker build')
    args = ['cargo', 'build', '--locked', '--offline', '--release',
            '--no-default-features', '--example', 'koblitz_pdp_export']
    started = time.monotonic_ns()
    with (out/'build.log').open('x') as log:
        process = subprocess.run(args, cwd=source_root, env=environment,
                                 stdout=log, stderr=subprocess.STDOUT, check=False)
    build_ns = time.monotonic_ns()-started
    write(out/'build-exit.json', dict(exit_code=process.returncode,
                                     diagnostic_wall_ns=build_ns), exclusive=True)
    require(process.returncode == 0, 'pinned exporter build failed; retain build log')
    check_source_checkout(source_root)
    require(source_manifest(source_root, metadata) == parent_source
            and digest(source_root/EXPORTER) == panel['exporter_source_sha256'],
            'source changed during exporter build')
    executable = Path(metadata['target_directory'])/'release/examples/koblitz_pdp_export'
    retained = out/'exporter'
    shutil.copy2(executable, retained)
    record_out = dict(schema_version=1, source_commit=panel['source_commit'],
                      parent_source_manifest_sha256=sha_bytes(json.dumps(parent_source,
                          sort_keys=True, separators=(',', ':')).encode()),
                      exporter_source_sha256=panel['exporter_source_sha256'],
                      exporter_executable_sha256=digest(retained),
                      build_command=args, build_log_sha256=digest(out/'build.log'),
                      rustc=command(['rustc', '-Vv'], source_root, environment),
                      cargo=command(['cargo', '-Vv'], source_root, environment),
                      platform=dict(system=platform.system(), machine=platform.machine()),
                      diagnostic_build_wall_ns=build_ns,
                      scope='controlled local source/dependency build; external CMS binary separate')
    write(out/'build-record.json', record_out, exclusive=True)
    return retained, record_out


def validate_export(manifest, item, panel, base, instance):
    require(manifest['kind'] == 'binary_koblitz_pdp_cross_solver_instance'
            and manifest['target_mode'] == 'explicit_affine'
            and manifest['n'] == 19 and manifest['ell'] == 6 and manifest['m'] == 3
            and manifest['curve_a'] == 0
            and manifest['irreducible_low_terms'] == panel['irreducible_low_terms']
            and manifest['representation'] == 'symmetrised_s4'
            and manifest['factor_base_basis_bitmasks'] == ['1', '2', '4', '8', '16', '32']
            and manifest['factor_base_geometry']['curve_points'] == len(base)
            and [int(manifest['target'][k]) for k in ('x', 'y')] == item['point']
            and manifest['seed'] == panel['export_nonce']
            and manifest['blind_instance_id'] == f"control-{item['trial']:02d}"
            and manifest['native_sat']['status'] == 'not_run_in_export_process'
            and manifest['mitm']['status'] == 'not_run_in_export_process',
            'exported instance differs from registered curve, base, point or source mode')
    exports = manifest['exports']
    require(set(exports) == {'wdsat_anf', 'cryptominisat_xor_dimacs', 'magma_boolean_f4'},
            'missing or extra source export')
    digests = {}
    for name, descriptor in exports.items():
        path = instance/descriptor['path']
        require(path.is_file() and path.parent == instance
                and path.stat().st_size == descriptor['bytes'],
                'missing, unsafe or changed source export')
        digests[name] = dict(bytes=path.stat().st_size, sha256=digest(path))
    return digests


def lift_source_assignment(model, manifest, curve, base, target):
    ell = manifest['ell']
    basis = [int(value) for value in manifest['factor_base_basis_bitmasks']]
    xs = [0, 0, 0]
    for i in range(3):
        for j in range(ell):
            if model[i*ell+j]:
                xs[i] ^= basis[j]
    allowed = set(base)
    choices = [tuple(point for point in lifts(curve, x) if point in allowed) for x in xs]
    for triple in product(*choices):
        if curve.add(curve.add(triple[0], triple[1]), triple[2]) == target:
            index_of = {p: i for i, p in enumerate(base)}
            return dict(x_coordinates=xs, point_indices=[index_of[p] for p in triple],
                        points=[list(p) for p in triple], group_replay=True)
    return dict(x_coordinates=xs, point_indices=None, points=None, group_replay=False)


def one_control(panel, item, exporter, cms, curve, base, out):
    trial = item['trial']
    directory = out/f'trial-{trial:02d}'
    directory.mkdir()
    instance = directory/'instance'
    export_command = [exporter, '19', '6', 'standard',
                      str(panel['export_nonce']), str(panel['cms_conflict_budget']),
                      instance, '0', '0', '--target-x', str(item['point'][0]),
                      '--target-y', str(item['point'][1]),
                      '--blind-instance-id', f'control-{trial:02d}', '--export-only']
    export = meter(export_command, directory, 'export', panel['export_timeout_seconds'])
    row = dict(trial=trial, probe_scalar=item['probe_scalar'],
               exact_relation_exists=item['exact_relation_exists'],
               public_point=item['point'], exporter=export,
               cms=None, source_model_valid=None, point_witness=None,
               status='EXPORT_FAILURE')
    if export['timed_out'] or export['returncode'] != 0:
        return row
    try:
        manifest = read(instance/'manifest.json')
        row['exports'] = validate_export(manifest, item, panel, base, instance)
        row['manifest_sha256'] = digest(instance/'manifest.json')
    except (OSError, ValueError, KeyError, TypeError) as error:
        row.update(status='INVALID_EXPORT', reason=f'{type(error).__name__}: {error}')
        return row
    cms_command = [cms, '--verb', '1', '--threads', '1', '--random', '1',
                   '--maxsol', '1',
                   '--maxconfl', str(panel['cms_conflict_budget']),
                   str(instance/'instance.xor.cnf')]
    run = meter(cms_command, directory, 'cms', panel['cms_timeout_seconds'])
    row['cms'] = run
    stdout = (directory/'cms.stdout').read_text()
    if run['timed_out']:
        row['status'] = 'TIMEOUT'
    elif run['returncode'] == 10 and 's SATISFIABLE' in stdout:
        max_var = manifest['exports']['cryptominisat_xor_dimacs']['variables']
        model = parse_cms_model(stdout, max_var)
        if model is None or not validate_xor_dimacs(instance/'instance.xor.cnf', model):
            row['status'] = 'INVALID_SOURCE_MODEL'
        else:
            row['source_model_valid'] = True
            row['source_model_sha256'] = sha_bytes(bytes(model))
            witness = lift_source_assignment(model, manifest, curve, base,
                                             tuple(item['point']))
            row['point_witness'] = witness
            row['status'] = ('CONTRADICTS_EXACT_NEGATIVE' if witness['group_replay']
                             and not item['exact_relation_exists'] else
                             'VALID_POINT_WITNESS' if witness['group_replay'] else
                             'SOURCE_MODEL_NONLIFTING')
    elif run['returncode'] == 20 and 's UNSATISFIABLE' in stdout:
        row['status'] = ('CONTRADICTS_EXACT_POSITIVE'
                         if item['exact_relation_exists'] else 'SOURCE_UNSAT')
    elif run['returncode'] == 0:
        row['status'] = 'UNKNOWN_INCONCLUSIVE'
    else:
        row['status'] = 'SOLVER_ERROR'
    return row


def run(panel, source_root, out):
    files, curve, base = preflight(panel)
    check_source_checkout(source_root)
    require(digest(source_root/EXPORTER) == panel['exporter_source_sha256'],
            'pinned exporter source changed')
    require(not out.exists(), 'output exists; frozen controls cannot be rerun or overwritten')
    out.mkdir(parents=True)
    shutil.copy2(PANEL, out/'registered-panel.json')
    shutil.copy2(REGISTRATION/'PROTOCOL.md', out/'PROTOCOL.md')
    shutil.copy2(Path(__file__), out/'registered-runner.py')
    shutil.copy2(SCRIPTS/'process_meter.py', out/'process_meter.py')
    shutil.copy2(SCRIPTS/'run_koblitz_pdp_matrix.py', out/'cms_parser.py')
    write(out/'host.json', dict(system=platform.system(), machine=platform.machine(),
                                cpu_count=os.cpu_count(),
                                scope='local binary-bound diagnostic host'), exclusive=True)
    parent_source = record(files, 'build/source-manifest.json')
    write(out/'parent-source-manifest.json', parent_source, exclusive=True)
    build_dir = out/'build'
    build_dir.mkdir()
    exporter, build_record = build_exporter(panel, source_root, build_dir, parent_source)
    cms = Path(panel['cms_path'])
    shutil.copy2(cms, out/'cms-executable')
    require(digest(out/'cms-executable') == panel['cms_executable_sha256'],
            'copied external SAT binary changed')
    rows = []
    for item in panel['schedule']:
        row = one_control(panel, item, exporter, out/'cms-executable', curve, base, out)
        rows.append(row)
        with (out/'progress.jsonl').open('a') as stream:
            stream.write(json.dumps(row, sort_keys=True)+'\n')
        print(json.dumps(dict(trial=row['trial'], status=row['status'])), flush=True)
    result = dict(schema_version=1, panel_sha256=PANEL_SHA256,
                  source_commit=panel['source_commit'],
                  parent_exact_result_sha256=panel['parent_exact_result_sha256'],
                  exporter_build=build_record,
                  cms_binary_sha256=panel['cms_executable_sha256'], rows=rows,
                  full_sat_ic_admission=False, natural_yield_estimate=None,
                  speedup=None, scope=panel['claim_boundary'])
    write(out/'summary.json', result, exclusive=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run(read(PANEL), args.source_root.resolve(), args.out.resolve())


if __name__ == '__main__':
    main()
