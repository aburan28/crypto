#!/usr/bin/env python3
"""One-shot source-receipted static CMS controls on disclosed n17a1 queries."""
import argparse
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

from generic_solver_feasibility import check_source_checkout
from oracle import Curve, require
from run_cms_s4_controls import (
    EXPORTER, build_exporter, digest, lift_source_assignment, meter,
    parse_cms_model, sha_bytes, validate_xor_dimacs,
)
from run_generic_exact_yield_audit import (
    PANEL as EXACT_PANEL, audit as audit_exact, load_evidence, record,
)
from tournament import read, write

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCRIPTS = ROOT/'scripts'
REGISTRATION = HERE/'goal_20260924/static-cms-s4-controls'
PANEL = REGISTRATION/'panel.json'
PANEL_SHA256 = '08290ed3c823bc878c123fdfdba98006e696b41a323db656c45ada658c4deec7'
BUNDLE = ROOT/'research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-20-phase-b-terminal-evidence-successor-04-20260910'
sys.path.insert(0, str(SCRIPTS))
from verify_stage20_phase_b_terminal_evidence import verify as verify_build_bundle  # noqa: E402


def admit(panel, *, require_local_binary=True):
    require(digest(PANEL) == PANEL_SHA256, 'static SAT panel changed')
    require(panel['schema_version'] == 1
            and panel['status'] == 'REGISTERED_BEFORE_EXECUTION'
            and panel['source_commit'] == '765c3c5f19032bd852163805f257c56babef2040'
            and panel['curve_cell'] == 'n17a1'
            and panel['field_degree'] == 17 and panel['curve_a'] == 1
            and panel['irreducible_low_terms'] == [0, 3]
            and panel['factor_base'] == {'kind': 'standard_subspace', 'dimension': 6}
            and panel['summands'] == 3
            and panel['source_encoding'] == 'wide-symmetrised-S4-circuit-XOR-DIMACS'
            and panel['export_only'] is True
            and panel['export_nonce'] == 2026092932
            and panel['export_timeout_seconds'] == 60
            and panel['cms_timeout_seconds'] == 120
            and panel['cms_conflict_budget'] == 1000000
            and panel['cms_threads'] == panel['cms_random_seed']
                == panel['cms_max_models_per_query'] == 1
            and [item['trial'] for item in panel['schedule']]
                == [0, 3, 4, 67, 71, 75],
            'static SAT controls or limits changed')
    require(digest(SCRIPTS/'process_meter.py') == panel['process_meter_sha256']
            and digest(SCRIPTS/'run_koblitz_pdp_matrix.py')
                == panel['cms_parser_sha256'],
            'process meter or model parser changed')
    require(digest(SCRIPTS/'verify_stage20_phase_b_terminal_evidence.py')
            == panel['cms_build_verifier_sha256']
            and digest(BUNDLE/'bundle-seal.json')
                == panel['cms_build_bundle_seal_sha256']
            and digest(BUNDLE/'builds/cryptominisat-receipt.json')
                == panel['cms_build_receipt_sha256'],
            'static solver build bundle changed')
    verified = verify_build_bundle(BUNDLE)
    require(verified['status'] == 'pass', 'archived static solver build did not verify')
    receipt = read(BUNDLE/'builds/cryptominisat-receipt.json')
    require(receipt['status'] == 'completed'
            and receipt['tool'] == 'cryptominisat'
            and receipt['source_commit'] == panel['cms_source_commit']
            and receipt['dependency_commits'] == panel['cms_dependency_commits']
            and receipt['binaries']['cryptominisat']['sha256']
                == panel['cms_executable_sha256'],
            'static solver receipt does not bind selected source and executable')
    if require_local_binary:
        cms = Path(panel['cms_path'])
        require(cms.is_file() and digest(cms) == panel['cms_executable_sha256'],
                'selected static SAT executable changed')
    require(digest(EXACT_PANEL.parent/'RESULT.json')
            == panel['parent_exact_result_sha256']
            and digest(HERE/panel['parent_recovery_archive_file'])
                == panel['parent_recovery_archive_sha256'],
            'parent exact or worker evidence changed')
    exact_panel = read(EXACT_PANEL)
    exact = read(EXACT_PANEL.parent/'RESULT.json')
    replay = audit_exact(exact_panel)
    replay['runner_sha256'] = exact['runner_sha256']
    require(replay == exact, 'parent exact-yield/source audit did not replay')
    files = load_evidence(exact_panel)
    report = record(files, 'jobs/n17a1/f5/stdout.json')
    fixture = report['fixture']
    curve = Curve(fixture)
    require(fixture['degree'] == 17 and fixture['curve_a'] == 1
            and fixture['irreducible']['low_terms'] == [0, 3],
            'n17a1 curve representation changed')
    base = tuple(curve.decode(point) for point in report['factor_base'])
    require(len(base) == len(set(base)) == 63, 'geometric base changed')
    source_row = next(row for row in exact['rows']
                      if row['cell'] == 'n17a1' and row['solver'] == 'f5')
    for item in panel['schedule']:
        label = source_row['queries'][item['trial']]
        require(label['trial'] == item['trial']
                and label['a'] == item['probe_scalar']
                and label['exact_feasible'] is item['exact_relation_exists']
                and curve.mul(curve.g, item['probe_scalar']) == tuple(item['point']),
                'registered scalar, public point or exact label changed')
    require([item['exact_relation_exists'] for item in panel['schedule']]
            == [False, False, True, True, True, True],
            'control classes changed')
    return files, curve, base, verified


def validate_export(manifest, item, panel, base, instance):
    require(manifest['kind'] == 'binary_koblitz_pdp_cross_solver_instance'
            and manifest['target_mode'] == 'explicit_affine'
            and manifest['n'] == 17 and manifest['ell'] == 6
            and manifest['m'] == 3 and manifest['curve_a'] == 1
            and manifest['irreducible_low_terms'] == [0, 3]
            and manifest['representation'] == 'symmetrised_s4'
            and manifest['factor_base_basis_bitmasks']
                == ['1', '2', '4', '8', '16', '32']
            and manifest['factor_base_geometry']['curve_points'] == len(base)
            and [int(manifest['target'][key]) for key in ('x', 'y')]
                == item['point']
            and manifest['seed'] == panel['export_nonce']
            and manifest['blind_instance_id'] == f"control-{item['trial']:02d}"
            and manifest['native_sat']['status'] == 'not_run_in_export_process'
            and manifest['direct_meet_in_the_middle']['status']
                == 'not_run_in_export_process',
            'exported instance differs from registered source mode or public point')
    exports = manifest['exports']
    require(set(exports) == {'wdsat_anf', 'cryptominisat_xor_dimacs',
                             'magma_boolean_f4'},
            'source export set changed')
    identities = {}
    for name, descriptor in exports.items():
        path = instance/descriptor['path']
        require(path.is_file() and path.parent == instance
                and path.stat().st_size == descriptor['bytes'],
                'missing, unsafe or changed source export')
        identities[name] = dict(bytes=path.stat().st_size, sha256=digest(path))
    return identities


def one_control(panel, item, exporter, cms, curve, base, out):
    trial = item['trial']
    directory = out/f'trial-{trial:02d}'
    directory.mkdir()
    instance = directory/'instance'
    export_command = [exporter, '17', '6', 'standard',
                      str(panel['export_nonce']),
                      str(panel['cms_conflict_budget']),
                      instance, '1', '0', '--target-x', str(item['point'][0]),
                      '--target-y', str(item['point'][1]),
                      '--blind-instance-id', f'control-{trial:02d}',
                      '--export-only']
    exported = meter(export_command, directory, 'export',
                     panel['export_timeout_seconds'])
    row = dict(trial=trial, probe_scalar=item['probe_scalar'],
               public_point=item['point'],
               exact_relation_exists=item['exact_relation_exists'],
               exporter=exported, cms=None, source_model_valid=None,
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
            row['status'] = ('CONTRADICTS_EXACT_NEGATIVE' if witness['group_replay']
                             and not item['exact_relation_exists'] else
                             'VALID_POINT_WITNESS' if witness['group_replay'] else
                             'SOURCE_MODEL_NONLIFTING')
    elif measured['returncode'] == 20 and 's UNSATISFIABLE' in stdout:
        row['status'] = ('CONTRADICTS_EXACT_POSITIVE'
                         if item['exact_relation_exists'] else 'SOURCE_UNSAT')
    elif measured['returncode'] == 0:
        row['status'] = 'UNKNOWN_INCONCLUSIVE'
    else:
        row['status'] = 'SOLVER_ERROR'
    return row


def run(panel, source_root, out):
    files, curve, base, build_verified = admit(panel)
    check_source_checkout(source_root)
    require(digest(source_root/EXPORTER) == panel['exporter_source_sha256'],
            'pinned exporter source changed')
    require(not out.exists(), 'output exists; no registered job may be retried')
    out.mkdir(parents=True)
    for source, target in ((PANEL, 'registered-panel.json'),
                           (REGISTRATION/'PROTOCOL.md', 'PROTOCOL.md'),
                           (Path(__file__), 'registered-runner.py'),
                           (SCRIPTS/'process_meter.py', 'process_meter.py'),
                           (SCRIPTS/'run_koblitz_pdp_matrix.py', 'cms_parser.py'),
                           (BUNDLE/'builds/cryptominisat-receipt.json',
                            'cms-build-receipt.json'),
                           (BUNDLE/'bundle-seal.json', 'cms-build-bundle-seal.json')):
        shutil.copy2(source, out/target)
    write(out/'cms-build-verification.json', build_verified, exclusive=True)
    write(out/'host.json', dict(system=platform.system(),
                                machine=platform.machine(),
                                cpu_count=os.cpu_count(),
                                scope='local stage correctness diagnostic'),
          exclusive=True)
    parent_source = record(files, 'build/source-manifest.json')
    write(out/'parent-source-manifest.json', parent_source, exclusive=True)
    cms = out/'cms-executable'
    shutil.copy2(panel['cms_path'], cms)
    require(digest(cms) == panel['cms_executable_sha256'],
            'copied static SAT executable changed')
    linkage = subprocess.check_output(['otool', '-L', str(cms)], text=True)
    (out/'cms-linkage.txt').write_text(linkage)
    require('@rpath' not in linkage and '/opt/homebrew' not in linkage,
            'copied SAT executable has a nonportable library path')
    version = meter([cms, '--version'], out, 'cms-preflight', 10)
    if (version['returncode'] != 0 or version['timed_out']
            or 'CryptoMiniSat version 5.14.7'
                not in (out/'cms-preflight.stdout').read_text()):
        write(out/'summary.json', dict(schema_version=1,
              status='PREFLIGHT_FAILURE', panel_sha256=PANEL_SHA256,
              rows=[], cms_preflight=version,
              full_sat_ic_admission=False, natural_yield_estimate=None,
              online_speedup=None), exclusive=True)
        return
    build_dir = out/'build'
    build_dir.mkdir()
    exporter, build_record = build_exporter(panel, source_root, build_dir,
                                             parent_source)
    rows = []
    for item in panel['schedule']:
        row = one_control(panel, item, exporter, cms, curve, base, out)
        rows.append(row)
        with (out/'progress.jsonl').open('a') as stream:
            stream.write(json.dumps(row, sort_keys=True)+'\n')
        print(json.dumps(dict(trial=row['trial'], status=row['status'])),
              flush=True)
    write(out/'summary.json', dict(schema_version=1, status='CONTROL_PANEL_COMPLETE',
          panel_sha256=PANEL_SHA256, source_commit=panel['source_commit'],
          parent_exact_result_sha256=panel['parent_exact_result_sha256'],
          exporter_build=build_record,
          cms_build_receipt_sha256=panel['cms_build_receipt_sha256'],
          cms_binary_sha256=panel['cms_executable_sha256'],
          cms_preflight=version, rows=rows,
          full_sat_ic_admission=False, natural_yield_estimate=None,
          online_speedup=None, scope=panel['claim_boundary']), exclusive=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run(read(PANEL), args.source_root.resolve(), args.out.resolve())


if __name__ == '__main__':
    main()
