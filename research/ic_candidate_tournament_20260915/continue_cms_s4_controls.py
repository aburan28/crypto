#!/usr/bin/env python3
"""Resume only the four preregistered SAT controls from sealed valid exports."""
import argparse
import io
import json
import os
from pathlib import Path
import shutil
import tarfile

from oracle import require
from run_cms_s4_controls import (
    PANEL, PANEL_SHA256, REGISTRATION, digest, lift_source_assignment,
    meter, parse_cms_model, preflight, sha_bytes, validate_export,
    validate_xor_dimacs,
)
from tournament import read, write

STAGE_A = REGISTRATION/'stage-a-evidence.tar.gz'
STAGE_A_SHA256 = '2b1fa7abb093c94e741cebecdce35ea8fd59dbc0fb301e9aac9f1b95e0b195ef'
STAGE_A_SUMMARY_SHA256 = 'b75c9e2c2f621bdef499e2f441d772a4435d93d1b1bfc49278ad8adafb68778c'
ORIGINAL_RUNNER_SHA256 = 'fc555329cbe4fb299d03fac68fecf22d5a6c29abd9743bf5fefff25ddf015f79'
EXPORTER_BINARY_SHA256 = '4da7f5781da1aba2f76d6b9ce8f6e4ea43845f2465cf8edf0612615931d0629a'


def stage_a_files():
    data = STAGE_A.read_bytes()
    require(sha_bytes(data) == STAGE_A_SHA256, 'sealed original stage-A archive changed')
    with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as archive:
        members = archive.getmembers()
        names = [m.name for m in members]
        require(len(names) == len(set(names))
                and all(m.isfile() and not m.name.startswith('/')
                        and '..' not in Path(m.name).parts for m in members),
                'unsafe or duplicate stage-A archive member')
        return {m.name: archive.extractfile(m).read() for m in members}


def verify_stage_a(files, panel, curve, base):
    require(files['registered-panel.json'] == PANEL.read_bytes()
            and sha_bytes(files['registered-panel.json']) == PANEL_SHA256
            and sha_bytes(files['registered-runner.py']) == ORIGINAL_RUNNER_SHA256
            and sha_bytes(files['summary.json']) == STAGE_A_SUMMARY_SHA256
            and sha_bytes(files['build/exporter']) == EXPORTER_BINARY_SHA256
            and sha_bytes(files['cms-executable']) == panel['cms_executable_sha256']
            and sha_bytes(files['process_meter.py']) == panel['process_meter_sha256']
            and sha_bytes(files['cms_parser.py']) == panel['cms_parser_sha256'],
            'original panel, code, build or binary changed')
    original = json.loads(files['summary.json'])
    build = json.loads(files['build/build-record.json'])
    parent_source = json.loads(files['parent-source-manifest.json'])
    require(build['exporter_executable_sha256'] == EXPORTER_BINARY_SHA256
            and build['exporter_source_sha256'] == panel['exporter_source_sha256']
            and build['source_commit'] == panel['source_commit']
            and build['parent_source_manifest_sha256'] == sha_bytes(json.dumps(
                parent_source, sort_keys=True, separators=(',', ':')).encode())
            and original['panel_sha256'] == PANEL_SHA256
            and original['cms_binary_sha256'] == panel['cms_executable_sha256']
            and len(original['rows']) == len(panel['schedule']) == 4,
            'original source/build or summary binding changed')
    validations = []
    for item, row in zip(panel['schedule'], original['rows']):
        trial = item['trial']
        prefix = f'trial-{trial:02d}/'
        require(row['trial'] == trial and row['public_point'] == item['point']
                and row['probe_scalar'] == item['probe_scalar']
                and row['exact_relation_exists'] is item['exact_relation_exists']
                and row['status'] == 'INVALID_EXPORT'
                and row['reason'] == "KeyError: 'mitm'"
                and row['exporter']['returncode'] == 0
                and row['exporter']['timed_out'] is False
                and row['cms'] is None
                and prefix+'cms.metrics.json' not in files
                and prefix+'cms.stdout' not in files,
                'not the exact registered export-only gate failure')
        manifest = json.loads(files[prefix+'instance/manifest.json'])
        require('mitm' not in manifest
                and manifest['direct_meet_in_the_middle']['status']
                    == 'not_run_in_export_process',
                'original manifest does not contain the documented key')
        # Change only the reader's local key. The raw manifest is not altered.
        compatible = dict(manifest, mitm=manifest['direct_meet_in_the_middle'])
        # `validate_export` inspects filesystem bytes, so run it after safe
        # extraction; this preliminary check freezes exact source byte counts.
        for descriptor in manifest['exports'].values():
            name = prefix+'instance/'+descriptor['path']
            require(name in files and len(files[name]) == descriptor['bytes'],
                    'original source export bytes are missing or changed')
        validations.append((item, row, compatible))
    return original, validations


def run_sat(item, original_row, manifest, panel, curve, base, source_root, output_root):
    trial = item['trial']
    instance = source_root/f'trial-{trial:02d}'/'instance'
    destination = output_root/f'trial-{trial:02d}'
    destination.mkdir()
    exports = validate_export(manifest, item, panel, base, instance)
    cms = source_root/'cms-executable'
    command = [cms, '--verb', '1', '--threads', '1', '--random', '1',
               '--maxsol', '1', '--maxconfl', str(panel['cms_conflict_budget']),
               str(instance/'instance.xor.cnf')]
    metrics = meter(command, destination, 'cms', panel['cms_timeout_seconds'])
    row = dict(trial=trial, probe_scalar=item['probe_scalar'],
               public_point=item['point'], exact_relation_exists=item['exact_relation_exists'],
               original_stage_a_status=original_row['status'],
               original_exporter=original_row['exporter'],
               exports=exports, raw_manifest_sha256=digest(instance/'manifest.json'),
               cms=metrics, source_model_valid=None, point_witness=None,
               status='SOLVER_ERROR')
    stdout = (destination/'cms.stdout').read_text()
    if metrics['timed_out']:
        row['status'] = 'TIMEOUT'
    elif metrics['returncode'] == 10 and 's SATISFIABLE' in stdout:
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
    elif metrics['returncode'] == 20 and 's UNSATISFIABLE' in stdout:
        row['status'] = ('CONTRADICTS_EXACT_POSITIVE'
                         if item['exact_relation_exists'] else 'SOURCE_UNSAT')
    elif metrics['returncode'] == 0:
        row['status'] = 'UNKNOWN_INCONCLUSIVE'
    return row


def run(panel, out):
    _, curve, base = preflight(panel)
    require(not out.exists(), 'continuation output already exists; no SAT job is retried')
    files = stage_a_files()
    original, registered = verify_stage_a(files, panel, curve, base)
    out.mkdir(parents=True)
    source_root = out/'original-stage-a'
    source_root.mkdir()
    for name, payload in files.items():
        path = source_root/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
    os.chmod(source_root/'cms-executable', 0o755)
    # Recheck the actual extracted export bytes for all four before solving
    # any of them. A failure here does not launch a partial solver schedule.
    for item, _, manifest in registered:
        validate_export(manifest, item, panel, base,
                        source_root/f"trial-{item['trial']:02d}"/'instance')
    shutil.copy2(Path(__file__), out/'registered-continuation.py')
    shutil.copy2(REGISTRATION/'GATE-REPAIR.md', out/'GATE-REPAIR.md')
    stage_b = out/'stage-b'
    stage_b.mkdir()
    rows = []
    for item, original_row, manifest in registered:
        row = run_sat(item, original_row, manifest, panel, curve, base,
                      source_root, stage_b)
        rows.append(row)
        with (out/'progress.jsonl').open('a') as stream:
            stream.write(json.dumps(row, sort_keys=True)+'\n')
        print(json.dumps(dict(trial=row['trial'], status=row['status'])), flush=True)
    result = dict(schema_version=1, status='SAT_CONTROLS_CONTINUED',
                  panel_sha256=PANEL_SHA256,
                  original_stage_a_archive_sha256=STAGE_A_SHA256,
                  original_stage_a_summary_sha256=STAGE_A_SUMMARY_SHA256,
                  original_stage_a_statuses=[r['status'] for r in original['rows']],
                  rows=rows, source_bound_complete_sat_ic=False,
                  natural_yield_estimate=None, online_speedup=None,
                  scope=panel['claim_boundary'])
    write(out/'summary.json', result, exclusive=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run(read(PANEL), args.out.resolve())


if __name__ == '__main__':
    main()
