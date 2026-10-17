#!/usr/bin/env python3
"""Bind four intermediate-base archives, provider runners and selector."""

import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


prior = module('n83_intermediate_prior_aggregate',
               ROOT / 'research/n83_e2e_final_holdout_20261016/aggregate.py')
runner = module('n83_intermediate_aggregate_runner', HERE / 'run.py')
selector = module('n83_intermediate_aggregate_select', HERE / 'select.py')
identity = module('n83_intermediate_provider_runner_identity',
                  HERE / 'provider_jobs.py')


def read(path):
    return json.loads(path.read_text())


def sha256(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def write_once(path, value):
    blob = json.dumps(value, sort_keys=True, indent=2) + '\n'
    if path.exists() and path.read_text() != blob:
        raise ValueError(f'sealed output differs: {path}')
    path.write_text(blob)


def bind_stratum(stratum, results, archive, provider_sha, run_id, attempt,
                 head_sha, reuse_local_replay):
    cloud_analysis, archive_sha = prior.finegrid.archive_analysis(
        archive, provider_sha)
    if (results / 'analysis.json').read_bytes() != cloud_analysis:
        raise ValueError(f'stratum {stratum}: extracted analysis differs from ZIP')
    member_count = prior.finegrid.bind_extraction(archive, results)
    log = results.parent
    stdout_path = log / 'local_replay.stdout.txt'
    stderr_path = log / 'local_replay.stderr.txt'
    status_path = log / 'local_replay.exit_status.txt'
    if reuse_local_replay:
        status = int(status_path.read_text().strip())
        stdout, stderr = stdout_path.read_text(), stderr_path.read_text()
    else:
        command = [sys.executable, str(HERE / 'verify.py'), str(results),
                   '--stratum', str(stratum)]
        process = subprocess.run(command, cwd=ROOT, capture_output=True,
                                 text=True, check=False)
        status, stdout, stderr = (process.returncode, process.stdout,
                                  process.stderr)
        stdout_path.write_text(stdout)
        stderr_path.write_text(stderr)
        status_path.write_text(f'{status}\n')
    analysis = json.loads(cloud_analysis)
    if (status not in (0, 2) or stderr or json.loads(stdout) != analysis or
            (results / 'analysis.json').read_bytes() != cloud_analysis):
        raise ValueError(f'stratum {stratum}: independent local replay differs')
    if (status == 0) != (analysis['status'] == 'HALF_ADMITTED'):
        raise ValueError(f'stratum {stratum}: replay exit and admission differ')
    original_verifier_log = log / 'verify.log'
    original_failure = (original_verifier_log.is_file() and
                        "AttributeError: module 'n83_intermediate_choice' has no attribute 'five_run'"
                        in original_verifier_log.read_text())
    host, source = read(results / 'host.json'), read(results / 'source.json')
    arms, choice = runner.choice_module.choose()
    if (analysis['schema'] != 'n83.intermediate-stratum-analysis/v1' or
            analysis['stratum'] != stratum or
            analysis['point_indices'] != list(range(stratum * 5, stratum * 5 + 5)) or
            source['source_sha256'] != runner.SOURCE_SHA256 or
            source['targets_sha256'] != runner.TARGETS_SHA256 or
            source['answers_sha256'] != runner.ANSWERS_SHA256 or
            source['runner_sha256'] != sha256(HERE / 'run.py') or
            source['choice'] != choice or
            source['arms'] != [list(arm) for arm in arms] or
            host['git_sha'] != head_sha or host['provider'] != 'github_actions' or
            host['parent_instance_id'] != f'gh-{run_id}-{attempt}-intermediate-stratum' or
            host['instance_id'] !=
            f'gh-{run_id}-{attempt}-intermediate-stratum-intermediate-stratum{stratum}' or
            host['stratum'] != stratum or not host['runner_name']):
        raise ValueError(f'stratum {stratum}: host, source, choice or coverage differs')
    return {'stratum': stratum, 'analysis': analysis, 'source': source, 'host': host,
            'results': str(results), 'raw_artifact_sha256': archive_sha,
            'raw_artifact_zip': str(archive),
            'bound_result_members': member_count,
            'analysis_sha256': hashlib.sha256(cloud_analysis).hexdigest(),
            'local_replay_exit': status,
            'archived_verifier_failure': original_failure,
            'archived_verifier_log_sha256': sha256(original_verifier_log)
            if original_failure else None,
            'local_replay_stdout_sha256': sha256(stdout_path)}


def aggregate(strata, provider_jobs, run_id, attempt, head_sha):
    first = strata[0]
    keys = ('hashes', 'source_sha256', 'binary_sha256', 'rho_binary_sha256',
            'runner_sha256', 'targets_sha256', 'answers_sha256',
            'fixture_manifest_sha256', 'arms', 'choice', 'rho_seed_rule',
            'calibration_seed')
    same_source = all(item['source'][key] == first['source'][key]
                      for item in strata for key in keys)
    proof = identity.verify(provider_jobs, strata, run_id, attempt,
                            'intermediate-stratum')
    points = sorted([point for item in strata for point in item['analysis']['points']],
                    key=lambda item: item['point_index'])
    ids = [item[3] for item in first['source']['arms']]
    same_base = all(len({arm['raw_replay'].get('base_hash') for point in points
                         for arm in point['arms']
                         if arm['candidate_id'] == candidate and
                         arm['status'] == 'PASS_INDEPENDENT_REPLAY'}) == 1
                    for candidate in ids)
    host_strata = [
        {'stratum': item['stratum'],
         'point_indices': item['analysis']['point_indices'],
         'tier_status': item['analysis']['tier_status'],
         'aa_status': item['analysis']['aa_status'],
         'instance_id': item['host']['instance_id'],
         'runner_name': item['host']['runner_name'],
         'raw_artifact_sha256': item['raw_artifact_sha256'],
         'analysis_sha256': item['analysis_sha256']}
        for item in strata]
    evidence = ('PASS_INDEPENDENT_REPLAY' if same_source and same_base and
                len(points) == 20 and
                {point['point_index'] for point in points} == set(range(20)) and
                all(item['analysis']['evidence_status'] ==
                    'PASS_INDEPENDENT_REPLAY' and
                    item['analysis']['status'] == 'HALF_ADMITTED' and
                    item['local_replay_exit'] == 0 for item in strata) else
                'VERIFICATION_FAILED')
    validated = {
        'schema': 'n83.intermediate-holdout-validated/v1',
        'source_sha256': runner.SOURCE_SHA256,
        'targets_sha256': runner.TARGETS_SHA256,
        'answers_sha256': runner.ANSWERS_SHA256,
        'candidate_ids': ids, 'choice': first['source']['choice'],
        'evidence_status': evidence,
        'run_id': run_id, 'run_attempt': attempt, 'head_sha': head_sha,
        'same_source_and_binary': same_source,
        'different_cloud_instances': proof['status'] ==
                                     'PASS_FOUR_DISTINCT_CONCURRENT_RUNNERS',
        'provider_runner_identity': proof,
        'same_candidate_bases_across_strata': same_base,
        'host_strata': host_strata, 'points': points,
    }
    return validated, selector.select(validated)


def main():
    parser = argparse.ArgumentParser()
    for stratum in range(4):
        parser.add_argument(f'--stratum{stratum}-results', type=Path, required=True)
        parser.add_argument(f'--stratum{stratum}-zip', type=Path, required=True)
        parser.add_argument(f'--stratum{stratum}-provider-sha256', required=True)
    parser.add_argument('--provider-jobs', type=Path, required=True)
    parser.add_argument('--run-id', type=int, required=True)
    parser.add_argument('--run-attempt', type=int, default=1)
    parser.add_argument('--head-sha', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reuse-local-replay', action='store_true')
    args = parser.parse_args()
    strata = []
    for stratum in range(4):
        strata.append(bind_stratum(
            stratum, getattr(args, f'stratum{stratum}_results').resolve(),
            getattr(args, f'stratum{stratum}_zip').resolve(),
            getattr(args, f'stratum{stratum}_provider_sha256'), args.run_id,
            args.run_attempt, args.head_sha, args.reuse_local_replay))
    validated, result = aggregate(strata, args.provider_jobs.resolve(),
                                  args.run_id, args.run_attempt, args.head_sha)
    args.output.mkdir(parents=True, exist_ok=True)
    write_once(args.output / 'validated_rows.json', validated)
    write_once(args.output / 'selection.json', result)
    print(json.dumps({'evidence_status': validated['evidence_status'],
                      'online_selected': result['online']['selected_candidate_id'],
                      'cold_selected': result['cold']['selected_candidate_id']},
                     sort_keys=True))


if __name__ == '__main__':
    main()
