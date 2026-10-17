#!/usr/bin/env python3
"""Check four-host identity against duplicate and nonoverlap evidence."""

import json
import tempfile
from pathlib import Path

import provider_jobs


with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    run_id, attempt = 12345, 1
    strata, jobs = [], []
    for number in range(4):
        results = root / f'stratum{number}'
        results.mkdir()
        host = {'runner_name': f'runner-{number}',
                'run_url': f'https://github.com/aburan28/crypto/actions/runs/{run_id}',
                'recorded_utc': '2026-10-17T12:30:00Z',
                'parent_instance_id': 'shared-matrix-label'}
        (results / 'host.json').write_text(json.dumps(host))
        strata.append({'stratum': number, 'host': host, 'results': str(results)})
        jobs.append({'name': f'intermediate-stratum ({number})',
                     'run_id': run_id, 'run_attempt': attempt,
                     'status': 'completed', 'conclusion': 'success',
                     'id': 100 + number, 'runner_id': 200 + number,
                     'runner_name': f'runner-{number}',
                     'started_at': '2026-10-17T12:00:00Z',
                     'completed_at': '2026-10-17T13:00:00Z'})
    path = root / 'jobs.json'

    def check():
        path.write_text(json.dumps({'total_count': 4, 'jobs': jobs}))
        return provider_jobs.verify(path, strata, run_id, attempt,
                                    'intermediate-stratum')

    assert check()['status'] == 'PASS_FOUR_DISTINCT_CONCURRENT_RUNNERS'
    jobs[0]['conclusion'] = 'failure'
    try:
        check()
        raise AssertionError('unexplained provider failure was admitted')
    except ValueError:
        pass
    strata[0].update({'archived_verifier_failure': True,
                      'archived_verifier_log_sha256': 'a' * 64,
                      'local_replay_exit': 0})
    assert check()['jobs'][0]['provider_job_conclusion'] == 'failure'
    jobs[3]['runner_id'] = jobs[2]['runner_id']
    try:
        check()
        raise AssertionError('duplicate runner ID was admitted')
    except ValueError:
        pass
    jobs[3]['runner_id'] = 203
    jobs[3]['started_at'] = '2026-10-17T13:01:00Z'
    jobs[3]['completed_at'] = '2026-10-17T14:00:00Z'
    strata[3]['host']['recorded_utc'] = '2026-10-17T13:30:00Z'
    try:
        check()
        raise AssertionError('serial runner allocation was admitted')
    except ValueError:
        pass

print(json.dumps({'status': 'PASS', 'cases': 5}))
