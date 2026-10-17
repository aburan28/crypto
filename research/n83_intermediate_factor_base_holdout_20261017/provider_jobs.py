#!/usr/bin/env python3
"""Bind four overlapping GitHub allocations to archived host receipts."""

import hashlib
import json
from datetime import datetime
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def instant(value):
    result = datetime.fromisoformat(value.replace('Z', '+00:00'))
    if result.tzinfo is None:
        raise ValueError('provider timestamp lacks a timezone')
    return result


def verify(jobs_path, strata, run_id, attempt, job_stem):
    provider = json.loads(jobs_path.read_text())
    jobs = provider.get('jobs')
    if (provider.get('total_count') != 4 or not isinstance(jobs, list) or
            len(jobs) != 4 or len(strata) != 4):
        raise ValueError('provider did not return four jobs and strata')
    names = {f'{job_stem} ({number})' for number in range(4)}
    if {job.get('name') for job in jobs} != names:
        raise ValueError('provider matrix job names differ')
    by_name = {job['name']: job for job in jobs}
    rows = []
    for number, item in enumerate(strata):
        job = by_name[f'{job_stem} ({number})']
        host = item['host']
        if (item['stratum'] != number or job.get('run_id') != run_id or
                job.get('run_attempt') != attempt or
                job.get('status') != 'completed' or
                job.get('conclusion') != 'success' or
                not isinstance(job.get('id'), int) or job['id'] <= 0 or
                not isinstance(job.get('runner_id'), int) or
                job['runner_id'] <= 0 or
                job.get('runner_name') != host.get('runner_name') or
                host.get('run_url') !=
                f'https://github.com/aburan28/crypto/actions/runs/{run_id}'):
            raise ValueError(f'stratum {number}: provider job or host differs')
        start, end = instant(job['started_at']), instant(job['completed_at'])
        if not start <= instant(host['recorded_utc']) <= end:
            raise ValueError(f'stratum {number}: host receipt outside job interval')
        rows.append({'stratum': number, 'job_id': job['id'],
                     'runner_id': job['runner_id'],
                     'runner_name': job['runner_name'],
                     'started_at': job['started_at'],
                     'completed_at': job['completed_at'],
                     'host_receipt_sha256': digest(
                         Path(item['results']) / 'host.json')})
    if (len({row['job_id'] for row in rows}) != 4 or
            len({row['runner_id'] for row in rows}) != 4 or
            len({row['runner_name'] for row in rows}) != 4):
        raise ValueError('four distinct provider allocations required')
    overlap_start = max(instant(row['started_at']) for row in rows)
    overlap_end = min(instant(row['completed_at']) for row in rows)
    if overlap_start >= overlap_end:
        raise ValueError('four provider jobs did not overlap')
    return {'schema': 'n83.intermediate-provider-runner-identity/v1',
            'status': 'PASS_FOUR_DISTINCT_CONCURRENT_RUNNERS',
            'provider_jobs_sha256': digest(jobs_path),
            'run_id': run_id, 'run_attempt': attempt,
            'shared_parent_instance_id':
                strata[0]['host']['parent_instance_id'],
            'overlap_start': overlap_start.isoformat(),
            'overlap_end': overlap_end.isoformat(),
            'jobs': rows}
