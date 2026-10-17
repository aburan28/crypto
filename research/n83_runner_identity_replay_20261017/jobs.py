#!/usr/bin/env python3
"""Bind two overlapping GitHub job allocations to archived runner receipts."""

import hashlib
import json
from datetime import datetime
from pathlib import Path


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def instant(value):
    parsed = datetime.fromisoformat(value.replace('Z', '+00:00'))
    if parsed.tzinfo is None:
        raise ValueError('provider job timestamp has no timezone')
    return parsed


def verify(jobs_path: Path, halves, run_id: int, attempt: int, job_stem: str):
    """Require distinct provider runner IDs on concurrent, successful jobs.

    The producer's parent_instance_id is a workflow job *label*, shared by
    both matrix rows. It is checked by each frozen bind_half, then ignored
    for cross-host identity. The provider job IDs and runner IDs are the
    allocation identifiers; host receipts must independently name the same
    runner assignments.
    """
    provider = json.loads(jobs_path.read_text())
    jobs = provider.get('jobs')
    if provider.get('total_count') != 2 or not isinstance(jobs, list) or len(jobs) != 2:
        raise ValueError('provider did not return exactly two jobs')
    expected_names = {f'{job_stem} ({half})' for half in (0, 1)}
    if {job.get('name') for job in jobs} != expected_names:
        raise ValueError('provider matrix job names differ')
    by_name = {job['name']: job for job in jobs}
    rows = []
    for half in (0, 1):
        item = halves[half]
        job = by_name[f'{job_stem} ({half})']
        host = item['host']
        if (item['half'] != half or job.get('run_id') != run_id or
                job.get('run_attempt') != attempt or
                job.get('status') != 'completed' or
                job.get('conclusion') != 'success' or
                not isinstance(job.get('id'), int) or job['id'] <= 0 or
                not isinstance(job.get('runner_id'), int) or
                job['runner_id'] <= 0 or
                job.get('runner_name') != host.get('runner_name') or
                host.get('run_url') !=
                f'https://github.com/aburan28/crypto/actions/runs/{run_id}'):
            raise ValueError(f'half {half}: provider job or runner receipt differs')
        start, end = instant(job['started_at']), instant(job['completed_at'])
        recorded = instant(host['recorded_utc'])
        if not start <= recorded <= end:
            raise ValueError(f'half {half}: host receipt is outside its provider job')
        rows.append({'half': half, 'job_id': job['id'],
                     'runner_id': job['runner_id'],
                     'runner_name': job['runner_name'],
                     'started_at': job['started_at'],
                     'completed_at': job['completed_at'],
                     'host_receipt_sha256': sha256(
                         Path(item['results']) / 'host.json')})
    if (rows[0]['job_id'] == rows[1]['job_id'] or
            rows[0]['runner_id'] == rows[1]['runner_id'] or
            rows[0]['runner_name'] == rows[1]['runner_name']):
        raise ValueError('provider jobs did not use distinct runner allocations')
    overlap_start = max(instant(row['started_at']) for row in rows)
    overlap_end = min(instant(row['completed_at']) for row in rows)
    if overlap_start >= overlap_end:
        raise ValueError('provider jobs did not overlap in time')
    return {'schema': 'n83.provider-runner-identity/v1',
            'status': 'PASS_DISTINCT_CONCURRENT_RUNNERS',
            'provider_jobs_sha256': sha256(jobs_path),
            'run_id': run_id, 'run_attempt': attempt,
            'shared_parent_instance_id': halves[0]['host']['parent_instance_id'],
            'overlap_start': overlap_start.isoformat(),
            'overlap_end': overlap_end.isoformat(),
            'jobs': rows}
