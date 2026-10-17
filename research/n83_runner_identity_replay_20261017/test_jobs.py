#!/usr/bin/env python3
"""Regression checks for provider-backed cloud runner identity."""

import json
import tempfile
import unittest
from pathlib import Path

import jobs


class ProviderIdentityTest(unittest.TestCase):
    def test_distinct_overlapping_jobs_and_mutations(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            jobs_path = root / 'jobs.json'
            halves = []
            rows = []
            for half in (0, 1):
                results = root / str(half)
                results.mkdir()
                (results / 'host.json').write_text('{}\n')
                name = f'GitHub Actions {100 + half}'
                halves.append({
                    'half': half, 'results': str(results),
                    'host': {'runner_name': name,
                             'recorded_utc': '2026-10-17T10:30:00Z',
                             'parent_instance_id': 'shared-workflow-label',
                             'run_url': 'https://github.com/aburan28/crypto/actions/runs/123'}})
                rows.append({'name': f'panel ({half})', 'id': 200 + half,
                             'run_id': 123, 'run_attempt': 1,
                             'status': 'completed', 'conclusion': 'success',
                             'runner_id': 100 + half, 'runner_name': name,
                             'started_at': '2026-10-17T10:00:00Z',
                             'completed_at': '2026-10-17T11:00:00Z'})

            def check(jobs):
                jobs_path.write_text(json.dumps({'total_count': 2, 'jobs': jobs}))
                return jobs_module_verify(jobs_path, halves)

            def jobs_module_verify(path, items):
                return jobs.verify(path, items, 123, 1, 'panel')

            self.assertEqual(check(rows)['status'],
                             'PASS_DISTINCT_CONCURRENT_RUNNERS')
            mutations = [
                {'runner_id': 100},
                {'runner_name': 'another runner'},
                {'completed_at': '2026-10-17T10:00:00Z'},
                {'conclusion': 'cancelled'},
            ]
            for mutation in mutations:
                altered = [dict(row) for row in rows]
                altered[1].update(mutation)
                with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                    check(altered)


if __name__ == '__main__':
    unittest.main()
