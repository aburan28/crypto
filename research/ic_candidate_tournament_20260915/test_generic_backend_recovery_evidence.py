"""Replay the committed, source-bound disclosed recovery evidence.

This rechecks every exited worker report, including bounded failures. It never
executes a worker or replaces the registered measurement with a fresh run.
"""
import hashlib
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest

from generic_build import verify_build_record
from run_generic_backend_recovery_pilot import (
    CELLS, PANEL, SOLVERS, audit_report, declared_job, fixtures_from_inventory,
    resources, validate_panel,
)
from tournament import read


EVIDENCE = (Path(__file__).resolve().parent / 'goal_20260924'
            / 'generic-backend-recovery-pilot' / 'evidence.tar.gz')


def payloads():
    expected_hash = (EVIDENCE.parent / 'EVIDENCE.sha256').read_text().split()[0]
    data = EVIDENCE.read_bytes()
    assert hashlib.sha256(data).hexdigest() == expected_hash
    with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as archive:
        members = archive.getmembers()
        names = [member.name for member in members]
        assert len(names) == len(set(names))
        assert all(member.isfile() and not member.name.startswith('/')
                   and '..' not in Path(member.name).parts for member in members)
        return {member.name: archive.extractfile(member).read() for member in members}


def record(files, name):
    return json.loads(files[name])


class RecoveryEvidenceTests(unittest.TestCase):
    def test_replay_all_registered_jobs_and_admission_receipts(self):
        files = payloads()
        panel = read(PANEL)
        validate_panel(panel)
        self.assertEqual(files['registered-panel.json'], PANEL.read_bytes())
        fixtures = fixtures_from_inventory(panel)
        summary = record(files, 'summary.json')
        self.assertEqual([(row['cell'], row['solver']) for row in summary['rows']],
                         [(cell, solver) for cell in CELLS for solver in SOLVERS])
        self.assertEqual(len(summary['rows']), 20)
        self.assertEqual([json.loads(line) for line in
                          files['progress.jsonl'].splitlines()], summary['rows'])
        self.assertEqual(summary['status'], 'DISCLOSED_RECOVERY_NOT_ESTABLISHED')
        self.assertTrue(summary['n17_f4_f5_complete'])
        self.assertFalse(summary['n17_sat_complete'])

        source = record(files, 'build/source-manifest.json')
        build = record(files, 'build/build-record.json')
        verify_build_record(build, source)
        self.assertEqual(summary['source_manifest_sha256'],
                         build['source_manifest_sha256'])
        self.assertEqual(summary['worker_sha256'], build['worker_sha256'])
        with tarfile.open(fileobj=io.BytesIO(files['build/root-source.tar.gz']),
                          mode='r:gz') as root_archive:
            self.assertEqual(set(root_archive.getnames()), set(source['root_files']))
            for name, expected_hash in source['root_files'].items():
                self.assertEqual(hashlib.sha256(root_archive.extractfile(name).read()).hexdigest(),
                                 expected_hash, name)

        with tempfile.TemporaryDirectory() as temp:
            worker = Path(temp) / 'worker'
            worker.write_bytes(files['build/worker'])
            self.assertEqual(hashlib.sha256(worker.read_bytes()).hexdigest(),
                             build['worker_sha256'])
            collection_queries = {}
            for number, row in enumerate(summary['rows'], 1):
                cell, solver = row['cell'], row['solver']
                prefix = f'jobs/{cell}/{solver}/'
                job = record(files, prefix + 'job.json')
                self.assertEqual(job, declared_job(panel, cell, solver, fixtures[cell]))
                self.assertEqual(record(files, prefix + 'result.json'), row)
                process = record(files, prefix + 'process.json')
                for key in ('cell', 'solver', 'exit_code', 'diagnostic_process_wall_ns',
                            'sampled_peak_rss_bytes', 'memory_policy'):
                    self.assertEqual(process[key], row[key])
                if row['disposition'] == 'TIMEOUT':
                    self.assertEqual(files[prefix + 'stdout.json'], b'')
                    self.assertNotIn(prefix + 'admission.json', files)
                    self.assertIsNone(row['online_wall_ns'])
                    self.assertIsNone(row['cold_wall_ns'])
                    self.assertEqual(row['audit_status'], 'NOT_AVAILABLE')
                    continue
                self.assertIn(row['disposition'],
                              ('VERIFIED_COMPLETE', 'BOUNDED_INCOMPLETE_REPORT'))
                report = record(files, prefix + 'stdout.json')
                queries = [(attempt['a'], attempt['b'])
                           for batch in report['collection_reports']
                           for attempt in batch['attempts']]
                prior = collection_queries.setdefault(cell, queries)
                self.assertEqual(queries, prior,
                                 'reported arms must share the same ordinary query sequence')
                receipt, diagnostic = audit_report(
                    report, fixtures[cell], job, build, source, worker,
                    process['diagnostic_process_wall_ns'], panel, cell, number)
                recorded = record(files, prefix + 'admission.json')
                # Audit runtime is external to the worker and varies on replay.
                receipt['run'].pop('independent_audit_wall_ns')
                recorded['run'].pop('independent_audit_wall_ns')
                self.assertEqual(receipt, recorded)
                for key, value in diagnostic.items():
                    self.assertEqual(row[key], value)
                self.assertEqual(row['audit_status'], 'PASS')
                self.assertIsNone(row['paired_rho_speedup'])
                self.assertFalse(row['promotion_eligible'])

        f5 = next(row for row in summary['rows'] if
                  (row['cell'], row['solver']) == ('n17a1', 'f5'))
        self.assertTrue(f5['verified_complete'])
        self.assertEqual(f5['actual_base_points'], 62)
        self.assertEqual(f5['folded_columns'], f5['final_rank'])
        self.assertEqual(f5['pdp_outcomes'], {'witness': 29, 'incomplete': 75})
        self.assertEqual(record(files, 'jobs/n17a1/f5/admission.json')
                         ['run']['certificate']['solutions'], ['40605'])


if __name__ == '__main__':
    unittest.main()
