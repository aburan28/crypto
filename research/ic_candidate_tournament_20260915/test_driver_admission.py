"""Real retained arithmetic reports; synthetic build tags only for identity tests."""
import copy
import json
from pathlib import Path
import unittest
import tempfile

from driver_admission import (make_admission, run_record, check_admission, online_table,
                             freeze_admission, distinct_candidates)
from identity import sha256
from oracle import InvalidEvidence
from tournament import comparison, gate

DATA = Path(__file__).parent/'producer/testdata'


class DriverAdmissionTests(unittest.TestCase):
    def setUp(self):
        vector = json.loads((DATA/'n13-public.json').read_text())['ic']
        self.job, self.report = vector['job'], vector['report']
        self.ns = vector['process']['process_wall_ns']
        self.fixture = self.report['fixture']
        self.manifest = {p: str(i)*64 for i, p in enumerate((
            '.cargo/config.toml', 'examples/ic_tournament_worker.rs',
            'src/cryptanalysis/koblitz_tiny_ic.rs', 'src/cryptanalysis/ic_phase.rs'), 1)}
        self.source = sha256(self.manifest)
        self.report['diagnostics']['source_manifest_sha256'] = self.source
        self.inventory = dict(status='inventory', fixture=self.fixture,
            phase_schema=3, target_input='supplied_public_point',
            factor_base_orbits=self.report['factor_base_orbits'], columns=self.report['columns'],
            column_convention='representative', source_manifest_sha256=self.source, field_kernel='portable')
        self.kwargs = dict(job=self.job, fixture=self.fixture, report=self.inventory, manifest=self.manifest,
            metadata=dict(preparation=dict(source_manifest_sha256=self.source, reference='scaled'),
                build=dict(compiler='test compiler', target='test target', cargo_config_sha256='1'*64)),
            resources=dict(worker_threads=1, cpu=None, memory_bytes=None, timeout_seconds='180'),
            worker_sha256='5'*64)

    def test_float_process_duration_is_saved_outside_canonical_identity(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)/'inventory'
            process = dict(exit_code=0, status='EXITED', whole_process_wall_seconds=.001)
            def execute(binary, job, directory):
                directory.mkdir()
                (directory/'stdout.json').write_text(json.dumps(self.inventory))
                return process
            kwargs = {k:v for k,v in self.kwargs.items() if k != 'report'}
            admitted = freeze_admission(directory, binary=Path('test worker'), execute=execute, **kwargs)
            self.assertEqual(json.loads((directory/'process.json').read_text()), process)
            self.assertEqual(json.loads((directory/'admission.json').read_text()), admitted)
            self.assertEqual(json.loads((directory/'job.json').read_text())['mode'], 'inventory')

    def complete(self, admitted):
        return run_record(admitted, number=0, host_id='test host', status='complete',
                          native=self.report, process_wall_ns=self.ns)

    def test_native_record_verifies_without_inventing_instruction_cost(self):
        admitted = make_admission(**self.kwargs)
        record = self.complete(admitted)
        self.assertIsNone(record['total_operations'])
        self.assertFalse(record['complete_cold_cost'])
        self.assertTrue(record['candidate_id'].startswith('IC1N13'))
        self.assertGreater(record['stage_audit']['queries']['final_rank'], 0)
        self.assertEqual(record['native_timing']['cold']['wall_ns'], self.ns)
        self.assertTrue(record['native_timing']['online']['scalar_replay_included'])

    def test_failures_keep_key_but_cannot_retain_verified_online_time(self):
        admitted = make_admission(**self.kwargs)
        complete = self.complete(admitted)
        for status in ('error', 'oom', 'timeout'):
            record = run_record(admitted, number=0, host_id='test host', status=status)
            self.assertEqual(record['run_id'], complete['run_id'])
            self.assertIsNone(record['certificate'])
            self.assertIsNone(record['total_operations'])
            self.assertIsNone(record['native_timing'])
        with self.assertRaises(InvalidEvidence):
            run_record(admitted, number=0, host_id='host', status='error', native=self.report)

    def test_public_point_and_actual_backend_required_before_run(self):
        for change in ('public_targets', 'source', 'backend'):
            kwargs = copy.deepcopy(self.kwargs)
            if change == 'public_targets':
                kwargs['job'].pop('public_targets')
            elif change == 'source':
                kwargs['report']['source_manifest_sha256'] = '9'*64
            else:
                kwargs['job']['config']['linear_algebra'] = 'sparse'
            with self.subTest(change=change), self.assertRaises(InvalidEvidence):
                make_admission(**kwargs)

    def test_unused_rho_flag_cannot_create_a_new_ic_competitor(self):
        a = make_admission(**self.kwargs)
        self.kwargs['job']['config']['rho_parallel_walks'] = 1
        b = make_admission(**self.kwargs)
        with self.assertRaises(InvalidEvidence):
            distinct_candidates([('incumbent', a), ('fake_competitor', b)])
        distinct_candidates([('incumbent', a), ('aa_control', b)])

    def test_changed_canonical_admission_and_timing_fail_reconstruction(self):
        admitted = make_admission(**self.kwargs)
        kwargs = {k:v for k,v in self.kwargs.items() if k != 'report'}
        check_admission(admitted, **kwargs)
        altered = copy.deepcopy(admitted)
        altered['candidate']['candidate_id'] += '0'
        with self.assertRaises(InvalidEvidence):
            check_admission(altered, **kwargs)
        self.report['diagnostics']['online_wall_ns'] += 1
        with self.assertRaises(InvalidEvidence):
            self.complete(admitted)

    def test_rho_uses_reference_identity_and_same_canonical_workload(self):
        ic = make_admission(**self.kwargs)
        vector = json.loads((DATA/'n13-rho-prepared.json').read_text())
        report = vector['report']
        report['source_manifest_sha256'] = self.source
        self.kwargs['job'] = dict(self.job, mode='rho')
        rho = make_admission(**self.kwargs)
        record = run_record(rho, number=1, host_id='test host', status='complete', native=report,
                            process_wall_ns=vector['process']['process_wall_ns'])
        self.assertNotIn('candidate_id', record)
        self.assertTrue(record['reference_id'].startswith('RHO1h'))
        self.assertEqual(record['workload_id'], ic['workload']['workload_id'])
        self.assertTrue(record['native_timing']['online']['scalar_replay_included'])

    def test_online_pair_preserves_failed_reference_instead_of_a_win(self):
        ic = self.complete(make_admission(**self.kwargs))
        self.kwargs['job'] = dict(self.job, mode='rho')
        admitted = make_admission(**self.kwargs)
        failed = run_record(admitted, number=1, host_id='test host', status='timeout')
        rows = [dict(case='one', arm='incumbent', repetition=0, status='VERIFIED', measurement=ic),
                dict(case='one', arm='rho', repetition=0, status='TIMEOUT', measurement=failed)]
        table = online_table(rows, [dict(id='one', fixture=self.fixture)],
                             [dict(id='incumbent'), dict(id='rho')], 1, ['rho'])
        self.assertFalse(table[0]['verified'])
        self.assertIsNone(table[0]['online_speedup'])
        self.assertEqual(len(table[0]['run_ids']), 2)

    def test_admission_alone_cannot_promote(self):
        from test_tournament import rows
        result = comparison(rows(.5), 'candidate', draws=40)
        self.assertFalse(gate(result, dict(scientific_admission=True, reference_qualification=None)))


if __name__ == '__main__':
    unittest.main()
