"""Real retained arithmetic evidence with synthetic binary tags for adapter tests."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

from driver_admission import make_admission, check_admission, run_record, INPUT_LAW
from generic_driver import ADAPTER, INVENTORY_POLICY
from generic_build import digest
from generic_phases import PHASES, parse_profiles
from identity import workload_manifest
from oracle import InvalidEvidence
from test_generic_admission import build_fixture, retained, row
from tournament import executed_job, qualification_references


class GenericDriverTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.worker = Path(self.temporary.name)/'worker'
        self.worker.write_bytes(b'synthetic binary tag; never executed')
        self.build, self.source = build_fixture()
        self.build['worker_sha256'] = digest(self.worker)
        self.metadata = dict(adapter=ADAPTER, build_record=self.build)
        self.item = row()
        self.report, self.job = self.item['report'], self.item['job']
        self.fixture = self.report['fixture']
        self.inventory = copy.deepcopy(next(r['report'] for r in retained()
            if r['name'] == 'inventory-n9-0-factor'))
        self.inventory.update(field_kernel=self.report['collector_dispatch']['field_kernel'],
                              inventory_policy=INVENTORY_POLICY)
        self.resources = dict(cpu=0, worker_threads=1, memory_bytes=8*1024**3, timeout_seconds='60')

    def admission(self, **changes):
        return make_admission(**dict(dict(job=self.job, fixture=self.fixture, report=self.inventory,
            manifest=self.source, metadata=self.metadata, resources=self.resources,
            worker_sha256=digest(self.worker), executable=self.worker), **changes))

    def run_pair(self, admitted=None, **changes):
        # Synthetic instruction counts exercise record semantics only. Actual
        # Callgrind parsing and the paired process runner have separate controls.
        costs = {p: None if p == 'rho_solve' else 100 for p in PHASES}
        return run_record(admitted or self.admission(), **dict(dict(number=10, host_id='test', status='complete',
            native=self.report, profile=copy.deepcopy(self.report), costs=costs,
            process_wall_ns=self.item['process_wall_ns'], profile_wall_ns=self.item['process_wall_ns'],
            native_status='EXITED', profile_status='EXITED', manifest=self.source, executable=self.worker), **changes))

    def test_inventory_and_observed_pipeline_have_one_candidate_and_common_workload(self):
        admitted = self.admission()
        record = self.run_pair(admitted)
        self.assertEqual(record['candidate_id'], admitted['candidate']['candidate_id'])
        self.assertEqual(record['total_operations'], 1300)
        self.assertNotEqual(record['run_id'], record['profile_execution']['run_id'])
        self.assertEqual(record['profile_execution']['total_operations'], record['total_operations'])
        self.assertEqual(admitted['workload'], workload_manifest(self.fixture, input_law=INPUT_LAW,
            algorithm_seed=self.job['algorithm_seed'], resource_envelope=self.resources))
        timing = record['native_timing']
        self.assertEqual(sum(timing['cold']['phase_wall_ns'].values()), self.item['process_wall_ns'])
        self.assertEqual(sum(timing['online']['phase_wall_ns'].values()), timing['online']['wall_ns'])
        self.assertFalse(record['promotion_eligible'])

    def test_different_kernel_or_actual_dispatch_cannot_use_inventory_candidate(self):
        wrong = copy.deepcopy(self.inventory)
        wrong['field_kernel'] = 'portable' if wrong['field_kernel'] != 'portable' else 'pmull'
        with self.assertRaisesRegex(InvalidEvidence, 'method differs'):
            self.run_pair(self.admission(report=wrong))
        wrong = copy.deepcopy(self.report)
        wrong['collector_dispatch']['pair_table'] = False
        with self.assertRaises(InvalidEvidence):
            self.run_pair(profile=wrong)

    def test_inventory_rejects_target_work_and_changed_binary(self):
        for changes in (dict(online_wall_ns=1), dict(scalar_replay_included=True),
                        dict(inventory_policy='unknown'), dict(solutions=[{}])):
            with self.assertRaises(InvalidEvidence):
                self.admission(report=dict(self.inventory, **changes))
        admitted = self.admission()
        self.worker.write_bytes(b'changed')
        with self.assertRaises(InvalidEvidence):
            check_admission(admitted, job=self.job, fixture=self.fixture, manifest=self.source,
                metadata=self.metadata, resources=self.resources,
                worker_sha256=self.build['worker_sha256'], executable=self.worker)

    def test_failure_preserves_identity_and_process_state_with_unknown_cost(self):
        admitted = self.admission()
        for status in ('timeout', 'oom', 'error'):
            record = self.run_pair(admitted, status=status, native=None, profile=None, costs=None,
                process_wall_ns=None, native_status='NOT_RUN', profile_status='TIMEOUT')
            self.assertEqual(record['candidate_id'], admitted['candidate']['candidate_id'])
            self.assertIsNone(record['total_operations'])
            self.assertIsNone(record['native_timing'])
            self.assertIsNone(record['certificate'])
            self.assertEqual(record['native_process_status'], 'NOT_RUN')
            self.assertEqual(record['profile_execution']['process_status'], 'TIMEOUT')

    def test_short_process_clock_or_missing_profile_never_yields_a_cost(self):
        for changes in (dict(profile_wall_ns=1), dict(process_wall_ns=1), dict(profile=None),
                        dict(costs=dict.fromkeys(PHASES)), dict(number=True)):
            with self.assertRaises(InvalidEvidence):
                self.run_pair(**changes)

    def test_rho_keeps_reference_identity_and_has_no_fictitious_ic_stages(self):
        self.item = copy.deepcopy(next(r for r in retained() if r['job']['mode'] == 'rho'))
        self.job, self.report = self.item['job'], self.item['report']
        self.fixture = self.report['fixture']
        self.inventory = {k:copy.deepcopy(self.report[k]) for k in ('mode', 'fixture',
            'generic_build', 'generic_runtime_policy', 'generic_admission_schema', 'effective_config', 'field_kernel')}
        self.inventory.update(status='inventory', inventory_policy=INVENTORY_POLICY,
                              online_wall_ns=None, scalar_replay_included=False)
        admitted = self.admission()
        costs = {p:100 if p in {'setup','precompute','rho_solve','recovery_check'} else None for p in PHASES}
        record = self.run_pair(admitted, costs=costs)
        self.assertNotIn('candidate_id', record)
        self.assertNotIn('candidate', admitted)
        self.assertEqual(record['reference_id'], admitted['reference']['reference_id'])
        self.assertEqual(record['total_operations'], 400)
        self.assertEqual(record['rho_context']['public_target'], self.fixture['targets'][0])
        self.assertEqual(record['rho_context']['worker_count'], 1)
        self.assertEqual(record['rho_context']['requested_walks'], self.report['rho_dispatch']['requested_walks'])
        self.assertEqual(record['rho_context']['effective_walks'], self.report['solutions'][0]['effective_walks'])
        self.assertIsNone(record['rho_context']['distinguished_point_peak_bytes'])
        self.assertEqual(record['native_timing']['online']['wall_ns'], self.report['online_wall_ns'])
        self.assertNotEqual(record['run_id'], record['profile_execution']['run_id'])
        malformed = copy.deepcopy(self.report)
        malformed['rho_dispatch']['requested_walks'] = 1
        with self.assertRaises(InvalidEvidence):
            self.run_pair(admitted, profile=malformed, costs=costs)

    def test_profile_compression_round_trips_without_weakening_checksum(self):
        import gzip
        root = Path(self.temporary.name)/'profile'
        root.mkdir()
        entered = [p for p in PHASES if self.report['generic_phase_timing']['phases_ns'][p] is not None]
        for ordinal, phase in enumerate(entered, 1):
            trigger = 'Program termination' if phase == 'setup' else 'Client Request: generic_ic_'+phase
            (root/f'callgrind.out.{ordinal}').write_text(
                f'events: Ir\npart: {ordinal}\nsummary: 100\ntotals: 100\ndesc: Trigger: {trigger}\n')
        (root/'stderr.txt').write_text(f'Collected : {100*len(entered)}\n')
        original = parse_profiles(root, self.report, self.job)
        for path in root.glob('callgrind.out*'):
            path.with_suffix(path.suffix+'.gz').write_bytes(gzip.compress(path.read_bytes(), mtime=0))
            path.unlink()
        self.assertEqual(parse_profiles(root, self.report, self.job, compressed=True), original)
        with self.assertRaises(InvalidEvidence):
            parse_profiles(root, self.report, self.job)
        (root/'stderr.txt').write_text('Collected : 1\n')
        with self.assertRaisesRegex(InvalidEvidence, 'do not close'):
            parse_profiles(root, self.report, self.job, compressed=True)

    def test_job_transformation_is_explicit_and_leaves_prepared_fixture_unchanged(self):
        case = dict(job=dict(self.job, exclusive_phases=False))
        frozen = copy.deepcopy(case)
        arm = dict(id='generic', adapter=ADAPTER, config=self.job['config'])
        self.assertTrue(executed_job(case, arm)['exclusive_phases'])
        self.assertEqual(case, frozen)
        refs = qualification_references([dict(arm, source_manifest_sha256='a'*64)], [1,4])
        self.assertEqual([executed_job(case,r)['mode'] for r in refs], ['rho','rho'])


if __name__ == '__main__':
    unittest.main()
