"""Scientific producer regressions use a retained real optimized IC report."""
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from identity import sha256
from oracle import InvalidEvidence, verify
from producer.evidence import audit_stages, check_build_identity, scientific_ledger
from producer.prepare import verify_source
from producer.timing import native_intervals, RAW_PHASES
from tournament import parse_profiles

VECTOR = json.loads((Path(__file__).parent/'producer/testdata/n13.json').read_text())


class StageEvidenceTests(unittest.TestCase):
    def audit(self, report):
        return audit_stages(report, VECTOR['report']['fixture'], VECTOR['job']['algorithm_seed'])

    def test_real_optimized_report_replays_queries_matrix_and_rank(self):
        result = self.audit(VECTOR['report'])
        self.assertEqual(result['queries']['ordinary_queries'], 7)
        self.assertEqual(result['queries']['final_rank'], 7)
        self.assertEqual(result['matrix_nonzeros'], 17)
        self.assertEqual(result['matrix_sha256'], '4c74b498fc0148e5cc1fcd831d884648965e04e34e3caf3e422d0e1f618cafbc')

    def test_query_rank_and_matrix_counter_corruption_rejected(self):
        for key in ('ordinary_queries', 'identity_queries', 'pdp_attempts', 'final_rank', 'matrix_nonzeros'):
            report = copy.deepcopy(VECTOR['report'])
            report['diagnostics'][key] += 1
            with self.subTest(key=key), self.assertRaises(InvalidEvidence):
                self.audit(report)
        report = copy.deepcopy(VECTOR['report'])
        report['diagnostics']['rank_events'][0][2] = 0
        with self.assertRaises(InvalidEvidence):
            self.audit(report)

    def test_rank_replay_reads_coefficients_after_previous_elimination(self):
        vector = json.loads((Path(__file__).parent/'producer/testdata/n23.json').read_text())
        report = vector['report']
        result = audit_stages(report, report['fixture'], vector['job']['algorithm_seed'])
        self.assertEqual(result['queries']['final_rank'], 8)
        self.assertEqual(result['rank_events'], [[i, i, i+1] for i in range(8)])
        # Row six needs a coefficient created by an earlier elimination. The
        # stale-list iterator incorrectly reported this row as dependent.
        report['diagnostics']['rank_events'][5][2] = 5
        with self.assertRaises(InvalidEvidence):
            audit_stages(report, report['fixture'], vector['job']['algorithm_seed'])

    def test_correct_relation_cannot_be_relabelled_as_another_query(self):
        report = copy.deepcopy(VECTOR['report'])
        # The group relation remains correct, but its claimed query position is
        # now wrong. The old complete-scalar checker alone would accept it.
        report['relations'][0]['trial'] = 6
        with self.assertRaises(InvalidEvidence):
            self.audit(report)

    def test_bounded_misses_cannot_be_relabelled_unsat(self):
        report = copy.deepcopy(VECTOR['report'])
        report['diagnostics']['pdp_outcomes']['proved_unsat'] = 1
        with self.assertRaises(InvalidEvidence):
            self.audit(report)

    def test_declared_backend_must_match_real_dispatch(self):
        report = copy.deepcopy(VECTOR['report'])
        report['executed_method']['relation_la'] = 'block_wiedemann'
        with self.assertRaises(InvalidEvidence):
            self.audit(report)

    def test_stale_build_and_changed_field_dispatch_are_rejected(self):
        report = {'diagnostics': {'source_manifest_sha256': '1'*64, 'field_kernel': 'portable'}}
        check_build_identity(report, '1'*64, 'portable')
        with self.assertRaises(InvalidEvidence):
            check_build_identity(report, '2'*64, 'portable')
        with self.assertRaises(InvalidEvidence):
            check_build_identity(report, '1'*64, 'pclmulqdq')
        with self.assertRaises(InvalidEvidence):
            check_build_identity({}, '1'*64)


class ScientificProfilesTests(unittest.TestCase):
    def profile(self, root, labels):
        for part, label in enumerate(labels, 1):
            trigger = 'Program termination' if label is None else 'Client Request: '+label
            (root/f'callgrind.out.{part}').write_text(
                f'part: {part}\nevents: Ir\nsummary: 10\ntotals: 10\ndesc: Trigger: {trigger}\n')
        (root/'stderr.txt').write_text(f'Collected : {10*len(labels)}\n')

    def test_repeated_exclusive_phases_and_termination_close(self):
        labels = ['ic_'+phase for phase in ('setup', 'factor_base', 'precompute', 'queries', 'pdp',
            'relation_check', 'matrix_build', 'relation_la', 'target_descent', 'recovery_check')]
        labels += ['ic_pdp', None]
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self.profile(root, labels)
            costs = parse_profiles(root, phase_schema=2)
        ledger = scientific_ledger(VECTOR['report'], costs)
        self.assertEqual(ledger['cold_operations'], 120)
        self.assertEqual(costs['setup'], 20)
        self.assertEqual(costs['pdp'], 20)
        self.assertEqual(ledger['operations']['isogeny'], 0)

    def test_legacy_labels_missing_phases_and_bad_checksum_fail(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self.profile(root, ['factor_base_and_tables', None])
            with self.assertRaises(InvalidEvidence):
                parse_profiles(root, phase_schema=2)
            self.profile(root, ['ic_setup', None])
            costs = parse_profiles(root, phase_schema=2)
            with self.assertRaises(InvalidEvidence):
                scientific_ledger(VECTOR['report'], costs)
            (root/'stderr.txt').write_text('Collected : 21\n')
            with self.assertRaises(InvalidEvidence):
                parse_profiles(root, phase_schema=2)

    def test_target_subphases_fold_into_cold_descent_once(self):
        report = copy.deepcopy(VECTOR['report'])
        report['phase_schema'] = 3
        from measurement import PHASES
        phases = sorted((set(PHASES)-{'isogeny'}) |
                        {'target_query', 'target_pdp', 'target_relation_check'})
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self.profile(root, ['ic_'+p for p in phases] + [None])
            costs = parse_profiles(root, phase_schema=3)
            with self.assertRaises(InvalidEvidence):
                parse_profiles(root, phase_schema=2)
        ledger = scientific_ledger(report, costs)
        self.assertEqual(ledger['cold_operations'], 140)
        self.assertEqual(ledger['operations']['target_descent'], 40)
        self.assertEqual(ledger['operations']['setup'], 20)
        del costs['target_pdp']
        with self.assertRaises(InvalidEvidence):
            scientific_ledger(report, costs)

    def test_source_manifest_and_contents_both_must_match(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root/'source').mkdir()
            code = root/'source/test.rs'
            code.write_bytes(b'original')
            manifest = {'test.rs': hashlib.sha256(b'original').hexdigest()}
            (root/'source-manifest.json').write_text(json.dumps(manifest))
            self.assertEqual(verify_source(root, sha256(manifest)), manifest)
            with self.assertRaises(InvalidEvidence):
                verify_source(root, '0'*64)
            code.write_bytes(b'changed')
            with self.assertRaises(InvalidEvidence):
                verify_source(root, sha256(manifest))


class NativeIntervalTests(unittest.TestCase):
    def test_real_public_point_ic_and_rho_vectors_replay_and_close(self):
        vector = json.loads((Path(__file__).parent/'producer/testdata/n13-public.json').read_text())
        fixture = vector['ic']['report']['fixture']
        for mode in ('ic', 'rho'):
            item = vector[mode]
            self.assertEqual(item['job']['public_targets'], fixture['targets'])
            verify(item['report'], fixture, expected_mode=mode)
            result = native_intervals(item['report'], item['process']['process_wall_ns'])
            self.assertEqual(result['online']['wall_ns'], sum(result['online']['phase_wall_ns'].values()))
            self.assertEqual(result['cold']['wall_ns'], sum(result['cold']['phase_wall_ns'].values()))
        audit_stages(vector['ic']['report'], fixture, vector['ic']['job']['algorithm_seed'])

    def report(self, mode='ic'):
        phases = dict.fromkeys(RAW_PHASES, 0)
        phases.update(setup=100, recovery_check=5)
        if mode == 'ic':
            phases.update(factor_base=30, precompute=40, queries=7, pdp=20,
                          relation_check=3, matrix_build=5, relation_la=15,
                          target_query=2, target_pdp=8, target_relation_check=3,
                          target_descent=7)
        else:
            phases['reference_solve'] = 20
        report = {'phase_schema': 3, 'status': 'complete', 'mode': mode,
                  'fixture': {'targets': [['1', '2']]}}
        intervals = {'target_input': 'supplied_public_point',
                     'phase_wall_ns': phases, 'online_wall_ns': 25}
        if mode == 'ic':
            report['diagnostics'] = intervals
        else:
            report.update(intervals)
        return report

    def test_online_excludes_setup_and_cold_accounts_for_external_tail(self):
        for mode in ('ic', 'rho'):
            with self.subTest(mode=mode):
                report = self.report(mode)
                result = native_intervals(report, 400)
                self.assertEqual(result['online']['wall_ns'], 25)
                self.assertEqual(sum(result['cold']['phase_wall_ns'].values()), 400)
                self.assertGreater(result['cold']['external_setup_remainder_ns'], 0)
                self.assertTrue(result['online']['scalar_replay_included'])
                self.assertFalse(result['online']['target_generation_included'])
                self.assertEqual(result['unit'], 'native_monotonic_ns')
                if mode == 'ic':
                    self.assertEqual(result['cold']['phase_wall_ns']['target_descent'], 20)

    def test_incomplete_multiple_generated_and_unknown_intervals_fail_closed(self):
        mutations = [
            lambda r: r.update(phase_schema=2),
            lambda r: r.update(status='incomplete'),
            lambda r: r['fixture']['targets'].append(['3', '4']),
            lambda r: r['diagnostics'].update(target_input='generated'),
            lambda r: r['diagnostics'].update(online_wall_ns=24),
            lambda r: r['diagnostics']['phase_wall_ns'].pop('target_query'),
            lambda r: r['diagnostics']['phase_wall_ns'].update(reference_solve=1),
            lambda r: r['diagnostics']['phase_wall_ns'].update(setup=401),
        ]
        for index, mutate in enumerate(mutations):
            report = self.report()
            mutate(report)
            with self.subTest(index=index), self.assertRaises(InvalidEvidence):
                native_intervals(report, 400)

    def test_interval_units_require_nonnegative_integer_nanoseconds(self):
        for value in (None, -1, True, 1.5, '10'):
            report = self.report()
            report['diagnostics']['phase_wall_ns']['target_query'] = value
            with self.subTest(value=value), self.assertRaises(InvalidEvidence):
                native_intervals(report, 400)
        report = self.report('rho')
        report['phase_wall_ns']['target_pdp'] = 1
        with self.assertRaises(InvalidEvidence):
            native_intervals(report, 400)


if __name__ == '__main__':
    unittest.main()
