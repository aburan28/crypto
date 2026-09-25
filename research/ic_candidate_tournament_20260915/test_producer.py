"""Scientific producer regressions use a retained real optimized IC report."""
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from identity import sha256
from oracle import InvalidEvidence
from producer.evidence import audit_stages, check_build_identity, scientific_ledger
from producer.prepare import verify_source
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


if __name__ == '__main__':
    unittest.main()
