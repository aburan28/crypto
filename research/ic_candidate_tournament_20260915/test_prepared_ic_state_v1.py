"""Actual accepted preparations and adversarial import controls; no native jobs."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

from identity import sha256
from oracle import InvalidEvidence
from prepared_ic_state_v1 import (accepted_files, create, ordinary_inputs,
                                  reconstruct, verify)


class PreparedICStateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.files = {family:accepted_files(family) for family in ('f5', 'sat')}
        cls.documents = {family:create(family, files) for family, files in cls.files.items()}

    def changed(self, family='f5'):
        return copy.deepcopy(self.documents[family])

    def test_accepted_families_independently_recover_identical_mathematical_state(self):
        f5, sat = (self.documents[family] for family in ('f5', 'sat'))
        self.assertEqual(f5['record'], sat['record'])
        self.assertEqual(f5['state_id'], 'ICP1hedbff76da644')
        self.assertEqual(f5['state_id'], sat['state_id'])
        self.assertNotEqual(sha256(f5), sha256(sat))
        for family, queries, rows in [('f5', 216, 61), ('sat', 149, 37)]:
            result = verify(self.documents[family], sha256(self.documents[family]))
            self.assertEqual((result['rank'], result['ordinary_queries'],
                              result['verified_relation_rows']), (29, queries, rows))
            self.assertFalse(result['target_input_present'] or result['native_execution']
                             or result['promotion_eligible'])
            self.assertIsNone(result['online_speedup'])

    def test_committed_certificates_equal_the_accepted_archive_derivations(self):
        root = Path(__file__).resolve().parent/'goal_20260924/prepared-ic-state-v1'
        for family, expected in [('f5', 'd8c5d6679fe89561606154785fdfba763835bf261dceb08fae193eff9f5636ae'),
                                 ('sat', '91856ab78550436d3f668367f9aebd9e2c0604bd64b1472d9d19ec318e2b144e')]:
            retained = json.loads((root/(family+'-preparation.json')).read_text())
            self.assertEqual(retained, self.documents[family])
            self.assertEqual(verify(retained, expected)['state_id'], 'ICP1hedbff76da644')

    def test_target_answers_workloads_and_timings_do_not_enter_preparation(self):
        for family in ('f5', 'sat'):
            files = copy.deepcopy(self.files[family])
            name = ('execution/entry-output/pipeline.stdout' if family == 'f5'
                    else 'execution/entry-output/summary.json')
            report = files[name]
            report.update(recovered_scalar='fake', target_attempts=['fake'],
                          solutions=['fake'], online_wall_ns=123)
            if family == 'f5':
                report['fixture'].update(targets=['fake'], target_seeds=['fake'])
            self.assertEqual(ordinary_inputs(family, files),
                             ordinary_inputs(family, self.files[family]))
        for document in self.documents.values():
            encoded = json.dumps(document['record'])
            for forbidden in ('targets', 'target_seeds', 'recovered_scalar', 'online_wall_ns',
                              'target_attempts', 'workload_id', 'run_id'):
                self.assertNotIn('"'+forbidden+'"', encoded)

    def test_target_dependent_native_relation_is_rejected_before_projection(self):
        files = copy.deepcopy(self.files['f5'])
        files['execution/entry-output/pipeline.stdout']['collection_reports'][0]['attempts'][0]['b'] = 1
        with self.assertRaisesRegex(InvalidEvidence, 'depends on a target'):
            ordinary_inputs('f5', files)

    def test_target_bearing_preparation_fixture_is_rejected(self):
        changed = self.changed()
        changed['certificate']['inputs']['fixture']['targets'] = [[52411, 72106]]
        with self.assertRaisesRegex(InvalidEvidence, 'contains target input'):
            verify(changed, sha256(changed))

    def test_wrong_log_or_column_order_is_rejected_even_with_new_external_seal(self):
        for family in ('f5', 'sat'):
            changed = self.changed(family)
            changed['certificate']['inputs']['logs'][0]['log'] = '63340'
            with self.assertRaisesRegex(InvalidEvidence, 'column or logarithm differs'):
                verify(changed, sha256(changed))
            changed = self.changed(family)
            logs = changed['certificate']['inputs']['logs']
            logs[0], logs[1] = logs[1], logs[0]
            with self.assertRaisesRegex(InvalidEvidence, 'column or logarithm differs'):
                verify(changed, sha256(changed))

    def test_geometry_cannot_drop_torsion_duplicate_or_reorder_points(self):
        for change in ('drop_torsion', 'duplicate', 'reorder'):
            changed = self.changed()
            base = changed['certificate']['inputs']['base']
            if change == 'drop_torsion':
                base.pop(0)
            elif change == 'duplicate':
                base[1] = base[0]
            else:
                base[1], base[3] = base[3], base[1]
            with self.subTest(change=change), self.assertRaises(InvalidEvidence):
                verify(changed, sha256(changed))

    def test_rank_loss_and_false_group_relation_cannot_certify_preparation(self):
        changed = self.changed()
        attempts = changed['certificate']['inputs']['attempts']
        witness = next(item for item in attempts if item['outcome'] == 'witness')
        witness['scalar'] = (witness['scalar']+1) % 65587
        with self.assertRaisesRegex(InvalidEvidence, 'does not re-add'):
            verify(changed, sha256(changed))
        inputs = copy.deepcopy(self.documents['sat']['certificate']['inputs'])
        inputs['attempts'] = inputs['attempts'][:1]
        with self.assertRaisesRegex(InvalidEvidence, 'rank deficient'):
            reconstruct(inputs)

    def test_inconclusive_queries_cannot_contribute_rows(self):
        changed = self.changed('sat')
        row = next(item for item in changed['certificate']['inputs']['attempts']
                   if item['outcome'] == 'CONFLICT_BUDGET_INCONCLUSIVE')
        row['indices'] = [0, 1, 2]
        with self.assertRaisesRegex(InvalidEvidence, 'failed ordinary query'):
            verify(changed, sha256(changed))

    def test_wrong_external_seal_provenance_or_promotion_is_rejected(self):
        with self.assertRaisesRegex(InvalidEvidence, 'external certificate seal'):
            verify(self.documents['f5'], '0'*64)
        for field, value in [('archive_sha256', '0'*64),
                             ('parent_candidate_record_sha256', '84b504d19844'+'0'*52)]:
            changed = self.changed()
            changed['provenance'][field] = value
            with self.assertRaisesRegex(InvalidEvidence, 'provenance differs'):
                verify(changed, sha256(changed))
        changed = self.changed()
        changed['online_speedup'] = '100'
        with self.assertRaisesRegex(InvalidEvidence, 'measurement or promotion'):
            verify(changed, sha256(changed))

    def test_valid_mathematical_preparation_cannot_swap_source_provenance(self):
        changed = self.changed()
        changed['certificate'] = copy.deepcopy(self.documents['sat']['certificate'])
        with self.assertRaisesRegex(InvalidEvidence, 'inputs differ from accepted source'):
            verify(changed, sha256(changed))

    def test_altered_bundle_receipt_cannot_change_accepted_archive_binding(self):
        # This rejects before reading/decompressing another archive or writing.
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            (bundle/'evidence.tar.gz').write_bytes(b'fake archive')
            (bundle/'receipt.json').write_text(json.dumps(dict(
                archive_bytes=12, archive_sha256='0'*64, execution_sha256='0'*64)))
            with self.assertRaisesRegex(InvalidEvidence, 'accepted external binding'):
                accepted_files('f5', bundle)


if __name__ == '__main__':
    unittest.main()
