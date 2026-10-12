"""Checks that future comparisons cannot silently change the benchmark target."""
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from compare import compare, load_run
from run import aggregate
from run_pilot import GF

HERE = Path(__file__).resolve().parent
BASELINE = HERE / 'results' / 'baseline_v2'


class RegressionTests(unittest.TestCase):
    def test_algebraic_only_trace_certificate(self):
        certificate = json.loads((HERE / 'negative_control.json').read_text())
        field = GF(certificate['field_degree'], certificate['field_modulus_integer'])
        self.assertEqual(field.s4(*certificate['x_witness'], certificate['target_x']), 0)
        for item in certificate['trace_certificates']:
            rhs = item['x'] ^ 1 ^ field.sq(field.inv(item['x']))
            trace = z = rhs
            for _ in range(1, certificate['field_degree']):
                z = field.sq(z)
                trace ^= z
            self.assertEqual(trace, 1)
            self.assertEqual(rhs, item['artin_schreier_rhs'])
            self.assertEqual(field.lift(item['x']), [])

    def candidate(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        path = Path(temp.name)
        for name in ('summary.json', 'raw.jsonl', 'contract.json', 'corpus_manifest.json'):
            shutil.copyfile(BASELINE / name, path / name)
        return path

    def rewrite(self, path, change):
        records = [json.loads(line) for line in (path / 'raw.jsonl').read_text().splitlines()]
        summary = json.loads((path / 'summary.json').read_text())
        change(records, summary)
        raw = ''.join(json.dumps(r) + '\n' for r in records)
        (path / 'raw.jsonl').write_text(raw)
        summary['raw_sha256'] = hashlib.sha256(raw.encode()).hexdigest()
        summary['cells'], summary['groups'] = aggregate(records, json.loads((path / 'contract.json').read_text()))
        (path / 'summary.json').write_text(json.dumps(summary))

    def test_baseline_is_a_complete_repeatable_target(self):
        report = compare(BASELINE, BASELINE)
        self.assertEqual(len(report['per_case']), 120)
        for variant in report['variants']:
            self.assertEqual(variant['conflict_ratio_candidate_over_baseline'], 1)
            self.assertFalse(variant['diagnostic_20_percent_counter_target_met'])
            self.assertIsNone(variant['full_dlp_S'])
        summary, rows = load_run(BASELINE)
        self.assertEqual(len(rows), 360)
        self.assertEqual(sum(r['point_relation_verified'] for r in rows), 180)
        self.assertEqual(sum(r['status'] == 'SAT_ALGEBRAIC_ONLY' for r in rows), 6)
        self.assertEqual(summary['exhaustive_s4_checks'], 929280)

    def test_missing_or_duplicate_trials_are_rejected(self):
        for action in ('missing', 'duplicate'):
            with self.subTest(action=action):
                path = self.candidate()
                self.rewrite(path, lambda rs, s: rs.pop() if action == 'missing' else rs.append(rs[0]))
                with self.assertRaisesRegex(ValueError, 'trial slots'):
                    compare(BASELINE, path)

    def test_modified_inputs_are_rejected(self):
        path = self.candidate()
        self.rewrite(path, lambda rs, s: rs[0].update(normalized_sha256='changed'))
        with self.assertRaisesRegex(ValueError, 'input fingerprint'):
            compare(BASELINE, path)

    def test_counter_unit_changes_are_rejected(self):
        path = self.candidate()
        self.rewrite(path, lambda rs, s: s['metadata'].update(counter_unit='different-counter'))
        with self.assertRaisesRegex(ValueError, 'counter_unit'):
            compare(BASELINE, path)

    def test_non_lifting_control_cannot_be_counted_as_relation(self):
        path = self.candidate()
        def change(rs, summary):
            next(r for r in rs if r['status'] == 'SAT_ALGEBRAIC_ONLY')['point_relation_verified'] = True
        self.rewrite(path, change)
        with self.assertRaisesRegex(ValueError, 'mislabeled as a relation'):
            compare(BASELINE, path)

    def test_counter_regressions_are_reported(self):
        path = self.candidate()
        def change(rs, summary):
            for row in rs:
                row['conflicts'] *= 2
                lines = row['stdout'].strip().splitlines()
                lines[-1] = str(row['conflicts'])
                row['stdout'] = '\n'.join(lines) + '\n'
        self.rewrite(path, change)
        result = compare(BASELINE, path)
        for variant in result['variants']:
            self.assertEqual(variant['conflict_ratio_candidate_over_baseline'], 2)
            self.assertTrue(variant['cases_regressing_over_10_percent'])
            self.assertFalse(variant['diagnostic_20_percent_counter_target_met'])

    def test_counter_must_match_raw_stdout(self):
        path = self.candidate()
        self.rewrite(path, lambda rs, s: rs[0].update(conflicts=rs[0]['conflicts'] + 1))
        with self.assertRaisesRegex(ValueError, 'does not match stdout'):
            compare(BASELINE, path)

    def test_corrupted_raw_evidence_is_rejected(self):
        path = self.candidate()
        with (path / 'raw.jsonl').open('a') as stream:
            stream.write('{}\n')
        with self.assertRaisesRegex(ValueError, 'integrity mismatch'):
            compare(BASELINE, path)


if __name__ == '__main__':
    unittest.main()
