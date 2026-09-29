"""Portable audit of the retained wide-S4 SAT pilot outcomes."""
import hashlib
import io
import json
import tarfile
import unittest

from continue_cms_s4_controls import STAGE_A_SHA256
from run_cms_s4_controls import PANEL, PANEL_SHA256, REGISTRATION
from tournament import read


STAGE_B = REGISTRATION/'stage-b-evidence.tar.gz'
STAGE_B_SHA256 = 'c9246c617ef5d60ecc5281631decb6721e985ca3604514b6546c6177f05b064b'
STAGE_B_SUMMARY_SHA256 = 'ffb5b334d62c6eeb8c807f7083e049534e1f5e0c70ce748f4f8b01a8e7402b07'
CONTINUATION_SHA256 = 'd97d794858ade8a81dfa45a66d20b9cb882f075a6244107688fa055c6b504933'


class CmsS4EvidenceTests(unittest.TestCase):
    def test_four_loader_failures_are_retained_without_sat_verdicts(self):
        blob = STAGE_B.read_bytes()
        self.assertEqual(hashlib.sha256(blob).hexdigest(), STAGE_B_SHA256)
        with tarfile.open(fileobj=io.BytesIO(blob), mode='r:gz') as archive:
            members = archive.getmembers()
            names = [member.name for member in members]
            self.assertEqual(len(names), len(set(names)))
            self.assertTrue(all(member.isfile() and not member.name.startswith('/')
                                and '..' not in member.name.split('/') for member in members))
            files = {member.name: archive.extractfile(member).read()
                     for member in members}

        def digest(name):
            return hashlib.sha256(files[name]).hexdigest()

        panel = read(PANEL)
        self.assertEqual(digest('original-stage-a/registered-panel.json'), PANEL_SHA256)
        self.assertEqual(digest('registered-continuation.py'), CONTINUATION_SHA256)
        self.assertEqual(digest('summary.json'), STAGE_B_SUMMARY_SHA256)
        self.assertEqual(digest('original-stage-a/cms-executable'),
                         panel['cms_executable_sha256'])
        self.assertEqual(hashlib.sha256((REGISTRATION/'stage-a-evidence.tar.gz')
                                        .read_bytes()).hexdigest(), STAGE_A_SHA256)
        original = json.loads(files['original-stage-a/summary.json'])
        result = json.loads(files['summary.json'])
        self.assertEqual(result['original_stage_a_archive_sha256'], STAGE_A_SHA256)
        self.assertEqual(result['original_stage_a_statuses'], ['INVALID_EXPORT']*4)
        self.assertEqual([row['status'] for row in original['rows']],
                         ['INVALID_EXPORT']*4)
        self.assertEqual(len(result['rows']), len(panel['schedule']))
        self.assertIs(result['source_bound_complete_sat_ic'], False)
        self.assertIsNone(result['natural_yield_estimate'])
        self.assertIsNone(result['online_speedup'])
        for item, row in zip(panel['schedule'], result['rows']):
            trial = item['trial']
            prefix = f'stage-b/trial-{trial:02d}/'
            self.assertEqual(row['trial'], trial)
            self.assertEqual(row['public_point'], item['point'])
            self.assertEqual(row['exact_relation_exists'], item['exact_relation_exists'])
            self.assertEqual(row['status'], 'SOLVER_ERROR')
            self.assertIsNone(row['source_model_valid'])
            self.assertIsNone(row['point_witness'])
            self.assertEqual(row['cms']['returncode'], -6)
            self.assertIs(row['cms']['timed_out'], False)
            self.assertEqual(row['cms'], json.loads(files[prefix+'cms.metrics.json']))
            self.assertEqual(files[prefix+'cms.stdout'], b'')
            self.assertIn(b'Library not loaded: @rpath/libcryptominisat5.5.14.dylib',
                          files[prefix+'cms.stderr'])
            self.assertIn(b'(no such file)', files[prefix+'cms.stderr'])


if __name__ == '__main__':
    unittest.main()
