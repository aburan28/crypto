"""Seal checks for the reconstructed censored-seed 2026092902 exposure census."""
import hashlib
import json
from pathlib import Path
import unittest

from target_history import validate_exposure_source

HERE = Path(__file__).resolve().parent
REG = HERE / 'goal_20260924/generic-backend-qualification-v2'
EXPORT = REG / 'lost-v2-campaign-exposures.json'
EXPORT_SHA256 = '0cc792cceb8c7190a533e6f4e665e8486911153ad57609ac8d1243af282b8e54'
PRIOR = REG / 'lost-campaign-exposures.json'
PRIOR_SHA256 = 'a728677b199eac02800d8338204d5306f391ec5da757c910bec1e51955fe7b41'
INDEPENDENT = REG / 'independent-v2-prepare-replay-fixtures.json'
INDEPENDENT_SHA256 = 'e7af64f1afcda01c2db4238375d7467adc651fd6513acd5d9c7e7434eff95f6d'
RECEIPT = REG / 'v2-exposure-reconstruction-receipt.json'


class V2ExposureReconstructionTests(unittest.TestCase):
    def test_sealed_export_hash_and_schedule(self):
        digest = hashlib.sha256(EXPORT.read_bytes()).hexdigest()
        self.assertEqual(digest, EXPORT_SHA256)
        validate_exposure_source(EXPORT.read_bytes())
        data = json.loads(EXPORT.read_text())
        self.assertEqual(data['schema_version'], 1)
        self.assertEqual(data['seed'], 2026092902)
        self.assertEqual(data['source_run_id'], 36580669479)
        self.assertEqual(data['source_commit'],
                         '765c3c5f19032bd852163805f257c56babef2040')
        self.assertEqual(data['prior_censored_exposures_sha256'], PRIOR_SHA256)
        accepted = [row for row in data['attempts'] if row['accepted']]
        self.assertEqual(len(accepted), 25)
        self.assertEqual(
            {stage: sum(1 for row in accepted if row['stage'] == stage)
             for stage in ('aa', 'smoke', 'development')},
            dict(aa=5, smoke=5, development=15))

    def test_no_overlap_with_first_run_accepted_points(self):
        self.assertEqual(hashlib.sha256(PRIOR.read_bytes()).hexdigest(), PRIOR_SHA256)
        prior = {(row['fixture']['degree'], row['fixture']['curve_a'],
                  tuple(map(int, row['fixture']['targets'][0])))
                 for row in json.loads(PRIOR.read_text())['attempts'] if row['accepted']}
        current = {(row['fixture']['degree'], row['fixture']['curve_a'],
                    tuple(map(int, row['fixture']['targets'][0])))
                   for row in json.loads(EXPORT.read_text())['attempts'] if row['accepted']}
        self.assertEqual(len(prior), 25)
        self.assertEqual(len(current), 25)
        self.assertFalse(prior & current)

    def test_prepare_replay_export_matches_sealed_points(self):
        self.assertEqual(hashlib.sha256(INDEPENDENT.read_bytes()).hexdigest(),
                         INDEPENDENT_SHA256)
        sealed = {(row['stage'], row['case']): row['fixture']
                  for row in json.loads(EXPORT.read_text())['attempts'] if row['accepted']}
        replay = json.loads(INDEPENDENT.read_text())
        matched = 0
        for stage, cases in replay.items():
            for case in cases:
                self.assertEqual(sealed[(stage, case['id'])], case['fixture'])
                matched += 1
        self.assertEqual(matched, 25)
        receipt = json.loads(RECEIPT.read_text())
        self.assertEqual(receipt['export_sha256'], EXPORT_SHA256)
        self.assertEqual(receipt['independent_export_sha256'], INDEPENDENT_SHA256)
        self.assertTrue(receipt['accepted_points_equal'])


if __name__ == '__main__':
    unittest.main()
