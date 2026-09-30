"""Seal the reconstructed public points from censored seed 2026092902."""
import hashlib
import json
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REGISTRATION = HERE / 'goal_20260924/generic-backend-qualification-v2'
THIS_RUN = REGISTRATION / 'this-run-exposures.json'
THIS_RUN_SHA256 = '64f30e19f0ef1a4c4b95d4168c977b13de9559cb8a934cfcc9c3bf8a63fa0a25'
PRIOR = REGISTRATION / 'lost-campaign-exposures.json'
PRIOR_SHA256 = 'a728677b199eac02800d8338204d5306f391ec5da757c910bec1e51955fe7b41'


class ThisRunExposures(unittest.TestCase):
    def test_sealed_corpus_identity_and_disjointness(self):
        digest = hashlib.sha256(THIS_RUN.read_bytes()).hexdigest()
        self.assertEqual(digest, THIS_RUN_SHA256)
        self.assertEqual(hashlib.sha256(PRIOR.read_bytes()).hexdigest(), PRIOR_SHA256)
        data = json.loads(THIS_RUN.read_text())
        prior = json.loads(PRIOR.read_text())
        self.assertEqual(data['schema_version'], 1)
        self.assertEqual(data['seed'], 2026092902)
        self.assertEqual(data['source_run_id'], 36580669479)
        self.assertEqual(data['source_commit'],
                         '765c3c5f19032bd852163805f257c56babef2040')
        self.assertEqual(data['prior_censored_exposures_sha256'], PRIOR_SHA256)
        accepted = [row for row in data['attempts'] if row['accepted']]
        self.assertEqual(len(accepted), 25)
        self.assertEqual(len(data['attempts']), 25)
        stages = {stage: sum(1 for row in accepted if row['stage'] == stage)
                  for stage in ('aa', 'smoke', 'development')}
        self.assertEqual(stages, dict(aa=5, smoke=5, development=15))
        points = {tuple(map(int, row['fixture']['targets'][0])) for row in accepted}
        self.assertEqual(len(points), 25)
        prior_points = {tuple(map(int, row['fixture']['targets'][0]))
                        for row in prior['attempts'] if row['accepted']}
        self.assertEqual(len(prior_points), 25)
        self.assertFalse(points & prior_points)


if __name__ == '__main__':
    unittest.main()
