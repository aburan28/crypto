"""Check the fresh F5/rho registration against retained source and target."""
import json
from pathlib import Path
import tempfile
import unittest

from register_paired_generic import identities
from run_paired_generic import run_worker
from static_sat_registration_v2 import REGISTRATION as SAT_REGISTRATION
from tournament import read


class FreshGenericRegistrationTests(unittest.TestCase):
    def test_worker_process_receipt_captures_exit_and_peak_rss(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            receipt = run_worker(Path('/bin/cat'), {'probe': 7}, root)
            self.assertEqual(receipt['exit_code'], 0)
            self.assertGreater(receipt['process_wall_ns'], 0)
            self.assertEqual(json.loads((root/'stdout.json').read_text()),
                             {'probe': 7})
            self.assertIsNotNone(receipt['memory_peak_bytes'])
            self.assertGreater(receipt['memory_peak_bytes'], 0)

    def test_source_bound_f5_and_rho_share_the_frozen_public_point(self):
        data = identities()
        sat = read(SAT_REGISTRATION/'panel.json')
        point = [str(value) for value in sat['target_input']['point']]
        self.assertEqual(data['f5_job']['public_targets'], [point])
        self.assertEqual(data['rho_job']['public_targets'], [point])
        self.assertEqual(data['f5_job']['target_seeds'],
                         [sat['target_input']['seed']])
        self.assertEqual(data['rho_job']['target_seeds'],
                         [sat['target_input']['seed']])
        self.assertEqual(data['f5_job']['algorithm_seed'],
                         sat['relation_query_seed'])
        self.assertEqual(data['candidate']['record']['factor_base']
                         ['inventory']['usable_point_count'], 62)
        self.assertEqual(data['candidate']['record']['factor_base']
                         ['inventory']['effective_columns'], 29)
        self.assertEqual(data['candidate']['record']['point_decomposition']
                         ['solver'], 'f5')
        self.assertEqual(data['panel']['worker_sha256'],
                         data['build']['worker_sha256'])
        self.assertIsNone(data['panel']['resource_envelope']
                          ['total_wall_limit_seconds'])


if __name__ == '__main__':
    unittest.main()
