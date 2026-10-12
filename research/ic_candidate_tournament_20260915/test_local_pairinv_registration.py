"""Guard the local incumbent's public target and accepted-source gate."""
import json
from pathlib import Path
import tempfile
import unittest

from oracle import InvalidEvidence
from register_local_pairinv import TARGET, build_inputs, jobs
from static_sat_registration_v2 import REGISTRATION as SAT_REGISTRATION
from tournament import read


class LocalPairinvRegistrationTests(unittest.TestCase):
    def test_incumbent_and_rho_receive_the_same_frozen_public_point(self):
        target = read(TARGET)['target_input']
        sat = read(SAT_REGISTRATION/'panel.json')['target_input']
        ic, rho = jobs()
        self.assertEqual(target['point'], sat['point'])
        self.assertFalse(target['known_scalar_supplied'])
        for job, mode in ((ic, 'ic'), (rho, 'rho')):
            self.assertEqual(job['mode'], mode)
            self.assertEqual(job['public_targets'],
                             [[str(v) for v in target['point']]])
            self.assertEqual(job['target_seeds'], [target['seed']])
            self.assertEqual(job['factor_base'],
                             {'kind': 'subgroup_orbits', 'seed': 43,
                              'points': 102})
            self.assertNotIn('target_scalar', job)
        self.assertEqual(rho['config']['rho_parallel_walks'], 4)

    def test_unaccepted_prepared_source_cannot_enter_local_build(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prepared, built = root/'prepared', root/'built'
            prepared.mkdir()
            built.mkdir()
            (prepared/'preparation.json').write_text(json.dumps({
                'reference': 'scaled',
                'source_manifest_sha256': '0'*64}))
            with self.assertRaisesRegex(InvalidEvidence,
                                        'prepared source is not accepted'):
                build_inputs(prepared, built)


if __name__ == '__main__':
    unittest.main()
