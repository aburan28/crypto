"""Development points cannot reappear as fresh confirmation after a round boundary."""
import copy
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

from oracle import InvalidEvidence
from target_history import (extend, freeze_exposures, history_sets, key_for,
                            validate_fresh, verify_exposures)
from test_certificate import FIXTURE
from tournament import prepare
from test_campaign_rules import contract
from campaign_rules import validate_contract


class TargetExposureTests(unittest.TestCase):
    def test_later_contract_cannot_drop_the_exposure_replay_gate(self):
        original = contract()
        validate_contract(original)
        later = dict(original, attempt_number=2, seed=2026092552, target_exposure_schema=1)
        validate_contract(later)
        for value in (None, True, 2):
            with self.assertRaisesRegex(InvalidEvidence, 'sealed supplemental'):
                validate_contract(dict(later, target_exposure_schema=value))

    def test_retained_source_replays_without_original_and_blocks_changed_seed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / 'development.json'
            source.write_text(json.dumps({'controls': [{'fixture': FIXTURE}]}))
            data = source.read_bytes()
            output = root / 'frozen'
            baseline = dict(schema_version=1, sources=[], curves=[])
            original = copy.deepcopy(baseline)
            result = freeze_exposures(baseline, [source], output)
            self.assertEqual(original, baseline)
            self.assertEqual((output/'0000/fixtures.json').read_bytes(), data)
            source.unlink()
            receipt = verify_exposures(output, result)
            self.assertEqual(receipt['fixture_sources'], 1)
            self.assertEqual(receipt['excluded_points'], len(FIXTURE['targets']))
            changed = copy.deepcopy(FIXTURE)
            changed['target_seeds'] = [987654321] * len(changed['targets'])
            with self.assertRaisesRegex(InvalidEvidence, 'reused public point'):
                validate_fresh({'confirmation': [{'fixture': changed}]}, result)

    def test_omitted_point_changed_source_and_unlisted_source_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / 'development.json'
            source.write_text(json.dumps(FIXTURE))
            output = root / 'frozen'
            result = freeze_exposures(dict(schema_version=1, curves=[]), [source], output)
            omitted = copy.deepcopy(result)
            omitted['curves'][0]['points'].pop()
            with self.assertRaisesRegex(InvalidEvidence, 'exclusions differ'):
                verify_exposures(output, omitted)
            path = output / '0000/fixtures.json'
            data = path.read_bytes()
            path.write_text('{}')
            with self.assertRaisesRegex(InvalidEvidence, 'changed supplemental fixture'):
                verify_exposures(output, result)
            path.write_bytes(data)
            extra = output / '0001/fixtures.json'
            extra.parent.mkdir()
            extra.write_bytes(data)
            with self.assertRaisesRegex(InvalidEvidence, 'unlisted supplemental'):
                verify_exposures(output, result)

    def test_later_round_inherits_points_absent_from_prior_round_workload(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first, second = (copy.deepcopy(FIXTURE) for _ in range(2))
            first['targets'], second['targets'] = FIXTURE['targets'][:1], FIXTURE['targets'][1:]
            source = root / 'development.json'
            source.write_text(json.dumps(first))
            prior_history = freeze_exposures(dict(schema_version=1, curves=[]), [source], root/'exposures')
            prior = root / 'prior-round'
            prior.mkdir()
            (prior/'target-history.json').write_text(json.dumps(prior_history))
            (prior/'fixtures.json').write_text(json.dumps({'development': [{'fixture': second}]}))
            inherited = extend(dict(schema_version=1, curves=[]), [prior])
            self.assertEqual(history_sets(inherited)[key_for(FIXTURE)],
                             {tuple(map(int, p)) for p in FIXTURE['targets']})
            self.assertIn('target-history.json', inherited['round_additions'][0]['files'])
            self.assertEqual(verify_exposures(root/'exposures', prior_history)['status'], 'VERIFIED')

    def test_empty_or_mistyped_evidence_fails_before_any_measurement(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for ordinal, value in enumerate(({}, {'fixtures': []}, dict(FIXTURE, targets=[]))):
                source = root / f'source-{ordinal}.json'
                source.write_text(json.dumps(value))
                with self.assertRaisesRegex(InvalidEvidence, 'no exposed public targets'):
                    freeze_exposures(dict(schema_version=1, curves=[]), [source], root/str(ordinal))
            with self.assertRaisesRegex(InvalidEvidence, 'bounded campaign'):
                prepare(SimpleNamespace(out=root/'unused', exposed_fixtures=[source], attempt_number=0))
            self.assertFalse((root/'unused').exists())

    def test_zero_supplemental_sources_preserve_all_parent_points(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root/'source.json'
            source.write_text(json.dumps(FIXTURE))
            parent = freeze_exposures(dict(schema_version=1, curves=[]), [source], root/'first')
            result = freeze_exposures(parent, [], root/'next')
            self.assertEqual(history_sets(parent), history_sets(result))
            self.assertEqual(verify_exposures(root/'next', result)['fixture_sources'], 0)

    def test_coordinate_coercions_and_invalid_points_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for ordinal, coordinate in enumerate((True, 1.5, '-1', 1 << FIXTURE['degree'])):
                fixture = copy.deepcopy(FIXTURE)
                fixture['targets'][0][0] = coordinate
                source = root / f'source-{ordinal}.json'
                source.write_text(json.dumps(fixture))
                with self.assertRaises(InvalidEvidence):
                    freeze_exposures(dict(schema_version=1, curves=[]), [source], root/str(ordinal))


if __name__ == '__main__':
    unittest.main()
