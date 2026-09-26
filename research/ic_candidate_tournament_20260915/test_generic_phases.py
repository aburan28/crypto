"""Adversarial closure controls. Real-worker records are checked by CI integration."""
from pathlib import Path
import gzip
import hashlib
import json
import tempfile
import unittest

from generic_phases import ONLINE, PHASES, REQUIRED, parse_profiles, verify_native
from oracle import InvalidEvidence


def fixture():
    cold = {phase: 10 if phase in REQUIRED['ic'] else None for phase in PHASES}
    online = {phase: cold[phase] for phase in ONLINE}
    job = dict(mode='ic', exclusive_phases=True, public_targets=[['1', '2']])
    report = dict(mode='ic', status='complete', fixture=dict(targets=job['public_targets']),
                  generic_phase_policy='exclusive-owner-thread-v1', online_timing_schema=2,
                  reusable_setup_excluded=True, target_input='supplied_public_point',
                  scalar_replay_included=True, outer_online_wall_ns=45, online_wall_ns=50,
                  generic_phase_timing=dict(schema_version=1, phases_ns=cold, observed_wall_ns=130,
                                            online_phases_ns=online, online_wall_ns=50))
    return job, report


class GenericPhaseTests(unittest.TestCase):
    def test_retained_worker_intervals_replay_without_trusting_receipts(self):
        directory = Path(__file__).parent / 'goal_20260924/generic-exclusive-phases'
        for path in (directory / 'worker-raw.jsonl', directory / 'final/worker-raw.jsonl'):
            checked = 0
            for line in path.read_text().splitlines():
                row = json.loads(line)
                if row['mode'] == 'exclusive':
                    receipt = verify_native(row['report'], row['job'], process_wall_ns=row['process_wall_ns'])
                    self.assertFalse(receipt['promotion_eligible'])
                    checked += 1
            self.assertEqual(checked, 147, str(path))

    def test_measured_sources_are_preserved_after_guard_and_test_changes(self):
        directory = Path(__file__).parent / 'goal_20260924/generic-exclusive-phases'
        expected = json.loads((directory / 'environment.json').read_text())['sources']
        retained = {
            'src/cryptanalysis/koblitz_index_calculus.rs': gzip.decompress((directory / 'initial-library-source.rs.gz').read_bytes()),
            'src/cryptanalysis/ic_measurement.rs': (directory / 'initial-measurement.rs').read_bytes(),
            'examples/ic_tournament_worker.rs': (directory / 'initial-worker.rs').read_bytes(),
            'research/ic_candidate_tournament_20260915/goal_20260924/generic-exclusive-phases/PROTOCOL.md':
                (directory / 'initial-protocol.md').read_bytes(),
        }
        for name, data in retained.items():
            self.assertEqual(hashlib.sha256(data).hexdigest(), expected[name], name)
        final = directory / 'final'
        final_expected = json.loads((final / 'environment.json').read_text())['sources']
        self.assertEqual(
            hashlib.sha256(gzip.decompress((final / 'measured-library-source.rs.gz').read_bytes())).hexdigest(),
            final_expected['src/cryptanalysis/koblitz_index_calculus.rs'])

    def test_external_remainder_is_charged_exactly_once_to_setup(self):
        job, report = fixture()
        receipt = verify_native(report, job, process_wall_ns=200)
        self.assertEqual(receipt['external_setup_remainder_ns'], 70)
        self.assertEqual(receipt['process_phases_ns']['setup'], 80)
        self.assertEqual(sum(v for v in receipt['process_phases_ns'].values() if v is not None), 200)
        self.assertEqual(receipt['cold_wall_ns'], 200)
        self.assertFalse(receipt['promotion_eligible'])

    def test_complete_result_cannot_hide_missing_or_misattributed_phases(self):
        mutations = (
            lambda r: r['generic_phase_timing']['phases_ns'].update(pdp=None),
            lambda r: r['generic_phase_timing']['phases_ns'].update(pdp=True),
            lambda r: r['generic_phase_timing']['online_phases_ns'].update(target_pdp=9),
            lambda r: r.update(online_wall_ns=49),
            lambda r: r.update(outer_online_wall_ns=51),
            lambda r: r.update(online_timing_schema=1),
            lambda r: r.update(scalar_replay_included=False),
        )
        for mutation in mutations:
            job, report = fixture()
            mutation(report)
            with self.assertRaises(InvalidEvidence):
                verify_native(report, job, process_wall_ns=200)

    def test_unstarted_target_and_unattempted_la_remain_unknown(self):
        job, report = fixture()
        trace = report['generic_phase_timing']
        for phase in (*ONLINE, 'relation_la'):
            trace['phases_ns'][phase] = None
        trace.update(online_phases_ns=dict.fromkeys(ONLINE), online_wall_ns=None, observed_wall_ns=70)
        report.update(status='incomplete', online_wall_ns=None, outer_online_wall_ns=None,
                      scalar_replay_included=False)
        receipt = verify_native(report, job, process_wall_ns=100)
        self.assertIsNone(receipt['cold_wall_ns'])
        self.assertIsNone(receipt['process_phases_ns']['relation_la'])
        self.assertIn('relation_la', receipt['missing_required_phases'])

    def write_profiles(self, root, report):
        entered = [p for p, v in report['generic_phase_timing']['phases_ns'].items() if v is not None]
        for part, phase in enumerate(entered, 1):
            (root / f'callgrind.out.{part}').write_text(
                f'part: {part}\ndesc: Trigger: Client Request: generic_ic_{phase}\n'
                'events: Ir\nsummary: 10\ntotals: 10\n')
        (root / 'callgrind.out').write_text(
            f'part: {len(entered) + 1}\ndesc: Trigger: Program termination\n'
            'events: Ir\nsummary: 20\ntotals: 20\n')
        (root / 'stderr.txt').write_text(f'==9== Collected : {len(entered) * 10 + 20}\n')

    def test_profiler_closure_includes_the_process_tail(self):
        job, report = fixture()
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self.write_profiles(root, report)
            receipt = parse_profiles(root, report, job)
            self.assertEqual(receipt['process_instructions'], 150)
            self.assertEqual(receipt['phases']['setup'], 30)
            self.assertEqual(receipt['online_instructions'], 50)
            self.assertIsNone(receipt['phases']['rho_solve'])

    def test_profiler_rejects_missing_duplicate_mixed_and_inconsistent_intervals(self):
        for mutation in ('missing', 'duplicate', 'legacy', 'checksum', 'coverage'):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as d:
                root = Path(d)
                job, report = fixture()
                self.write_profiles(root, report)
                p = root / 'callgrind.out.1'
                if mutation == 'missing':
                    p.unlink()
                elif mutation == 'duplicate':
                    (root / 'callgrind.out.duplicate').write_text(p.read_text())
                elif mutation == 'legacy':
                    p.write_text(p.read_text().replace('generic_ic_setup', 'ic_setup'))
                elif mutation == 'checksum':
                    (root / 'stderr.txt').write_text('==9== Collected : 151\n')
                else:
                    p.write_text(p.read_text().replace('generic_ic_setup', 'generic_ic_rho_solve'))
                with self.assertRaises(InvalidEvidence):
                    parse_profiles(root, report, job)


if __name__ == '__main__':
    unittest.main()
