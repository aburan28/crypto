"""Qualification budgets, baseline binding and observer diagnostic boundaries."""
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import generic_reference_qualification as qualification
from generic_observer import checked_observation, effects, entries, number_base
from identity import sha256
from oracle import InvalidEvidence
from test_generic_admission import build_fixture, row


class GenericQualificationTests(unittest.TestCase):
    def test_observer_execution_numbers_cannot_overlap_parent_range(self):
        parent = dict(run_aliases=list(range(23)), repetitions=3,
                      execution_number_stride=2, run_number_base=17 << 16)
        self.assertEqual(number_base(parent), (17 << 16) + (1 << 15))
        parent['run_aliases'] = list(range(1000))
        with self.assertRaises(InvalidEvidence):
            number_base(parent)

    def test_panel_preserves_qualified_sources_and_configuration_independence(self):
        result = qualification.registry({'both': Path('/prepared/both')}, Path('/generic'))
        self.assertEqual([a['id'] for a in result],
                         ['incumbent', 'prepared_both', 'generic_dense', 'generic_sparse'])
        self.assertEqual(result[1]['source_root'], '/prepared/both/source')
        self.assertEqual([a['config']['linear_algebra'] for a in result], ['tiny_gauss', 'tiny_gauss', 'dense', 'sparse'])
        result[0]['config']['max_trials'] = 1
        self.assertEqual(result[1]['config']['max_trials'], 65536)
        self.assertEqual(qualification.CONFIG['max_trials'], 65536)

    def test_bound_prepared_artifact_rejects_changed_source_bytes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = {'Cargo.toml': hashlib.sha256(b'original').hexdigest()}
            expected = sha256(manifest)
            for name in ('both', 'pairinv'):
                directory = root / f'ic-producer-{name}-123-1/ic-producer'
                (directory / 'source').mkdir(parents=True)
                (directory / 'source/Cargo.toml').write_bytes(b'original')
                (directory / 'source-manifest.json').write_text(json.dumps(manifest))
                (directory / 'preparation.json').write_text(json.dumps(dict(reference=name,
                    instrumented=True, source_manifest_sha256=expected)))
            with patch.object(qualification, 'IC_SOURCE', expected), patch.object(qualification, 'COLD_RHO_SOURCE', expected):
                sources = qualification.prepared_sources(root)
                (sources['both'] / 'source/Cargo.toml').write_bytes(b'changed')
                with self.assertRaises(InvalidEvidence):
                    qualification.prepared_sources(root)

    def test_schedule_is_bounded_paired_and_counterbalanced(self):
        fixtures = dict(development=[dict(id=f'n{cell}-{i}', cell=f'n{cell}')
                                     for cell in range(5) for i in range(3)])
        arms = [dict(id=f'arm{i}') for i in range(8)]
        schedule = entries(fixtures, arms)
        self.assertEqual(len(schedule), 360)
        self.assertEqual(schedule, entries(fixtures, arms))
        self.assertEqual(len({r['name'] for r in schedule}), 360)
        for case in fixtures['development']:
            for arm in arms:
                group = sorted((r for r in schedule if r['case'] == case['id'] and r['arm'] == arm['id']),
                               key=lambda r:r['repetition'])
                self.assertEqual([r['repetition'] for r in group], [0, 1, 2])
                self.assertEqual(group[0]['order'], group[1]['order'][::-1])
                self.assertEqual(group[0]['order'], group[2]['order'])

    def test_uncertainty_uses_points_not_processes_and_preserves_failures(self):
        rows = []
        for cell in ('a', 'b'):
            for case in range(3):
                for repetition, value in enumerate((80, 100, 120)):
                    rows.append(dict(cell=cell, case=f'{cell}-{case}', repetition=repetition,
                        status='MATCH_COMPLETE', observations=dict(
                            enabled=dict(outer_online_wall_ns=2*value, process_wall_ns=3*value),
                            legacy=dict(outer_online_wall_ns=value, process_wall_ns=value))))
        result = effects(rows)
        self.assertEqual(result['target_count'], 6)
        self.assertEqual(result['cell_count'], 2)
        self.assertAlmostEqual(result['metrics']['outer_online']['enabled_over_legacy'], 2)
        for value in result['metrics']['process_wall']['descriptive_percentile_95']:
            self.assertAlmostEqual(value, 3)
        with self.assertRaises(InvalidEvidence):
            effects(rows[:-1])
        rows[0]['status'] = 'FAILED'
        self.assertIsNone(effects(rows))

    def test_legacy_format_cannot_acquire_exclusive_admission_or_invalid_clock(self):
        # Retained real group evidence, synthetic executable tag, and a derived
        # legacy-format record exercise checking only. CI runs actual paired modes.
        item = row()
        report, job = copy.deepcopy(item['report']), copy.deepcopy(item['job'])
        build, source = build_fixture()
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            worker = directory / 'worker'
            worker.write_bytes(b'synthetic tag, never executed')
            build['worker_sha256'] = hashlib.sha256(worker.read_bytes()).hexdigest()
            process = dict(process_status='EXITED', exit_code=0, process_wall_ns=item['process_wall_ns'])
            (directory / 'process.json').write_text(json.dumps(process))

            def check(mode):
                (directory / 'stdout.json').write_text(json.dumps(report))
                return checked_observation(directory, mode=mode, job=job,
                    fixture=report['fixture'], build=build, source=source, worker=worker)

            enabled = check('enabled')
            self.assertIsNotNone(enabled['stage_audit'])
            self.assertIsNotNone(enabled['phase_audit'])
            job['exclusive_phases'] = False
            with self.assertRaises(InvalidEvidence):
                check('legacy')
            report['online_wall_ns'] = report['outer_online_wall_ns']
            report.pop('generic_phase_timing')
            report['generic_admission_schema'] = None
            legacy = check('legacy')
            self.assertEqual(sha256(legacy['signature']), sha256(enabled['signature']))
            self.assertIsNone(legacy['stage_audit'])
            self.assertIsNone(legacy['phase_audit'])
            self.assertIsNone(legacy['enabled_native_timing'])
            report['online_wall_ns'] = True
            with self.assertRaises(InvalidEvidence):
                check('legacy')


if __name__ == '__main__':
    unittest.main()
