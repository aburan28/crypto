"""Adversarial checks for fresh targets, accepted references and final inference."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace

import campaign_rules as rules
from oracle import InvalidEvidence
from portfolio import retain
from target_history import (extend, file_hash, key_for, reserve_fresh,
                            validate_fresh, verify_sources)
from test_certificate import FIXTURE


def contract():
    config = dict(rules.CONFIG)
    bindings = {}
    for name, selected, source in (
            ('incumbent', 'pairinv', rules.IC_SOURCE),
            ('rho', 'rho_incumbent_4', rules.COLD_RHO_SOURCE),
            ('rho_online', 'rho_pairinv_4', rules.IC_SOURCE)):
        bindings[name] = dict(selected_alias=selected, source_manifest_sha256=source,
                             configuration=config if name == 'incumbent' else dict(config, rho_parallel_walks=4))
    return dict(purpose=rules.PURPOSE, familywise_rule=copy.deepcopy(rules.RULE),
        attempt_number=1, seed=2026092551, confirmation_cases=72, repetitions=3,
        cells=rules.CELLS, holdout_cells=rules.HOLDOUT_CELLS,
        confirmation_cases_per_cell={cell:12 for cell in rules.CELLS + rules.HOLDOUT_CELLS},
        confirmation_ratio=.8, max_cell_ratio=1.1, require_native_progress=True,
        objective='incumbent', selection_width=6, exploration_slots=1, target_count=1,
        limits=dict(max_profiled_jobs=3500, memory_bytes=8*1024**3, timeout_seconds=180, worker_threads=1),
        comparison_kind='factor-base-policy', reference_qualification=dict(
            qualification_sha256=rules.QUALIFICATION_SHA256,
            archive_sha256=rules.REFERENCE_ARCHIVE_SHA256, bindings=bindings))


def paired_rows(ratio):
    rows = []
    for cell in rules.CELLS + rules.HOLDOUT_CELLS:
        for target in range(12):
            for rep in range(3):
                for alias, scale in (('incumbent', 1), ('challenger', ratio)):
                    cost = int(1000000 * scale)
                    rows.append(dict(arm=alias, case=f'{cell}-{target}', cell=cell,
                        repetition=rep, case_sha256=f'{cell}-{target}', status='VERIFIED',
                        total_operations=cost, mode='ic',
                        certificate={'factor_base_sha256':alias},
                        native_process={'process_wall_seconds':cost/1e9},
                        measurement={'native_timing':{'cold':{'wall_ns':cost}, 'online':{'wall_ns':cost}}}))
    return rows


class CampaignRulesTests(unittest.TestCase):
    def test_preflight_rejects_budget_and_panel_changes_before_any_fresh_target(self):
        args = SimpleNamespace(attempt_number=1, seed=2026092551, qualification=False,
            profile='pilot', confirmation_cases='', cells='17a1,19a0,23a0,23a1,31a0',
            holdout_cells='29a1', selection_width=6, exploration_slots=1, timeout=180,
            max_processes=3500, targets=1, require_native_progress=True,
            objective='incumbent', comparison_kind='factor-base-policy')
        rules.validate_preparation(args)
        for key, value in (('profile','standard'), ('max_processes',3501), ('targets',2),
                           ('attempt_number',4), ('cells','17a1,19a0,23a0,23a1')):
            changed = copy.copy(args); setattr(changed, key, value)
            with self.assertRaises(InvalidEvidence): rules.validate_preparation(changed)

    def test_registered_round_has_real_combinations_and_distinct_executed_policies(self):
        from producer.evidence import executed_policy
        from run_improvement import registry, PANEL
        panel = json.loads(PANEL.read_text())
        rows = registry(panel, Path('/candidate/source'))
        self.assertEqual(len(rows), 16)
        self.assertEqual(len({json.dumps([row.get('source_root'),row['config']], sort_keys=True)
                              for row in rows}), 16)
        policies = {row['id']:executed_policy(row['config'], panel['candidate_panel']) for row in rows[1:]}
        self.assertEqual(policies['stop3_word_full'], dict(policies['stop3'], row_kernel='word', full_pair_table=True))
        with self.assertRaises(InvalidEvidence): executed_policy(rows[1]['config'], None)
        for key, value in (('orbit_batch',True), ('orbit_target',0), ('row_kernel','ignored'),
                           ('full_pair_table',1)):
            with self.assertRaises(InvalidEvidence):
                executed_policy(dict(rows[1]['config'], **{key:value}), panel['candidate_panel'])

    def test_accepted_qualification_binds_exact_sources_and_both_rho_settings(self):
        path = Path(__file__).parent/'goal_20260924/reference-qualification/RESULTS.json'
        report = json.loads(path.read_text())['qualification']
        c = contract()
        arms = [dict(id=name, source_manifest_sha256=value['source_manifest_sha256'],
                     config=value['configuration']) for name,value in c['reference_qualification']['bindings'].items()]
        self.assertEqual(rules.qualified_binding(report, arms[0], arms[1:]), c['reference_qualification'])
        arms[1]['config'] = dict(arms[1]['config'], rho_parallel_walks=1)
        with self.assertRaises(InvalidEvidence):
            rules.qualified_binding(report, arms[0], arms[1:])
        with self.assertRaises(InvalidEvidence):
            rules.qualified_binding(report, arms[0], [])

    def test_round_budget_and_reference_flags_cannot_relax_confirmation(self):
        for key, value in (('attempt_number',4), ('confirmation_cases',59), ('repetitions',1),
                           ('familywise_rule',dict(rules.RULE, comparisons=1)),
                           ('reference_qualification',{'complete':True})):
            c = contract(); c[key] = value
            with self.assertRaises(InvalidEvidence, msg=key): rules.validate_contract(c)
        c = contract(); c['reference_qualification']['bindings']['rho']['source_manifest_sha256'] = 'f'*64
        with self.assertRaises(InvalidEvidence): rules.validate_contract(c)

    def test_known_uniform_gain_has_correct_ratio_and_target_sample_count(self):
        c = contract(); result = rules.final_comparison(paired_rows(.7), 'challenger', c)
        self.assertEqual(result['familywise']['paired_targets'], 72)
        self.assertEqual(result['familywise']['monte_carlo_tail_index'], 137)
        for metric in rules.METRICS:
            self.assertAlmostEqual(result['familywise']['metrics'][metric]['ratio'], .7)
            self.assertAlmostEqual(result['familywise']['metrics'][metric]['upper'], .7)
        self.assertTrue(rules.promotion_passes(result, c))
        # A fast aggregate cannot hide a regressing cell, even in one metric.
        bad = copy.deepcopy(result)
        bad['familywise']['metrics']['online_ns']['per_cell'][rules.CELLS[0]] = 1.10001
        self.assertFalse(rules.promotion_passes(bad, c))
        bad = copy.deepcopy(result); bad['familywise']['metrics']['cold_ns']['upper'] = 1
        self.assertFalse(rules.promotion_passes(bad, c))

    def test_missing_target_or_timeout_never_becomes_an_inference_win(self):
        c = contract(); rows = paired_rows(.5)
        self.assertFalse(rules.final_comparison(rows[:-1], 'challenger', c)['eligible'])
        rows[-1]['status'] = 'TIMEOUT'
        self.assertFalse(rules.final_comparison(rows, 'challenger', c)['eligible'])

    def test_all_repetitions_of_a_point_cannot_be_dropped_as_if_workload_complete(self):
        c = contract(); rows = paired_rows(.5)
        rows = [row for row in rows if row['case'] != rules.CELLS[0]+'-0']
        with self.assertRaisesRegex(InvalidEvidence, 'allocation'):
            rules.final_comparison(rows, 'challenger', c)

    def test_online_leader_survives_cold_ranking_with_exploration(self):
        arms = [dict(id=str(i), config={}) for i in range(8)]
        comparisons = []
        for i in range(8):
            cost = .5 + i/10
            online = .1 if i == 7 else 1
            comparisons.append(dict(candidate=str(i), eligible=True, candidate_over_baseline=cost,
                native_wall_candidate_over_baseline=cost, per_cell={'c':cost},
                native_wall_per_cell={'c':cost}, online=dict(candidate_over_baseline=online, per_cell={'c':online})))
        selected = retain(comparisons, arms, width=6, exploration=1, seed=17)
        self.assertEqual(selected[0]['candidate'], '7')
        self.assertIn('0', [item['candidate'] for item in selected])
        self.assertEqual(len(selected), 6)
        self.assertEqual(selected[-1]['reason'], 'predeclared exploration slot')


class TargetHistoryTests(unittest.TestCase):
    def test_new_seed_cannot_reuse_an_old_public_point(self):
        fixture = copy.deepcopy(FIXTURE)
        used = {}
        self.assertTrue(reserve_fresh(fixture, used))
        fixture['target_seeds'] = [999999] * len(fixture['targets'])
        self.assertFalse(reserve_fresh(fixture, used))

    def test_cross_stage_reuse_is_rejected_and_replay_is_intentional(self):
        case = {'fixture':copy.deepcopy(FIXTURE)}
        empty = dict(schema_version=1, curves=[])
        fixtures = dict(confirmation=[case], replay=[copy.deepcopy(case)])
        self.assertEqual(validate_fresh(fixtures, empty), len(FIXTURE['targets']))
        fixtures['development'] = [copy.deepcopy(case)]
        with self.assertRaises(InvalidEvidence): validate_fresh(fixtures, empty)

    def test_interrupted_preparation_points_are_retained(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root/'fixture_generation/confirmation/case/attempt-0/stdout.json'
            path.parent.mkdir(parents=True); path.write_text(json.dumps({'fixture':FIXTURE}))
            extended = extend(dict(schema_version=1, curves=[], sources=[]), [root])
            self.assertEqual(extended['curves'][0]['curve_id'], key_for(FIXTURE))
            self.assertEqual(extended['curves'][0]['points'], sorted([list(map(int,p)) for p in FIXTURE['targets']]))
            self.assertIsNone(extended['round_additions'][0]['contract_sha256'])

    def test_source_reconstruction_detects_omitted_points(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); path = root/'fixture.json'; path.write_text(json.dumps(FIXTURE))
            history = dict(schema_version=1, sources=[dict(path='fixture.json',sha256=file_hash(path))],
                curves=[dict(curve_id=key_for(FIXTURE),points=[list(map(int,p)) for p in FIXTURE['targets']])])
            self.assertEqual(verify_sources(history, root)['status'], 'VERIFIED')
            history['curves'][0]['points'].pop()
            with self.assertRaises(InvalidEvidence): verify_sources(history, root)


if __name__ == '__main__':
    unittest.main()
