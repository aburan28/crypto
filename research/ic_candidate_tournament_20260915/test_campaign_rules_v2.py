"""Reference substitutions and cross-metric regressions cannot earn promotion."""
import copy
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import campaign_rules as legacy
import campaign_rules_v2 as rules
from oracle import InvalidEvidence
from test_campaign_rules import contract as old_contract, paired_rows
from tournament import bounded_protocol, gate, is_reference
import tournament

HERE = Path(__file__).parent
EVIDENCE = HERE / 'goal_20260924/generic-reference-qualification/RESULTS.json'


def contract():
    result = old_contract()
    result.update(purpose=rules.PURPOSE, attempt_number=2, seed=2026092552,
        scientific_admission=True, schema_version=2, unit='valgrind-3.22-amd64-Ir',
        stages=['aa', 'smoke', 'development', 'selection', 'confirmation', 'replay'],
        evidence_scope='synthetic statistics test; no measured worker',
        target_exposure_schema=1, reference_qualification=rules.expected_binding(),
        metric_references=copy.deepcopy(rules.METRIC_REFERENCES), execution_number_stride=2,
        candidate_count=11, scheduled_slot_bound=3480, bootstrap_draws=2000)
    return result


def rows(online_reference=.8):
    result = paired_rows(.7)
    for row in list(result):
        if row['arm'] == 'incumbent':
            item = copy.deepcopy(row)
            item['arm'] = 'ic_online'
            item['certificate']['factor_base_sha256'] = 'ic_online'
            item['measurement']['native_timing']['online']['wall_ns'] = int(1000000 * online_reference)
            result.append(item)
    return result


def references():
    binding = rules.expected_binding()['bindings']
    result = []
    for name, value in binding.items():
        arm = dict(id=name, source_manifest_sha256=value['source_manifest_sha256'],
                   config=copy.deepcopy(value['configuration']))
        if value['adapter']:
            arm.update(adapter=value['adapter'], build_sha256=value['build_sha256'])
        if name != 'incumbent':
            arm['kind'] = 'ic-reference' if name == 'ic_online' else 'rho-reference'
        result.append(arm)
    return result


class VersionedCampaignTests(unittest.TestCase):
    def test_binds_accepted_evidence_and_every_separate_reference(self):
        data = json.loads(EVIDENCE.read_text())
        arms = references()
        self.assertEqual(rules.qualified_binding(data['qualification'], data['observer'],
                         arms[0], arms[1:]), rules.expected_binding())
        self.assertEqual(rules.RULE, legacy.RULE)
        self.assertIs(bounded_protocol(contract()), rules)
        self.assertIs(bounded_protocol(old_contract()), legacy)
        for index, key, value in ((1, 'source_manifest_sha256', rules.IC_SOURCE),
                                  (2, 'config', dict(rules.CONFIG, rho_parallel_walks=4)),
                                  (3, 'build_sha256', 'f' * 64),
                                  (1, 'kind', 'rho-reference'), (0, 'kind', 'rho-reference')):
            changed = copy.deepcopy(arms)
            changed[index][key] = value
            with self.assertRaises(InvalidEvidence):
                rules.qualified_binding(data['qualification'], data['observer'], changed[0], changed[1:])
        changed = copy.deepcopy(data['observer'])
        changed['all_observer_costs_charged'] = False
        with self.assertRaises(InvalidEvidence):
            rules.qualified_binding(data['qualification'], changed, arms[0], arms[1:])
        self.assertTrue(is_reference(arms[1]))
        self.assertFalse(is_reference(arms[0]))

    def test_full_schedule_fits_without_reducing_confirmation_or_portfolio(self):
        self.assertEqual(rules.scheduled_slot_bound(11), 3480)
        with self.assertRaises(InvalidEvidence):
            rules.scheduled_slot_bound(12)
        for key, value in (('attempt_number', 1), ('candidate_count', 12),
                           ('scheduled_slot_bound', 3300), ('confirmation_cases', 60),
                           ('metric_references', dict(rules.METRIC_REFERENCES, online_ns='incumbent')),
                           ('familywise_rule', dict(rules.RULE, comparisons=12)),
                           ('scientific_admission', False), ('unit', 'curve-additions')):
            changed = contract()
            changed[key] = value
            with self.assertRaises(InvalidEvidence, msg=key):
                rules.validate_contract(changed)

    def test_invalid_references_fail_before_host_preflight_or_any_worker(self):
        args = SimpleNamespace(out=Path('/unused'), exposed_fixtures=[], attempt_number=2,
            seed=2026092552, qualification=False, profile='pilot', confirmation_cases='',
            cells='17a1,19a0,23a0,23a1,31a0', holdout_cells='29a1', selection_width=6,
            exploration_slots=1, timeout=180, max_processes=3500, targets=1,
            require_native_progress=True, objective='incumbent', comparison_kind='factor-base-policy',
            candidates=Path('/candidate-registry'), campaign_version=2,
            reference_registry=Path('/reference-registry'), qualified_observer=Path('/observer'),
            qualified_report=Path('/qualification'), target_history=Path('/history'),
            rho_source_root=None, rho_config=None)
        arms = [dict(id='incumbent', config=rules.CONFIG), dict(id='challenger', config=rules.CONFIG)]
        refs = references()[1:]
        refs[-1]['config']['rho_parallel_walks'] = 32
        with patch.object(tournament, 'read', side_effect=[arms, refs]), \
             patch.object(tournament.platform, 'machine') as host, \
             patch.object(tournament, 'snapshot_build') as build:
            with self.assertRaisesRegex(InvalidEvidence, 'declared reference'):
                tournament.prepare(args)
            host.assert_not_called()
            build.assert_not_called()

    def test_actual_stage_schedule_retains_references_within_the_pair_cap(self):
        c = contract()
        c['reference_arms'] = references()[1:]
        arms = [references()[0]] + [dict(id=f'challenger{i}', config=rules.CONFIG) for i in range(10)]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / 'summaries').mkdir()
            tournament.write(root / 'contract.json', c)
            tournament.write(root / 'summaries/smoke.json', dict(failures=[]))
            tournament.write(root / 'summaries/development.json', dict(
                retained_portfolio=[dict(candidate=a['id']) for a in arms[1:7]]))
            tournament.write(root / 'summaries/selection.json', dict(provisional_challenger=arms[1]['id']))
            count = 0
            for stage, targets, expected_arms in (
                    ('aa', 5, 2), ('smoke', 5, 14), ('development', 15, 14),
                    ('selection', 15, 10), ('confirmation', 72, 5), ('replay', 72, 5)):
                active = tournament.stage_arms(root, stage, arms)
                self.assertEqual(len(active), expected_arms, stage)
                if stage != 'aa':
                    self.assertEqual([a['id'] for a in active if is_reference(a)], list(rules.REFERENCE_ROLES))
                count += targets * len(active) * c['repetitions']
            self.assertEqual(count, rules.scheduled_slot_bound(len(arms)))

    def test_final_inference_uses_stronger_online_reference_and_original_family_budget(self):
        c = contract()
        result = rules.final_comparison(rows(.8), 'challenger', c)
        self.assertTrue(gate(result, c))
        self.assertEqual(result['familywise']['paired_targets'], 72)
        self.assertEqual(result['familywise']['monte_carlo_tail_index'], 137)
        self.assertAlmostEqual(result['familywise']['metrics']['instructions']['ratio'], .7)
        self.assertAlmostEqual(result['familywise']['metrics']['cold_ns']['ratio'], .7)
        self.assertAlmostEqual(result['familywise']['metrics']['online_ns']['ratio'], .875)
        self.assertEqual(result['familywise']['metric_references'], rules.METRIC_REFERENCES)
        regression = rules.final_comparison(rows(.6), 'challenger', c)
        self.assertAlmostEqual(regression['online']['candidate_over_baseline'], 7/6)
        self.assertFalse(gate(regression, c))
        tampered = copy.deepcopy(result)
        tampered['familywise']['metric_references']['online_ns'] = 'incumbent'
        self.assertFalse(gate(tampered, c))

    def test_missing_or_failed_online_reference_is_never_a_successful_subset(self):
        c = contract()
        measured = rows()
        removed = [r for r in measured if not (r['arm'] == 'ic_online' and r['case'] == 'n17a1-0')]
        self.assertFalse(rules.final_comparison(removed, 'challenger', c)['eligible'])
        measured[-1]['status'] = 'TIMEOUT'
        self.assertFalse(rules.final_comparison(measured, 'challenger', c)['eligible'])

    def test_reference_is_reported_but_cannot_enter_the_candidate_portfolio(self):
        c = contract()
        c['bootstrap_draws'] = 20
        arms = references()[:2] + [dict(id='challenger', config=rules.CONFIG)]
        with patch.object(tournament, 'load_stage', return_value=rows()), \
             patch.object(tournament, 'online_table', return_value=[]):
            # The stored raw process clocks are irrelevant to this selection control.
            with patch.object(tournament, 'process_seconds', return_value=0):
                summary = tournament.summarize(Path('/unused'), c, 'development', [], arms, save=False)
        self.assertEqual([r['candidate'] for r in summary['comparisons']], ['challenger'])
        self.assertEqual([r['candidate'] for r in summary['retained_portfolio']], ['challenger'])
        self.assertTrue(summary['ic_reference_comparisons']['ic_online']['eligible'])

    def test_failed_reference_blocks_the_final_decision_even_if_candidate_gates_pass(self):
        c = contract()
        summary = dict(comparisons=[dict(candidate='challenger', eligible=True)],
            rho_comparisons={name: dict(eligible=True) for name in ('rho', 'rho_online')},
            ic_reference_comparisons=dict(ic_online=dict(eligible=False)), rho_over_incumbent=None)
        def read(path):
            return dict(provisional_challenger='challenger') if path.name == 'selection.json' else summary
        with patch.object(tournament, 'read', side_effect=read), patch.object(tournament, 'gate', return_value=True):
            result = tournament.decision(Path('/unused'), c,
                dict(confirmation=[], replay=[]), [], save=False)
        self.assertEqual(result['status'], 'inconclusive')
        self.assertIsNone(result['winner'])
        self.assertFalse(result['promotion_eligible'])

    def test_recomputes_real_development_comparison_against_online_leader(self):
        data = json.loads(EVIDENCE.read_text())
        aliases = dict(incumbent='incumbent', prepared_both='ic_online', generic_dense='challenger')
        measured = []
        for run in data['development_runs']:
            if run['alias'] not in aliases:
                continue
            m = run['measurement']
            measured.append(dict(arm=aliases[run['alias']], case=run['case'], cell=run['cell'],
                repetition=run['repetition'], case_sha256=run['workload_id'], status=run['status'],
                total_operations=m['total_operations'], mode='ic', measurement=m,
                certificate=m['certificate'], native_process=dict(
                    process_wall_seconds=m['native_timing']['cold']['wall_ns']/1e9)))
        result = rules.comparison(measured, 'challenger', contract())
        table = {row['alias']: row['comparison_to_archived_incumbent'] for row in data['qualification']['table']}
        self.assertTrue(result['eligible'])
        self.assertEqual(result['paired_cases'], 15)
        self.assertAlmostEqual(result['candidate_over_baseline'], table['generic_dense']['candidate_over_baseline'], places=11)
        self.assertAlmostEqual(result['online']['candidate_over_baseline'],
            table['generic_dense']['online']['candidate_over_baseline'] /
            table['prepared_both']['online']['candidate_over_baseline'], places=12)

    def test_registered_round_two_panel_uses_new_source_and_pair_table_modes(self):
        from producer.evidence import executed_policy
        from run_improvement_v2 import registry, PANEL, EXPOSED
        panel = json.loads(PANEL.read_text())
        rows = registry(panel, Path('/candidate/source'))
        self.assertEqual(panel['round'], 2)
        self.assertEqual(panel['candidate_panel'], 'round2-v1')
        self.assertEqual(len(rows), 11)
        self.assertNotEqual(panel['candidate_source_sha256'],
                            json.loads((HERE/'goal_20260924/improvement/round1.json').read_text())
                            ['candidate_source_sha256'])
        self.assertEqual(len({json.dumps([row.get('source_root'), row['config']], sort_keys=True)
                              for row in rows}), 11)
        policies = {row['id']: executed_policy(row['config'], panel['candidate_panel'])
                    for row in rows[1:]}
        self.assertEqual(policies['compatibility']['pair_table'], 'heuristic')
        self.assertEqual(policies['half_table']['pair_table'], 'half')
        self.assertEqual(policies['cover_table']['pair_table'], 'cover')
        self.assertEqual(policies['stop6_word_half'],
                         dict(orbit_batch=1, orbit_target=6, row_kernel='word', pair_table='half'))
        with self.assertRaises(InvalidEvidence):
            executed_policy(dict(rows[1]['config'], full_pair_table=False), panel['candidate_panel'])
        with self.assertRaises(InvalidEvidence):
            executed_policy(dict(rows[1]['config'], pair_table='quarter'), panel['candidate_panel'])
        for path in EXPOSED:
            self.assertTrue(path.is_file(), path.name)
        self.assertEqual(2026092550 + 2, 2026092552)

    def test_final_registered_panel_uses_distinct_mechanism_combinations(self):
        from producer.evidence import executed_policy
        from run_improvement_v3 import registry, PANEL, PANEL_SHA256, EXPOSED
        from tournament import digest
        panel = json.loads(PANEL.read_text())
        rows = registry(panel, Path('/candidate/source'))
        self.assertEqual(digest(PANEL), PANEL_SHA256)
        self.assertEqual(panel['round'], 3)
        self.assertEqual(panel['candidate_panel'], 'round3-v1')
        self.assertEqual(len(rows), 11)
        self.assertEqual(rules.scheduled_slot_bound(len(rows)), 3480)
        self.assertEqual(len({json.dumps(row['config'], sort_keys=True) for row in rows}), 11)
        for prior in (HERE/'goal_20260924/improvement/round1.json',
                      HERE/'goal_20260924/improvement-v2/round2.json'):
            earlier = json.loads(prior.read_text())
            exact = {json.dumps(a['config'], sort_keys=True) for a in earlier['candidates'][1:]}
            self.assertFalse(exact & {json.dumps(a['config'], sort_keys=True) for a in rows[1:]})
        policies = {row['id']: executed_policy(row['config'], panel['candidate_panel'])
                    for row in rows[1:]}
        self.assertEqual(policies['stop4_word_half'],
                         dict(orbit_batch=1, orbit_target=4, row_kernel='word', pair_table='half'))
        self.assertEqual(policies['stop6_bounded_half']['row_kernel'], 'bounded')
        self.assertEqual(policies['word_cover']['pair_table'], 'cover')
        self.assertEqual(rows[-1]['config']['batch_trials'], 2)
        for path in EXPOSED:
            self.assertTrue(path.is_file(), path.name)


if __name__ == '__main__':
    unittest.main()
