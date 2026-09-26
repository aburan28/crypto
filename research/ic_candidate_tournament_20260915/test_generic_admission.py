"""Adversarial scientific controls against retained real worker executions."""
import copy
from functools import lru_cache
import gzip
import json
from pathlib import Path
import tempfile
import unittest

from generic_admission import admit, method_record, scientific_ledger
from execution_ids import BLOCK_SIZE, allocation, audit_runs, ci_start
from generic_bases import complete_factors, poly_product, verify_base
from generic_build import digest, verify_binding, verify_build_record
from generic_phases import verify_native
from generic_queries import verify_queries
from generic_stages import effective_config, verify_stages
from identity import candidate_manifest, sha256
from oracle import InvalidEvidence, require, verify

ROOT = Path(__file__).parent / 'goal_20260924/generic-scientific-admission/final'


@lru_cache(maxsize=1)
def retained():
    return [json.loads(line) for line in gzip.decompress((ROOT / 'worker-raw.jsonl.gz').read_bytes()).splitlines()]


def row(name='a0-pair_table-dense'):
    return copy.deepcopy(next(item for item in retained() if item['name'] == name))


def check(item):
    return verify_stages(item['report'], item['report']['fixture'], item['job'])


def build_fixture():
    return (json.loads((ROOT / 'build-record.json').read_text()),
            json.loads(gzip.decompress((ROOT / 'source-manifest.json.gz').read_bytes())))


class GenericAdmissionTests(unittest.TestCase):
    def test_run_number_is_required_and_rejects_boolean_or_negative_values(self):
        kwargs = dict(executable=None, process_wall_ns=None, resources=None)
        with self.assertRaises(TypeError):
            admit(None, None, None, None, None, **kwargs)
        for value in (True, -1, 1.0, '1'):
            with self.assertRaisesRegex(InvalidEvidence, 'run number'):
                admit(None, None, None, None, None, number=value, **kwargs)

    def test_execution_allocation_separates_lanes_attempts_and_workflow_runs(self):
        occupied = set()
        for workflow, attempt in ((1, 1), (1, 2), (2, 1), (36215166278, 1),
                                  (36215166278, 2), (36215166279, 1), (1, 1000001)):
            for lane in ('controls', 'integration'):
                start = ci_start(workflow, attempt, lane)
                plan = allocation(start, [f'job-{i}' for i in range(BLOCK_SIZE)])
                numbers = {item['number'] for item in plan['executions']}
                self.assertFalse(occupied & numbers)
                occupied |= numbers
        for labels in ([], ['duplicate', 'duplicate'], ['']):
            with self.assertRaises(InvalidEvidence):
                allocation(0, labels)
        with self.assertRaises(InvalidEvidence):
            allocation(0, [str(i) for i in range(BLOCK_SIZE + 1)])
        with self.assertRaises(InvalidEvidence):
            allocation(True, ['job'])

    def test_corrected_runs_preserve_measurements_and_expose_original_collisions(self):
        original = list(retained())
        for panel in ('failed-la', 'sparse-core'):
            original += [json.loads(line) for line in gzip.decompress(
                (ROOT.parent / panel / 'worker-raw.jsonl.gz').read_bytes()).splitlines()]
        old = {item['name']: item['receipt']['admission'] for item in original if item['job']['mode'] == 'ic'}
        with self.assertRaisesRegex(InvalidEvidence, 'duplicate canonical run key'):
            audit_runs([item['run'] for item in old.values()])
        export = json.loads(gzip.decompress((ROOT.parent / 'run-records-v2.json.gz').read_bytes()))
        records = export['records']
        self.assertEqual(audit_runs([item['run'] for item in records])['unique_keys'], 53)
        by_name = {item['name']: item for item in records}
        self.assertEqual(set(by_name), set(old))
        for name, item in by_name.items():
            self.assertEqual(item['candidate'], old[name]['candidate'])
            self.assertEqual(item['workload'], old[name]['workload'])
            restored = dict(item['run'], run_id=item['original_run_id'])
            self.assertEqual(restored, old[name]['run'])
        for n in (9, 13):
            for la in ('dense', 'sparse'):
                left, right = (by_name[f'n{n}-window{w}-{la}'] for w in (0, 1000))
                self.assertEqual(left['candidate'], right['candidate'])
                self.assertEqual(left['workload'], right['workload'])
                self.assertNotEqual(left['run']['run_id'], right['run']['run_id'])
        with self.assertRaisesRegex(InvalidEvidence, 'duplicate canonical run key'):
            audit_runs([records[0]['run'], records[0]['run']])

    def test_boolean_integer_aliases_cannot_change_dispatch_counters_or_matrix(self):
        def column(r):
            entry = next(e for row in r['relation_matrix']['rows'] for e in row['entries'] if e[0] in (0, 1))
            entry[0] = bool(entry[0])
        def sparse(r):
            r['matrix_batches'][-1]['sparse_report']['attempts'] = False
            r['log_table_report']['sparse_report']['attempts'] = False
        for name, mutate in [('a0-pair_table-dense', lambda r: r['collector_dispatch'].update(pair_table=1)),
                             ('a0-pair_table-dense', lambda r: r['descent_dispatch'].update(direct_collision=0)),
                             ('a0-pair_table-dense', column), ('a0-pair_table-sparse', sparse),
                             ('a0-pair_table-dense', lambda r: r.update(generic_admission_schema=True))]:
            item = row(name)
            mutate(item['report'])
            verify(item['report'], item['report']['fixture'], summands=2)
            with self.assertRaises(InvalidEvidence):
                check(item)

    def test_actual_block_wiedemann_controls_and_parameter_substitution(self):
        rows = [json.loads(line) for line in gzip.decompress(
            (ROOT.parent / 'sparse-core/worker-raw.jsonl.gz').read_bytes()).splitlines()]
        self.assertEqual(len(rows), 2)
        for item in rows:
            receipt = check(item)
            self.assertTrue(receipt['matrix']['certified_logs'])
            core = item['report']['matrix_batches'][-1]['sparse_report']['wiedemann']
            self.assertGreater(core['products'], 0)
            core['block_m'] += 1
            with self.assertRaisesRegex(InvalidEvidence, 'Wiedemann dispatch'):
                check(item)

    def test_failed_la_calls_are_preserved_before_failure_and_success(self):
        rows = [json.loads(line) for line in gzip.decompress(
            (ROOT.parent / 'failed-la/worker-raw.jsonl.gz').read_bytes()).splitlines()]
        self.assertEqual(len(rows), 4)
        for item in rows:
            check(item)
            failed = [b for b in item['report']['matrix_batches'] if b['solve_attempts'] and not b['verified']]
            self.assertTrue(failed)
            failed[0]['solve_attempts'] = 0
            with self.assertRaisesRegex(InvalidEvidence, 'batch accounting'):
                check(item)

    def test_runtime_overrides_were_rejected_before_measurement(self):
        summary = json.loads((ROOT / 'summary.json').read_text())
        self.assertEqual(len(summary['environment_overrides_rejected']), 11)
        self.assertIn('IC_REDIS_URL', summary['environment_overrides_rejected'])

    def test_every_frozen_stage_case_and_failure_replays(self):
        rows = retained()
        self.assertEqual(len(rows), 63)
        self.assertEqual(sum(item['report']['status'] == 'inventory' for item in rows), 14)
        self.assertEqual(sum(item['report']['status'] == 'incomplete' for item in rows), 7)
        for item in rows:
            with self.subTest(case=item['name']):
                self.assertEqual(item['receipt']['status'], 'PASS')
                if item['job']['mode'] != 'rho':
                    receipt = check(item)
                    self.assertFalse(receipt['promotion_eligible'])
                    if item['report']['status'] == 'complete':
                        self.assertEqual(receipt['matrix']['rank'], item['report']['columns'])

    def test_full_factorization_including_factors_omitted_by_producer(self):
        for n in range(5, 32, 2):
            product = 1
            for factor in complete_factors(n):
                product = poly_product(product, factor)
            self.assertEqual(product, (1 << n) | 1)
        self.assertEqual(complete_factors(29), (3, 536870911))

    def test_valid_group_points_cannot_claim_another_recipe(self):
        item = row()
        verify_queries(item['report'], item['report']['fixture'], 2)
        recipe = dict(kind='divisor', indices=[0, 2])
        item['job']['factor_base'] = recipe
        item['report']['effective_factor_base'] = recipe
        with self.assertRaisesRegex(InvalidEvidence, 'construction or ordering'):
            verify_base(item['report'], item['report']['fixture'], item['job'])

    def test_base_point_order_is_part_of_identity(self):
        item = row()
        item['report']['factor_base'].reverse()
        with self.assertRaisesRegex(InvalidEvidence, 'construction or ordering'):
            verify_base(item['report'], item['report']['fixture'], item['job'])

    def test_true_group_equations_cannot_substitute_frontend_tags(self):
        for name, substitute in [('a0-f4-dense', 'f5'), ('a0-sat_xor-dense', 'sat_cnf')]:
            item = row(name)
            item['job']['config']['solver'] = substitute
            item['report']['effective_config']['solver'] = substitute
            verify_queries(item['report'], item['report']['fixture'], 2)
            with self.assertRaisesRegex(InvalidEvidence, 'substitution'):
                check(item)

    def test_actual_dispatch_cannot_be_removed_or_relabelled(self):
        mutations = [lambda r: r['collector_dispatch'].update(pair_table=False),
                     lambda r: r['collector_dispatch'].update(field_kernel=None),
                     lambda r: r['descent_dispatch'].update(direct_collision=True),
                     lambda r: r.update(generic_runtime_policy=None),
                     lambda r: r.update(generic_admission_schema=None)]
        for mutate in mutations:
            item = row()
            mutate(item['report'])
            with self.assertRaises(InvalidEvidence):
                check(item)

    def test_stored_matrix_tampering_with_valid_original_certificate(self):
        def coefficient(r):
            entries = next(value['entries'] for value in r['relation_matrix']['rows'] if value['entries'])
            entries[0][1] = str(int(entries[0][1])+1)
        mutations = [coefficient,
                     lambda r: r['relation_matrix']['rows'][0].update(rhs='0'),
                     lambda r: r['relation_matrix']['column_points'].reverse(),
                     lambda r: r['relation_matrix']['rows'].pop(),
                     lambda r: r['relation_matrix'].update(solver='sparse-filter-block-wiedemann')]
        for mutate in mutations:
            item = row()
            mutate(item['report'])
            verify(item['report'], item['report']['fixture'], summands=2)
            with self.assertRaises(InvalidEvidence):
                check(item)

    def test_batch_attempt_counts_and_stop_cannot_be_forged(self):
        mutations = [lambda r: r['matrix_batches'].pop(),
                     lambda r: r['matrix_batches'][-1].update(solve_attempts=0),
                     lambda r: r['matrix_batches'][-1].update(accepted_rows=0),
                     lambda r: r['matrix_batches'][-1].update(verified=False),
                     lambda r: r['log_table_report'].update(solve_attempts=0)]
        for mutate in mutations:
            item = row()
            mutate(item['report'])
            with self.assertRaises(InvalidEvidence):
                check(item)

    def test_incomplete_preparation_retains_rank_attempts_and_unknown_cost(self):
        item = row('incomplete-pair_table')
        receipt = check(item)
        self.assertFalse(receipt['matrix']['certified_logs'])
        phases = verify_native(item['report'], item['job'], process_wall_ns=item['process_wall_ns'])
        ledger = scientific_ledger(phases['process_phases_ns'], unit='native_wall_ns',
                                   process_total=item['process_wall_ns'])
        self.assertIsNone(ledger['cold_operations'])
        self.assertIsNone(phases['online_wall_ns'])
        item['job']['config']['max_trials'] = 2
        item['report']['effective_config']['max_trials'] = 2
        with self.assertRaisesRegex(InvalidEvidence, 'stopped early'):
            check(item)

    def test_phase_fold_has_no_double_charge_or_missing_to_zero(self):
        item = row()
        phases = verify_native(item['report'], item['job'], process_wall_ns=item['process_wall_ns'])
        values = phases['process_phases_ns']
        ledger = scientific_ledger(values, unit='native_wall_ns', process_total=item['process_wall_ns'])
        self.assertEqual(ledger['cold_operations'], item['process_wall_ns'])
        self.assertEqual(ledger['operations']['target_descent'], sum(values[p] for p in
            ('target_query', 'target_pdp', 'target_relation_check', 'target_descent')))
        values['target_pdp'] = None
        ledger = scientific_ledger(values, unit='native_wall_ns', process_total=item['process_wall_ns'])
        self.assertIsNone(ledger['operations']['target_descent'])
        self.assertIsNone(ledger['cold_operations'])

    def test_source_build_binary_and_report_binding(self):
        build, source = build_fixture()
        self.assertIn('docs/ic/calibration.json', source['root_files'])
        identity = verify_build_record(build, source)
        self.assertEqual(row()['report']['generic_build'], identity)
        for mutate in [lambda b, s: s['root_files'].update({'Cargo.toml': '0'*64}),
                       lambda b, s: b['build']['flags'].update(rustflags='-C opt-level=0'),
                       lambda b, s: b['identity'].update(schema_version=True),
                       lambda b, s: b['identity'].update(source_manifest_sha256='0'*64)]:
            b, s = copy.deepcopy(build), copy.deepcopy(source)
            mutate(b, s)
            with self.assertRaises(InvalidEvidence):
                verify_build_record(b, s)
        # Exercise file integrity without executing a platform-specific artifact.
        with tempfile.TemporaryDirectory() as directory:
            executable = Path(directory) / 'worker'
            executable.write_bytes(b'fixture executable bytes')
            build['worker_sha256'] = digest(executable)
            report = row()['report']
            verify_binding(report, build, source, executable=executable)
            executable.write_bytes(b'changed')
            with self.assertRaisesRegex(InvalidEvidence, 'executable digest'):
                verify_binding(report, build, source, executable=executable)
            build['worker_sha256'] = digest(executable)
            report['generic_build']['build_sha256'] = '0'*64
            with self.assertRaisesRegex(InvalidEvidence, 'different build identity'):
                verify_binding(report, build, source, executable=executable)

    def test_source_bound_candidate_ignores_run_data_and_inactive_flags(self):
        item, (build, _) = row(), build_fixture()
        stages = check(item)
        method = method_record(item['job'], item['report'], stages, build)
        candidate = candidate_manifest(item['report']['fixture'], item['report'], method)
        changed = copy.deepcopy(item)
        changed['job']['algorithm_seed'] += 1
        changed['job']['target_seeds'] = [99]
        changed['job']['config']['groebner_degree'] = 4  # ignored by pair-table backend
        changed['report']['effective_config']['groebner_degree'] = 4
        self.assertEqual(method_record(changed['job'], changed['report'], stages, build), method)
        self.assertIn('fb54PDP2pair', candidate['candidate_id'])
        build['build_sha256'] = '0'*64
        other = method_record(item['job'], item['report'], stages, build)
        self.assertNotEqual(sha256(other), sha256(method))

    def test_recursive_sparse_defaults_match_serde_without_boolean_integer_alias(self):
        item = row('a0-pair_table-sparse')
        item['job']['config']['sparse'] = {'filter': {'target_excess': 32}}
        effective_config(item['job'], item['report'])
        item['job']['config']['sparse']['filter']['remove_duplicates'] = 1
        with self.assertRaises(InvalidEvidence):
            effective_config(item['job'], item['report'])


if __name__ == '__main__':
    unittest.main()
