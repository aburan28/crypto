"""Adversarial controls for the registered two-family qualification rule."""
import unittest

from generic_backend_gate import ALIASES, CELLS, FAMILIES, GENERIC, REPETITIONS, STAGES, evaluate
from oracle import InvalidEvidence


def fixture(qualified=('generic_f4_dense', 'generic_sat_xor_dense')):
    failures = set(FAMILIES['f4_f5'] + FAMILIES['sat']) - set(qualified)
    summary = dict(trial_slots=750, verified_native_profile_pairs=750-len(failures),
                   retained_failures=len(failures), distinct_record_ids=1230,
                   exposed_points=25, status='AUDITED_WITH_FAILURES',
                   promotion_eligible=False)
    qualification = dict(schema_version=1, status='DEVELOPMENT_REFERENCES_SELECTED',
                         cells=list(CELLS), cases=15, repetitions=REPETITIONS,
                         measured_stages=['aa', 'smoke', 'development'],
                         heldout_data_used=False, promotion_eligible=False, table=[])
    natural = dict(schema_version=1, status='AUDITED',
                   planted_decompositions_used_for_yield=False,
                   observations=[], rows=[])
    for alias in ALIASES:
        if alias in GENERIC:
            failed = alias in failures
            qualification['table'].append(dict(alias=alias, mode='ic',
                scheduled_runs=45, verified_runs=44 if failed else 45,
                failures=[dict(status='INVALID_OR_INCOMPLETE')] if failed else [],
                smoke_failures=[], qualified=not failed,
                comparison_to_archived_incumbent=dict(eligible=not failed),
                source_manifest_sha256='0'*64))
        else:
            qualification['table'].append(dict(alias=alias))
    for stage, cases_per_cell in STAGES.items():
        for arm in GENERIC:
            for cell in CELLS:
                observations = []
                for case_number in range(cases_per_cell):
                    for rep in range(REPETITIONS):
                        failed = (arm in failures and stage == 'development'
                                  and cell == CELLS[0] and case_number == 0 and rep == 0)
                        observations.append(dict(stage=stage, arm=arm, cell=cell,
                            case=f'{stage}-{cell}-{case_number}', repetition=rep,
                            audited=True, execution_status='INVALID_OR_INCOMPLETE' if failed else 'VERIFIED',
                            query_count=2, witness_count=1,
                            outcome_mix=dict(witness=1, unresolved=1)))
                natural['observations'].extend(observations)
                count = len(observations)
                natural['rows'].append(dict(stage=stage, arm=arm, cell=cell,
                    scheduled_runs=count, audited_bounded_runs=count, censored_runs=0,
                    verified_complete_runs=sum(o['execution_status']=='VERIFIED'
                                               for o in observations),
                    observed_ordinary_queries=2*count, verified_witness_queries=count,
                    independently_checked_outcomes=dict(unresolved=count, witness=count),
                    actual_base_and_columns=[54, 3]))
    return summary, qualification, natural


class FamilyGateTests(unittest.TestCase):
    def test_reference_selection_does_not_qualify_failed_new_families(self):
        result = evaluate(*fixture(qualified=()))
        self.assertEqual(result['status'], 'NEGATIVE_FAMILY_QUALIFICATION')
        self.assertEqual(result['qualified_families'], dict(f4_f5=[], sat=[]))
        self.assertEqual(result['arms']['generic_f4_dense']['development_verified'], 44)
        self.assertEqual(result['arms']['generic_f4_dense']['ordinary_non_witness_queries_observed_across_processes'], 60)

    def test_one_complete_arm_from_each_family_passes_without_selecting_local_winner(self):
        result = evaluate(*fixture())
        self.assertEqual(result['status'], 'F4_F5_AND_SAT_QUALIFIED')
        self.assertEqual(result['qualified_families'],
                         dict(f4_f5=['generic_f4_dense'], sat=['generic_sat_xor_dense']))
        self.assertEqual(result['arms']['generic_f4_dense']['actual_base_and_columns_by_cell']['n31a0'], [54, 3])
        self.assertEqual(result['arms']['generic_sat_xor_dense']['censored_processes'], 0)

    def test_censored_profile_cannot_support_a_complete_arm(self):
        summary, qualification, natural = fixture()
        item = next(o for o in natural['observations'] if o['stage']=='smoke'
                    and o['arm']=='generic_f4_dense' and o['cell']==CELLS[0])
        item.update(audited=False, query_count=None, witness_count=None, outcome_mix=None)
        row = next(r for r in natural['rows'] if r['stage']=='smoke'
                   and r['arm']=='generic_f4_dense' and r['cell']==CELLS[0])
        row['audited_bounded_runs'] -= 1
        row['censored_runs'] += 1
        row['observed_ordinary_queries'] -= 2
        row['verified_witness_queries'] -= 1
        row['independently_checked_outcomes']['unresolved'] -= 1
        row['independently_checked_outcomes']['witness'] -= 1
        with self.assertRaises(InvalidEvidence):
            evaluate(summary, qualification, natural)

    def test_missing_repetition_is_not_a_negative_solver_result(self):
        summary, qualification, natural = fixture()
        natural['observations'].pop()
        with self.assertRaises(InvalidEvidence):
            evaluate(summary, qualification, natural)


if __name__ == '__main__':
    unittest.main()
