"""Unit checks for the F4/F5-only smoke gate."""
import unittest

from generic_backend_gate_v3 import CELLS, F4_F5, GENERIC, evaluate_v3
from oracle import InvalidEvidence


def _observation(stage, arm, cell, case, *, verified=True, audited=True):
    return dict(stage=stage, arm=arm, cell=cell, case=case, repetition=0,
                execution_status='VERIFIED' if verified else 'TIMEOUT',
                audited=audited, query_count=2 if audited else None,
                witness_count=0 if audited else None,
                outcome_mix={'unsat': 2} if audited else None,
                accepted_rows=0 if audited else None, rank=0 if audited else None,
                base_points=6 if audited else None, folded_columns=6 if audited else None,
                target_descent_queries=0 if audited else None,
                target_outcome_mix={} if audited else None,
                complete_instruction_phases=None if not verified else
                dict(queries=1, pdp=1, relation_check=1),
                censor_reason=None if audited else 'no complete profiler report')


def _row(stage, arm, cell, *, verified=1, audited=1, censored=0):
    return dict(stage=stage, arm=arm, cell=cell, scheduled_runs=1,
                verified_complete_runs=verified, audited_bounded_runs=audited,
                censored_runs=censored, independently_checked_outcomes={'unsat': 2},
                actual_base_and_columns=[6, 6] if audited else None)


def _bundle(*, complete=True):
    observations = []
    rows = []
    for arm in GENERIC:
        for cell in CELLS:
            case = f'{cell}-000'
            verified = complete or arm == 'generic_pair_subspace_dense'
            audited = verified
            observations.append(_observation('smoke', arm, cell, case,
                                             verified=verified, audited=audited))
            rows.append(_row('smoke', arm, cell, verified=int(verified),
                             audited=int(audited), censored=int(not audited)))
    summary = dict(registration_panel_sha256=
                   'df92d5507785446a2a5b333bd7a04776781de4f54c63c5ea3ffa921bc99ad2e4',
                   prior_censored_exposures_sha256=
                   'a728677b199eac02800d8338204d5306f391ec5da757c910bec1e51955fe7b41',
                   prior_v2_censored_exposures_sha256=
                   '0cc792cceb8c7190a533e6f4e665e8486911153ad57609ac8d1243af282b8e54',
                   qualification_schedule='smoke', trial_slots=45, repetitions=1,
                   exposed_points=10, verified_native_profile_pairs=45 if complete else 35,
                   retained_failures=0 if complete else 10,
                   status='AUDITED_COMPLETE' if complete else 'AUDITED_WITH_FAILURES',
                   promotion_eligible=False, sat_in_scope=False)
    natural = dict(schema_version=1, status='AUDITED', observations=observations,
                   rows=rows, qualification_schedule='smoke',
                   planted_decompositions_used_for_yield=False)
    return summary, natural


class GateV3Test(unittest.TestCase):
    def test_complete_f4_f5_smoke_qualifies(self):
        result = evaluate_v3(*_bundle(complete=True))
        self.assertEqual(result['status'], 'F4_F5_SMOKE_QUALIFIED')
        self.assertEqual(result['qualified_f4_f5_arms'], list(F4_F5))
        self.assertFalse(result['sat_in_scope'])

    def test_incomplete_algebraic_smoke_is_negative(self):
        result = evaluate_v3(*_bundle(complete=False))
        self.assertEqual(result['status'], 'NEGATIVE_FAMILY_QUALIFICATION')
        self.assertEqual(result['qualified_f4_f5_arms'], [])
        self.assertTrue(result['arms']['generic_pair_subspace_dense']['complete'])

    def test_wrong_schedule_is_rejected(self):
        summary, natural = _bundle(complete=True)
        summary['qualification_schedule'] = 'full'
        with self.assertRaises(InvalidEvidence):
            evaluate_v3(summary, natural)


if __name__ == '__main__':
    unittest.main()
