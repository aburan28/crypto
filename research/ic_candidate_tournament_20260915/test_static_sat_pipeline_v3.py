"""Pure mathematical controls exercise the full controller, never a SAT claim.

Only the external PDP callback/preflight and source admission are replaced.
The real curve, real 63-point base, relation matrix, LA and scalar replay run.
"""
import contextlib
import io
import json
from pathlib import Path
import platform
import tempfile
import unittest
from unittest.mock import patch

from generic_bases import construct
from oracle import Curve
from run_generic_exact_yield_audit import exact_three_sum,pair_index
from static_sat_pipeline_v3 import run,save_progress

FIXTURE=dict(cofactor='2',curve_a=1,degree=17,generator=['43693','23339'],
             group_order='131174',irreducible=dict(degree=17,low_terms=[0,3]),
             **{'lambda':'17184'},subgroup_order='65587',target_scalar_constructed=False,
             target_seeds=[],targets=[])


class StaticSatPipelineV3Tests(unittest.TestCase):
    def controlled(self,temporary,*,failure=None):
        root=Path(temporary)
        out=root/'entry-output'
        out.mkdir()
        curve=Curve(FIXTURE)
        base,_=construct(curve,dict(kind='standard_subspace',dimension=6))
        pairs=pair_index(curve,base)
        panel=dict(relation_query_seed=2026093001,descent_query_seed=2026093002,
                   max_relation_queries=512,max_descent_queries=64,
                   target_input=dict(point=[52411,72106]))
        if failure=='rank':
            panel['max_relation_queries']=3
        if failure=='target':
            panel['max_descent_queries']=3
        args=dict(panel=panel,source={'scope':'mocked-source-admission; exact-oracle control'},
                  method={'control':True},candidate={'control':True},workload={'control':True},
                  seal=dict(candidate_id='CONTROL',workload_id='CONTROL',run_id='CONTROL',
                            panel_sha256='CONTROL'))
        (root/'execution.json').write_text(json.dumps(dict(arguments=args,asset_manifest={})))
        native=dict(platform=dict(system=platform.system(),machine=platform.machine()))
        timestamps=[]
        def oracle_query(panel,item,execution,c,points,directory):
            target=None if item['point'] is None else tuple(item['point'])
            indices=exact_three_sum(c,points,pairs,target)
            target_index=item['trial']-panel['max_relation_queries']
            status='VALID_POINT_WITNESS' if indices is not None else 'SOURCE_UNSAT'
            if failure=='rank' or (failure=='target' and target_index>=0):
                status='TIMEOUT'
            elif target_index==0:
                status='CONFLICT_BUDGET_INCONCLUSIVE'
            elif target_index==1:
                status='TIMEOUT'
            return dict(trial=item['trial'],probe_scalar=item['probe_scalar'],
                        public_point=item['point'],status=status,verification_wall_ns=0,
                        point_witness={'point_indices':indices} if status=='VALID_POINT_WITNESS' else None)
        def progress(directory,name,row):
            if name=='descent.progress.jsonl' and row.get('scalar_replay_verified'):
                # Deliberately expensive final logging must follow the endpoint.
                import time
                timestamps.append(time.monotonic_ns())
                time.sleep(0.03)
            save_progress(directory,name,row)
        def preflight(*unused):
            (out/'cms_preflight.stdout').write_text('CryptoMiniSat version 5.14.7\n')
            return dict(returncode=0,timed_out=False)
        with contextlib.ExitStack() as stack:
            stack.enter_context(patch('static_sat_pipeline_v3.check_extracted_assets',
                                     return_value={'exporter/build-record.json':b'{}'}))
            stack.enter_context(patch('static_sat_pipeline_v3.native_admission',
                                     return_value=(FIXTURE,{},curve,base,native)))
            stack.enter_context(patch('static_sat_pipeline_v3.mathematical_registration',return_value=args))
            stack.enter_context(patch('static_sat_pipeline_v3.meter',side_effect=preflight))
            stack.enter_context(patch('static_sat_pipeline_v3.one_query',side_effect=oracle_query))
            stack.enter_context(patch('static_sat_pipeline_v3.save_progress',side_effect=progress))
            stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
            outcome=run(args,out)
        return outcome,json.loads((out/'summary.json').read_text()),timestamps,curve

    def test_full_rows_and_descent_include_failures_and_stop_before_final_logging(self):
        with tempfile.TemporaryDirectory() as temporary:
            outcome,result,timestamps,curve=self.controlled(temporary)
            self.assertEqual(outcome['status'],'COMPLETE')
            self.assertEqual(result['matrix']['rank'],29)
            self.assertGreater(result['matrix']['dependent_relations'],0)
            self.assertEqual([row['status'] for row in result['target_attempts'][:2]],
                             ['CONFLICT_BUDGET_INCONCLUSIVE','TIMEOUT'])
            self.assertEqual(curve.mul(curve.g,int(result['recovered_scalar'])),(52411,72106))
            self.assertTrue(result['scalar_verified'])
            self.assertEqual(sum(result['online_phases_ns'].values()),result['online_wall_ns'])
            self.assertEqual(result['online_attempt_wall_ns'],result['online_wall_ns'])
            self.assertEqual(result['online_stop_event'],'independent-scalar-replay')
            self.assertLess(result['online_stop_monotonic_ns'],timestamps[0])
            self.assertIsNone(result['online_speedup'])

    def test_rank_cap_preserves_failed_attempts_and_never_starts_target(self):
        with tempfile.TemporaryDirectory() as temporary:
            outcome,result,_,_=self.controlled(temporary,failure='rank')
            self.assertEqual(outcome['status'],'INCOMPLETE_RELATION_RANK')
            self.assertEqual(len(result['collection']),3)
            self.assertTrue(all(row['status']=='TIMEOUT' for row in result['collection']))
            self.assertEqual(result['matrix']['rank'],0)
            self.assertEqual(result['target_attempts'],[])
            self.assertFalse(result['scalar_verified'])
            self.assertIsNone(result['online_wall_ns'])

    def test_failed_target_keeps_charged_interval_without_verified_cost(self):
        with tempfile.TemporaryDirectory() as temporary:
            outcome,result,_,_=self.controlled(temporary,failure='target')
            self.assertEqual(outcome['status'],'INCOMPLETE_TARGET')
            self.assertEqual(len(result['target_attempts']),3)
            self.assertFalse(result['scalar_verified'])
            self.assertIsNone(result['recovered_scalar'])
            self.assertIsNone(result['online_wall_ns'])
            self.assertEqual(sum(result['online_phases_ns'].values()),result['online_attempt_wall_ns'])
            self.assertEqual(result['online_stop_event'],'frozen-target-attempt-cap')


if __name__=='__main__':
    unittest.main()
