import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

from autolab import check_arms, run_native, summarize, verify_trial
from oracle import InvalidEvidence
from portfolio import factorial_candidates, recombine, retain
from tournament import BASE_CONFIG, digest, objhash, write
from test_certificate import FIXTURE, REPORT


class PortfolioTests(unittest.TestCase):
    def test_separate_rho_reference_survives_every_tournament_stage(self):
        from tournament import stage_arms
        arms=[dict(id='incumbent',config=BASE_CONFIG),dict(id='candidate',config=dict(BASE_CONFIG,batch_trials=1))]
        rho=dict(id='rho',config=dict(BASE_CONFIG,rho_parallel_walks=8),binary_relative='rho_reference/worker')
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            write(root/'contract.json',{'rho_reference':rho})
            write(root/'summaries/smoke.json',{'failures':[]})
            write(root/'summaries/development.json',{'retained_portfolio':[{'candidate':'candidate'}]})
            write(root/'summaries/selection.json',{'provisional_challenger':'candidate'})
            for stage in ('smoke','development','selection','confirmation','replay'):
                active=stage_arms(root,stage,arms)
                self.assertEqual(next(a for a in active if a['id']=='rho'),rho)
            self.assertNotIn('rho',[a['id'] for a in stage_arms(root,'aa',arms)])

    def test_confirmation_points_cannot_repeat_development_under_new_seeds(self):
        from tournament import reserve_targets
        seen = set()
        self.assertTrue(reserve_targets({'targets':[['3','5'],['7','9']]},seen))
        before = set(seen)
        self.assertFalse(reserve_targets({'targets':[[3,5],[11,13]]},seen))
        self.assertEqual(seen,before)
        self.assertFalse(reserve_targets({'targets':[[11,13],[11,13]]},seen))
        self.assertTrue(reserve_targets({'targets':[[11,13],[15,17]]},seen))

    def test_nonfinite_and_nonpositive_costs_cannot_promote(self):
        from tournament import comparison
        from test_tournament import rows
        for cost in (0,-1,float('nan'),float('inf'),True):
            data=rows(.5)
            data[0]['total_operations']=cost
            self.assertFalse(comparison(data,'candidate',draws=20)['eligible'])

    def row(self, name, cost, cells, wall=None):
        return {'candidate': name, 'eligible': True, 'candidate_over_baseline': cost,
                'native_wall_candidate_over_baseline': wall or cost, 'per_cell': cells,
                'native_wall_per_cell': cells}

    def test_cell_specialist_and_exploration_survive_aggregate_ranking(self):
        arms = [{'id': n, 'config': dict(BASE_CONFIG, solver=s)} for n,s in
                [('fast','pair_table'), ('similar','pair_table'), ('specialist','sat_xor'),
                 ('slow','f5'), ('other','f4'), ('timeout','sat_cnf')]]
        rows = [self.row('fast',.7,{'small':.7,'large':.7}),
                self.row('similar',.71,{'small':.71,'large':.71}),
                self.row('specialist',.9,{'small':1.5,'large':.3}),
                self.row('slow',2,{'small':2,'large':2}),
                self.row('other',3,{'small':3,'large':3}),
                {'candidate':'timeout','eligible':False}]
        chosen = retain(rows, arms, width=3, exploration=1, seed=12)
        self.assertEqual([r['candidate'] for r in chosen[:2]], ['fast','specialist'])
        self.assertEqual(chosen[-1]['reason'], 'predeclared exploration slot')
        self.assertNotIn('timeout', [r['candidate'] for r in chosen])
        self.assertEqual(chosen, retain(rows, arms, width=3, exploration=1, seed=12))

    def test_combinations_include_individually_losing_parents(self):
        arms = [dict(id='incumbent',config=BASE_CONFIG),
                dict(id='batch',config=dict(BASE_CONFIG,batch_trials=1)),
                dict(id='dense',config=dict(BASE_CONFIG,linear_algebra='dense'))]
        out = recombine(arms,['batch','dense'])
        self.assertEqual(len(out),2)
        self.assertEqual(out[1]['parents'], ['batch','dense'])
        self.assertEqual(out[1]['config']['batch_trials'],1)
        self.assertEqual(out[1]['config']['linear_algebra'],'dense')
        arms[2]['source_root'] = '/different/source'
        self.assertEqual(len(recombine(arms,['batch','dense'])),1)

    def test_factorial_grid_covers_interactions_and_all_solver_families(self):
        arms = factorial_candidates(BASE_CONFIG)
        check_arms(arms)
        self.assertTrue(any(a['config']['batch_trials']==1 and a['config']['linear_algebra']=='dense'
                            and a['config'].get('collection_window')==8 for a in arms))
        algebra = factorial_candidates(BASE_CONFIG,panel='algebra')
        check_arms(algebra)
        self.assertEqual({a['config']['solver'] for a in algebra},
                         {'pair_table','enumerate','f4','f5','inherited_f4','sat_xor','sat_cnf'})
        check_arms(factorial_candidates(BASE_CONFIG,panel='factor-base'))

    def test_invalid_cost_and_unmatched_cells_are_rejected(self):
        arms = [dict(id='a',config=BASE_CONFIG),dict(id='b',config=BASE_CONFIG)]
        with self.assertRaises(InvalidEvidence):
            retain([self.row('a',float('nan'),{'c':1})],arms)
        with self.assertRaises(InvalidEvidence):
            retain([self.row('a',1,{'c':1}),self.row('b',1,{'other':1})],arms)


class NativeEvidenceTests(unittest.TestCase):
    def test_native_ratio_pairs_cases_and_weights_cells_equally(self):
        from autolab import paired_wall_ratio
        cases=[{'id':'a','job':{'degree':9,'curve_a':0}},
               {'id':'b','job':{'degree':13,'curve_a':0}},
               {'id':'c','job':{'degree':13,'curve_a':0}}]
        contract={'cases':cases,'repetitions':1}
        rows=[]
        for case,base,ratio in [('a',100,2),('b',1,.5),('c',2,.5)]:
            for arm,cost in [('incumbent',base),('candidate',base*ratio)]:
                rows.append(dict(case=case,arm=arm,repetition=0,status='VERIFIED',
                                 process={'whole_process_wall_seconds':cost}))
        self.assertAlmostEqual(paired_wall_ratio(contract,rows,'candidate'),1)
        rows[-1]['status']='TIMEOUT'
        self.assertIsNone(paired_wall_ratio(contract,rows,'candidate'))

    def test_watchdog_retains_timeout_and_raw_job(self):
        with tempfile.TemporaryDirectory() as tmp:
            # An executable script is enough to exercise the real process path.
            script = Path(tmp)/'sleeper'
            script.write_text('#!'+sys.executable+'\nimport time; time.sleep(5)\n')
            script.chmod(0o755)
            directory = Path(tmp)/'trial'
            p = run_native(script, {'test':1}, directory, .1)
            self.assertEqual(p['status'],'TIMEOUT')
            self.assertNotEqual(p['exit_code'],0)
            self.assertLess(p['whole_process_wall_seconds'],3)
            self.assertEqual(json.loads((directory/'job.json').read_text()), {'test':1})
            with self.assertRaises(FileExistsError):
                run_native(script, {}, directory, .1)

    def test_incomplete_fast_arm_cannot_receive_a_cost_or_promotion(self):
        contract = {'cases':[{'id':'c'}], 'repetitions':2, 'comparison_kind':'fixed-support',
                    'arms':[{'id':'incumbent','mode':'ic'},{'id':'failed','mode':'ic'}]}
        rows = [dict(arm=arm,case='c',repetition=i,status=status,
                     process={'whole_process_wall_seconds':cost},
                     certificate={'rank':3,'signed_base_size':55,'factor_base_sha256':'base'})
                for arm,status,cost in [('incumbent','VERIFIED',1),('failed','TIMEOUT',.01)] for i in range(2)]
        summary = summarize(contract,rows)
        self.assertFalse(summary['promotion_eligible'])
        self.assertIsNone(summary['table'][1]['median_cold_wall_seconds'])
        self.assertIsNone(summary['table'][0]['total_operations'])
        self.assertEqual(summary,summarize(contract,list(reversed(rows))))
        rows[1]['certificate']['factor_base_sha256'] = 'changed'
        with self.assertRaises(InvalidEvidence):
            summarize(contract, rows)

    def test_trial_replay_rejects_corruption_and_wrong_job(self):
        from oracle import verify
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case = {'id':'c','fixture':FIXTURE,'job':{'mode':'fixture'}}
            arm = {'id':'incumbent','mode':'ic','config':BASE_CONFIG}
            directory = root/'trials/c/incumbent/rep-0'
            directory.mkdir(parents=True)
            write(directory/'job.json',dict(case['job'],mode='ic',config=BASE_CONFIG))
            write(directory/'stdout.json',REPORT)
            (directory/'stderr.txt').write_text('')
            process = dict(status='EXITED',exit_code=0,whole_process_wall_seconds=.1)
            write(directory/'process.json',process)
            receipt = {'arm':'incumbent','case':'c','repetition':0,'status':'VERIFIED',
                       'case_sha256':objhash(case),'arm_sha256':objhash(arm),'process':process,
                       'certificate':verify(REPORT,FIXTURE,expected_mode='ic',summands=3),
                       'artifacts':{p.name:digest(p) for p in directory.iterdir()}}
            write(directory/'receipt.json',receipt)
            verify_trial(root,{},case,arm,0)
            bad = copy.deepcopy(REPORT)
            bad['solutions'][0]['recovered']='1622'
            write(directory/'stdout.json',bad)
            with self.assertRaisesRegex(InvalidEvidence,'changed trial artifact'):
                verify_trial(root,{},case,arm,0)
            # Rehashing a forged answer cannot bypass independent arithmetic.
            receipt['artifacts']['stdout.json']=digest(directory/'stdout.json')
            write(directory/'receipt.json',receipt)
            with self.assertRaisesRegex(InvalidEvidence,'incorrect scalar'):
                verify_trial(root,{},case,arm,0)

    def test_rho_null_rank_is_not_an_ic_rank_or_zero_cost(self):
        contract = {'cases':[{'id':'c'}], 'repetitions':1, 'comparison_kind':'fixed-support',
                    'arms':[{'id':'incumbent','mode':'ic'},{'id':'rho_w32','mode':'rho'}]}
        rows = [dict(arm='incumbent',case='c',repetition=0,status='VERIFIED',
                     process={'whole_process_wall_seconds':1},
                     certificate={'rank':3,'signed_base_size':55,'factor_base_sha256':'base'}),
                dict(arm='rho_w32',case='c',repetition=0,status='VERIFIED',
                     process={'whole_process_wall_seconds':.5},
                     certificate={'rank':None,'signed_base_size':None,'factor_base_sha256':None})]
        result = summarize(contract,rows)
        self.assertIsNone(result['table'][1]['rank_range'])
        self.assertEqual(result['table'][1]['median_cold_wall_seconds'],.5)


if __name__ == '__main__':
    unittest.main()
