import contextlib
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import shutil
import tempfile
import unittest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('solve_analysis',HERE/'analyze.py')
analysis=importlib.util.module_from_spec(spec);spec.loader.exec_module(analysis)

class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='solve-evidence-');self.root=Path(self.temp.name)/'run'
        shutil.copytree(HERE/'run_01',self.root)
    def tearDown(self):self.temp.cleanup()
    def rejected(self):
        with contextlib.redirect_stdout(io.StringIO()),self.assertRaises(AssertionError):analysis.main(self.root)
    def mutate(self,predicate,change):
        path=self.root/'raw-n24.jsonl';rows=[json.loads(line) for line in path.read_text().splitlines()]
        row=next(r for r in rows if predicate(r));cell=row['cell'];change(row)
        row['raw_line']=json.dumps({k:v for k,v in row.items() if k not in ['cell','split','raw_line']})+'\n'
        path.write_text(''.join(json.dumps(r)+'\n' for r in rows))
        path=self.root/'receipts.json';receipts=json.loads(path.read_text());receipt=next(r for r in receipts if r['cell']==cell)
        receipt['stdout_sha256']=hashlib.sha256(''.join(r['raw_line'] for r in rows if r['cell']==cell).encode()).hexdigest();path.write_text(json.dumps(receipts))
    def test_frozen_run_and_manifest_replay(self):
        for name,digest in json.loads((self.root/'manifest.json').read_text())['files'].items():
            self.assertEqual(analysis.sha(self.root/name),digest)
        expected=(self.root/'results.json').read_bytes();report=(self.root/'RESULT.md').read_bytes()
        with contextlib.redirect_stdout(io.StringIO()):analysis.main(self.root)
        self.assertEqual(expected,(self.root/'results.json').read_bytes())
        self.assertEqual(report,(self.root/'RESULT.md').read_bytes())

    def test_source_lineage_and_retained_failed_grid(self):
        predecessor=HERE.parent/'boolean_quadratic_frontier_20260923'
        current=json.loads((HERE/'run_01/protocol.json').read_text())
        prior=json.loads((predecessor/'run_01/protocol.json').read_text())
        self.assertEqual(current['regression_seeds'],prior['regression_seeds']+prior['holdout_seeds'])
        self.assertEqual(current['discovery_seeds'],prior['discovery_seeds'])
        for dirname in ['boolean_quadratic_frontier_20260923','boolean_quadratic_state_20260923','boolean_solve_trace_20260923']:
            for path in (HERE.parent/dirname).glob('run_*/protocol.json'):
                earlier=json.loads(path.read_text())
                self.assertFalse(set(current['holdout_seeds'])&set(earlier['holdout_seeds']))
        lineage=json.loads((HERE/'SOURCE_LINEAGE.json').read_text())
        self.assertEqual(current['predecessor_source_hashes'],lineage['source_files'])
        for name,digest in lineage['source_files'].items():self.assertEqual(analysis.sha(HERE.parents[1]/name),digest)
        for name in ['kernel.rs','quadratic.rs','packed.rs']:
            self.assertEqual((HERE/name).read_bytes(),(predecessor/name).read_bytes())
        for name in ['worker.rs','kernel.rs','quadratic.rs','packed.rs','basis.rs','basis_tests.rs','protocol.json','run.py','analyze.py']:
            self.assertEqual((HERE/name).read_bytes(),(HERE/'run_01'/name).read_bytes())

    def test_discovery_probes_are_separate_and_hash_bound(self):
        for run,arms in [('resource_probe_01',11),('resource_probe_02',13)]:
            root=HERE/run
            for name,digest in json.loads((root/'manifest.json').read_text())['files'].items():self.assertEqual(analysis.sha(root/name),digest)
            receipt=json.loads((root/'receipt.json').read_text());self.assertTrue(receipt['complete'])
            self.assertTrue(all(r['exit_code']==0 for r in receipt['receipts']))
            for path in root.glob('*.jsonl'):
                records=[json.loads(line) for line in path.read_text().splitlines()]
                self.assertEqual(len(records),arms+1)
                fixture=records[0];self.assertEqual((fixture['n'],fixture['seed']),(24,17))
                for row in records[1:]:
                    self.assertTrue(row['verified'])
                    if row['outcome']=='SAT':self.assertTrue(analysis.satisfies(fixture['polys'],row['model']))
                    else:self.assertEqual((row['outcome'],fixture['search_reference']),('UNSAT','UNSAT'))

    def test_fastest_control_cannot_be_replaced_by_flat(self):
        # Synthetic timing mutation in a disposable copy: flat looks 5x slower,
        # but another retained method is twice as fast as the candidate.
        path=self.root/'raw-n16.jsonl'
        rows=[json.loads(line) for line in path.read_text().splitlines()]
        for row in rows:
            if row['type']!='sample':continue
            total={'packed_state':100,'basis_wide':200}.get(row['variant'],1000)
            row.update(total_ns=total,solve_ns=total*9//10,validation_ns=total//10,
                       kernel_ns=total//2 if row['kernel_calls'] else 0)
            row['raw_line']=json.dumps({k:v for k,v in row.items() if k not in ['cell','split','raw_line']})+'\n'
        path.write_text(''.join(json.dumps(r)+'\n' for r in rows))
        path=self.root/'receipts.json';receipts=json.loads(path.read_text())
        for receipt in receipts:
            if receipt['n']==16:
                receipt['stdout_sha256']=hashlib.sha256(''.join(r['raw_line'] for r in rows if r['cell']==receipt['cell']).encode()).hexdigest()
        path.write_text(json.dumps(receipts))
        with contextlib.redirect_stdout(io.StringIO()):analysis.main(self.root)
        result=json.loads((self.root/'results.json').read_text())
        for gate in result['basis_details']:
            if gate['candidate']!='basis_wide' or gate['reference']!='strongest':continue
            self.assertIn('packed_state',gate['controls'])
            self.assertEqual(len(gate['controls']),12)
            if gate['n']==16:
                self.assertEqual(gate['paired_ratio_median'],0.5)
                self.assertEqual(gate['ci95_paired_median'],[0.5,0.5])
                self.assertFalse(gate['pass'])
        self.assertEqual(result['basis_wide_dramatic_gate'],'REJECTED')

    def test_changed_source_rejected(self):
        with (self.root/'kernel.rs').open('a') as out:out.write('// changed\n')
        self.rejected()
    def test_missing_sample_rejected(self):
        path=self.root/'raw-n24.jsonl';path.write_text(''.join(path.read_text().splitlines(keepends=True)[:-1]));self.rejected()
    def test_incomplete_run_rejected(self):
        path=self.root/'metadata.json';data=json.loads(path.read_text());data['complete']=False;path.write_text(json.dumps(data));self.rejected()
    def test_bad_model_rejected_by_original_equations(self):
        records=[json.loads(line) for line in (self.root/'raw-n24.jsonl').read_text().splitlines()]
        sample=next(r for r in records if r['type']=='sample' and r['outcome']=='SAT')
        fixture=next(r for r in records if r['type']=='fixture' and r['cell']==sample['cell'])
        bad=next(m for m in range(1<<16) if not analysis.satisfies(fixture['polys'],m))
        self.mutate(lambda r:r['type']=='sample' and r['cell']==sample['cell'] and r['variant']==sample['variant'] and r['rep']==sample['rep'],lambda r:r.update(model=bad))
        self.rejected()
    def test_policy_trace_difference_rejected(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='word_tail',lambda r:r.update(trace=r['trace']^1));self.rejected()
    def test_impossible_kernel_time_rejected(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='word_tail',lambda r:r.update(kernel_ns=r['total_ns']+1));self.rejected()
    def test_wrong_active_width_profile_rejected(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='word_tail',lambda r:r['kernel_calls_by_active'].__setitem__(24,1));self.rejected()
    def test_false_unsat_for_planted_fixture_rejected(self):
        self.mutate(lambda r:r['type']=='sample' and '-planted' in r['cell'],lambda r:r.update(outcome='UNSAT',model=None));self.rejected()

    def test_search_trace_change_is_rejected(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='merge_search',lambda r:r.update(trace=r['trace']^1));self.rejected()

    def test_resource_censor_blocks_completion_and_promotion(self):
        path=self.root/'raw-n24.jsonl';rows=[json.loads(line) for line in path.read_text().splitlines()]
        cell=next(r['cell'] for r in rows if r['type']=='fixture' and r['split']=='holdout' and r['family']=='planted')
        limit=json.loads((self.root/'protocol.json').read_text())['limits']['nodes']
        for row in rows:
            if row['cell']==cell and row['type']=='sample' and row['variant'] in ['small_flat','word_tail']:
                row.update(outcome='UNKNOWN',model=None,reason='NODE_CAP',verified=False,nodes=limit)
                row['raw_line']=json.dumps({k:v for k,v in row.items() if k not in ['cell','split','raw_line']})+'\n'
        path.write_text(''.join(json.dumps(r)+'\n' for r in rows))
        path=self.root/'receipts.json';receipts=json.loads(path.read_text())
        receipt=next(r for r in receipts if r['cell']==cell)
        receipt['stdout_sha256']=hashlib.sha256(''.join(r['raw_line'] for r in rows if r['cell']==cell).encode()).hexdigest()
        path.write_text(json.dumps(receipts))
        with contextlib.redirect_stdout(io.StringIO()):analysis.main(self.root)
        result=json.loads((self.root/'results.json').read_text());self.assertFalse(result['all_results_completed_and_verified'])
        self.assertEqual(result['word_solve_gate'],'REJECTED');self.assertEqual(result['merge_search_gate'],'REJECTED')
        self.assertEqual(result['basis_wide_dramatic_gate'],'REJECTED')
        affected=next(c for c in result['cells_detail'] if c['cell']==cell)
        self.assertIsNone(affected['arms']['word_tail']['completion_ns'])
        self.assertGreater(affected['arms']['word_tail']['observed_total_ns'],0)

    def test_censored_basis_search_has_no_completion_cost(self):
        path=self.root/'raw-n24.jsonl'
        rows=[json.loads(line) for line in path.read_text().splitlines()]
        cell=next(r['cell'] for r in rows if r['type']=='fixture' and r['split']=='holdout' and r['family']=='planted')
        limit=json.loads((self.root/'protocol.json').read_text())['limits']['nodes']
        for row in rows:
            if row['cell']==cell and row['type']=='sample' and row['variant'] in ['basis_list','basis_wide']:
                row.update(outcome='UNKNOWN',model=None,reason='NODE_CAP',verified=False,nodes=limit)
                row['raw_line']=json.dumps({k:v for k,v in row.items() if k not in ['cell','split','raw_line']})+'\n'
        path.write_text(''.join(json.dumps(r)+'\n' for r in rows))
        path=self.root/'receipts.json';receipts=json.loads(path.read_text())
        receipt=next(r for r in receipts if r['cell']==cell)
        receipt['stdout_sha256']=hashlib.sha256(''.join(r['raw_line'] for r in rows if r['cell']==cell).encode()).hexdigest()
        path.write_text(json.dumps(receipts))
        with contextlib.redirect_stdout(io.StringIO()):analysis.main(self.root)
        result=json.loads((self.root/'results.json').read_text())
        self.assertEqual(result['basis_wide_dramatic_gate'],'REJECTED')
        affected=next(c for c in result['cells_detail'] if c['cell']==cell)
        for arm in ['basis_list','basis_wide']:
            self.assertIsNone(affected['arms'][arm]['completion_ns'])
            self.assertGreater(affected['arms'][arm]['observed_total_ns'],0)

    def test_basis_trace_change_rejected(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='basis_wide',lambda r:r.update(trace=r['trace']^1))
        self.rejected()

    def test_basis_work_change_rejected(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='basis_wide',lambda r:r.update(specialized_terms=r['specialized_terms']+1))
        self.rejected()

    def test_changed_basis_source_rejected(self):
        with (self.root/'basis.rs').open('a') as out:out.write('// changed\n')
        self.rejected()

    def test_policy_groups_are_preserved_without_cross_policy_trace_claims(self):
        result=json.loads((self.root/'results.json').read_text())
        keys=['outcome','model','reason','nodes','kernel_calls','decisions','forced','specialized_terms','source_rows','source_columns','max_depth','trace']
        for cell in result['cells_detail']:
            arms=cell['arms']
            for left,right in [('basis_list','basis_wide'),('tail_list','tail_wide')]:
                self.assertEqual({k:arms[left][k] for k in keys},{k:arms[right][k] for k in keys})
                self.assertGreater(arms[right]['kernel_calls'],0)
            self.assertEqual(arms['packed_state']['kernel_calls'],0)
        self.assertIsNone(result['production_solver_cost']);self.assertIsNone(result['full_ic_cost']);self.assertIsNone(result['rho_ratio'])

    def test_tail_trace_difference_is_rejected(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='tail_wide',lambda r:r.update(trace=r['trace']^1));self.rejected()

    def test_affine_opportunities_match_the_uninstrumented_discovery_traces(self):
        root=HERE/'affine_opportunity_01'
        for name,digest in json.loads((root/'manifest.json').read_text())['files'].items():self.assertEqual(analysis.sha(root/name),digest)
        receipt=json.loads((root/'receipt.json').read_text())
        self.assertTrue(receipt['complete']);self.assertEqual(receipt['exit_code'],0)
        self.assertEqual(receipt['primary_manifest_sha256'],analysis.sha(HERE/'run_01/manifest.json'))
        expected={}
        keys=['outcome','model','trace','nodes','decisions','forced','kernel_calls','specialized_terms','source_rows','source_columns','max_depth']
        for n in [16,24]:
            for line in (HERE/'run_01'/f'raw-n{n}.jsonl').read_text().splitlines():
                row=json.loads(line)
                if row['type']=='sample' and row['split']=='discovery' and row['variant']=='tail_wide' and row['rep']==0:
                    expected[row['cell']]={k:row[k] for k in keys}
        rows=[json.loads(line) for line in (root/'raw.jsonl').read_text().splitlines()]
        self.assertEqual(len(rows),12)
        for row in rows:
            self.assertIn(row['seed'],[17,937])
            cell=f"n{row['n']}-discovery-{row['seed']}-{row['family']}"
            self.assertEqual({k:row[k] for k in keys},expected[cell])
            self.assertLessEqual(row['affine_opportunity_calls'],row['decisions'])
            self.assertLessEqual(row['affine_opportunity_calls'],row['affine_opportunity_rank_sum'])
            self.assertLessEqual(row['affine_opportunity_rank_sum'],row['n']*row['affine_opportunity_calls'])
        summary=json.loads((root/'summary.json').read_text())
        self.assertEqual(summary['total_decisions'],sum(r['decisions'] for r in rows))
        self.assertEqual(summary['opportunity_calls'],sum(r['affine_opportunity_calls'] for r in rows))
        self.assertEqual(summary['rank_sum'],sum(r['affine_opportunity_rank_sum'] for r in rows))
        self.assertIsNone(summary['substitution_solver_cost'])
        replay=HERE/'affine_opportunity_02'
        for name,digest in json.loads((replay/'manifest.json').read_text())['files'].items():self.assertEqual(analysis.sha(replay/name),digest)
        for name in ['raw.jsonl','summary.json','worker.rs','kernel.rs','quadratic.rs','packed.rs','basis.rs','basis_tests.rs','SOURCE_LINEAGE.json']:
            self.assertEqual((root/name).read_bytes(),(replay/name).read_bytes())

if __name__=='__main__':unittest.main()
