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
    def test_frozen_runs_and_manifests_replay(self):
        for run in ['run_01','run_02']:
            target=Path(self.temp.name)/run
            shutil.copytree(HERE/run,target)
            for name,digest in json.loads((target/'manifest.json').read_text())['files'].items():
                self.assertEqual(analysis.sha(target/name),digest)
            expected=(target/'results.json').read_bytes();report=(target/'RESULT.md').read_bytes()
            with contextlib.redirect_stdout(io.StringIO()):analysis.main(target)
            self.assertEqual(expected,(target/'results.json').read_bytes())
            self.assertEqual(report,(target/'RESULT.md').read_bytes())

    def test_confirmation_changes_only_held_out_protocol_inputs(self):
        primary=json.loads((HERE/'run_01/protocol.json').read_text())
        confirm=json.loads((HERE/'run_02/protocol.json').read_text())
        self.assertFalse(set(primary['holdout_seeds'])&set(confirm['holdout_seeds']))
        self.assertEqual(confirm.pop('primary_manifest_sha256'),analysis.sha(HERE/'run_01/manifest.json'))
        confirm.pop('confirmation')
        confirm['holdout_seeds']=primary['holdout_seeds']
        self.assertEqual(primary,confirm)
        for name in ['worker.rs','kernel.rs','quadratic.rs','run.py','analyze.py']:
            self.assertEqual((HERE/'run_01'/name).read_bytes(),(HERE/'run_02'/name).read_bytes())
        self.assertEqual((HERE/'confirmation_protocol.json').read_bytes(),(HERE/'run_02/protocol.json').read_bytes())

    def test_lineage_and_fresh_holdouts(self):
        predecessor=HERE.parent/'boolean_solve_trace_20260923'
        frozen=json.loads((HERE/'run_01/protocol.json').read_text())
        for run in ['run_01','run_02','run_03']:
            earlier=json.loads((predecessor/run/'protocol.json').read_text())
            self.assertFalse(set(frozen['holdout_seeds'])&set(earlier['holdout_seeds']))
        repo=HERE.parents[1]
        lineage=json.loads((HERE/'SOURCE_LINEAGE.json').read_text())
        self.assertEqual(lineage['source_files'],frozen['predecessor_source_hashes'])
        for name,digest in lineage['source_files'].items():
            self.assertEqual(analysis.sha(repo/name),digest)
        self.assertEqual((HERE/'kernel.rs').read_bytes(),(predecessor/'kernel.rs').read_bytes())
        for name in ['worker.rs','kernel.rs','quadratic.rs','protocol.json','run.py','analyze.py']:
            self.assertEqual((HERE/name).read_bytes(),(HERE/'run_01'/name).read_bytes())

    def test_fastest_control_cannot_be_replaced_by_flat(self):
        # Synthetic timing mutation in a disposable copy: flat looks 5x slower,
        # but another retained method is twice as fast as the candidate.
        path=self.root/'raw-n16.jsonl'
        rows=[json.loads(line) for line in path.read_text().splitlines()]
        for row in rows:
            if row['type']!='sample':continue
            total={'bucket':100,'quadratic_state':200}.get(row['variant'],1000)
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
        for gate in result['quadratic_state_details']:
            self.assertEqual(set(gate['controls']),{'search','flat','bucket','hybrid','small_flat','word_tail','merge_search'})
            if gate['n']==16:
                self.assertEqual(gate['paired_ratio_median'],0.5)
                self.assertEqual(gate['ci95_paired_median'],[0.5,0.5])
                self.assertFalse(gate['pass'])
        self.assertEqual(result['quadratic_state_gate'],'REJECTED')

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
        self.assertEqual(result['quadratic_state_gate'],'REJECTED')
        affected=next(c for c in result['cells_detail'] if c['cell']==cell)
        self.assertIsNone(affected['arms']['word_tail']['completion_ns'])
        self.assertGreater(affected['arms']['word_tail']['observed_total_ns'],0)

    def test_censored_quadratic_search_has_no_completion_cost(self):
        path=self.root/'raw-n24.jsonl'
        rows=[json.loads(line) for line in path.read_text().splitlines()]
        cell=next(r['cell'] for r in rows if r['type']=='fixture' and r['split']=='holdout' and r['family']=='planted')
        limit=json.loads((self.root/'protocol.json').read_text())['limits']['nodes']
        for row in rows:
            if row['cell']==cell and row['type']=='sample' and row['variant'] in ['search','merge_search','quadratic_state']:
                row.update(outcome='UNKNOWN',model=None,reason='NODE_CAP',verified=False,nodes=limit)
                row['raw_line']=json.dumps({k:v for k,v in row.items() if k not in ['cell','split','raw_line']})+'\n'
        path.write_text(''.join(json.dumps(r)+'\n' for r in rows))
        path=self.root/'receipts.json';receipts=json.loads(path.read_text())
        receipt=next(r for r in receipts if r['cell']==cell)
        receipt['stdout_sha256']=hashlib.sha256(''.join(r['raw_line'] for r in rows if r['cell']==cell).encode()).hexdigest()
        path.write_text(json.dumps(receipts))
        with contextlib.redirect_stdout(io.StringIO()):analysis.main(self.root)
        result=json.loads((self.root/'results.json').read_text())
        self.assertEqual(result['quadratic_state_gate'],'REJECTED')
        affected=next(c for c in result['cells_detail'] if c['cell']==cell)
        for arm in ['search','merge_search','quadratic_state']:
            self.assertIsNone(affected['arms'][arm]['completion_ns'])
            self.assertGreater(affected['arms'][arm]['observed_total_ns'],0)

    def test_quadratic_trace_change_rejected(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='quadratic_state',lambda r:r.update(trace=r['trace']^1))
        self.rejected()

    def test_quadratic_work_change_rejected(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='quadratic_state',lambda r:r.update(specialized_terms=r['specialized_terms']+1))
        self.rejected()

    def test_changed_quadratic_source_rejected(self):
        with (self.root/'quadratic.rs').open('a') as out:out.write('// changed\n')
        self.rejected()

    def test_quadratic_state_does_no_kernel_work(self):
        result=json.loads((self.root/'results.json').read_text())
        for cell in result['cells_detail']:
            arm=cell['arms']['quadratic_state']
            self.assertEqual(arm['kernel_calls'],0)
            self.assertEqual(arm['median_kernel_ns'],0)
            self.assertFalse(any(arm['kernel_calls_by_active']))
        self.assertIsNone(result['production_solver_cost'])
        self.assertIsNone(result['full_ic_cost'])
        self.assertIsNone(result['rho_ratio'])

if __name__=='__main__':unittest.main()
