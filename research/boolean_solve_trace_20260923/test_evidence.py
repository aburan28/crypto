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
        shutil.copytree(HERE/'run_03',self.root)
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
    def test_both_runs_and_manifests_replay(self):
        for run in ['run_01','run_02','run_03']:
            target=Path(self.temp.name)/run;shutil.copytree(HERE/run,target)
            for name,digest in json.loads((target/'manifest.json').read_text())['files'].items():self.assertEqual(analysis.sha(target/name),digest)
            expected=(target/'results.json').read_bytes();report=(target/'RESULT.md').read_bytes()
            with contextlib.redirect_stdout(io.StringIO()):analysis.main(target)
            self.assertEqual(expected,(target/'results.json').read_bytes());self.assertEqual(report,(target/'RESULT.md').read_bytes())
    def test_lineage_and_fresh_holdouts(self):
        first=json.loads((HERE/'run_01/protocol.json').read_text());second=json.loads((HERE/'run_02/protocol.json').read_text())
        self.assertFalse(set(first['holdout_seeds'])&set(second['holdout_seeds']))
        for name,digest in second['predecessor_source_hashes'].items():self.assertEqual(analysis.sha(HERE/'run_01'/name),digest)
        third=json.loads((HERE/'run_03/protocol.json').read_text())
        self.assertFalse(set(third['holdout_seeds'])&(set(first['holdout_seeds'])|set(second['holdout_seeds'])))
        for name,digest in third['predecessor_source_hashes'].items():self.assertEqual(analysis.sha(HERE/'run_02'/name),digest)
        for name in ['worker.rs','kernel.rs']:self.assertEqual((HERE/name).read_bytes(),(HERE/'run_03'/name).read_bytes())
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
        affected=next(c for c in result['cells_detail'] if c['cell']==cell)
        self.assertIsNone(affected['arms']['word_tail']['completion_ns'])
        self.assertGreater(affected['arms']['word_tail']['observed_total_ns'],0)

if __name__=='__main__':unittest.main()
