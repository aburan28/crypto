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
spec=importlib.util.spec_from_file_location('tail_analysis',HERE/'analyze.py')
analysis=importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='tail-evidence-')
        self.root=Path(self.temp.name)/'run'
        shutil.copytree(HERE/'run_09',self.root)

    def tearDown(self):self.temp.cleanup()

    def rejected(self):
        with contextlib.redirect_stdout(io.StringIO()),self.assertRaises(AssertionError):analysis.main(self.root)

    def mutate_record(self,n,predicate,mutate):
        path=self.root/f'raw-n{n}.jsonl'
        rows=[json.loads(line) for line in path.read_text().splitlines()]
        row=next(r for r in rows if predicate(r));cell=row['cell'];mutate(row)
        row['raw_line']=json.dumps({k:v for k,v in row.items() if k not in ('cell','split','raw_line')})+'\n'
        path.write_text(''.join(json.dumps(r)+'\n' for r in rows))
        receipts_path=self.root/'receipts.json';receipts=json.loads(receipts_path.read_text())
        receipt=next(r for r in receipts if r['cell']==cell)
        receipt['stdout_sha256']=hashlib.sha256(''.join(r['raw_line'] for r in rows if r['cell']==cell).encode()).hexdigest()
        receipts_path.write_text(json.dumps(receipts))

    def test_all_complete_runs_replay_exactly(self):
        for run in ['run_01','run_02','run_04','run_05','run_06','run_07','run_08','run_09']:
            with self.subTest(run=run):
                target=Path(self.temp.name)/run;shutil.copytree(HERE/run,target)
                for name,expected in json.loads((target/'manifest.json').read_text())['files'].items():
                    self.assertEqual(analysis.sha(target/name),expected,name)
                expected=(target/'results.json').read_bytes();report=(target/'RESULT.md').read_bytes()
                with contextlib.redirect_stdout(io.StringIO()):analysis.main(target)
                self.assertEqual(expected,(target/'results.json').read_bytes())
                self.assertEqual(report,(target/'RESULT.md').read_bytes())

    def test_source_lineage_and_fresh_holdouts(self):
        seen=set()
        for run in ['run_01','run_02','run_04','run_05','run_06','run_07','run_08','run_09']:
            protocol=json.loads((HERE/run/'protocol.json').read_text())
            self.assertFalse(seen & set(protocol['holdout_seeds']));seen.update(protocol['holdout_seeds'])
        for child,parent in [('run_02','run_01'),('run_04','run_02'),('run_05','run_04'),('run_06','run_05'),('run_07','run_06'),('run_08','run_07'),('run_09','run_08')]:
            p=json.loads((HERE/child/'protocol.json').read_text())
            self.assertEqual(p['predecessor_worker_sha256'],analysis.sha(HERE/parent/'worker.rs'))
        self.assertEqual((HERE/'worker.rs').read_bytes(),(HERE/'run_09/worker.rs').read_bytes())
        four=json.loads((HERE/'run_04/protocol.json').read_text());five=json.loads((HERE/'run_05/protocol.json').read_text())
        for field in ['performance_gate','portfolio_gate','variables','families','batches','repetitions']:
            self.assertEqual(four[field],five[field])

    def test_failed_launch_is_preserved_not_counted_as_a_run(self):
        failure=json.loads((HERE/'LAUNCH_FAILURE.json').read_text())
        self.assertEqual(failure['completed_cells'],0);self.assertEqual(failure['measured_samples'],0)
        for name,expected in failure['files'].items():self.assertEqual(analysis.sha(HERE/name),expected,name)
        metadata=json.loads((HERE/'run_03/metadata.json').read_text());self.assertFalse(metadata['complete'])
        receipts=json.loads((HERE/'run_03/receipts.json').read_text())
        self.assertEqual(len(receipts),1);self.assertNotEqual(receipts[0]['exit_code'],0)
        self.assertEqual(receipts[0]['command'].count('with-census'),2)
        self.assertEqual((HERE/'run_03/failed.stdout').read_bytes(),b'')
        for name in ['worker.rs','protocol.json','analyze.py']:
            self.assertEqual((HERE/'run_03'/name).read_bytes(),(HERE/'run_04'/name).read_bytes())

    def test_changed_source_rejected(self):
        with (self.root/'worker.rs').open('a') as out:out.write('// post-measurement change\n')
        self.rejected()

    def test_missing_sample_rejected(self):
        path=self.root/'raw-n36.jsonl';path.write_text(''.join(path.read_text().splitlines(keepends=True)[:-1]));self.rejected()

    def test_incomplete_campaign_rejected(self):
        path=self.root/'metadata.json';row=json.loads(path.read_text());row['complete']=False;path.write_text(json.dumps(row));self.rejected()

    def test_wrong_exact_column_census_rejected(self):
        self.mutate_record(36,lambda r:r['type']=='sample' and r['variant']=='stream_census',lambda r:r.update(source_columns=r['source_columns']-1));self.rejected()

    def test_false_cache_hit_rejected(self):
        self.mutate_record(36,lambda r:r['type']=='sample' and r['variant']=='flat_cached',lambda r:r.update(layout_hits=r['layout_hits']+1));self.rejected()

    def test_wrong_dispatch_rejected(self):
        self.mutate_record(36,lambda r:r['type']=='sample' and r['variant']=='hybrid_census' and 'restricted_cycle' in r['cell'],lambda r:r.update(flat_dispatches=0));self.rejected()

    def test_missing_census_storage_rejected(self):
        self.mutate_record(36,lambda r:r['type']=='sample' and r['variant']=='stream_census',lambda r:r.update(retained_census_entries=0));self.rejected()

    def test_dropped_affine_consequences_rejected(self):
        self.mutate_record(36,lambda r:r['type']=='sample' and r['variant']=='hybrid_census' and 'cross_cancel' in r['cell'],lambda r:r.update(tail_rank=0));self.rejected()

    def test_wrong_occurring_mask_rejected(self):
        self.mutate_record(36,lambda r:r['type']=='fixture' and r['family']=='restricted_cycle',lambda r:r['inputs'][0].update(active=(1<<36)-1));self.rejected()

    def test_inner_batch_normalization_is_explicit(self):
        result=json.loads((self.root/'results.json').read_text())
        self.assertEqual(result['inner_batches'],16)
        records=[json.loads(line) for line in (self.root/'raw-n36.jsonl').read_text().splitlines()]
        sample=next(r for r in records if r['type']=='sample' and r['variant']=='flat_cached' and 'quadratic-b8' in r['cell'])
        cell=next(c for c in result['cells_detail'] if c['cell']==sample['cell'])
        import statistics
        values=[r['total_ns'] for r in records if r['cell']==sample['cell'] and r['type']=='sample' and r['variant']=='flat_cached']
        self.assertEqual(cell['arms']['flat_cached']['median_total_ns'],statistics.median(values)/16)
        self.assertEqual(sample['verified_outputs'],128)

    def test_caches_cannot_be_warmed_across_inner_batches(self):
        self.mutate_record(36,lambda r:r['type']=='sample' and r['variant']=='flat_cached' and 'restricted_cycle-b8' in r['cell'],
                           lambda r:r.update(layout_hits=r['layout_hits']+15))
        self.rejected()

    def test_false_inner_batch_count_rejected(self):
        self.mutate_record(36,lambda r:r['type']=='sample' and r['variant']=='hybrid_census',lambda r:r.update(inner_batches=1));self.rejected()

    def test_predecessor_pairing_cannot_be_changed(self):
        self.mutate_record(36,lambda r:r['type']=='sample',lambda r:r.update(rep=(r['rep']+1)%42));self.rejected()


if __name__=='__main__':unittest.main()
