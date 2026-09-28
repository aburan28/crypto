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
spec=importlib.util.spec_from_file_location('initial_analysis',HERE/'analyze.py')
analysis=importlib.util.module_from_spec(spec);spec.loader.exec_module(analysis)

class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.assertTrue((HERE/'run_01/manifest.json').exists(),'Primary evidence must be complete')
        self.temp=tempfile.TemporaryDirectory(prefix='initial-evidence-')
        self.root=Path(self.temp.name)/'run'
        shutil.copytree(HERE/'run_01',self.root)
    def tearDown(self):self.temp.cleanup()
    def replay(self):
        with contextlib.redirect_stdout(io.StringIO()):analysis.main(self.root)
    def rejected(self):
        with self.assertRaises(AssertionError):self.replay()
    def mutate(self,predicate,change,all_matches=False):
        path=self.root/'raw-n24.jsonl';rows=[json.loads(line) for line in path.read_text().splitlines()]
        affected=set()
        for row in rows:
            if predicate(row):
                change(row);affected.add(row['cell'])
                row['raw_line']=json.dumps({k:v for k,v in row.items() if k not in ['cell','split','raw_line']})+'\n'
                if not all_matches:break
        self.assertTrue(affected)
        path.write_text(''.join(json.dumps(r)+'\n' for r in rows))
        path=self.root/'receipts.json';receipts=json.loads(path.read_text())
        for receipt in receipts:
            if receipt['cell'] in affected:
                receipt['stdout_sha256']=hashlib.sha256(''.join(r['raw_line'] for r in rows if r['cell']==receipt['cell']).encode()).hexdigest()
        path.write_text(json.dumps(receipts))
    def test_frozen_manifest_and_exact_analysis_replay(self):
        for name,digest in json.loads((self.root/'manifest.json').read_text())['files'].items():
            self.assertEqual(analysis.sha(self.root/name),digest)
        old=(self.root/'results.json').read_bytes();report=(self.root/'RESULT.md').read_bytes()
        self.replay()
        self.assertEqual((self.root/'results.json').read_bytes(),old)
        self.assertEqual((self.root/'RESULT.md').read_bytes(),report)
    def test_changed_source_rejected(self):
        with (self.root/'initial.rs').open('a') as f:f.write('// changed\n')
        self.rejected()
    def test_incomplete_producer_rejected(self):
        p=self.root/'metadata.json';m=json.loads(p.read_text());m['complete']=False;p.write_text(json.dumps(m));self.rejected()
    def test_missing_observation_rejected(self):
        p=self.root/'raw-n24.jsonl';p.write_text('\n'.join(p.read_text().splitlines()[:-1])+'\n');self.rejected()
    def test_bad_model_rejected_against_original_equations(self):
        def change(r):
            rows=[json.loads(s) for s in (self.root/'raw-n24.jsonl').read_text().splitlines()]
            fixture=next(s for s in rows if s['cell']==r['cell'] and s['type']=='fixture')
            r['model']=next(p for p in range(1<<12) if not analysis.satisfies(fixture['polys'],p))
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='initial_simd' and r['outcome']=='SAT',change)
        self.rejected()
    def test_transport_checksum_mismatch_rejected(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='gray_delta_simd',lambda r:r.update(trace=r['trace']^1));self.rejected()
    def test_missing_initial_phase_rejected(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='initial_simd',lambda r:r.update(substitution_ns=None));self.rejected()
    def test_overlapping_phase_costs_rejected(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='initial_simd',lambda r:r.update(enumeration_ns=r['solve_ns']));self.rejected()
    def test_rank_cannot_be_claimed_after_full_quadratic_rank_certificate(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='initial_simd' and r['projection_certified']==1,lambda r:r.update(affine_eliminated=1,derived_rows=1));self.rejected()
    def test_projection_work_must_match_reference(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='initial_simd',lambda r:r.update(projection_xors=r['projection_xors']+1));self.rejected()
    def test_unknown_retains_work_and_blocks_all_gates(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant'] in ['initial_list','initial_simd'],
            lambda r:r.update(outcome='UNKNOWN',model=None,reason='ENUM_CAP',verified=False),True)
        self.replay();result=json.loads((self.root/'results.json').read_text())
        self.assertFalse(result['all_complete_verified'])
        for verdict in result['gates'].values():self.assertTrue(all(v=='REJECTED' for v in verdict.values()))
        for cell in result['cells_detail']:
            if cell['n']==24:
                for arm in ['initial_list','initial_simd']:
                    self.assertIsNone(cell['arms'][arm]['completion_ns'])
                    self.assertGreater(cell['arms'][arm]['observed_total_ns'],0)
    def test_retained_simd_control_cannot_be_omitted(self):
        p=self.root/'protocol.json';v=json.loads(p.read_text());v['reference_arms'].remove('gray_simd');p.write_text(json.dumps(v))
        p=self.root/'metadata.json';m=json.loads(p.read_text());m['source_hashes']['protocol.json']=analysis.sha(self.root/'protocol.json');p.write_text(json.dumps(m))
        self.rejected()
    def test_incorrect_arm_position_rejected(self):
        self.mutate(lambda r:r['type']=='sample',lambda r:r.update(order=(r['order']+1)%28));self.rejected()
    def test_discovery_bundles_remain_complete_and_unpromoted(self):
        for name in ['resource_probe_01','resource_probe_02']:
            root=HERE/name
            for file,digest in json.loads((root/'manifest.json').read_text())['files'].items():self.assertEqual(analysis.sha(root/file),digest)
            p=json.loads((root/'protocol.json').read_text());self.assertEqual(p['seeds'],[17,937])
            s=json.loads((root/'summary.json').read_text());self.assertEqual(s['observations'],1008);self.assertIsNone(s['promotion'])
            self.assertTrue(all(c['verified'] for c in s['cells']))
    def test_complete_predecessor_grid_and_fresh_holdouts_are_retained(self):
        def fixtures(root):
            out={}
            for p in root.glob('raw-n*.jsonl'):
                for line in p.read_text().splitlines():
                    r=json.loads(line)
                    if r['type']=='fixture':out[r['n'],r['seed'],r['family']]=(r['polys'],r['planted_witness'])
            return out
        old=fixtures(HERE.parent/'boolean_affine_substitution_20260923/run_01');new=fixtures(self.root)
        self.assertEqual(len(old),120);self.assertEqual(len(new),144)
        for key,value in old.items():self.assertEqual(new[key],value)
        self.assertEqual(len(set(new)-set(old)),24)
        prior=json.loads((HERE.parent/'boolean_affine_substitution_20260923/run_01/results.json').read_text())
        current=json.loads((self.root/'results.json').read_text())
        current={(c['n'],c['seed'],c['family']):c for c in current['cells_detail']}
        fields=['outcome','model','reason','nodes','kernel_calls','decisions','forced','specialized_terms','source_rows','source_columns','max_depth','affine_eliminated','derived_rows','enumeration_points','enumeration_batches','enumeration_leaves','trace']
        for cell in prior['cells_detail']:
            got=current[cell['n'],cell['seed'],cell['family']]
            for arm,reference in cell['arms'].items():
                self.assertEqual({k:got['arms'][arm][k] for k in fields},{k:reference[k] for k in fields})
    def test_all_positions_are_balanced_over_each_size_cohort(self):
        p=json.loads((self.root/'protocol.json').read_text())
        for n in p['variables']:
            counts={a:[0]*28 for a in p['variants']}
            for line in (self.root/f'raw-n{n}.jsonl').read_text().splitlines():
                r=json.loads(line)
                if r['type']=='sample':counts[r['variant']][r['order']]+=1
            self.assertTrue(all(c==[9]*28 for c in counts.values()))
    def test_source_lineage_is_bound_and_unchanged_modules_match(self):
        lineage=json.loads((HERE/'SOURCE_LINEAGE.json').read_text())
        repo=HERE.parent.parent
        for name,digest in lineage['files'].items():self.assertEqual(analysis.sha(repo/name),digest)
        for name in ['kernel.rs','quadratic.rs','packed.rs','basis.rs','affine.rs','compact.rs','linearization.rs','syndrome.rs','untraced.rs']:
            self.assertEqual((HERE/name).read_bytes(),(HERE.parent/'boolean_affine_substitution_20260923'/name).read_bytes())

if __name__=='__main__':unittest.main()
