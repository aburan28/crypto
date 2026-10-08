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
spec=importlib.util.spec_from_file_location('construction_analysis',HERE/'corrected_analysis.py')
analysis=importlib.util.module_from_spec(spec);spec.loader.exec_module(analysis)

class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.assertTrue((HERE/'run_01/manifest.json').exists())
        self.temp=tempfile.TemporaryDirectory(prefix='construction-scratch-',dir=HERE)
        self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)/'run'
        shutil.copytree(HERE/'run_01',self.root)
        self.output=Path(self.temp.name)/'analysis';shutil.copytree(HERE/'analysis_01',self.output)
    def replay(self):
        with contextlib.redirect_stdout(io.StringIO()):analysis.main(self.root,self.output)
    def rejected(self):
        with self.assertRaises(AssertionError):self.replay()
    def mutate(self,predicate,change,all_matches=False):
        p=self.root/'n24-17-unplanted.jsonl';rows=[json.loads(line) for line in p.read_text().splitlines()];changed=0
        for row in rows:
            if predicate(row):
                change(row);changed+=1
                if not all_matches:break
        self.assertGreater(changed,0);p.write_text(''.join(json.dumps(r)+'\n' for r in rows))
        q=self.root/'receipts.json';receipts=json.loads(q.read_text())
        next(r for r in receipts if r['cell']=='n24-17-unplanted')['stdout_sha256']=analysis.sha(p)
        q.write_text(json.dumps(receipts))
    def test_frozen_hashes_and_exact_analysis_replay(self):
        for name,digest in json.loads((self.root/'manifest.json').read_text())['files'].items():self.assertEqual(analysis.sha(self.root/name),digest)
        old=(self.output/'results.json').read_bytes();note=(self.output/'RESULT.md').read_bytes();self.replay()
        self.assertEqual((self.output/'results.json').read_bytes(),old);self.assertEqual((self.output/'RESULT.md').read_bytes(),note)
    def test_original_sources_are_retained_without_algorithm_changes(self):
        source=HERE.parent/'boolean_byte_sieve_20260923'
        for name,digest in json.loads((HERE/'SOURCE_LINEAGE.json').read_text())['files'].items():
            self.assertEqual(analysis.sha(HERE.parent.parent/name),digest)
            actual=HERE/Path(name).name
            if actual.name=='worker.rs':self.assertTrue(actual.read_bytes().startswith((source/actual.name).read_bytes()))
            else:self.assertEqual(actual.read_bytes(),(source/actual.name).read_bytes())
    def test_changed_scan_source_rejected(self):
        with (self.root/'construction.rs').open('a') as f:f.write('// changed\n')
        self.rejected()
    def test_missing_sample_rejected(self):
        p=self.root/'n24-17-unplanted.jsonl';p.write_text('\n'.join(p.read_text().splitlines()[:-1])+'\n');self.rejected()
    def test_profile_cannot_drop_the_scan_phase(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='construction64_profile',lambda r:r['phases'].pop('scan_ns'));self.rejected()
    def test_unprofiled_costs_are_not_zero_filled(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='construction64_rebound',lambda r:r.update(phases={k:0 for k in analysis.PHASES}));self.rejected()
    def test_phase_sum_cannot_exceed_total(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='construction64_profile',lambda r:r['phases'].update(scan_ns=r['solve_ns']+1));self.rejected()
    def test_changed_scan_work_rejected_against_original_policy(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='construction64_profile',lambda r:r['logical'].update(enumeration_points=r['logical']['enumeration_points']-64));self.rejected()
    def test_zero_scan_has_no_finite_ceiling(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='construction64_profile',lambda r:r['phases'].update(scan_ns=0),True)
        self.replay();result=json.loads((self.output/'results.json').read_text())
        group=next(v for v in result['bounds'] if (v['n'],v['family'],v['kind'])==(24,'unplanted','64'))
        self.assertIsNone(group['optimistic_ceiling_ci95']);self.assertEqual(group['decision'],'INCONCLUSIVE')
    def test_large_observer_cost_blocks_representativeness(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='construction64_profile',lambda r:r.update(total_ns=r['total_ns']*10),True)
        self.replay();result=json.loads((self.output/'results.json').read_text())
        group=next(v for v in result['bounds'] if (v['n'],v['family'],v['kind'])==(24,'unplanted','64'))
        self.assertFalse(group['comparable']);self.assertEqual(group['decision'],'INCONCLUSIVE')
    def test_adjustment_only_makes_the_ceiling_more_optimistic(self):
        self.assertEqual(analysis.adjusted_scan(100,180,160,10),70)
        self.assertEqual(analysis.adjusted_scan(100,150,160,10),90)
        self.assertEqual(analysis.adjusted_scan(100,300,160,10),0)
        for observer in range(0,200,7):
            self.assertLessEqual(analysis.adjusted_scan(100,160+observer,160,10),100)
    def test_reference_roster_cannot_omit_a_fast_retained_method(self):
        p=self.root/'protocol.json';v=json.loads(p.read_text());v['reference_arms'].remove('wide64_unrolled');p.write_text(json.dumps(v))
        p=self.root/'metadata.json';m=json.loads(p.read_text());m['source_hashes']['protocol.json']=analysis.sha(self.root/'protocol.json');p.write_text(json.dumps(m));self.rejected()
    def test_censored_work_stays_observed_without_a_constructor_claim(self):
        names=['construction16_profile','construction16_rebound','word16_unrolled']
        self.mutate(lambda r:r['type']=='sample' and r['variant'] in names,lambda r:r.update(outcome='UNKNOWN',model=None,reason='ENUM_CAP',verified=False),True)
        self.replay();result=json.loads((self.output/'results.json').read_text());self.assertFalse(result['all_complete_verified'])
        group=next(v for v in result['bounds'] if (v['n'],v['family'],v['kind'])==(24,'unplanted','16'))
        self.assertIsNone(group['optimistic_ceiling_ci95']);self.assertEqual(group['decision'],'INCONCLUSIVE')
        rows=[r for r in result['costs'] if r['n']==24 and r['family']=='unplanted' and r['variant'] in names]
        self.assertTrue(all(r['completion_ns'] is None and r['observed_ns']>0 for r in rows))
        self.assertIsNone(result['performance_promotion']);self.assertFalse(result['constructor_optimization_implemented'])
    def test_random_orders_are_exact_reverse_pairs(self):
        p=json.loads((self.root/'protocol.json').read_text())
        for path in self.root.glob('n*-*.jsonl'):
            rows=[json.loads(s) for s in path.read_text().splitlines()][1:]
            for rep in [0,2,4,6]:
                a=[r['variant'] for r in rows if r['rep']==rep];b=[r['variant'] for r in rows if r['rep']==rep+1]
                self.assertEqual(a[::-1],b);self.assertEqual(sorted(a),sorted(p['variants']))
    def test_analysis_correction_is_separate_and_hash_bound(self):
        source=HERE/'analysis_01';meta=json.loads((source/'metadata.json').read_text())
        self.assertEqual(meta['input_manifest_sha256'],analysis.sha(self.root/'manifest.json'))
        self.assertEqual(meta['original_analyzer_sha256'],analysis.sha(self.root/'analyze.py'))
        self.assertEqual(meta['corrected_analyzer_sha256'],analysis.sha(source/'corrected_analysis.py'))
        for name,digest in json.loads((source/'manifest.json').read_text())['files'].items():self.assertEqual(analysis.sha(source/name),digest)
        self.assertFalse((self.root/'results.json').exists())
        self.assertFalse(json.loads((self.root/'EXECUTION_STATUS.json').read_text())['analysis_complete'])
    def test_original_analysis_failure_is_preserved(self):
        spec=importlib.util.spec_from_file_location('original_failed_analysis',self.root/'analyze.py')
        original=importlib.util.module_from_spec(spec);spec.loader.exec_module(original)
        with self.assertRaises(AssertionError):original.main(self.root)
        self.assertEqual(original.adjusted_scan(100,180,160,10),analysis.adjusted_scan(100,180,160,10))

if __name__=='__main__':unittest.main()
