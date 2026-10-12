import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

HERE=Path(__file__).resolve().parent
RUN=HERE/'run_01'
spec=importlib.util.spec_from_file_location('minor_analysis',HERE/'analyze.py')
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)


class EvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=a.read(RUN/'results.json')
        cls.cell='n12-17-planted'
        cls.rows=[json.loads(line) for line in (RUN/(cls.cell+'.jsonl')).read_text().splitlines()]
        cls.sample=next(s for s in cls.rows[1:] if s['variant']=='symbolic_minor4')

    def overlay(self):
        tmp=tempfile.TemporaryDirectory(dir=HERE,prefix='minor-scratch-');self.addCleanup(tmp.cleanup)
        root=Path(tmp.name)
        for path in RUN.iterdir():
            if path.is_file():(root/path.name).symlink_to(path)
        return root

    def replace(self,root,name,data):
        path=root/name;self.assertTrue(path.is_symlink());path.unlink();path.write_text(data)

    def replace_rows(self,root,rows):
        self.replace(root,self.cell+'.jsonl',''.join(json.dumps(r)+'\n' for r in rows))
        receipts=a.read(RUN/'receipts.json');receipts[0]['stdout_sha256']=a.sha(root/(self.cell+'.jsonl'))
        self.replace(root,'receipts.json',json.dumps(receipts))

    def test_manifest_sources_and_lineage_are_exact(self):
        manifest=a.read(RUN/'manifest.json')
        self.assertEqual(set(manifest['files']),{p.name for p in RUN.iterdir() if p.is_file() and p.name!='manifest.json'})
        for name,digest in manifest['files'].items():self.assertEqual(a.sha(RUN/name),digest,name)
        for name in a.read(RUN/'metadata.json')['source_hashes']:self.assertEqual(a.sha(HERE/name),a.sha(RUN/name),name)
        lineage=a.read(HERE/'SOURCE_LINEAGE.json')
        for name,digest in lineage['files'].items():self.assertEqual(a.sha(HERE/lineage['predecessor']/name),digest)
        self.assertEqual(a.sha(HERE/lineage['reference_python']['path']),lineage['reference_python']['sha256'])

    def test_independent_exact_replay(self):
        self.assertEqual(a.validate(RUN),self.result)
        self.assertEqual((self.result['cells'],self.result['observations']),(12,5280))
        self.assertFalse(self.result['complete_candidate_implemented'])

    def test_truth_transform_matches_direct_evaluation_and_inverts(self):
        for n in range(1,5):
            for code in range(min(1<<(1<<n),256)):
                result=a.mobius(code,n)
                expected=sum((sum((code>>m)&1 for m in range(1<<n) if x&m==m)%2)<<x for x in range(1<<n))
                self.assertEqual(result,expected)
                self.assertEqual(a.mobius(result,n),code)

    def test_rank_deficient_minors_do_not_establish_consistency(self):
        self.assertEqual(a.determinant([16,0,0,0,0]),0)
        self.assertFalse(a.consistent([0,0,0,0],1))
        self.assertTrue(a.consistent([0,0,0,0],0))

    def test_overlap_zero_and_duplicate_minors_are_retained(self):
        data={'accepted':[1],'minors':[[10],[10],[12],[0]]}
        stats=a.filter_statistics(data,2)
        self.assertEqual(stats['individual_rejected'],[2,2,2,0])
        self.assertEqual(stats['joint_rejected'],3)
        self.assertEqual(stats['accepted_prefixes'],1)
        self.assertEqual((stats['zero_minors'],stats['duplicate_minors']),(1,1))
        self.assertEqual(stats['overlap'][0][1],2)

    def test_corrupted_compiler_source_is_rejected(self):
        root=self.overlay();self.replace(root,'minors.rs',(RUN/'minors.rs').read_text()+'\n// corrupt\n')
        with self.assertRaises(AssertionError):a.validate(root)

    def test_missing_receipt_is_rejected(self):
        root=self.overlay();receipts=a.read(RUN/'receipts.json');self.replace(root,'receipts.json',json.dumps(receipts[:-1]))
        with self.assertRaises(AssertionError):a.validate(root)

    def test_worker_failure_is_rejected(self):
        root=self.overlay();receipts=a.read(RUN/'receipts.json');receipts[0]['exit_code']=1
        self.replace(root,'receipts.json',json.dumps(receipts))
        with self.assertRaises(AssertionError):a.validate(root)

    def test_altered_mask_is_rejected_after_receipt_rebind(self):
        root=self.overlay();rows=copy.deepcopy(self.rows)
        sample=next(s for s in rows[1:] if s['variant']=='symbolic_minor4')
        sample['filter']['accepted'][0]^=1
        self.replace_rows(root,rows)
        with self.assertRaises(AssertionError):a.validate(root)

    def test_false_original_solution_count_is_rejected(self):
        root=self.overlay();rows=copy.deepcopy(self.rows);rows[0]['original_solutions']+=1
        self.replace_rows(root,rows)
        with self.assertRaises(AssertionError):a.validate(root)

    def test_false_complete_solver_model_is_rejected(self):
        root=self.overlay();rows=copy.deepcopy(self.rows)
        bad=next(x for x in range(1<<12) if not a.ref.satisfies(rows[0]['polys'],x))
        next(s for s in rows[1:] if s['variant']=='word_dispatch')['model']=bad
        self.replace_rows(root,rows)
        with self.assertRaises(AssertionError):a.validate(root)

    def test_product_and_degree_corruptions_are_rejected(self):
        exact=a.oracle(tuple(tuple(p) for p in self.rows[0]['polys']),12)
        stats=a.filter_statistics(exact['oracle_minor'],8)
        for key in ['product_toggles','cancelled_toggles','support_sum','dp_states','truth_word_xors','dp_payload_peak_bytes']:
            work=copy.deepcopy(self.sample['work']);work[key]+=1
            with self.assertRaises(AssertionError,msg=key):a.check_work(work,'symbolic_minor4',exact,8,stats,True)
        work=copy.deepcopy(self.sample['work']);work['degrees'][0]=7
        with self.assertRaises(AssertionError):a.check_work(work,'symbolic_minor4',exact,8,stats,True)

    def test_censoring_never_creates_an_opportunity(self):
        root=self.overlay();rows=copy.deepcopy(self.rows)
        for sample in rows[1:]:
            if sample['variant']=='symbolic_minor4':sample.update(outcome='CAPPED',reason='PRODUCT_CAP',verified=False,filter=None)
        self.replace_rows(root,rows)
        result=a.validate(root)
        self.assertFalse(result['all_complete'])
        self.assertFalse(result['complete_candidate_justified'])
        self.assertTrue(all(g['decision']=='INCONCLUSIVE' for g in result['gates']))

    def test_missing_costs_remain_unmeasured(self):
        for key in ['production_solver_cost','full_ic_cost','rho_ratio','calibrated_operation_ratio']:
            self.assertIsNone(self.result[key]);self.assertIsNone(a.read(HERE/'protocol.json')[key])

    def test_report_and_scoreboard_bind_frozen_evidence(self):
        ledger=a.read(HERE/'RUN_LEDGER.json')
        for key,path in [('manifest_sha256',RUN/'manifest.json'),('results_sha256',RUN/'results.json'),
                         ('report_sha256',HERE/'report.py'),('summary_sha256',HERE/'SUMMARY.json'),('conclusion_sha256',HERE/'CONCLUSION.md')]:
            self.assertEqual(ledger[key],a.sha(path))
        summary=a.read(HERE/'SUMMARY.json')
        self.assertEqual(summary['gates'],self.result['gates'])
        self.assertEqual(summary['n16_ms'],self.result['n16_ms'])
        page=(HERE.parent.parent/'docs/index-calculus-scoreboard.html').read_text()
        marker='<section class="panel" id="boolean-implicit-minors-20260923">'
        self.assertEqual(page.count(marker),1)
        section=page.split(marker,1)[1].split('</section>',1)[0]
        for row in summary['n16_ms']:self.assertIn('<td>'+row['variant']+'</td>',section)


if __name__=='__main__':unittest.main()
