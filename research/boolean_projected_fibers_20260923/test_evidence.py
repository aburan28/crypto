import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
RUN = HERE/'run_01'
spec = importlib.util.spec_from_file_location('projected_analysis', HERE/'analyze.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


class EvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = a.read(RUN/'results.json')
        cls.first = 'n12-discovery-17-planted'
        cls.rows = [json.loads(x) for x in (RUN/(cls.first+'.jsonl')).read_text().splitlines()]
        cls.sample = next(s for s in cls.rows[1:] if s['variant']=='projected4')

    def overlay(self):
        tmp = tempfile.TemporaryDirectory(dir=HERE, prefix='projected-scratch-')
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        for path in RUN.iterdir():
            if path.is_file():
                (root/path.name).symlink_to(path)
        return root

    def replace(self, root, name, data):
        path = root/name
        self.assertTrue(path.is_symlink())
        path.unlink()
        path.write_text(data)

    def test_manifest_sources_and_immutable_predecessor_lineage(self):
        manifest = a.read(RUN/'manifest.json')
        self.assertEqual(set(manifest['files']), {p.name for p in RUN.iterdir() if p.is_file() and p.name!='manifest.json'})
        for name,digest in manifest['files'].items():
            self.assertEqual(a.sha(RUN/name),digest,name)
        for name in a.read(RUN/'metadata.json')['source_hashes']:
            self.assertEqual(a.sha(HERE/name),a.sha(RUN/name),name)
        lineage=a.read(HERE/'SOURCE_LINEAGE.json')
        predecessor=HERE/lineage['predecessor']
        for name,digest in lineage['files'].items():
            self.assertEqual(a.sha(predecessor/name),digest,name)

    def test_exact_complete_replay(self):
        self.assertEqual(a.validate(RUN),self.result)
        self.assertEqual(self.result['cells'],216)
        self.assertEqual(self.result['samples'],89856)

    def test_independent_generator_covers_every_prior_fixture(self):
        for r in a.read(HERE/'REFERENCE_FIXTURES.json'):
            polys,witness=a.fixture(r['n'],r['seed'],r['family'])
            self.assertEqual((polys,witness),(r['polys'],r['planted_witness']))

    def test_mutated_timed_source_is_rejected(self):
        root=self.overlay()
        self.replace(root,'projected.rs',(RUN/'projected.rs').read_text()+'\n// mutation\n')
        with self.assertRaises(AssertionError):a.validate(root)

    def test_missing_worker_receipt_is_rejected(self):
        root=self.overlay()
        receipts=a.read(RUN/'receipts.json')
        self.replace(root,'receipts.json',json.dumps(receipts[:-1]))
        with self.assertRaises(AssertionError):a.validate(root)

    def test_worker_failure_is_rejected(self):
        root=self.overlay()
        receipts=a.read(RUN/'receipts.json')
        receipts[0]['exit_code']=1
        self.replace(root,'receipts.json',json.dumps(receipts))
        with self.assertRaises(AssertionError):a.validate(root)

    def test_altered_raw_stream_is_rejected(self):
        root=self.overlay()
        self.replace(root,self.first+'.jsonl',(RUN/(self.first+'.jsonl')).read_text()+'\n')
        with self.assertRaises(AssertionError):a.validate(root)

    def test_false_verified_original_model_is_rejected_after_receipt_rebind(self):
        root=self.overlay()
        rows=copy.deepcopy(self.rows)
        wrong=next(x for x in range(1<<12) if not a.satisfies(rows[0]['polys'],x))
        next(s for s in rows[1:] if s['variant']=='projected4')['model']=wrong
        self.replace(root,self.first+'.jsonl',''.join(json.dumps(r)+'\n' for r in rows))
        receipts=a.read(RUN/'receipts.json')
        receipts[0]['stdout_sha256']=a.sha(root/(self.first+'.jsonl'))
        self.replace(root,'receipts.json',json.dumps(receipts))
        with self.assertRaises(AssertionError):a.validate(root)

    def test_screen_only_consistency_cannot_count_as_verified_model(self):
        polys=[[3,0]]
        self.assertFalse(a.satisfies(polys,0))
        self.assertTrue(a.satisfies(polys,3))
        self.assertEqual(a.quotient_rank(polys,2),1)

    def test_rank_count_and_free_space_corruption_are_rejected(self):
        original=self.sample['projected']
        for key in ['low_variables','annihilated_rank','quotient_dimension','extension_space','original_rejected','prefixes']:
            work=copy.deepcopy(original);work[key]+=1
            with self.assertRaises(AssertionError,msg=key):a.check_work(work,self.rows[0]['polys'],12,4,'SAT')

    def test_unsat_requires_all_prefixes_and_extensions(self):
        work=copy.deepcopy(self.sample['projected'])
        with self.assertRaises(AssertionError):a.check_work(work,self.rows[0]['polys'],12,4,'UNSAT')

    def test_cap_cannot_promote_as_complete(self):
        root=self.overlay()
        rows=copy.deepcopy(self.rows)
        for sample in rows[1:]:
            if sample['variant']=='projected4':
                sample.update(outcome='UNKNOWN',model=None,reason='PROJECTED_EXTENSION_CAP',verified=False)
                sample['projected']['extensions_checked']=sample['projected']['original_rejected']
        self.replace(root,self.first+'.jsonl',''.join(json.dumps(r)+'\n' for r in rows))
        receipts=a.read(RUN/'receipts.json')
        receipts[0]['stdout_sha256']=a.sha(root/(self.first+'.jsonl'))
        self.replace(root,'receipts.json',json.dumps(receipts))
        result=a.validate(root)
        self.assertFalse(result['all_complete'])
        self.assertFalse(any(d['dramatic_pass'] or d['incremental_pass'] for d in result['decisions'].values()))

    def test_orders_are_exact_reverse_pairs(self):
        for count in [3,49,52]:
            for rep in range(0,8,2):
                order=a.paired_order(count,19437,rep)
                self.assertEqual(sorted(order),list(range(count)))
                self.assertEqual(order,a.paired_order(count,19437,rep+1)[::-1])

    def test_missing_costs_remain_null(self):
        for key in ['production_solver_cost','full_ic_cost','rho_ratio','calibrated_operation_ratio']:
            self.assertIsNone(self.result[key])
            self.assertIsNone(a.read(HERE/'protocol.json')[key])

    def test_published_report_binds_the_frozen_evidence(self):
        ledger=a.read(HERE/'RUN_LEDGER.json')
        for key,path in [('run_manifest_sha256',RUN/'manifest.json'),('results_sha256',RUN/'results.json'),
                         ('report_sha256',HERE/'report.py'),('summary_sha256',HERE/'SUMMARY.json'),
                         ('conclusion_sha256',HERE/'CONCLUSION.md')]:
            self.assertEqual(ledger[key],a.sha(path))
        summary=a.read(HERE/'SUMMARY.json')
        self.assertEqual(summary['decisions'],self.result['decisions'])
        self.assertEqual(summary['n24_holdout_ms'],self.result['n24_holdout_ms'])
        page=(HERE.parent.parent/'docs/index-calculus-scoreboard.html').read_text()
        marker='<section class="panel" id="boolean-projected-fibers-20260923">'
        self.assertEqual(page.count(marker),1)
        section=page.split(marker,1)[1].split('</section>',1)[0]
        for row in summary['n24_holdout_ms']:
            self.assertIn('<td>'+row['variant']+'</td>',section)


if __name__ == '__main__':unittest.main()
