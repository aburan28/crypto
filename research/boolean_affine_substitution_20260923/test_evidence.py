import contextlib
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock

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
        predecessor=HERE.parent/'boolean_basis_transport_20260923'
        current=json.loads((HERE/'run_01/protocol.json').read_text());prior=json.loads((predecessor/'run_01/protocol.json').read_text())
        self.assertEqual(current['regression_seeds'],prior['regression_seeds']+prior['holdout_seeds'])
        self.assertEqual(current['discovery_seeds'],prior['discovery_seeds'])
        for dirname in ['boolean_basis_transport_20260923','boolean_quadratic_frontier_20260923','boolean_quadratic_state_20260923','boolean_solve_trace_20260923']:
            for path in (HERE.parent/dirname).glob('run_*/protocol.json'):
                earlier=json.loads(path.read_text());self.assertFalse(set(current['holdout_seeds'])&set(earlier['holdout_seeds']))
        lineage=json.loads((HERE/'SOURCE_LINEAGE.json').read_text());self.assertEqual(current['predecessor_source_hashes'],lineage['source_files'])
        for name,digest in lineage['source_files'].items():self.assertEqual(analysis.sha(HERE.parents[1]/name),digest)
        for name in ['kernel.rs','quadratic.rs','packed.rs','basis.rs','basis_tests.rs']:
            self.assertEqual((HERE/name).read_bytes(),(predecessor/name).read_bytes())
        for name in json.loads((HERE/'run_01/metadata.json').read_text())['source_hashes']:
            self.assertEqual((HERE/name).read_bytes(),(HERE/'run_01'/name).read_bytes())

    def test_all_discovery_probes_preserve_their_own_sources_and_results(self):
        for index,arms in [(1,15),(2,16),(3,18),(4,20),(5,24),(6,30),(7,32)]:
            root=HERE/f'resource_probe_{index:02d}'
            for name,digest in json.loads((root/'manifest.json').read_text())['files'].items():self.assertEqual(analysis.sha(root/name),digest)
            receipt=json.loads((root/'receipt.json').read_text());self.assertTrue(receipt['complete'])
            self.assertTrue(all(r['exit_code']==0 for r in receipt['receipts']))
            for path in root.glob('*.jsonl'):
                rows=[json.loads(s) for s in path.read_text().splitlines()];self.assertEqual(len(rows),arms+1)
                fixture=rows[0];self.assertEqual((fixture['n'],fixture['seed']),(24,17))
                for row in rows[1:]:
                    self.assertTrue(row['verified'])
                    if row['outcome']=='SAT':self.assertTrue(analysis.satisfies(fixture['polys'],row['model']))
                    else:self.assertEqual((row['outcome'],fixture['search_reference']),('UNSAT','UNSAT'))

    def test_native_sat_probe_is_separate_and_its_models_verify(self):
        root=HERE/'native_resource_probe_01'
        for name,digest in json.loads((root/'manifest.json').read_text())['files'].items():self.assertEqual(analysis.sha(root/name),digest)
        metadata=json.loads((root/'metadata.json').read_text());self.assertTrue(metadata['complete']);self.assertEqual(metadata['threads'],1)
        for family,digest in metadata['source_fixtures'].items():self.assertEqual(analysis.sha(HERE/'resource_probe_04'/(family+'.jsonl')),digest)
        for row in json.loads((root/'results.json').read_text()):
            fixture=json.loads((HERE/'resource_probe_04'/(row['family']+'.jsonl')).read_text().splitlines()[0])
            if row['status']=='SAT':self.assertTrue(analysis.satisfies(fixture['polys'],row['model']))
            elif row['status']=='UNSAT':self.assertEqual(fixture['search_reference'],'UNSAT')
            else:self.assertIsNone(row['completion_ns'])

    def test_untraced_reference_cannot_be_replaced_by_a_slower_control(self):
        # Disposable synthetic timing mutation: a slow traced reference looks
        # favorable, but the mandatory untraced control beats the candidate.
        path=self.root/'raw-n16.jsonl';rows=[json.loads(s) for s in path.read_text().splitlines()]
        for row in rows:
            if row['type']!='sample':continue
            total={'packed_untraced':100,'gray_simd':200}.get(row['variant'],1000)
            row.update(total_ns=total,solve_ns=total*9//10,validation_ns=total//10,kernel_ns=total//4 if row['kernel_calls'] else 0)
            for field in ['substitution_ns','enumeration_ns']:
                if row[field] is not None:row[field]=total//4
            row['raw_line']=json.dumps({k:v for k,v in row.items() if k not in ['cell','split','raw_line']})+'\n'
        path.write_text(''.join(json.dumps(r)+'\n' for r in rows))
        path=self.root/'receipts.json';receipts=json.loads(path.read_text())
        for receipt in receipts:
            if receipt['n']==16:receipt['stdout_sha256']=hashlib.sha256(''.join(r['raw_line'] for r in rows if r['cell']==receipt['cell']).encode()).hexdigest()
        path.write_text(json.dumps(receipts))
        with contextlib.redirect_stdout(io.StringIO()):analysis.main(self.root)
        result=json.loads((self.root/'results.json').read_text())
        for gate in result['new_details']:
            if gate['candidate']=='gray_simd' and gate['reference']=='retained_frontier':
                self.assertIn('packed_untraced',gate['controls'])
                self.assertIn('gray_scalar',gate['controls'])
                if gate['n']==16:
                    self.assertEqual(gate['paired_ratio_median'],0.5);self.assertEqual(gate['ci95_paired_median'],[0.5,0.5]);self.assertFalse(gate['pass'])
        self.assertEqual(result['new_gates']['gray_simd']['retained_frontier'],'REJECTED')

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
        keys=['outcome','model','reason','nodes','kernel_calls','decisions','forced','specialized_terms','source_rows','source_columns','max_depth','affine_eliminated','derived_rows','enumeration_points','enumeration_batches','enumeration_leaves','trace']
        for cell in result['cells_detail']:
            arms=cell['arms']
            for left,right in [('basis_list','basis_wide'),('tail_list','tail_wide'),('affine_sl_basis_list','affine_sl_basis_fast'),('gray_scalar','gray_simd'),('packed_gray12_scalar','packed_gray12_simd'),('packed_gray16_scalar','packed_gray16_simd')]:
                self.assertEqual({k:arms[left][k] for k in keys},{k:arms[right][k] for k in keys})
                if right in ['basis_wide','tail_wide','affine_sl_basis_fast']:self.assertGreater(arms[right]['kernel_calls'],0)
            self.assertEqual(arms['packed_state']['kernel_calls'],0)
        self.assertIsNone(result['production_solver_cost']);self.assertIsNone(result['full_ic_cost']);self.assertIsNone(result['rho_ratio'])

    def test_tail_trace_difference_is_rejected(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='tail_wide',lambda r:r.update(trace=r['trace']^1));self.rejected()

    def test_enumeration_checksum_difference_is_rejected(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='gray_simd',lambda r:r.update(trace=r['trace']^1));self.rejected()

    def test_impossible_enumeration_count_is_rejected(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='gray_simd',lambda r:r.update(enumeration_points=(1<<24)+1));self.rejected()

    def test_missing_phase_cost_is_not_zero_filled(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='packed_state',lambda r:r.update(substitution_ns=0));self.rejected()

    def test_enumeration_phase_cannot_exceed_complete_solve(self):
        self.mutate(lambda r:r['type']=='sample' and r['variant']=='gray_simd',lambda r:r.update(enumeration_ns=r['solve_ns']+1));self.rejected()

    def test_enumeration_censor_retains_observed_work_but_blocks_every_new_gate(self):
        path=self.root/'raw-n24.jsonl';rows=[json.loads(line) for line in path.read_text().splitlines()]
        cell=next(r['cell'] for r in rows if r['type']=='fixture' and r['split']=='holdout' and r['family']=='planted')
        for row in rows:
            if row['cell']==cell and row['type']=='sample' and row['variant'] in ['gray_scalar','gray_simd']:
                row.update(outcome='UNKNOWN',model=None,reason='ENUM_CAP',verified=False)
                row['raw_line']=json.dumps({k:v for k,v in row.items() if k not in ['cell','split','raw_line']})+'\n'
        path.write_text(''.join(json.dumps(r)+'\n' for r in rows))
        path=self.root/'receipts.json';receipts=json.loads(path.read_text())
        next(r for r in receipts if r['cell']==cell)['stdout_sha256']=hashlib.sha256(''.join(r['raw_line'] for r in rows if r['cell']==cell).encode()).hexdigest()
        path.write_text(json.dumps(receipts))
        with contextlib.redirect_stdout(io.StringIO()):analysis.main(self.root)
        result=json.loads((self.root/'results.json').read_text())
        self.assertFalse(result['all_results_completed_and_verified'])
        for verdict in result['new_gates'].values():self.assertEqual(verdict['retained_frontier'],'REJECTED')
        affected=next(c for c in result['cells_detail'] if c['cell']==cell)
        for arm in ['gray_scalar','gray_simd']:
            self.assertIsNone(affected['arms'][arm]['completion_ns'])
            self.assertGreater(affected['arms'][arm]['observed_total_ns'],0)
            self.assertGreater(affected['arms'][arm]['enumeration_points'],0)

    def test_untraced_control_preserves_work_and_has_no_trace_claim(self):
        result=json.loads((self.root/'results.json').read_text())
        keys=['outcome','model','nodes','decisions','forced','specialized_terms','max_depth']
        for cell in result['cells_detail']:
            a=cell['arms']['packed_state'];b=cell['arms']['packed_untraced']
            self.assertEqual({k:a[k] for k in keys},{k:b[k] for k in keys});self.assertIsNone(b['trace'])

    def test_primary_rejection_prevents_confirmation_side_effects(self):
        spec = importlib.util.spec_from_file_location('confirmation_guard', HERE / 'confirm.py')
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        self.root.rename(self.root.with_name('run_01'))
        with mock.patch.object(module, 'HERE', Path(self.temp.name)), mock.patch.object(module.subprocess, 'run') as launch:
            with self.assertRaisesRegex(AssertionError, 'No primary treatment passed'):
                module.main()
            launch.assert_not_called()
        self.assertFalse((Path(self.temp.name) / 'confirmation_protocol.json').exists())
        self.assertFalse((Path(self.temp.name) / 'run_02').exists())

    def test_guard_receipt_binds_primary_and_unstarted_confirmation(self):
        receipt = json.loads((HERE / 'CONFIRMATION_GUARD.json').read_text())
        self.assertEqual(receipt['primary_manifest_sha256'], analysis.sha(self.root / 'manifest.json'))
        self.assertEqual(receipt['primary_results_sha256'], analysis.sha(self.root / 'results.json'))
        self.assertEqual(receipt['launcher_sha256'], analysis.sha(HERE / 'confirm.py'))
        self.assertEqual(receipt['exit_code'], 1)
        self.assertFalse(receipt['confirmation_started'])
        self.assertFalse(receipt['execution_protocol_exists'])
        self.assertFalse(receipt['run_02_exists'])
        self.assertFalse((HERE / 'confirmation_protocol.json').exists())
        self.assertFalse((HERE / 'run_02').exists())

    def test_repeating_primary_cannot_masquerade_as_confirmation(self):
        spec = importlib.util.spec_from_file_location('run_summary', HERE / 'summarize.py')
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        module.validate_run_names(['run_01'])
        for names in [['run_01', 'run_01'], ['run_02'], ['run_01', 'arbitrary_run'], []]:
            with self.assertRaises(AssertionError):
                module.validate_run_names(names)

if __name__=='__main__':unittest.main()
