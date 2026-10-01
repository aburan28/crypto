import contextlib
import copy
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

HERE=Path(__file__).resolve().parent
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
a=module('half_analysis',HERE/'analyze.py')
runner=module('half_isolated_runner',HERE/'isolated_run.py')


class EvidenceTests(unittest.TestCase):
    def temporary(self):
        tmp=tempfile.TemporaryDirectory(dir=HERE,prefix='half-scratch-');self.addCleanup(tmp.cleanup);return Path(tmp.name)

    def test_frozen_probes_and_sources_are_intact(self):
        for name in ['probe_01','probe_02']:
            root=HERE/name;manifest=a.read(root/'manifest.json')
            self.assertEqual(set(manifest['files']),{p.name for p in root.iterdir() if p.is_file() and p.name!='manifest.json'})
            for path,digest in manifest['files'].items():self.assertEqual(a.sha(root/path),digest,path)

    def test_historical_probe_replay_is_exact(self):
        for name in ['probe_01','probe_02']:
            self.assertEqual(a.validate(HERE/name),a.read(HERE/name/'results.json'))

    def test_retained_sources_and_baseline_are_exact(self):
        lineage=a.read(HERE/'SOURCE_LINEAGE.json')
        for name,digest in lineage['files'].items():
            data=(HERE/name).read_bytes()
            if name=='worker.rs':
                suffix=lineage['worker_suffix'].encode();self.assertTrue(data.endswith(suffix));data=data[:-len(suffix)]
            self.assertEqual(hashlib.sha256(data).hexdigest(),digest,name)
        self.assertEqual(a.sha(HERE/'baseline_worker.rs'),lineage['files']['worker.rs'])
        self.assertEqual(a.sha(HERE/'reference.py'),lineage['reference_python']['sha256'])

    def test_full_protocol_keeps_all_previous_inputs_and_unused_holdouts(self):
        p=a.read(HERE/'protocol_full.json');prior=a.read(HERE/'REFERENCE_PROTOCOL.json')
        self.assertEqual(p['regression_seeds'],prior['regression_seeds']+prior['holdout_seeds'])
        self.assertFalse(set(p['holdout_seeds'])&set(p['regression_seeds']+p['discovery_seeds']))
        self.assertEqual(sum(len(p[s+'_seeds']) for s in p['splits'])*12,240)
        self.assertEqual(len(p['variants']),61)
        for old in a.read(HERE/'REFERENCE_FIXTURES.json')['fixtures']:
            polys,witness=a.ref.fixture(old['n'],old['seed'],old['family'])
            self.assertEqual((polys,witness),(old['polys'],old['planted_witness']))

    def valid_conditions(self,root):
        cell='n12-discovery-17-planted';name=cell+'-ab-conditions.jsonl'
        state={'loadavg':[0,0,0],'psi_cpu':{'some':{'avg10':0.0}},'psi_memory':{'some':{'avg10':0.0}}}
        worker=['/retained/worker','12','17','planted','8','200000','1','ab']
        record={'schema':'isolated-bench/1','mode':'run','reserved_cpus':[2],
            'host':{'logical_cpus':4,'machine':'aarch64','kernel':'test-kernel'},
            'preflight':{'settle':{'seconds':2,'other_cpu_seconds':0.0},'conditions':copy.deepcopy(state)},
            'run':{'cpus':[2],'command':worker,'exit_status':0,'contended':False,'wall_seconds':1.0,'other_cpu_seconds':0.0,
                   'before':copy.deepcopy(state),'after':copy.deepcopy(state),
                   'voluntary_switches':0,'involuntary_switches':0,'minor_faults':1,'major_faults':0,'max_rss_kib':1024}}
        (root/name).write_text(json.dumps(record)+'\n')
        stage={'conditions_file':name,'conditions_sha256':a.sha(root/name),'qualified':True,'exit_code':0,'timed_out':False,
               'worker_command':worker,'worker_peak_rss_bytes':1048576,
               'command':['python','isolated_bench.py','run','--settle','2','--max-other-cpu','0.10','--max-psi','5.0']}
        return cell,stage,record

    def test_complete_isolation_receipt_is_required(self):
        root=self.temporary();cell,stage,_=self.valid_conditions(root)
        a.verify_isolation(root,stage,[2],cell,'ab')
        self.assertEqual(runner.qualified_conditions(root/stage['conditions_file'],[2])[0],True)

    def test_contended_result_is_rejected_even_with_valid_masks(self):
        root=self.temporary();cell,stage,record=self.valid_conditions(root)
        record['run']['contended']=True
        path=root/stage['conditions_file'];path.write_text(json.dumps(record)+'\n');stage['conditions_sha256']=a.sha(path)
        with self.assertRaises(AssertionError):a.verify_isolation(root,stage,[2],cell,'ab')
        self.assertFalse(runner.qualified_conditions(path,[2])[0])
        self.assertFalse(runner.qualified_conditions(path,[2])[0])

    def test_false_uncontended_flag_cannot_hide_other_cpu_use(self):
        root=self.temporary();cell,stage,record=self.valid_conditions(root)
        record['run']['other_cpu_seconds']=.2
        path=root/stage['conditions_file'];path.write_text(json.dumps(record)+'\n');stage['conditions_sha256']=a.sha(path)
        with self.assertRaises(AssertionError):a.verify_isolation(root,stage,[2],cell,'ab')

    def test_missing_pressure_or_wrong_cpu_is_rejected(self):
        for mutation in ['psi','cpus','host','switches','preflight']:
            root=self.temporary();cell,stage,record=self.valid_conditions(root)
            if mutation=='psi':record['preflight']['conditions']['psi_cpu']=None
            if mutation=='cpus':record['run']['cpus']=[1]
            if mutation=='host':record['host']['machine']=''
            if mutation=='switches':record['run']['involuntary_switches']=-1
            if mutation=='preflight':record['preflight']['settle']['other_cpu_seconds']=.3
            path=root/stage['conditions_file'];path.write_text(json.dumps(record)+'\n');stage['conditions_sha256']=a.sha(path)
            with self.assertRaises((AssertionError,TypeError),msg=mutation):a.verify_isolation(root,stage,[2],cell,'ab')

    def test_missing_record_and_timeout_are_not_qualified(self):
        root=self.temporary();self.assertEqual(runner.qualified_conditions(root/'missing.jsonl',[2]),(False,None))
        cell,stage,_=self.valid_conditions(root);stage['timed_out']=True
        with self.assertRaises(AssertionError):a.verify_isolation(root,stage,[2],cell,'ab')

    def test_unsupported_platform_refuses_before_creating_run(self):
        root=self.temporary();destination=root/'must-not-exist'
        with patch.object(runner.platform,'system',return_value='Darwin'),patch.object(runner.sys,'argv',['isolated_run.py','--phase','discovery','--out',str(destination)]),contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as caught:runner.main()
            self.assertEqual(caught.exception.code,2);self.assertFalse(destination.exists())

    def quiet_stub(self,outcomes):
        calls=[]
        def preflight(args,exclude):
            calls.append((args.settle,args.max_other_cpu,args.max_psi,exclude))
            value=outcomes[len(calls)-1]
            if isinstance(value,BaseException):raise value
            return value
        tool=SimpleNamespace(preflight=preflight)
        spec=SimpleNamespace(loader=SimpleNamespace(exec_module=lambda module:None))
        return calls,patch.object(runner.importlib.util,'spec_from_file_location',return_value=spec),patch.object(runner.importlib.util,'module_from_spec',return_value=tool)

    def test_quiet_wait_retains_refusals_before_first_accepted_observation(self):
        root=self.temporary();_,_,condition=self.valid_conditions(root)
        calls,spec,module=self.quiet_stub([SystemExit('pressure too high'),condition['preflight']])
        path=root/'quiet.json'
        with spec,module:runner.wait_for_quiet(root/'tool.py',path)
        records=a.read(path)
        self.assertEqual([r['accepted'] for r in records],[False,True])
        self.assertEqual([c[:3] for c in calls],[(2.0,.10,5.0)]*2)
        contract=a.read(HERE/'protocol_discovery.json')['quiet_wait']
        a.verify_quiet_wait(path,contract)
        records[0]['accepted']=True;a.dump(path,records)
        with self.assertRaises(AssertionError):a.verify_quiet_wait(path,contract)

    def test_quiet_wait_is_bounded_and_failed_wait_cannot_be_admitted(self):
        root=self.temporary();calls,spec,module=self.quiet_stub([SystemExit('busy')]*30)
        path=root/'quiet.json'
        with spec,module,self.assertRaisesRegex(RuntimeError,'no timed worker was launched'):
            runner.wait_for_quiet(root/'tool.py',path)
        self.assertEqual(len(calls),30)
        self.assertEqual(len(a.read(path)),30)
        with self.assertRaises(AssertionError):a.verify_quiet_wait(path,a.read(HERE/'protocol_discovery.json')['quiet_wait'])

    def test_quiet_wait_cannot_override_pressure_or_cpu_threshold(self):
        root=self.temporary();_,_,condition=self.valid_conditions(root)
        contract=a.read(HERE/'protocol_discovery.json')['quiet_wait']
        for mutation in ['pressure','cpu']:
            pre=copy.deepcopy(condition['preflight'])
            if mutation=='pressure':pre['conditions']['psi_cpu']['some']['avg10']=5.01
            else:pre['settle']['other_cpu_seconds']=.21
            path=root/'quiet.json'
            a.dump(path,[{'attempt':1,'recorded_utc':'test','accepted':True,'preflight':pre}])
            with self.assertRaises(AssertionError):a.verify_quiet_wait(path,contract)

    def test_changed_projection_counters_are_rejected(self):
        root=HERE/'probe_01';rows=[json.loads(line) for line in (root/'n12-discovery-17-planted.jsonl').read_text().splitlines()]
        sample=next(s for s in rows[1:] if s['variant']=='half16_native')
        control=next(s for s in rows[1:] if s['variant']=='word16_unrolled' and s['rep']==sample['rep'])
        a.check_half(sample,control,12,'half16_native')
        for key in ['projected_hits','original_checks','rejected_hits','verified_hits','points','block_width']:
            bad=copy.deepcopy(sample);bad['half_work'][key]+=1
            with self.assertRaises(AssertionError,msg=key):a.check_half(bad,control,12,'half16_native')

    def test_model_and_logical_work_must_match_full_word_control(self):
        rows=[json.loads(line) for line in (HERE/'probe_01/n12-discovery-17-planted.jsonl').read_text().splitlines()]
        sample=next(s for s in rows[1:] if s['variant']=='half64_native')
        control=next(s for s in rows[1:] if s['variant']=='wide64_unrolled' and s['rep']==sample['rep'])
        for key in ['model','logical']:
            bad=copy.deepcopy(sample)
            if key=='model':bad['model']^=1
            else:bad['logical']['enumeration_points']+=64
            with self.assertRaises(AssertionError):a.check_half(bad,control,12,'half64_native')

    def test_missing_cryptanalytic_costs_stay_null(self):
        for name in ['protocol_discovery.json','protocol_full.json']:
            for field in ['full_ic_cost','production_solver_cost','rho_ratio','calibrated_operation_ratio']:
                self.assertIsNone(a.read(HERE/name)[field])

    def test_frozen_reference_membership_does_not_change(self):
        p=a.read(HERE/'protocol_full.json');old=a.read(HERE/'REFERENCE_PROTOCOL.json')
        self.assertEqual(p['reference_arms'],old['variants'])
        self.assertEqual(p['matched_controls']['half16_eor3'],'eor3_word16')
        self.assertEqual(p['matched_controls']['half64_eor3'],'eor3_word64')

    def aa_fixture(self):
        root=self.temporary();cell='n12-discovery-17-planted'
        rows=[json.loads(line) for line in (HERE/'probe_01'/(cell+'.jsonl')).read_text().splitlines()]
        source={**rows[0],'mode':'ab','architecture':'aarch64','eor3_available':True}
        pairs={(s['rep'],s['variant']):s for s in rows[1:]}
        seed=a.read(HERE/'probe_01/receipts.json')[0]['order_seed']
        samples=[]
        for rep in range(8):
            for order,index in enumerate(a.ref.paired_order(2,seed,rep)):
                sample=copy.deepcopy(pairs[rep,'word_dispatch']);sample.update(variant=['aa_a','aa_b'][index],order=order)
                samples.append(sample)
        path=root/(cell+'-aa.jsonl');path.write_text(''.join(json.dumps(r)+'\n' for r in [{**source,'mode':'aa'},*samples]))
        path.with_suffix('.stderr').write_text('')
        stage={'conditions_file':cell+'-aa-conditions.jsonl','stdout_sha256':a.sha(path),'stderr_sha256':a.sha(path.with_suffix('.stderr')),
               'worker_command':['/test/worker','12','17','planted','8','200000',str(seed),'aa']}
        return root,stage,source,pairs,seed,samples,path

    def test_aa_compares_identical_complete_work(self):
        root,stage,source,pairs,seed,_,_=self.aa_fixture()
        result=a.verify_aa(root,stage,source,pairs,seed)
        self.assertEqual(result['paired_ratios'],[1.0]*8)
        self.assertEqual(result['symmetric_ratios'],[1.0]*8)

    def test_paired_aa_uses_phase_slice_with_shared_resource_receipt(self):
        root,stage,source,pairs,seed,_,path=self.aa_fixture()
        stage['conditions_file']=path.name.removesuffix('-aa.jsonl')+'-paired-conditions.jsonl'
        stage['resource_mode']='paired';stage['worker_command'][-1]='paired'
        result=a.verify_aa(root,stage,source,pairs,seed)
        self.assertEqual(result['paired_ratios'],[1.0]*8)
        stage['conditions_file']=stage['conditions_file'].replace('-paired-','-ab-')
        with self.assertRaises(AssertionError):a.verify_aa(root,stage,source,pairs,seed)

    def test_aa_rejects_work_drift_after_hash_rebind(self):
        root,stage,source,pairs,seed,samples,path=self.aa_fixture()
        samples[0]['logical']['enumeration_points']+=16
        path.write_text(''.join(json.dumps(r)+'\n' for r in [{**source,'mode':'aa'},*samples]));stage['stdout_sha256']=a.sha(path)
        with self.assertRaises(AssertionError):a.verify_aa(root,stage,source,pairs,seed)

    def test_aa_rejects_missing_pairs(self):
        root,stage,source,pairs,seed,samples,path=self.aa_fixture()
        path.write_text(''.join(json.dumps(r)+'\n' for r in [{**source,'mode':'aa'},*samples[:-1]]));stage['stdout_sha256']=a.sha(path)
        with self.assertRaises(AssertionError):a.verify_aa(root,stage,source,pairs,seed)

    def test_failed_attempts_are_preserved_and_not_admitted(self):
        registry=a.read(HERE/'QUALIFIED_RUNS.json')
        for entry in registry['failed_attempts']:
            root=HERE/entry['path']
            self.assertEqual(a.sha(root/'manifest.json'),entry['manifest_sha256'])
            for name,digest in a.read(root/'manifest.json')['files'].items():self.assertEqual(a.sha(root/name),digest)
            self.assertTrue((root/'EXECUTION_STATUS.json').exists())
            self.assertNotEqual(registry['discovery'],entry)
            self.assertNotEqual(registry['full'],entry)

    def test_complete_paired_archive_replays_after_filename_fix_without_admission(self):
        replay=HERE/'analysis_replay_01'
        for name,digest in a.read(replay/'manifest.json')['files'].items():self.assertEqual(a.sha(replay/name),digest)
        receipt=a.read(replay/'receipt.json');root=HERE/receipt['source_bundle']
        self.assertEqual(a.sha(root/'manifest.json'),receipt['source_manifest_sha256'])
        self.assertTrue((root/'EXECUTION_STATUS.json').exists())
        self.assertFalse(receipt['performance_admitted'] or receipt['full_discovery_binding_eligible'])
        result=a.validate(root);self.assertEqual(result,a.read(replay/'results.json'))
        self.assertEqual((result['cells'],result['observations']),(24,14640))
        self.assertTrue(result['all_complete'] and result['qualified'])

    def test_accepted_qualified_runs_bind_current_source_and_exact_replay(self):
        registry=a.read(HERE/'QUALIFIED_RUNS.json')
        for phase in ['discovery','full']:
            entry=registry[phase]
            if entry is None:continue
            root=HERE/entry['path'];self.assertEqual(a.sha(root/'manifest.json'),entry['manifest_sha256'])
            for name,digest in a.read(root/'manifest.json')['files'].items():self.assertEqual(a.sha(root/name),digest)
            self.assertFalse((root/'EXECUTION_STATUS.json').exists())
            for name,digest in a.read(root/'metadata.json')['source_hashes'].items():
                if name.endswith('.rs'):self.assertEqual(a.sha(HERE/name),digest,name)
            self.assertEqual(a.sha(root/'protocol.json'),a.sha(HERE/('protocol_'+phase+'.json')))
            result=a.validate(root);self.assertEqual(result,a.read(root/'results.json'))
            self.assertTrue(result['qualified'] and result['all_complete'])
            self.assertEqual(result['phase'],phase)

    def test_qualified_predecessors_remain_exact_for_their_own_source(self):
        for entry in a.read(HERE/'QUALIFIED_RUNS.json').get('historical_qualified_runs',[]):
            root=HERE/entry['path'];self.assertEqual(a.sha(root/'manifest.json'),entry['manifest_sha256'])
            for name,digest in a.read(root/'manifest.json')['files'].items():self.assertEqual(a.sha(root/name),digest)
            result=a.validate(root);self.assertEqual(result,a.read(root/'results.json'))
            self.assertTrue(result['qualified'] and result['all_complete'])

    def test_paired_stream_split_is_lossless_and_rejects_offset_drift(self):
        root=self.temporary();cell='n12-discovery-17-planted'
        rows=[{'type':'fixture','mode':'aa'},{'type':'sample','value':1},
              {'type':'fixture','mode':'ab'},{'type':'sample','value':2}]
        path=root/(cell+'-paired.jsonl');path.write_text(''.join(json.dumps(r)+'\n' for r in rows))
        path.with_suffix('.stderr').write_text('')
        receipt={'stdout_sha256':a.sha(path),'stderr_sha256':a.sha(path.with_suffix('.stderr')),
                 'worker_command':['worker','12','17','planted','16','200000','1','paired'],
                 'conditions_file':cell+'-paired-conditions.jsonl','conditions_sha256':'unit-test-placeholder'}
        stages=runner.split_paired(path,cell,root,receipt)
        entry={'paired':receipt,'stages':stages};a.verify_paired_stream(root,entry,cell)
        self.assertEqual((root/(cell+'-aa.jsonl')).read_bytes()+(root/(cell+'.jsonl')).read_bytes(),path.read_bytes())
        stages['ab']['derived_from']['byte_start']+=1
        with self.assertRaises(AssertionError):a.verify_paired_stream(root,entry,cell)

    def test_paired_stream_requires_aa_before_ab(self):
        root=self.temporary();path=root/'bad.jsonl'
        path.write_text(json.dumps({'type':'fixture','mode':'ab'})+'\n'+json.dumps({'type':'fixture','mode':'aa'})+'\n')
        with self.assertRaises(AssertionError):runner.split_paired(path,'bad',root,{})


if __name__=='__main__':unittest.main()
