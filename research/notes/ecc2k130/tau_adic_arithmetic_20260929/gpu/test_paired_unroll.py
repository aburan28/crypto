"""Study wiring, paired statistics, isolation policy and process cleanup; no GPU."""
import contextlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch,Mock
import launch_gpu
import run_gpu
from test_readiness import cupy_stub,panel


class PairedTests(unittest.TestCase):
    def test_interruption_also_terminates_process_group(self):
        process=Mock(pid=12345)
        process.communicate.side_effect=[KeyboardInterrupt(),('','')]
        with patch.object(launch_gpu.subprocess,'Popen',return_value=process), \
             patch.object(launch_gpu.os,'killpg') as kill:
            with self.assertRaises(KeyboardInterrupt):launch_gpu.bounded(['unused'],{})
        kill.assert_called_once_with(12345,launch_gpu.signal.SIGTERM)

    def test_full_loop_retains_two_controls_and_alternates_round_order(self):
        cp,launches=cupy_stub()
        cp.cuda.Event=lambda:SimpleNamespace(record=lambda:None,synchronize=lambda:None)
        cp.cuda.get_elapsed_time=lambda begin,end:100.0  # Test constant, never a device measurement.
        cp.get_default_memory_pool=lambda:SimpleNamespace(free_all_blocks=lambda:None)
        fixture=panel();fixture.update(field={},curve={},source_receipt_sha256='test-only')
        with tempfile.TemporaryDirectory() as folder, \
             patch.dict('sys.modules',{'cupy':cp}), \
             patch.object(run_gpu.fx,'panels',return_value=iter([fixture])), \
             patch.object(run_gpu.fx.ref,'host',return_value={'test_backend':'host_stub'}), \
             patch.object(run_gpu,'hardware',return_value={'test_backend':'host_stub'}), \
             contextlib.redirect_stdout(io.StringIO()):
            output=Path(folder)/'test.json'
            self.assertEqual(run_gpu.run(output,'benchmark','square-unroll'),0)
            receipt=json.loads(output.read_text())
        row=receipt['panels'][0]
        self.assertEqual([len(pairs) for pairs in row['aa_by_control'].values()],[5,5])
        self.assertEqual(len(row['rounds']),7)
        for i,round_ in enumerate(row['rounds']):
            expected=receipt['variant_order'] if i%2==0 else receipt['variant_order'][::-1]
            self.assertEqual(round_['order'],expected)
        self.assertEqual(len(launches),4+4+1+20+28)
        self.assertFalse(any(receipt['screen_by_candidate'].values()))

    def test_compiler_cache_keeps_default_and_candidate_separate(self):
        options=[]
        def module(**kwargs):
            options.append(kwargs['options'])
            kernel=SimpleNamespace(attributes={})
            return SimpleNamespace(get_function=lambda name:kernel)
        cp=SimpleNamespace(RawModule=module);cache={}
        a,meta_a=run_gpu.load_kernel(cp,cache,83,1,0)
        b,meta_b=run_gpu.load_kernel(cp,cache,83,1,1)
        again,meta_again=run_gpu.load_kernel(cp,cache,83,1,0)
        self.assertIs(a,again);self.assertIsNot(a,b)
        self.assertEqual(len(options),2)
        self.assertEqual(options[1],options[0]+('-DLINEAR_SQUARE_NOUNROLL=1',))
        self.assertTrue(meta_again['reused_module'])

    def test_candidate_is_in_device_smoke_without_timing(self):
        cp,launches=cupy_stub();compiled=[];original=cp.RawModule
        def module(**kwargs):
            compiled.append(kwargs['options']);return original(**kwargs)
        cp.RawModule=module
        with tempfile.TemporaryDirectory() as folder, \
             patch.dict('sys.modules',{'cupy':cp}), \
             patch.object(run_gpu.fx,'panels',return_value=iter([panel()])), \
             patch.object(run_gpu.fx.ref,'host',return_value={'test_backend':'host_stub'}), \
             patch.object(run_gpu,'hardware',return_value={'test_backend':'host_stub'}), \
             contextlib.redirect_stdout(io.StringIO()):
            output=Path(folder)/'test.json'
            self.assertEqual(run_gpu.run(output,'smoke','square-unroll'),0)
            receipt=json.loads(output.read_text())
        self.assertEqual(launches,[1]*4)
        self.assertEqual(len(compiled),2)
        self.assertIn('-DLINEAR_SQUARE_NOUNROLL=1',compiled[1])
        self.assertEqual(receipt['panels'],[])
        self.assertNotIn('screen_by_candidate',receipt)

    def summary(self,times,noise=(0,0)):
        variants=run_gpu.variants_for('square-unroll')
        sample=lambda ms:{'event_ms':ms,'ms_per_launch':ms,'launches':1}
        aa={run_gpu.label(variants[i]):[[sample(100*(1+n)),sample(100)]]*5 for i,n in zip((0,2),noise)}
        row={'aa_by_control':aa,'rounds':[{'samples':{run_gpu.label(v):sample(t) for v,t in zip(variants,times)}} for _ in range(7)]}
        return run_gpu.summarize(row,variants,'square-unroll'),variants

    def test_each_recoding_uses_its_own_control(self):
        summary,v=self.summary([1000,800,100,120])
        binary=summary[run_gpu.label(v[1])];tau=summary[run_gpu.label(v[3])]
        self.assertEqual(binary['median_paired_kernel_cost_ratio'],.8)
        self.assertTrue(binary['screen_passed'])
        self.assertEqual(tau['median_paired_kernel_cost_ratio'],1.2)
        self.assertFalse(tau['screen_passed'])
        self.assertEqual(tau['control'],run_gpu.label(v[2]))

    def test_noise_and_short_samples_block_screen(self):
        summary,v=self.summary([100,80,100,80],noise=(.3,0))
        self.assertFalse(summary[run_gpu.label(v[1])]['screen_passed'])
        self.assertTrue(summary[run_gpu.label(v[3])]['screen_passed'])
        summary,v=self.summary([100,40,100,80])
        self.assertTrue(summary[run_gpu.label(v[1])]['short_samples'])
        self.assertFalse(summary[run_gpu.label(v[1])]['screen_passed'])

    def test_isolation_requires_full_siblings_and_spare_cpu(self):
        siblings=lambda cpu:{0,1} if cpu<2 else {2,3}
        self.assertEqual(launch_gpu.select_cpus({0,1,2,3},siblings),'0,1')
        with self.assertRaises(RuntimeError):launch_gpu.select_cpus({0,1},siblings)
        with self.assertRaises(RuntimeError):launch_gpu.select_cpus({0,2},siblings)

    def test_contended_and_incomplete_conditions_are_ineligible(self):
        record={'schema':'isolated-bench/1','mode':'run','run':{'exit_status':0,'contended':False},
                'left_on_reserved':{'user_threads':[]}}
        self.assertTrue(launch_gpu.timing_eligible(record))
        record['run']['contended']=True
        self.assertFalse(launch_gpu.timing_eligible(record))
        record['run']['contended']=False;record['left_on_reserved']['user_threads']=['other']
        self.assertFalse(launch_gpu.timing_eligible(record))
        self.assertFalse(launch_gpu.timing_eligible({}))

    def test_timed_launch_is_wrapped_and_smoke_is_not(self):
        command=launch_gpu.command_for('benchmark','square-unroll','result.json','conditions.jsonl','0,1')
        self.assertIn('isolated_entry.py',Path(command[1]).name)
        self.assertEqual(command[command.index('--study')+1],'square-unroll')
        smoke=launch_gpu.command_for('smoke','square-unroll','result.json','unused')
        self.assertEqual(Path(smoke[1]).name,'run_gpu.py')

    def test_timeout_terminates_child_group_and_unwinds_finally(self):
        code='''import subprocess,sys,time,signal
def stop(signum,frame): raise SystemExit(143)
signal.signal(signal.SIGTERM,stop)
child=subprocess.Popen([sys.executable,'-c','import time; time.sleep(60)'])
print(child.pid,flush=True)
try: time.sleep(60)
finally:
 child.wait(timeout=5)
 print('cleanup complete',flush=True)
'''
        result=launch_gpu.bounded([sys.executable,'-c',code],os.environ.copy(),seconds=.5)
        self.assertTrue(result['timed_out'])
        self.assertIn('cleanup complete',result['stdout'])
        child=int(result['stdout'].splitlines()[0])
        self.assertFalse(Path(f'/proc/{child}').exists())


if __name__=='__main__':unittest.main()
