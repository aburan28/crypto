"""Accounting/process tests. Simulated workers are never GPU speed measurements."""
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchmark_instance import Work, parse_client, run_cohort, run_ids, select_devices


class InstanceBenchmarkTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.devices = [dict(index=0, uuid='GPU-a'), dict(index=1, uuid='GPU-b')]

    def tearDown(self):
        self.temp.cleanup()

    def worker(self, mode='normal'):
        path = self.root / ('worker-' + mode)
        path.write_text('''#!/usr/bin/env python3
import os,sys,time
from pathlib import Path
args=sys.argv[1:]
def arg(name,default=0):return int(args[args.index(name)+1]) if name in args else default
workers=arg('--threads');steps=arg('--steps');launches=arg('--launches');verify=arg('--verify')
weight=arg('--dp-weight') if '--bench' not in args else 0
print(f'backend cuda-packed131: {workers} threads x 16 slots x 1 lanes = {workers*16} walks, dp weight {weight}, {steps} steps per launch',flush=True)
mode=''' + repr(mode) + '''
if mode=='exit':sys.exit(7)
time.sleep(1 if mode=='timeout' else (0.035 if os.environ['CUDA_VISIBLE_DEVICES']=='GPU-a' else 0.09))
points=max(verify,4) if weight else 0
if '--dp-file' in args:
 path=Path(args[args.index('--dp-file')+1]);path.write_bytes(bytes(32*(points-(1 if mode=='short-dp' else 0))))
print(f' {workers*16*steps*launches} iterations  {points} dp  {points} stored  0 dropped')
print(f' finished: 1000.000 M it/s, {points} distinguished points ({verify} verified against the reference, 0 dropped)')
''')
        path.chmod(0o755)
        return path

    def observer(self, binary, extra=None):
        def observe(owned):
            rows = [dict(uuid=uuid, pid=pid, name=str(binary)) for pid, uuid in owned.items()]
            if owned and extra:
                rows.append(extra)
            return rows
        return observe

    def test_unique_gpu_selection_and_run_ids(self):
        self.assertEqual(select_devices(self.devices, '1,0'), self.devices[::-1])
        self.assertEqual(select_devices(self.devices, 'all'), self.devices)
        with self.assertRaises(ValueError):select_devices(self.devices, '0,GPU-a')
        with self.assertRaises(ValueError):select_devices(self.devices, '2')
        self.assertEqual(run_ids(65532, 2, 1), [[65532, 65533], [65534, 65535]])
        with self.assertRaises(ValueError):run_ids(65533, 2, 1)
        values=sum(run_ids(52000, 8, 3), [])
        self.assertEqual(len(values), len(set(values)))

    def test_budget_geometry_and_drops(self):
        work = Work(2, 16, 3, 4)
        text = ('backend cuda-packed131: 2 threads x 16 slots x 1 lanes = 32 walks, dp weight 0, 3 steps per launch\n'
                '384 iterations\n finished: 1000.000 M it/s, 0 distinguished points (0 verified against the reference, 0 dropped)\n')
        self.assertEqual(parse_client(text, work)['iterations'], 384)
        for wrong in [text.replace('384 iterations','383 iterations'), text.replace('x 16 slots','x 32 slots'), text.replace('0 dropped','1 dropped'), text+'MISMATCH', text.replace('finished:', 'partial:')]:
            with self.assertRaises(ValueError):parse_client(wrong, work)
        with self.assertRaises(ValueError):Work(1,16,1,0).check()
        with self.assertRaises(ValueError):Work(2**30,16,2**30,2**30).check()
        Work(2**28,16,1,1).check()
        with self.assertRaises(ValueError):Work(2**28+1,16,1,1).check()

    def test_common_window_counts_work_and_owns_processes(self):
        binary=self.worker();work=Work(2,16,3,4);out=self.root/'cohort'
        result=run_cohort(binary,self.devices,[10,11],work,out,observer=self.observer(binary),poll_seconds=0.01)
        self.assertTrue(result['passed'])
        self.assertEqual(result['completed_scalar_iterations'],768)
        self.assertEqual(result['cohort_iterations_per_second'],768/result['elapsed_seconds'])
        self.assertGreater(result['elapsed_seconds'],0.09)
        self.assertLess(result['cohort_iterations_per_second'],100000)
        self.assertEqual([w['run_id'] for w in result['workers']],[10,11])
        for row in result['workers']:
            with self.assertRaises(ChildProcessError):os.waitpid(row['pid'],os.WNOHANG)
        with self.assertRaises(FileExistsError):run_cohort(binary,self.devices,[10,11],work,out,observer=self.observer(binary))

    def test_failed_or_timed_out_workers_have_no_throughput(self):
        for mode,timeout in [('exit',2),('timeout',0.025)]:
            binary=self.worker(mode)
            r=run_cohort(binary,self.devices,[20,21],Work(2,16,3,4),self.root/mode,timeout=timeout,observer=self.observer(binary),poll_seconds=0.01)
            self.assertFalse(r['passed']);self.assertIsNone(r['cohort_iterations_per_second']);self.assertIsNone(r['completed_scalar_iterations'])
            if mode=='timeout':self.assertTrue(r['timed_out'])
            for row in r['workers']:
                with self.assertRaises(ChildProcessError):os.waitpid(row['pid'],os.WNOHANG)

    def test_foreign_process_or_missing_gpu_observation_invalidates(self):
        binary=self.worker();work=Work(2,16,3,4)
        foreign=dict(pid=99999999,uuid='GPU-a',name='foreign-worker')
        r=run_cohort(binary,self.devices,[30,31],work,self.root/'foreign',observer=self.observer(binary,foreign),poll_seconds=0.01)
        self.assertFalse(r['passed']);self.assertTrue(r['observation_errors'])
        r=run_cohort(binary,self.devices,[30,31],work,self.root/'unobserved',observer=lambda _:[],poll_seconds=0.01)
        self.assertFalse(r['passed']);self.assertEqual(r['observed_worker_pids'],[])

    def test_duplicate_assignment_rejected_before_launch(self):
        binary=self.worker()
        for devices,ids in [(self.devices,[1,1]),([self.devices[0],self.devices[0]],[1,2]),(self.devices,[65535,65536])]:
            with self.assertRaises(ValueError):run_cohort(binary,devices,ids,Work(),self.root/'invalid',observer=self.observer(binary))
        self.assertFalse((self.root/'invalid').exists())

    def test_preflight_failures_are_saved_without_launching(self):
        binary=self.worker()
        def failed_observer(_):
            raise RuntimeError('GPU observer unavailable')
        busy_observer=lambda _: [dict(pid=99999999,uuid='GPU-a',name='foreign-worker')]
        for name,observer in [('observer-error',failed_observer),('busy',busy_observer)]:
            out=self.root/name
            with self.assertRaises(RuntimeError):
                run_cohort(binary,self.devices,[50,51],Work(),out,observer=observer)
            self.assertFalse(json.loads((out/'preflight.json').read_text())['passed'])
            self.assertTrue((out/'plan.json').exists())
            self.assertFalse((out/'owned-processes.json').exists())
            self.assertFalse(list(out.glob('gpu*.log')))

    def test_collection_record_count_and_replay_admission(self):
        work=Work(2,16,3,4,True,16,52)
        for mode,passed in [('normal',True),('short-dp',False)]:
            binary=self.worker(mode)
            r=run_cohort(binary,self.devices,[40,41],work,self.root/('collect-'+mode),observer=self.observer(binary),poll_seconds=0.01)
            self.assertEqual(r['passed'],passed)
            if passed:self.assertEqual(sum(w['cpu_replays'] for w in r['workers']),32)


if __name__=='__main__':unittest.main()
