"""Receipt durability and launch-gate tests with a CPU stub, never a GPU test."""
import json
import contextlib
import io
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
import receipt_io
import run_gpu


class HostArray(np.ndarray):
    def get(self):
        return np.asarray(self).copy()


def cupy_stub(mismatch=False):
    launches=[]

    def kernel(grid,block,args):
        # Fixtures in this orchestration test have scalar 1. No CUDA is executed.
        points,_,_,n,_,out,_=args
        out[:]=0 if mismatch else points
        launches.append(int(n))

    kernel.attributes={'test_backend':'host_stub'}
    module=SimpleNamespace(get_function=lambda name:kernel)
    runtime=SimpleNamespace(getDeviceProperties=lambda device:{'test_backend':'host_stub'},
                            runtimeGetVersion=lambda:0,driverGetVersion=lambda:0)
    fake=SimpleNamespace(__version__='host-stub-no-device',
                         asarray=lambda a:np.asarray(a).view(HostArray),
                         empty_like=lambda a:np.empty_like(a).view(HostArray),
                         RawModule=lambda **kwargs:module,
                         cuda=SimpleNamespace(runtime=runtime,Device=lambda:SimpleNamespace(id=0),
                                              Stream=SimpleNamespace(null=SimpleNamespace(synchronize=lambda:None))))
    # Deliberately no Event/timing API: a smoke run must not enter timing code.
    return fake,launches


def panel():
    return {'m':83,'holdout':False,'seed':7,'identity':{'test_fixture':True},
            'input_sha256':'test-only','sage_output_sha256':'test-only',
            'cases':[((1,0),1)],'expected':[(1,0)],'oracle_seconds':0}


class ReadinessTests(unittest.TestCase):
    def setUp(self):
        self.folder=tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.path=Path(self.folder.name)/'receipt.json'

    def test_existing_receipt_cannot_be_overwritten(self):
        receipt_io.reserve(self.path,{'status':'original'})
        with self.assertRaises(FileExistsError):
            receipt_io.reserve(self.path,{'status':'replacement'})
        self.assertEqual(receipt_io.read_partial(self.path)['status'],'original')

    def test_interrupted_replace_preserves_previous_valid_checkpoint(self):
        receipt_io.reserve(self.path,{'status':'previous'})
        with patch.object(receipt_io.os,'replace',side_effect=OSError('injected interruption')):
            with self.assertRaises(OSError):
                receipt_io.save(self.path,{'status':'new'})
        self.assertEqual(receipt_io.read_partial(self.path)['status'],'previous')
        self.assertEqual(list(self.path.parent.iterdir()),[self.path])

    def test_missing_or_malformed_partial_is_unknown_execution(self):
        self.assertIsNone(receipt_io.read_partial(self.path)['gpu_executed'])
        for text in ('{','[]'):
            self.path.write_text(text)
            result=receipt_io.read_partial(self.path)
            self.assertEqual(result['status'],'receipt_unavailable')
            self.assertIsNone(result['gpu_executed'])

    def invoke(self,mode,mismatch):
        fake,launches=cupy_stub(mismatch)
        with patch.dict('sys.modules',{'cupy':fake}), \
             patch.object(run_gpu.fx,'panels',return_value=iter([panel()])), \
             patch.object(run_gpu.fx.ref,'host',return_value={'test_backend':'host_stub'}), \
             patch.object(run_gpu,'hardware',return_value={'test_backend':'host_stub'}), \
             contextlib.redirect_stdout(io.StringIO()):
            status=run_gpu.run(self.path,mode)
        return status,json.loads(self.path.read_text()),launches

    def test_smoke_checks_four_arms_and_produces_no_timing_summary(self):
        status,result,launches=self.invoke('smoke',False)
        self.assertEqual(status,0)
        self.assertEqual(launches,[1]*4)
        self.assertEqual(result['smoke']['verified_outputs'],4)
        self.assertEqual(result['smoke']['status'],'passed')
        self.assertEqual(result['panels'],[])
        self.assertIsNone(result['rho_speedup'])
        self.assertIsNone(result['walk_iterations_per_second'])

    def test_device_mismatch_aborts_benchmark_and_preserves_launch_evidence(self):
        status,result,launches=self.invoke('benchmark',True)
        self.assertEqual(status,1)
        self.assertEqual(launches,[1])
        self.assertEqual(result['status'],'failed')
        self.assertTrue(result['gpu_executed'])
        self.assertTrue(result['gpu_launch_submitted'])
        self.assertEqual(result['smoke']['verified_outputs'],0)
        self.assertEqual(result['smoke']['status'],'failed')
        self.assertEqual(result['panels'],[])
        arm=next(iter(result['smoke']['panels'][0]['invocations'].values()))
        self.assertEqual(arm['status'],'mismatch')


if __name__=='__main__':
    unittest.main()
