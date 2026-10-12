"""Offline regression tests for instance attribution and benchmark separation."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import benchmark_hardware as hardware
spec=importlib.util.spec_from_file_location('hardware_report',ROOT/'benchmarks/by-hardware/report.py')
report=importlib.util.module_from_spec(spec);spec.loader.exec_module(report)


class HardwareBenchmarks(unittest.TestCase):
    def identity(self):
        return dict(provider='aws',instance_type='g7.2xlarge',hardware_type=None,
                    gpu_name='NVIDIA RTX PRO 4500 Blackwell Server Edition',driver_version='595.91.07',
                    power_limit_watts=165.0,memory_mib=32623,compute_capability='12.0',
                    measured_gpu_count=1,measurement_scope='single_gpu',cpu_affinity_count=8)

    def test_same_gpu_different_instance_sizes_cannot_be_compared(self):
        a=self.identity();b=dict(a,instance_type='g7.4xlarge')
        with self.assertRaisesRegex(ValueError,'different hardware or instance'):
            hardware.require_matched_hardware([{'hardware':a},{'hardware':b}])

    def test_provider_gpu_power_and_driver_changes_separate_results(self):
        a=self.identity()
        for change in ({'provider':'modal'},{'gpu_name':'RTX PRO 6000'},
                       {'power_limit_watts':300.0},{'driver_version':'580.95.05'}):
            with self.subTest(change=change),self.assertRaises(ValueError):
                hardware.require_matched_hardware([{'hardware':a},{'hardware':dict(a,**change)}])

    def test_unknown_instance_is_not_an_implicit_match(self):
        a=self.identity();a['instance_type']=None
        with self.assertRaisesRegex(ValueError,'explicit provider'):
            hardware.require_matched_hardware([{'hardware':a},{'hardware':copy.deepcopy(a)}])

    def test_same_type_and_limits_can_match_despite_uuid_and_capture_time(self):
        a=dict(self.identity(),gpu_uuid='GPU-A',captured_at_utc='one')
        b=dict(self.identity(),gpu_uuid='GPU-B',captured_at_utc='two')
        self.assertEqual(hardware.require_matched_hardware([{'hardware':a},{'hardware':b}]),hardware.hardware_key(a))

    def test_reused_labels_cannot_overwrite_another_instance_or_legacy_result(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'candidate.json'
            path.write_text(json.dumps({'hardware':self.identity(),'rate':1}))
            hardware.check_destination(path,self.identity())
            with self.assertRaisesRegex(ValueError,'different hardware'):
                hardware.check_destination(path,dict(self.identity(),instance_type='g7.4xlarge'))
            path.write_text(json.dumps({'rate':14.6}))
            with self.assertRaisesRegex(ValueError,'no hardware identity'):
                hardware.check_destination(path,self.identity())
            self.assertEqual(json.loads(path.read_text()),{'rate':14.6})

    def test_ec2_dmi_and_explicit_provider_do_not_leak_across_platforms(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/'sys_vendor').write_text('Amazon EC2\n')
            (root/'product_name').write_text('g7.2xlarge\n')
            self.assertEqual(hardware.platform_identity({},root)['instance_type'],'g7.2xlarge')
            modal=hardware.platform_identity({'ECC_BENCH_PROVIDER':'modal','ECC_BENCH_HARDWARE_TYPE':'RTX-PRO-6000'},root)
            self.assertIsNone(modal['instance_type'])
            self.assertEqual(modal['hardware_type'],'RTX-PRO-6000')

    def test_cuda_visibility_mask_selects_the_measured_gpu(self):
        result=type('Result',(),{'stdout':'RTX PRO 4500, GPU-test, 595.91.07, 165.0, 32623, 12.0\n'})()
        with patch.dict('os.environ',{'CUDA_VISIBLE_DEVICES':'3,5'},clear=True),patch.object(hardware,'platform_identity',return_value=self.identity()),patch.object(hardware.subprocess,'run',return_value=result) as run:
            captured=hardware.capture_hardware()
            command=run.call_args.args[0]
            self.assertEqual(command[command.index('-i')+1],'3')
            self.assertEqual(captured['measured_gpu_count'],1)
            self.assertEqual(captured['gpu_selector'],'3')

    def test_frozen_values_are_attributed_to_their_own_measured_hosts(self):
        rows=report.build_report()['rows']
        aws=next(r for r in rows if r['instance_type']=='g7.2xlarge')
        modal=next(r for r in rows if r['provider']=='modal')
        evidence=aws['evidence'][0]
        data=(ROOT.parent/evidence['path']).read_bytes()
        self.assertEqual(hashlib.sha256(data).hexdigest(),evidence['sha256'])
        local=json.loads(data)
        self.assertTrue(local['passed'])
        previous=json.loads((ROOT/'benchmarks/local-10b/final-inline-audit.json').read_text())
        native=json.loads((ROOT/'benchmarks/shared-sigma/native-audit.json').read_text())
        current=next(r for r in local['rows'] if r['family']=='final benchmark' and r['variant']=='candidate')
        self.assertEqual(aws['benchmark_billion_updates_per_second'],current['median_billion_updates_per_second'])
        self.assertEqual(aws['instance_type'],local['hardware']['instance_type'])
        self.assertEqual(aws['binary_sha256'],local['selected_binary_sha256'])
        self.assertEqual(aws['configuration'],local['configuration'])
        collection=next(r for r in local['rows'] if r['family']=='final collection' and r['variant']=='candidate')
        self.assertEqual(aws['collection_billion_updates_per_second'],collection['median_billion_updates_per_second'])
        self.assertIn(previous['summary']['benchmark']['median_million_updates_per_second']['candidate']/1000,
                      [r['benchmark_billion_updates_per_second'] for r in aws['historical_measurements']])
        for historical in aws['historical_measurements']:
            for source in historical['evidence']:
                self.assertEqual(hashlib.sha256((ROOT.parent/source['path']).read_bytes()).hexdigest(),source['sha256'])
        self.assertEqual(modal['benchmark_billion_updates_per_second'],native['benchmark']['rate']/1000)
        self.assertIsNone(modal['instance_type'])
        self.assertNotEqual(aws['configuration'],modal['configuration'])

    def test_untested_instances_do_not_inherit_or_multiply_modal_rates(self):
        rows=report.build_report()['rows']
        unknown=[r for r in rows if r['benchmark_status']=='unmeasured']
        self.assertEqual(len(unknown),11)
        for row in unknown:
            self.assertIsNone(row['benchmark_billion_updates_per_second'])
            self.assertIsNone(row['collection_billion_updates_per_second'])
        self.assertTrue(all(r['instance_aggregate_billion_updates_per_second'] is None for r in rows))


if __name__=='__main__':
    unittest.main()
