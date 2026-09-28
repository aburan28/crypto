import importlib.util
from pathlib import Path
import unittest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('tail_runner',HERE/'run.py')
runner=importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class CommandTests(unittest.TestCase):
    def test_mode_is_present_exactly_once(self):
        protocol={'repetitions':12,'enable_census':True}
        command=runner.worker_command(Path('/tmp/worker'),36,17,8,'quadratic',protocol)
        self.assertEqual(command,['/tmp/worker','36','17','8','quadratic','12','with-census'])
        self.assertEqual(command.count('with-census'),1)

    def test_legacy_command_has_no_optional_mode(self):
        for protocol in [{'repetitions':12},{'repetitions':12,'enable_census':False}]:
            self.assertEqual(runner.worker_command('/tmp/worker',12,937,1,'linear_drop',protocol),
                             ['/tmp/worker','12','937','1','linear_drop','12'])

    def test_balanced_mode_has_two_distinct_flags(self):
        protocol={'repetitions':42,'enable_census':True,'balanced_order':True}
        self.assertEqual(runner.worker_command('/tmp/worker',36,17,8,'quadratic',protocol),
                         ['/tmp/worker','36','17','8','quadratic','42','with-census','balanced-order'])

    def test_balancing_requires_the_census_arm_set(self):
        with self.assertRaises(ValueError):
            runner.worker_command('/tmp/worker',36,17,8,'quadratic',{'repetitions':42,'balanced_order':True})

    def test_matched_mode_uses_one_order_flag(self):
        p={'repetitions':42,'enable_census':True,'balanced_order':True,'matched_predecessors':True}
        self.assertEqual(runner.worker_command('/tmp/worker',36,17,8,'quadratic',p)[-2:],['with-census','matched-order'])

    def test_matched_mode_requires_balancing(self):
        with self.assertRaises(ValueError):
            runner.worker_command('/tmp/worker',36,17,8,'quadratic',{'repetitions':42,'matched_predecessors':True})

    def test_long_observation_retains_one_mode_and_explicit_batch_count(self):
        p={'repetitions':42,'enable_census':True,'balanced_order':True,'matched_predecessors':True,'inner_batches':16}
        self.assertEqual(runner.worker_command('/tmp/worker',36,17,8,'quadratic',p),
                         ['/tmp/worker','36','17','8','quadratic','42','with-census','matched-order','16'])


if __name__=='__main__':unittest.main()
