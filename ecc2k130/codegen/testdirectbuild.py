"""Exercise direct-reducer build and result helpers without Modal or a GPU."""
import ast
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import tempfile
import time
from types import SimpleNamespace
import unittest

from benchreport import benchResult, bestResult, summarizeSamples

ROOT = Path(__file__).resolve().parents[1]


def nodes(path):
    return ast.parse((ROOT / path).read_text()).body


def execute(body, env):
    exec(compile(ast.Module(body=body, type_ignores=[]), 'actual-wrapper', 'exec'), env)


def function(name, env, path='modal_app.py'):
    node = next(n for n in nodes(path) if isinstance(n, ast.FunctionDef) and n.name == name)
    node.decorator_list = []
    execute([node], env)
    return env[name]


def assignment(name, env):
    node = next(n for n in nodes('modal_app.py') if isinstance(n, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == name for t in n.targets))
    execute([node], env)


def environment(mode='0'):
    env = dict(re=re, hashlib=hashlib, pathlib=SimpleNamespace(Path=Path),
               subprocess=subprocess, time=time, json=json, benchResult=benchResult,
               summarizeSamples=summarizeSamples, bestResult=bestResult,
               CUDA_VERSION='13.0.0', DEFAULT_GPU='RTX-PRO-6000', BAKED_ARCHES=('120',), REMOTE='/unused',
               print=lambda *a, **k: None,
               bakedIntact=[True], computeCapability=lambda: '120', gpuName=lambda: 'fixture GPU')
    for key in ('SINGLE_PRODUCT', 'CACHE_DENOM', 'BY_VALUE', 'POLY_CHAIN',
                'UNROLL_INV', 'PAIR_PRODUCTS', 'POLY_STATE'):
        env['PACKED_' + key] = '1'
    env.update(PACKED_PERM_SIGMA='3', PACKED_DIRECT_REDUCE=mode)
    assignment('BAKED', env)
    return env


def raw(mode='0'):
    return (f'packed direct reduction: {mode}\n'
            'finished: 6000.000 M it/s, 0 distinguished points (0 verified against the reference, 0 dropped)\n')


class DirectBuildTests(unittest.TestCase):
    def test_environment_defaults_and_rejects_invalid_values(self):
        body = nodes('modal_app.py')
        index = next(i for i, n in enumerate(body) if isinstance(n, ast.Assign)
                     and any(isinstance(t, ast.Name) and t.id == 'PACKED_DIRECT_REDUCE' for t in n.targets))
        for value in (None, '0', '1', '', '2', '-1', 'true'):
            env = dict(os=SimpleNamespace(environ={} if value is None else {'ECC_PACKED_DIRECT_REDUCE': value}))
            with self.subTest(value=value):
                if value in (None, '0', '1'):
                    execute(body[index:index + 2], env)
                    self.assertEqual(env['PACKED_DIRECT_REDUCE'], value or '0')
                else:
                    with self.assertRaises(ValueError):
                        execute(body[index:index + 2], env)

    def test_image_preserves_flag_and_bakes_matching_binary(self):
        class Image:
            def __init__(self):
                self.calls = {}
            def __getattr__(self, name):
                def record(*args, **kwargs):
                    self.calls[name] = args
                    return self
                return record
        for mode in ('0', '1'):
            env = environment(mode)
            image = Image()
            env.update(modal=SimpleNamespace(Image=SimpleNamespace(from_registry=lambda *a, **k: image)),
                       GENCODE='fixture', LOCAL=ROOT)
            assignment('image', env)
            self.assertEqual(image.calls['env'][0]['ECC_PACKED_DIRECT_REDUCE'], mode)
            builds = [c for c in image.calls['run_commands'] if 'make gpu ' in c]
            self.assertEqual(len(builds), 1)
            self.assertIn('PACKED_DIRECT_REDUCE=' + mode, builds[0])
            self.assertEqual(env['BAKED']['packedDirectReduction'], mode == '1')

    def test_cache_mismatch_rebuilds_and_matching_bake_is_not_reused_afterward(self):
        env = environment('0')
        commands = []
        env.update(sh=lambda *a, **k: (0, ''),
                   shStream=lambda c, **k: (commands.append(c) or 0, 'built'))
        build = function('buildFor', env)
        self.assertTrue(build(32, 128, 0)[0])
        self.assertEqual(commands, [])
        env['PACKED_DIRECT_REDUCE'] = '1'
        self.assertTrue(build(32, 128, 0)[0])
        self.assertIn('PACKED_DIRECT_REDUCE=1', commands[-1])
        env['PACKED_DIRECT_REDUCE'] = '0'
        self.assertTrue(build(32, 128, 0)[0])
        self.assertIn('PACKED_DIRECT_REDUCE=0', commands[-1])
        self.assertEqual(len(commands), 2)

    def test_rebuild_false_rejects_a_different_reducer(self):
        env = environment('0')
        env['PACKED_DIRECT_REDUCE'] = '1'
        run = function('runBench', env)
        with self.assertRaisesRegex(ValueError, 'matching baked binary'):
            run(rebuild=False, packed=True)

    def test_result_gate_never_ranks_wrong_missing_duplicate_or_failed_identity(self):
        for mode in ('0', '1'):
            check = function('checkPackedReduction', environment(mode))
            good = benchResult('fixture', 0, raw(mode))
            self.assertTrue(check(good))
            self.assertEqual(good['packedDirectReduction'], mode == '1')
            for rc, output in [(0, raw(str(1 - int(mode)))), (0, raw(mode).split('\n', 1)[1]),
                               (0, raw(mode) + f'packed direct reduction: {mode}\n'),
                               (0, raw('2')), (1, raw(mode))]:
                sample = benchResult('fixture', rc, output)
                self.assertFalse(check(sample))
                self.assertEqual(sample['rate'], 0)
                self.assertFalse(summarizeSamples([good, sample])['valid'])

    def test_measurement_checks_packed_mode_but_preserves_other_backends(self):
        env = environment('1')
        env['sh'] = lambda *a, **k: (0, raw('0'))
        function('checkPackedReduction', env)
        measure = function('measureBench', env)
        self.assertFalse(measure(1, 1, 0, False, 1, packed=True)['valid'])
        self.assertTrue(measure(1, 1, 0, False, 1, packed=False)['valid'])

    def test_identity_and_compile_check_record_actual_requested_flag(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, text in [('Makefile', ''), ('modal_app.py', ''), ('ecc2k130', 'fixture'),
                               ('generated/eccF131.h', 'LEAF = 66')]:
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text)
            for mode in ('0', '1'):
                env = environment(mode)
                env.update(REMOTE=directory, sh=lambda *a, **k: (0, 'fixture'),
                           buildFor=lambda *a, **k: (True, 'built'))
                identity = function('benchmarkIdentity', env)
                self.assertEqual(identity(True)['packedDirectReduction'], mode == '1')
                self.assertIsNone(identity(False)['packedDirectReduction'])
                self.assertEqual(function('runCompileCheck', env)()['packedDirectReduction'], mode == '1')

    def test_bench_and_autotune_rows_retain_flag_without_changing_config_format(self):
        for mode in ('0', '1'):
            env = environment(mode)
            env.update(buildFor=lambda *a, **k: (True, 'built'),
                       benchmarkIdentity=lambda *a, **k: {},
                       measureBench=lambda *a, **k: dict(valid=True, rate=6000.0),
                       os=SimpleNamespace(makedirs=lambda *a, **k: None),
                       open=lambda *a, **k: io.StringIO(),
                       volume=SimpleNamespace(commit=lambda: None))
            bench = function('runBench', env)(rebuild=False, packed=True)
            self.assertEqual(bench['packedDirectReduction'], mode == '1')
            function('autotuneConfigs', env)
            tune = function('runAutotune', env)(configs='0:32:128:2', packed=True)
            self.assertEqual(tune['results'][0]['packedDirectReduction'], mode == '1')
            self.assertEqual(len(env['autotuneConfigs']('', '', '', '', '0:32:128:2')[0]), 7)

    def test_audit_forwards_flag_and_rejects_an_arithmetic_binary_with_wrong_mode(self):
        for mode in ('0', '1'):
            with tempfile.TemporaryDirectory() as directory:
                commands = []
                client = SimpleNamespace(**environment(mode))
                client.buildFor = lambda *a, **k: (True, 'built')
                client.benchmarkIdentity = lambda **k: {'packedDirectReduction': mode == '1'}
                client.volume = SimpleNamespace(commit=lambda: None)
                def run(command, **kwargs):
                    commands.append(command)
                    return SimpleNamespace(returncode=0, stdout='packed arithmetic direct reduction: '
                                           + str(1 - int(mode)) + '\n', stderr='')
                env = dict(client=client, re=re, json=json, time=time,
                           print=lambda *a, **k: None,
                           subprocess=SimpleNamespace(run=run),
                           Path=lambda path: Path(directory) / 'artifact')
                report = function('runAudit', env, 'packed_audit.py')()
                self.assertFalse(report['valid'])
                self.assertIn('reducer identity', report['error'])
                self.assertEqual(report['packedDirectReduction'], mode == '1')
                self.assertEqual(len(commands), 1)
                self.assertIn('PACKED_DIRECT_REDUCE=' + mode, commands[0])


if __name__ == '__main__':
    unittest.main()
