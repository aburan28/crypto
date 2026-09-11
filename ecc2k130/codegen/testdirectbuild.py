"""Exercise packed reducer/product build and result helpers without Modal or a GPU."""
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


def environment(mode='0', generated='0', tile='0', clmad='0'):
    env = dict(re=re, hashlib=hashlib, pathlib=SimpleNamespace(Path=Path),
               subprocess=subprocess, time=time, json=json, benchResult=benchResult,
               summarizeSamples=summarizeSamples, bestResult=bestResult,
               CUDA_VERSION='13.0.0', DEFAULT_GPU='RTX-PRO-6000', BAKED_ARCHES=('120',), REMOTE='/unused',
               print=lambda *a, **k: None,
               bakedIntact=[True], computeCapability=lambda: '120', gpuName=lambda: 'fixture GPU')
    for key in ('SINGLE_PRODUCT', 'CACHE_DENOM', 'BY_VALUE', 'POLY_CHAIN',
                'UNROLL_INV', 'PAIR_PRODUCTS', 'POLY_STATE'):
        env['PACKED_' + key] = '1'
    env.update(PACKED_PERM_SIGMA='3', PACKED_DIRECT_REDUCE=mode, PACKED_GENERATED_PRODUCT=generated, PACKED_STATE_TILE=tile, PACKED_CLMAD=clmad)
    assignment('BAKED', env)
    return env


def raw(mode='0', generated='0', tile='0', clmad='0'):
    return (f'packed direct reduction: {mode}\n'
            f'packed generated product: {generated}\n'
            f'packed state tile: {tile}\n'
            f'packed native carryless multiply: {clmad}\n'
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


class GeneratedProductBuildTests(unittest.TestCase):
    def test_environment_defaults_invalid_values_and_reducer_prerequisite(self):
        body = nodes('modal_app.py')
        index = next(i for i, n in enumerate(body) if isinstance(n, ast.Assign)
                     and any(isinstance(t, ast.Name) and t.id == 'PACKED_GENERATED_PRODUCT' for t in n.targets))
        for direct in ('0', '1'):
            for value in (None, '0', '1', '', '2', '-1', 'true'):
                env = dict(PACKED_DIRECT_REDUCE=direct,
                           os=SimpleNamespace(environ={} if value is None else {'ECC_PACKED_GENERATED_PRODUCT': value}))
                with self.subTest(value=value, direct=direct):
                    if value in (None, '0') or value == '1' and direct == '1':
                        execute(body[index:index + 3], env)
                        self.assertEqual(env['PACKED_GENERATED_PRODUCT'], value or '0')
                    else:
                        with self.assertRaises(ValueError):
                            execute(body[index:index + 3], env)

    def test_image_preserves_generated_mode_and_default_cuda_version(self):
        class Image:
            def __init__(self):
                self.calls = {}
            def __getattr__(self, name):
                def record(*args, **kwargs):
                    self.calls[name] = args
                    return self
                return record
        for mode in ('0', '1'):
            env = environment('1', mode)
            image = Image()
            env.update(modal=SimpleNamespace(Image=SimpleNamespace(from_registry=lambda *a, **k: image)),
                       GENCODE='fixture', LOCAL=ROOT)
            assignment('image', env)
            self.assertEqual(image.calls['env'][0]['ECC_PACKED_GENERATED_PRODUCT'], mode)
            self.assertEqual(image.calls['env'][0]['ECC_CUDA_VERSION'], env['CUDA_VERSION'])
            builds = [c for c in image.calls['run_commands'] if 'make gpu ' in c]
            self.assertEqual(len(builds), 1)
            self.assertIn('PACKED_GENERATED_PRODUCT=' + mode, builds[0])
            self.assertEqual(env['BAKED']['packedGeneratedProduct'], mode == '1')
        env = dict(os=SimpleNamespace(environ={}))
        assignment('CUDA_VERSION', env)
        self.assertEqual(env['CUDA_VERSION'], '12.8.1')

    def test_mode_change_rebuilds_and_cannot_reuse_a_stale_bake(self):
        env = environment('1', '0')
        commands = []
        env.update(sh=lambda *a, **k: (0, ''),
                   shStream=lambda c, **k: (commands.append(c) or 0, 'built'))
        build = function('buildFor', env)
        self.assertTrue(build(32, 128, 0)[0])
        self.assertEqual(commands, [])
        env['PACKED_GENERATED_PRODUCT'] = '1'
        with self.assertRaisesRegex(ValueError, 'matching baked binary'):
            function('runBench', env)(rebuild=False, packed=True)
        self.assertTrue(build(32, 128, 0)[0])
        self.assertIn('PACKED_GENERATED_PRODUCT=1', commands[-1])
        env['PACKED_GENERATED_PRODUCT'] = '0'
        self.assertTrue(build(32, 128, 0)[0])
        self.assertIn('PACKED_GENERATED_PRODUCT=0', commands[-1])
        self.assertEqual(len(commands), 2)
        with self.assertRaisesRegex(ValueError, 'matching baked binary'):
            function('runBench', env)(rebuild=False, packed=True)

    def test_actual_runtime_flag_is_required_for_every_packed_sample(self):
        for mode in ('0', '1'):
            env = environment('1', mode)
            check = function('checkPackedReduction', env)
            good = benchResult('fixture', 0, raw('1', mode))
            self.assertTrue(check(good))
            self.assertEqual(good['expectedPackedGeneratedProduct'], mode == '1')
            self.assertEqual(good['packedGeneratedProduct'], mode == '1')
            bad = [raw('1', str(1 - int(mode))),
                   raw('1', mode).replace(f'packed generated product: {mode}\n', ''),
                   raw('1', mode) + f'packed generated product: {mode}\n', raw('1', '2')]
            for text in bad:
                sample = benchResult('fixture', 0, text)
                self.assertFalse(check(sample))
                self.assertEqual(sample['rate'], 0)
                failed = summarizeSamples([good, sample])
                self.assertFalse(failed['valid'])
                self.assertEqual(failed['rate'], 0)
                self.assertIsNone(bestResult([failed]))
            env['sh'] = lambda *a, **k: (0, bad[0])
            measure = function('measureBench', env)
            self.assertFalse(measure(1, 1, 0, False, 1, packed=True)['valid'])
            self.assertTrue(measure(1, 1, 0, False, 1, packed=False)['valid'])

    def test_source_digest_covers_generator_and_generated_header(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, text in [('Makefile', ''), ('modal_app.py', ''), ('ecc2k130', 'fixture'),
                               ('generated/eccF131.h', 'LEAF = 66'),
                               ('codegen/genpackedproduct.py', 'generator'),
                               ('include/packedgeneratedproduct131.h', 'header')]:
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text)
            for mode in ('0', '1'):
                env = environment('1', mode)
                env.update(REMOTE=directory, sh=lambda *a, **k: (0, 'fixture'),
                           buildFor=lambda *a, **k: (True, 'built'))
                identity = function('benchmarkIdentity', env)
                self.assertEqual(identity(True)['packedGeneratedProduct'], mode == '1')
                self.assertIsNone(identity(False)['packedGeneratedProduct'])
                self.assertEqual(function('runCompileCheck', env)()['packedGeneratedProduct'], mode == '1')
                for name in ('codegen/genpackedproduct.py', 'include/packedgeneratedproduct131.h'):
                    before = identity(True)['sourceSha256']
                    path = root / name
                    path.write_text(path.read_text() + '\nchanged')
                    self.assertNotEqual(identity(True)['sourceSha256'], before)

    def test_benchmark_and_autotune_metadata_preserve_generated_mode(self):
        for mode in ('0', '1'):
            env = environment('1', mode)
            env.update(buildFor=lambda *a, **k: (True, 'built'),
                       benchmarkIdentity=lambda *a, **k: {},
                       measureBench=lambda *a, **k: dict(valid=True, rate=6000.0),
                       os=SimpleNamespace(makedirs=lambda *a, **k: None),
                       open=lambda *a, **k: io.StringIO(),
                       volume=SimpleNamespace(commit=lambda: None))
            self.assertEqual(function('runBench', env)(rebuild=False, packed=True)['packedGeneratedProduct'], mode == '1')
            function('autotuneConfigs', env)
            tune = function('runAutotune', env)(configs='0:32:128:2', packed=True)
            self.assertEqual(tune['results'][0]['packedGeneratedProduct'], mode == '1')

    def audit_fixture(self, mode, arithmetic=None, failed_phase=None, failed_marker=None,
                      clmad='0', clmad_arithmetic=None):
        with tempfile.TemporaryDirectory() as directory:
            commands = []
            env = environment('1', mode, clmad=clmad)
            workers = 192512
            def output(weight, phase):
                marker = mode if failed_phase != phase else failed_marker
                text = (f'packed direct reduction: 1\n'
                        f'packed generated product: {marker}\n'
                        'packed state tile: 0\n'
                        f'packed native carryless multiply: {clmad}\n'
                        f'backend cuda-packed131: {workers} threads x 32 slots x 1 lanes = {workers * 32} walks, '
                        f'dp weight {weight}, 1024 steps per launch\n'
                        '1.0 s 6000.000 M it/s 201863462912 iterations 1 dp 1 stored 0 dropped\n'
                        'finished: 6000.000 M it/s, 1 distinguished points (0 verified against the reference, 0 dropped)\n')
                return text.replace('packed generated product: None\n', '')
            env['sh'] = lambda *a, **k: (0, output(0, 'benchmark'))
            function('checkPackedReduction', env)
            function('measureBench', env)
            client = SimpleNamespace(**env)
            client.REMOTE = directory
            client.buildFor = lambda *a, **k: (True, 'built')
            client.benchmarkIdentity = lambda **k: {'packedGeneratedProduct': mode == '1'}
            client.volume = SimpleNamespace(commit=lambda: None)
            def run(command, **kwargs):
                commands.append(command)
                if command[0] == 'make':
                    text = ('packed arithmetic direct reduction: 1\n'
                             + (f'packed arithmetic native carryless multiply: {clmad}\n' if clmad_arithmetic is None else clmad_arithmetic)
                            + (f'packed arithmetic generated product: {mode}\n' if arithmetic is None else arithmetic))
                elif command[0] == 'python3':
                    text = 'integration passed\n'
                else:
                    Path(command[command.index('--dp-file') + 1]).write_bytes(bytes(32))
                    text = output(34, 'collection')
                return SimpleNamespace(returncode=0, stdout=text, stderr='')
            def path(value):
                return Path(directory) / 'artifact' if str(value).startswith('/data/') else Path(value)
            audit_env = dict(client=client, re=re, json=json, time=time, tempfile=tempfile,
                             print=lambda *a, **k: None, subprocess=SimpleNamespace(run=run), Path=path)
            function('checkScalarCounts', audit_env, 'packed_audit.py')
            result = function('runAudit', audit_env, 'packed_audit.py')(repeats=1, workers=workers)
            return result, commands

    def test_native_audit_forwards_flag_and_records_expected_and_actual_modes(self):
        for mode in ('0', '1'):
            result, commands = self.audit_fixture(mode)
            self.assertTrue(result['valid'], result.get('error'))
            self.assertIn('PACKED_GENERATED_PRODUCT=' + mode, commands[0])
            for row in [result, result['deviceArithmetic'], *result['benchmark']['samples'], *result['collection']]:
                self.assertEqual(row['expectedPackedGeneratedProduct'], mode == '1')
                self.assertEqual(row['packedGeneratedProduct'], mode == '1')
            self.assertEqual(result['collection'][0]['corpusBytes'], 32)
            self.assertEqual(result['benchmark']['samples'][0]['reportedIterations'], 201863462912)

    def test_arithmetic_identity_failure_stops_before_integration_or_timing(self):
        for mode in ('0', '1'):
            for marker in ('', f'packed arithmetic generated product: {1 - int(mode)}\n',
                           f'packed arithmetic generated product: {mode}\n' * 2,
                           'packed arithmetic generated product: 2\n'):
                result, commands = self.audit_fixture(mode, arithmetic=marker)
                self.assertFalse(result['valid'])
                self.assertIn('generated product identity', result['error'])
                self.assertEqual(len(commands), 1)
                self.assertNotIn('integration', result)

    def test_native_benchmark_and_collection_reject_stale_or_missing_markers(self):
        for mode in ('0', '1'):
            for phase in ('benchmark', 'collection'):
                for marker in (str(1 - int(mode)), None):
                    with self.subTest(mode=mode, phase=phase, marker=marker):
                        result, _ = self.audit_fixture(mode, failed_phase=phase, failed_marker=marker)
                        self.assertFalse(result['valid'])
                        row = (result['benchmark']['samples'] if phase == 'benchmark' else result['collection'])[0]
                        self.assertFalse(row['valid'])
                        self.assertEqual(row['rate'], 0)


class StateTileBuildTests(unittest.TestCase):
    def test_environment_and_dependencies(self):
        body = nodes('modal_app.py')
        index = next(i for i, node in enumerate(body) if isinstance(node, ast.Assign)
                     and any(isinstance(t, ast.Name) and t.id == 'PACKED_STATE_TILE' for t in node.targets))
        for value in (None, '0', '256', '', '1', '128', '0256', 'true'):
            for missing in (None, 'POLY_STATE', 'CACHE_DENOM', 'POLY_CHAIN'):
                values = {} if value is None else {'ECC_PACKED_STATE_TILE': value}
                env = dict(os=SimpleNamespace(environ=values), PACKED_POLY_STATE='1',
                           PACKED_CACHE_DENOM='1', PACKED_POLY_CHAIN='1')
                if missing:
                    env['PACKED_' + missing] = '0'
                valid = value in (None, '0') or (value == '256' and missing is None)
                with self.subTest(value=value, missing=missing):
                    if valid:
                        execute(body[index:index + 3], env)
                        self.assertEqual(env['PACKED_STATE_TILE'], value or '0')
                    else:
                        with self.assertRaises(ValueError):
                            execute(body[index:index + 3], env)

    def test_baked_image_geometry_matches_layout(self):
        class Image:
            def __init__(self):
                self.calls = {}
            def __getattr__(self, name):
                def record(*args, **kwargs):
                    self.calls[name] = args
                    return self
                return record
        for tile, threads in (('0', 128), ('256', 256)):
            env = environment('1', '1', tile)
            image = Image()
            env.update(modal=SimpleNamespace(Image=SimpleNamespace(from_registry=lambda *a, **k: image)),
                       GENCODE='fixture', LOCAL=ROOT)
            assignment('image', env)
            builds = [s for s in image.calls['run_commands'] if 'make gpu ' in s]
            self.assertEqual(len(builds), 1)
            self.assertEqual(image.calls['env'][0]['ECC_PACKED_STATE_TILE'], tile)
            self.assertIn('THREADS=' + str(threads), builds[0])
            self.assertIn('PACKED_STATE_TILE=' + tile, builds[0])
            self.assertEqual(env['BAKED']['threads'], threads)
            self.assertEqual(env['BAKED']['packedStateTile'], int(tile))

    def test_rebuild_cannot_reuse_another_layout(self):
        env = environment('1', '1', '0')
        commands = []
        env.update(sh=lambda *a, **k: (0, ''),
                   shStream=lambda command, **kwargs: (commands.append(command) or 0, ''))
        build = function('buildFor', env)
        self.assertTrue(build(32, 128, 0, arch='120')[0])
        self.assertEqual(commands, [])
        env['PACKED_STATE_TILE'] = '256'
        self.assertTrue(build(32, 256, 0, arch='120')[0])
        self.assertIn('PACKED_STATE_TILE=256', commands[-1])
        self.assertFalse(env['bakedIntact'][0])

    def test_incompatible_block_size_stops_before_build(self):
        env = environment('1', '1', '256')
        env['sh'] = env['shStream'] = lambda *a, **k: self.fail('invalid layout must not invoke compiler')
        ok, error = function('buildFor', env)(32, 128, 0, arch='120')
        self.assertFalse(ok)
        self.assertIn('256 threads per block', error)

    def test_missing_wrong_and_duplicate_tile_cannot_rank(self):
        for tile in ('0', '256'):
            env = environment('1', '1', tile)
            check = function('checkPackedReduction', env)
            text = raw('1', '1', tile)
            good = benchResult(['fixture'], 0, text)
            self.assertTrue(check(good))
            self.assertEqual(good['packedStateTile'], int(tile))
            for bad in (text.replace(f'packed state tile: {tile}\n', ''),
                        raw('1', '1', '256' if tile == '0' else '0'),
                        text + f'packed state tile: {tile}\n',
                        text + 'packed state tile: true\n',
                        text.replace(f'packed state tile: {tile}', 'packed state tile: 1')):
                sample = benchResult(['fixture'], 0, bad)
                self.assertFalse(check(sample))
                self.assertEqual(sample['rate'], 0)
                self.assertFalse(summarizeSamples([good, sample])['valid'])


class ClmadBuildTests(unittest.TestCase):
    def test_environment_defaults_and_rejects_invalid_values(self):
        body = nodes('modal_app.py')
        index = next(i for i, node in enumerate(body) if isinstance(node, ast.Assign)
                     and any(isinstance(t, ast.Name) and t.id == 'PACKED_CLMAD' for t in node.targets))
        for value in (None, '0', '1', '', '2', '-1', 'true'):
            env = dict(os=SimpleNamespace(environ={} if value is None else {'ECC_PACKED_CLMAD': value}))
            with self.subTest(value=value):
                if value in (None, '0', '1'):
                    execute(body[index:index + 2], env)
                    self.assertEqual(env['PACKED_CLMAD'], value or '0')
                else:
                    with self.assertRaises(ValueError):
                        execute(body[index:index + 2], env)

    def test_image_and_cache_preserve_native_mode(self):
        class Image:
            def __init__(self): self.calls = {}
            def __getattr__(self, name):
                def record(*args, **kwargs):
                    self.calls[name] = args
                    return self
                return record
        for mode in ('0', '1'):
            env = environment('1', '1', clmad=mode)
            image = Image()
            env.update(CUDA_VERSION='13.3.1', GENCODE='fixture', LOCAL=ROOT,
                       modal=SimpleNamespace(Image=SimpleNamespace(from_registry=lambda *a, **k: image)))
            assignment('image', env)
            self.assertEqual(image.calls['env'][0]['ECC_PACKED_CLMAD'], mode)
            builds = [c for c in image.calls['run_commands'] if 'make gpu ' in c]
            self.assertEqual(len(builds), 1)
            self.assertIn('PACKED_CLMAD=' + mode, builds[0])
            self.assertEqual(env['BAKED']['packedClmad'], mode == '1')
        env = environment('1', '1')
        commands = []
        env.update(sh=lambda *a, **k: (0, ''),
                   shStream=lambda c, **k: (commands.append(c) or 0, 'built'))
        build = function('buildFor', env)
        self.assertTrue(build(32, 128, 0)[0])
        self.assertEqual(commands, [])
        env['PACKED_CLMAD'] = '1'
        with self.assertRaisesRegex(ValueError, 'matching baked binary'):
            function('runBench', env)(rebuild=False, packed=True)
        self.assertTrue(build(32, 128, 0)[0])
        self.assertIn('PACKED_CLMAD=1', commands[-1])
        env['PACKED_CLMAD'] = '0'
        self.assertTrue(build(32, 128, 0)[0])
        self.assertIn('PACKED_CLMAD=0', commands[-1])

    def test_missing_wrong_or_duplicate_identity_cannot_rank(self):
        for mode in ('0', '1'):
            check = function('checkPackedReduction', environment('1', '1', clmad=mode))
            text = raw('1', '1', clmad=mode)
            good = benchResult('fixture', 0, text)
            self.assertTrue(check(good))
            self.assertEqual(good['packedClmad'], mode == '1')
            marker = f'packed native carryless multiply: {mode}\n'
            for bad in (text.replace(marker, ''), text + marker,
                        raw('1', '1', clmad=str(1 - int(mode))), raw('1', '1', clmad='2')):
                sample = benchResult('fixture', 0, bad)
                self.assertFalse(check(sample))
                self.assertEqual(sample['rate'], 0)
                self.assertFalse(summarizeSamples([good, sample])['valid'])

    def test_benchmark_and_autotune_record_native_mode(self):
        for mode in ('0', '1'):
            env = environment('1', '1', clmad=mode)
            env.update(buildFor=lambda *a, **k: (True, 'built'),
                       benchmarkIdentity=lambda *a, **k: {},
                       measureBench=lambda *a, **k: dict(valid=True, rate=6000.0),
                       os=SimpleNamespace(makedirs=lambda *a, **k: None),
                       open=lambda *a, **k: io.StringIO(), volume=SimpleNamespace(commit=lambda: None))
            self.assertEqual(function('runBench', env)(rebuild=False, packed=True)['packedClmad'], mode == '1')
            function('autotuneConfigs', env)
            tune = function('runAutotune', env)(configs='0:32:128:2', packed=True)
            self.assertEqual(tune['results'][0]['packedClmad'], mode == '1')

    def test_audit_binds_arithmetic_and_timed_modes(self):
        fixture = GeneratedProductBuildTests()
        for mode in ('0', '1'):
            result, commands = fixture.audit_fixture('1', clmad=mode)
            self.assertTrue(result['valid'], result.get('error'))
            self.assertIn('PACKED_CLMAD=' + mode, commands[0])
            for row in [result, result['deviceArithmetic'], *result['benchmark']['samples'], *result['collection']]:
                self.assertEqual(row['expectedPackedClmad'], mode == '1')
                self.assertEqual(row['packedClmad'], mode == '1')
            for marker in ('', f'packed arithmetic native carryless multiply: {1 - int(mode)}\n',
                           f'packed arithmetic native carryless multiply: {mode}\n' * 2):
                failed, commands = fixture.audit_fixture('1', clmad=mode, clmad_arithmetic=marker)
                self.assertFalse(failed['valid'])
                self.assertIn('CLMAD identity', failed['error'])
                self.assertEqual(len(commands), 1)
                self.assertNotIn('integration', failed)


if __name__ == '__main__':
    unittest.main()
