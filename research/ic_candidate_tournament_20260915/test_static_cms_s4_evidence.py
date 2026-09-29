"""Independent archive, formula, model and group-witness audit of six controls."""
import hashlib
import io
from itertools import product
import json
from pathlib import Path
import tarfile
import unittest

from generic_bases import lifts
from run_static_cms_s4_controls import PANEL, PANEL_SHA256, REGISTRATION, admit
from tournament import read


EVIDENCE = REGISTRATION/'evidence.tar.gz'
EVIDENCE_SHA256 = '3eb9c693bcc00519b4f94b39cab4b96adc5e8e67711a34e670e90f27e3106e97'
SUMMARY_SHA256 = '13a4b0203115756547236cea57b3abf760720109423af63d98ad07b89eff60a6'
RUNNER_SHA256 = '9b5cee5b4616f3ae7e7e57c88889e7b6512b7d14cc6798c4a6c71372a923da52'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def parse_complete_model(stdout, count):
    values = {}
    for line in stdout.decode().splitlines():
        if not line.startswith('v '):
            continue
        for token in line.split()[1:]:
            literal = int(token)
            if literal == 0:
                continue
            index = abs(literal)
            assert 1 <= index <= count
            bit = literal > 0
            assert index not in values or values[index] == bit
            values[index] = bit
    assert set(values) == set(range(1, count+1))
    return [values[index] for index in range(1, count+1)]


def independent_formula_check(data, assignment):
    lines = data.decode().splitlines()
    body = []
    header = None
    for line in lines:
        line = line.strip()
        if not line or line.startswith('c'):
            continue
        if line.startswith('p '):
            assert header is None
            parts = line.split()
            assert parts[:2] == ['p', 'cnf']
            header = int(parts[2]), int(parts[3])
            continue
        body.append(line)
    assert header == (len(assignment), len(body))
    for line in body:
        xor = line.startswith('x ')
        literals = [int(token) for token in line.split()[1 if xor else 0:]]
        assert literals[-1] == 0
        literals = literals[:-1]
        assert literals and all(1 <= abs(literal) <= header[0]
                                for literal in literals)
        truth = [(assignment[abs(literal)-1] if literal > 0
                  else not assignment[abs(literal)-1]) for literal in literals]
        assert (sum(truth) % 2 == 1) if xor else any(truth)


def independently_lift(curve, base, model, target):
    allowed = set(base)
    xs = [sum(1 << bit for bit in range(6) if model[6*slot+bit])
          for slot in range(3)]
    choices = [tuple(point for point in lifts(curve, x) if point in allowed)
               for x in xs]
    for triple in product(*choices):
        if curve.add(curve.add(triple[0], triple[1]), triple[2]) == target:
            return xs, [list(point) for point in triple]
    return xs, None


class StaticCmsS4EvidenceTests(unittest.TestCase):
    def test_six_source_systems_and_four_group_witnesses(self):
        blob = EVIDENCE.read_bytes()
        self.assertEqual(sha(blob), EVIDENCE_SHA256)
        with tarfile.open(fileobj=io.BytesIO(blob), mode='r:gz') as archive:
            members = archive.getmembers()
            names = [member.name for member in members]
            self.assertEqual(len(names), len(set(names)))
            self.assertTrue(all(member.isfile() and not member.name.startswith('/')
                                and '..' not in Path(member.name).parts
                                for member in members))
            files = {member.name: archive.extractfile(member).read()
                     for member in members}
        panel = read(PANEL)
        _, curve, base, _ = admit(panel, require_local_binary=False)
        self.assertEqual(sha(files['registered-panel.json']), PANEL_SHA256)
        self.assertEqual(sha(files['registered-runner.py']), RUNNER_SHA256)
        self.assertEqual(sha(files['summary.json']), SUMMARY_SHA256)
        self.assertEqual(sha(files['cms-executable']),
                         panel['cms_executable_sha256'])
        self.assertEqual(sha(files['cms-build-receipt.json']),
                         panel['cms_build_receipt_sha256'])
        self.assertEqual(sha(files['cms-build-bundle-seal.json']),
                         panel['cms_build_bundle_seal_sha256'])
        self.assertNotIn(b'@rpath', files['cms-linkage.txt'])
        preflight = json.loads(files['cms-preflight.metrics.json'])
        self.assertEqual(preflight['returncode'], 0)
        self.assertFalse(preflight['timed_out'])
        self.assertIn(b'CryptoMiniSat version 5.14.7',
                      files['cms-preflight.stdout'])
        result = json.loads(files['summary.json'])
        self.assertEqual(result['status'], 'CONTROL_PANEL_COMPLETE')
        self.assertEqual(result['panel_sha256'], PANEL_SHA256)
        self.assertEqual(result['cms_preflight'], preflight)
        self.assertEqual(result['cms_binary_sha256'],
                         panel['cms_executable_sha256'])
        self.assertFalse(result['full_sat_ic_admission'])
        self.assertIsNone(result['natural_yield_estimate'])
        self.assertIsNone(result['online_speedup'])
        self.assertEqual(len(result['rows']), len(panel['schedule']))
        self.assertEqual([json.loads(line) for line in
                          files['progress.jsonl'].decode().splitlines()],
                         result['rows'])
        for item, row in zip(panel['schedule'], result['rows']):
            trial = item['trial']
            prefix = f'trial-{trial:02d}/'
            instance = prefix+'instance/'
            self.assertEqual(row['trial'], trial)
            self.assertEqual(row['public_point'], item['point'])
            self.assertEqual(row['exact_relation_exists'],
                             item['exact_relation_exists'])
            self.assertEqual(row['exporter'],
                             json.loads(files[prefix+'export.metrics.json']))
            self.assertEqual(row['exporter']['returncode'], 0)
            self.assertFalse(row['exporter']['timed_out'])
            manifest = json.loads(files[instance+'manifest.json'])
            self.assertEqual(sha(files[instance+'manifest.json']),
                             row['manifest_sha256'])
            self.assertEqual(manifest['representation'], 'symmetrised_s4')
            self.assertEqual([int(manifest['target'][key])
                              for key in ('x', 'y')], item['point'])
            self.assertEqual(manifest['factor_base_geometry']['curve_points'],
                             len(base))
            self.assertEqual(manifest['direct_meet_in_the_middle']['status'],
                             'not_run_in_export_process')
            for name, descriptor in manifest['exports'].items():
                data = files[instance+descriptor['path']]
                self.assertEqual(len(data), descriptor['bytes'])
                self.assertEqual(row['exports'][name],
                                 {'bytes': len(data), 'sha256': sha(data)})
            cms = json.loads(files[prefix+'cms.metrics.json'])
            self.assertEqual(row['cms'], cms)
            self.assertFalse(cms['timed_out'])
            stdout = files[prefix+'cms.stdout']
            if item['exact_relation_exists']:
                self.assertEqual(row['status'], 'VALID_POINT_WITNESS')
                self.assertEqual(cms['returncode'], 10)
                self.assertIn(b's SATISFIABLE', stdout)
                count = manifest['exports']['cryptominisat_xor_dimacs']['variables']
                model = parse_complete_model(stdout, count)
                independent_formula_check(files[instance+'instance.xor.cnf'],
                                          model)
                self.assertEqual(row['source_model_sha256'], sha(bytes(model)))
                xs, witness = independently_lift(curve, base, model,
                                                  tuple(item['point']))
                self.assertIsNotNone(witness)
                self.assertEqual(row['point_witness']['x_coordinates'], xs)
                self.assertEqual(row['point_witness']['points'], witness)
                self.assertTrue(row['point_witness']['group_replay'])
            else:
                self.assertEqual(row['status'], 'SOURCE_UNSAT')
                self.assertEqual(cms['returncode'], 20)
                self.assertIn(b's UNSATISFIABLE', stdout)
                self.assertIsNone(row['point_witness'])
                self.assertIsNone(row['source_model_valid'])


if __name__ == '__main__':
    unittest.main()
