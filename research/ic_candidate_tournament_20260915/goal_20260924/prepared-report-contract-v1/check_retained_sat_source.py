#!/usr/bin/env python3
"""Substitute the one retained disclosed group witness into its original sources.

No SAT search, native execution, new target or modification of historical output.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[3]/'scripts'))
from identity import sha256, write_immutable  # noqa: E402
from oracle import Curve, require  # noqa: E402
from run_cms_s4_controls import lift_source_assignment  # noqa: E402
from run_koblitz_pdp_matrix import validate_wdsat_anf, validate_xor_dimacs  # noqa: E402

ARCHIVE_SHA256 = 'c94aba5c67afbe60b5109169d2d57c1c44ce89b218c63d7b023ed023bdb64b24'
EXECUTION_SHA256 = '43539f7d440289dae1ad4867ba1ca951bd4664bf01147b9070eaa68f3bad07cb'
PREFIX = 'sat-execution/entry-output/target/trial-01/instance/'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def check(bundle, out):
    bundle, out = Path(bundle), Path(out)
    require(not out.exists(), 'source-witness output already exists')
    data = (bundle/'evidence.tar.gz').read_bytes()
    receipt = json.loads((bundle/'receipt.json').read_text())
    require(digest(data) == receipt['archive_sha256'] == ARCHIVE_SHA256
            and len(data) == receipt['archive_bytes'], 'original published archive changed')
    inventory = {row['role']: row for row in receipt['inventory']}
    roles = ['sat-execution/execution.json', 'sat-execution/entry-output/summary.json',
             'sat-audit/admission.json', 'diagnosis/diagnosis.json']
    roles += [PREFIX+name for name in ('manifest.json', 'instance.anf', 'instance.xor.cnf', 'instance.magma')]
    files = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as archive:
        for role in roles:
            member = archive.getmember(role)
            require(member.isfile(), 'source member is not a regular file')
            value = archive.extractfile(member).read()
            require(len(value) == inventory[role]['bytes']
                    and digest(value) == inventory[role]['sha256'], 'source member changed')
            files[role] = value
    spec = json.loads(files['sat-execution/execution.json'])
    require(sha256(spec) == EXECUTION_SHA256, 'source control differs from sealed SAT invocation')
    summary = json.loads(files['sat-execution/entry-output/summary.json'])
    row = summary['target_attempts'][1]
    admission = json.loads(files['sat-audit/admission.json'])
    require(row['trial'] == 1 and row['status'] == 'CONFLICT_BUDGET_INCONCLUSIVE'
            and admission['status'] == 'ADMITTED_INCOMPLETE_PREPARED_SAT_CONTROL'
            and not admission['scalar_verified'], 'original incomplete outcome changed')
    inputs = spec['arguments']['preparation']['certificate']['inputs']
    curve = Curve(inputs['fixture'])
    base = tuple(curve.decode(p) for p in inputs['base'])
    diagnosis = json.loads(files['diagnosis/diagnosis.json'])
    indices = diagnosis['arms']['sat']['attempts'][1]['independently_readded_witness']
    require(curve.n == 17 and len(base) == 63 and indices == [29, 51, 2],
            'analysis must use the exact previously retained toy witness')
    target = curve.add(curve.mul(curve.g, row['a']),
                       curve.mul(curve.decode(summary['target_input']), row['b']))
    points = [base[index] for index in indices]
    require(target == (62577, 27783) and curve.add(curve.add(points[0], points[1]), points[2]) == target,
            'retained group witness no longer readds')
    manifest = json.loads(files[PREFIX+'manifest.json'])
    require(manifest['n'] == 17 and manifest['ell'] == 6 and manifest['m'] == 3
            and manifest['representation'] == 'symmetrised_s4'
            and manifest['source_variables'] == 51 and manifest['source_equations'] == 50
            and manifest['factor_base_basis_bitmasks'] == ['1', '2', '4', '8', '16', '32']
            and [int(manifest['target'][k]) for k in ('x', 'y')] == list(target)
            and manifest['exports']['cryptominisat_xor_dimacs']['variables'] == 767,
            'original symmetrised source layout differs')
    for name, info in manifest['exports'].items():
        raw = files[PREFIX+info['path']]
        require(len(raw) == info['bytes'] == row['exports'][name]['bytes']
                and digest(raw) == row['exports'][name]['sha256'], 'original source export differs')
    xs = [p[0] for p in points]
    symmetric = [xs[0]^xs[1]^xs[2],
                 curve.fm(xs[0], xs[1])^curve.fm(xs[0], xs[2])^curve.fm(xs[1], xs[2]),
                 curve.fm(curve.fm(xs[0], xs[1]), xs[2])]
    source_bits = []
    for value, width in zip(xs+symmetric, (6, 6, 6, 6, 11, 16)):
        require(0 <= value < 1 << width, 'symmetric coefficient overflow')
        source_bits.extend(bool(value & (1 << bit)) for bit in range(width))
    require(len(source_bits) == 51, 'source assignment width differs')
    model = source_bits+[False]*(767-51)
    definitions = {}
    cnf = files[PREFIX+'instance.xor.cnf'].decode().splitlines()
    require(cnf[0] == 'p cnf 767 2414', 'CNF header differs')
    ordinary, xor = 0, 0
    for line in cnf[1:]:
        if not line.strip() or line.startswith('c'):
            continue
        if line.startswith('x'):
            xor += 1
            continue
        ordinary += 1
        terms = [int(v) for v in line.split()]
        require(terms[-1] == 0 and all(1 <= abs(v) <= 767 for v in terms[:-1]),
                'malformed source clause')
        terms = terms[:-1]
        if len(terms) >= 3 and terms[-1] > 51 and all(-51 <= v < 0 for v in terms[:-1]):
            output, arguments = terms[-1], [-v for v in terms[:-1]]
            require(output not in definitions, 'duplicate auxiliary AND definition')
            definitions[output] = arguments
            model[output-1] = all(source_bits[v-1] for v in arguments)
    require(set(definitions) == set(range(52, 768)) and (ordinary, xor) == (2364, 50),
            'source monomial AND gates or row counts differ')
    out.mkdir(parents=True)
    for name in ('manifest.json', 'instance.anf', 'instance.xor.cnf', 'instance.magma'):
        with (out/name).open('xb') as stream:
            stream.write(files[PREFIX+name])
    require(files[PREFIX+'instance.anf'].decode().splitlines()[0] == 'p cnf 51 50',
            'ANF header differs')
    anf_valid = validate_wdsat_anf(out/'instance.anf', source_bits)
    cnf_valid = validate_xor_dimacs(out/'instance.xor.cnf', model)
    lifted = lift_source_assignment(model, manifest, curve, base, target)
    result = dict(schema_version=1, status=('PASS_RETAINED_LIFTING_SOURCE_MODEL'
                  if anf_valid and cnf_valid and lifted['group_replay'] else 'REJECTED_RETAINED_SOURCE_SUBSTITUTION'),
                  archive_sha256=ARCHIVE_SHA256, execution_sha256=EXECUTION_SHA256,
                  script_sha256=digest(Path(__file__).read_bytes()), trial=1,
                  query_point=list(target), group_witness_indices=indices,
                  source_variable_count=51, auxiliary_and_count=len(definitions),
                  anf_equations=50, cnf_clauses=ordinary, xor_rows=xor,
                  source_anf_valid=anf_valid, cnf_xor_valid=cnf_valid, lifted_witness=lifted,
                  source_assignment=source_bits, cnf_assignment=model,
                  source_assignment_sha256=digest(bytes(source_bits)),
                  cnf_assignment_sha256=digest(bytes(model)),
                  selected_original_files={role:dict(bytes=len(value), sha256=digest(value))
                                           for role, value in files.items()},
                  original_native_outcome=row['status'], original_ic_target_complete=False,
                  native_solvers_executed=0, fresh_targets_generated=0,
                  source_bound_scientific_runtime_admitted=False, online_speedup=None,
                  promotion_eligible=False,
                  scope='postexecution substitution of one disclosed witness; no solver search, new query or rate estimate')
    write_immutable(out/'result.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, default=HERE.parent/'prepared-one-target-controls-v1/publication')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = check(args.bundle, args.out)
    print(json.dumps({key:result[key] for key in ('status', 'source_anf_valid', 'cnf_xor_valid', 'native_solvers_executed')}))
