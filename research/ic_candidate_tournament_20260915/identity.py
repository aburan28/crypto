"""Canonical IC identities and independent factor-base inventory.

Labels are derived from records, never accepted as evidence of their contents.
This adapter currently admits the binary Koblitz fixtures checked by oracle.py.
Other curve families and verified isogeny routes need their own checked adapters.
"""
import copy
import hashlib
import json
import re
from pathlib import Path

from oracle import Curve, require


def canonical(value):
    """Stable, lossless JSON: in particular bool is not an integer measurement."""
    def check(x):
        require(x is None or type(x) in (bool, int, str, list, dict),
                'canonical records cannot contain floats or non-JSON values')
        if type(x) is dict:
            require(all(type(k) is str for k in x), 'non-string JSON key')
            for v in x.values():
                check(v)
        elif type(x) is list:
            for v in x:
                check(v)
    check(value)
    return json.dumps(value, sort_keys=True, separators=(',', ':'),
                      ensure_ascii=False, allow_nan=False).encode('utf-8')


def sha256(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def natural(value, name, *, positive=False):
    require(type(value) is int and value >= int(positive), 'invalid ' + name)
    return value


def digest(value, name):
    require(type(value) is str and re.fullmatch('[0-9a-f]{64}', value), 'invalid ' + name)
    return value


def fields(value, names, name):
    require(type(value) is dict and set(value) == set(names.split()),
            'missing or extra fields in ' + name)


def identity_only(value):
    """Reject run/provenance fields instead of silently dropping them from hashes."""
    forbidden = {'candidate_id', 'workload_id', 'run_id', 'algorithm_seed', 'run_seed',
                 'target_seeds', 'timestamp', 'created_at', 'measurements', 'wall_ns',
                 'total_operations', 'source_directory', 'source_root', 'path', 'paths'}
    if type(value) is dict:
        for k, v in value.items():
            require(k not in forbidden and not k.endswith(('_path', '_paths')),
                    'run data or path in candidate identity: ' + k)
            identity_only(v)
    elif type(value) is list:
        for v in value:
            identity_only(v)
    elif type(value) is str:
        require(not value.startswith(('/', 'file://')),
                'absolute path in candidate identity')


def curve_record(fixture):
    c = Curve(fixture)  # Independent group, subgroup, generator and Frobenius checks.
    record = {
        'field': {'p': 2, 'n': c.n, 'basis': 'polynomial',
                  'modulus': c.modulus, 'element_encoding': 'unsigned-integer-polynomial-bits'},
        'curve': {'tag': 'kb' + str(c.a), 'model': 'y^2+x*y=x^3+a2*x^2+a6',
                  'coefficients': {'a1': 1, 'a2': c.a, 'a3': 0, 'a4': 0, 'a6': 1},
                  'order': c.h*c.r, 'trace': (1 << c.n) + 1 - c.h*c.r,
                  'r': c.r, 'cofactor': c.h, 'generator': list(c.g),
                  'target_group': 'prime-order-subgroup-generated-by-G'}}
    record['curve']['curve_id'] = f'EC1N{c.n}Ckb{c.a}h{sha256(record)[:12]}'
    return record


def base_points(report, c):
    require(('factor_base' in report) != ('factor_base_orbits' in report),
            'supply exactly one factor-base representation')
    if 'factor_base' in report:
        points = [c.decode(p) for p in report['factor_base']]
    else:
        points = []
        for encoded in report['factor_base_orbits']:
            p = start = c.decode(encoded)
            require(p is not None, 'identity orbit representative')
            for _ in range(c.n):
                points.extend((p, c.neg(p)))
                p = c.frob(p)
            require(p == start, 'Frobenius orbit does not close')
    require(points and None not in points and len(set(points)) == len(points),
            'empty, duplicate or identity geometric base')
    return points


def subgroup_orbits(points, c):
    """Partition by the actual group action, including partially supplied orbits.

    On these F2-defined curves Frobenius and negation preserve r-torsion.
    Check one point of each regenerated orbit instead of repeating a scalar
    multiplication and orbit expansion for all 2*n images of the same point.
    """
    remaining = set(points)
    representatives = set()
    while remaining:
        p = min(remaining)
        require(c.mul(p, c.r) is None, 'non-subgroup factor-base point')
        orbit = set()
        for _ in range(c.n):
            orbit.update((p, c.neg(p)))
            p = c.frob(p)
        representatives.add(min(orbit))
        remaining.difference_update(orbit)
    return representatives


def factor_base_inventory(report, fixture):
    """Count distinct usable points BEFORE folding; keep raw geometry separately.

    The cofactor convention explicitly defines the usable base as the image
    [h]B, not as the raw lifted points. Image collisions and killed torsion
    are retained, so neither geometric size nor nominal size can become fb<B>.
    """
    c = Curve(fixture)
    points = base_points(report, c)
    convention = report.get('column_convention', 'cofactor')
    require(convention in ('representative', 'cofactor'), 'unknown column convention')
    projected = [c.mul(p, c.h) if convention == 'cofactor' else p for p in points]
    usable = set(projected) - {None}
    require(usable, 'base has no usable points')
    representatives = subgroup_orbits(usable, c)
    columns = natural(report['columns'], 'column count', positive=True)
    require(columns == len(representatives), 'claimed folded column count differs from actual orbits')
    if 'column_logs' in report:
        column_points = [c.decode(entry['point']) for entry in report['column_logs']]
        require(None not in column_points, 'identity column representative')
        column_reps = subgroup_orbits(column_points, c)
        require(len(report['column_logs']) == columns and column_reps == representatives,
                'column representatives do not cover the usable base exactly')
    encode = lambda ps: [list(p) for p in ps]
    return {'geometric_point_count': len(points),
            'geometric_set_sha256': sha256(encode(sorted(points))),
            'geometric_order_sha256': sha256(encode(points)),
            'usable_point_count': len(usable),
            'usable_set_sha256': sha256(encode(sorted(usable))),
            'subgroup_map': 'cofactor-multiplication' if convention == 'cofactor' else 'identity',
            'identity_images': projected.count(None),
            'duplicate_nonidentity_images': len(projected)-projected.count(None)-len(usable),
            'quotient': 'sign-and-Frobenius', 'effective_columns': columns}


STAGE_FIELDS = {
    'point_decomposition': 'summands solver summation_polynomial encoding equation_order '
                           'monomial_order internal_matrix_kernel limits cache_policy source_sha256',
    'relation_collection': 'collector query_distribution query_rule filtering verification '
                           'duplicates dependencies stop_rule source_sha256',
    'relation_linear_algebra': 'solver modulus matrix_construction orbit_quotient rank_criterion '
                               'block_parameters preconditioner source_sha256',
    'target_descent': 'method policy recursive_solvers success_rule stop_rule source_sha256',
}


def candidate_record(fixture, inventory, method):
    """Resolve an explicit method against an independently audited base inventory.

    A source digest is necessary but does not excuse unknown stage wiring.
    Unresolved catalog proposals must stay Q<number>, outside this constructor.
    """
    method = copy.deepcopy(method)
    fields(method, 'isogeny endomorphism factor_base point_decomposition relation_collection '
                   'relation_linear_algebra target_descent implementation', 'method')
    require(method['isogeny'] == 'none', 'ISO1 needs a verified route adapter')
    fields(method['endomorphism'], 'order_conductor frobenius_order_conductor volcano_levels', 'endomorphism')
    # Unknown conductor is null; a value/proof is bound as part of the method.
    for key in ('order_conductor', 'frobenius_order_conductor'):
        value = method['endomorphism'][key]
        if value is not None:
            fields(value, 'value proof_sha256', key)
            natural(value['value'], key, positive=True)
            digest(value['proof_sha256'], 'conductor proof')
    require(type(method['endomorphism']['volcano_levels']) is list, 'invalid volcano levels')
    for level in method['endomorphism']['volcano_levels']:
        fields(level, 'ell level proof_sha256', 'volcano level')
        natural(level['ell'], 'ell', positive=True)
        natural(level['level'], 'volcano level')
        digest(level['proof_sha256'], 'volcano proof')
    fields(method['factor_base'], 'construction nominal_bound', 'factor-base recipe')
    require(type(method['factor_base']['construction']) is dict and method['factor_base']['construction'],
            'missing exact factor-base construction')
    def resolved(value):
        if value is None or value == '':
            return False
        if type(value) is dict:
            return all(resolved(v) for v in value.values())
        if type(value) is list:
            return all(resolved(v) for v in value)
        return True

    for stage, names in STAGE_FIELDS.items():
        fields(method[stage], names, stage)
        require(resolved(method[stage]), 'unresolved ' + stage)
        digest(method[stage]['source_sha256'], stage + ' source')
    fields(method['implementation'], 'source_manifest_sha256 components flags', 'implementation')
    digest(method['implementation']['source_manifest_sha256'], 'implementation source manifest')
    require(type(method['implementation']['flags']) is dict, 'invalid implementation flags')
    components = method['implementation']['components']
    require(type(components) is list and components, 'missing executed components')
    roles = set()
    for component in components:
        fields(component, 'role sha256', 'component')
        require(type(component['role']) is str and component['role'] and component['role'] not in roles,
                'missing or duplicate component role')
        digest(component['sha256'], 'component digest')
        roles.add(component['role'])
    source_hashes = {c['sha256'] for c in components} | {method['implementation']['source_manifest_sha256']}
    require(all(method[stage]['source_sha256'] in source_hashes for stage in STAGE_FIELDS),
            'stage source is absent from the implementation manifest/components')
    pdp = method['point_decomposition']
    natural(pdp['summands'], 'summand count', positive=True)
    require(method['relation_linear_algebra']['modulus'] == int(fixture['subgroup_order']),
            'relation LA must use the subgroup scalar field')
    codes = [pdp['solver'], method['relation_collection']['collector'],
             method['relation_linear_algebra']['solver'], method['target_descent']['method']]
    require(all(type(v) is str and re.fullmatch('[a-z][a-z0-9]*', v) for v in codes),
            'invalid compact stage code')
    method['factor_base']['inventory'] = copy.deepcopy(inventory)
    identity_only(method)
    record = dict(schema_version=1, **curve_record(fixture), **method)
    canonical(record)
    return record


def candidate_manifest(fixture, report, method):
    inventory = factor_base_inventory(report, fixture)
    record = candidate_record(fixture, inventory, method)
    c = record['curve']
    label = (f"IC1N{record['field']['n']}C{c['tag']}fb{inventory['usable_point_count']}"
             f"PDP{record['point_decomposition']['summands']}{record['point_decomposition']['solver']}"
             f"RC{record['relation_collection']['collector']}"
             f"LA{record['relation_linear_algebra']['solver']}"
             f"TD{record['target_descent']['method']}ISO0h{sha256(record)[:12]}")
    return {'candidate_id': label, 'record_sha256': sha256(record), 'record': record}


def workload_manifest(fixture, *, input_law, algorithm_seed, resource_envelope,
                      cache_policy='cold'):
    require(cache_policy in ('cold', 'warm'), 'unknown workload cache policy')
    require(type(input_law) is str and input_law, 'missing target input law')
    natural(algorithm_seed, 'algorithm seed')
    targets = fixture['targets']
    require(type(targets) is list and targets, 'empty target workload')
    require(len(fixture['target_seeds']) == len(targets), 'missing target seeds')
    require(fixture['target_scalar_constructed'] is False, 'planted target workload')
    require(type(resource_envelope) is dict and resource_envelope, 'missing resources')
    c = Curve(fixture)
    for encoded in targets:
        p = c.decode(encoded)
        require(p is not None and c.mul(p, c.r) is None, 'target outside subgroup')
    record = {'curve_id': curve_record(fixture)['curve']['curve_id'],
              'targets': [[int(x) for x in p] for p in targets], 'input_law': input_law,
              'target_seeds': fixture['target_seeds'], 'algorithm_seed': algorithm_seed,
              'cache_policy': cache_policy, 'target_count': len(targets),
              'resource_envelope': copy.deepcopy(resource_envelope)}
    h = sha256(record)
    return {'workload_id': h[:12], 'record_sha256': h, 'record': record}


def run_id(candidate_id, workload_id, number):
    require(type(candidate_id) is str and re.fullmatch(
        r'IC1N[1-9][0-9]*C[a-z][a-z0-9]*fb[1-9][0-9]*PDP[1-9][0-9]*[a-z][a-z0-9]*'
        r'RC[a-z][a-z0-9]*LA[a-z][a-z0-9]*TD[a-z][a-z0-9]*ISO[01]h[0-9a-f]{12,64}', candidate_id),
        'invalid candidate ID')
    require(type(workload_id) is str and re.fullmatch('[0-9a-f]{12,64}', workload_id),
            'invalid workload ID')
    natural(number, 'run number')
    return f'{candidate_id}W{workload_id}R{number}'


def write_immutable(path, value):
    """A reused truncated ID must match the entire record; never overwrite it."""
    path = Path(path)
    data = canonical(value) + b'\n'
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with path.open('xb') as f:
            f.write(data)
    except FileExistsError:
        require(path.read_bytes() == data, 'immutable record differs; extend a colliding digest')
