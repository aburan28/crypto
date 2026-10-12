"""Independently reconstruct disclosed S3 ANFs and bounded Macaulay row spaces.

No solver search or target recovery. Field-valued Boolean polynomials use
square-free masks; multiplication cancels equal masks over characteristic two.
"""
import itertools

from identity import sha256
from oracle import Curve, require


def add(*polys):
    out = {}
    for poly in polys:
        for mask, coefficient in poly.items():
            out[mask] = out.get(mask, 0) ^ coefficient
    return {mask: coefficient for mask, coefficient in out.items() if coefficient}


def multiply(curve, left, right):
    out = {}
    for a, x in left.items():
        for b, y in right.items():
            mask = a | b
            out[mask] = out.get(mask, 0) ^ curve.fm(x, y)
    return {mask: coefficient for mask, coefficient in out.items() if coefficient}


def square(curve, poly):
    # Cross terms cancel in characteristic two; x_i^2=x_i on the cube.
    return {mask: curve.fm(coefficient, coefficient) for mask, coefficient in poly.items()}


def s3_polynomial(curve, x, y, z):
    product = multiply(curve, x, y)
    return add(multiply(curve, square(curve, add(x, y)), square(curve, z)),
               multiply(curve, product, z), square(curve, product), {0: 1})


def equations(curve, basis, target_x):
    ell = len(basis)
    summands = [{1 << (i*ell+j): value for j, value in enumerate(basis)} for i in range(3)]
    intermediate = {1 << (3*ell+j): 1 << j for j in range(curve.n)}
    links = [s3_polynomial(curve, summands[0], summands[1], intermediate),
             s3_polynomial(curve, intermediate, summands[2], {0: target_x})]
    return [sorted(mask for mask, value in link.items() if value & (1 << bit))
            for link in links for bit in range(curve.n)]


def checked_rows(rows, n_vars):
    require(type(rows) is list and type(n_vars) is int and 0 < n_vars <= 64,
            'malformed Boolean matrix or width')
    require(all(type(row) is list and len(row) == len(set(row))
                and all(type(mask) is int and 0 <= mask < 1 << n_vars for mask in row)
                for row in rows), 'malformed square-free Boolean row')
    return [sorted(row) for row in rows]


def evaluate(row, assignment):
    return sum(mask & assignment == mask for mask in row) % 2


def rename_mask(mask, permutation):
    return sum(1 << new for old, new in enumerate(permutation) if mask & (1 << old))


def rename_rows(rows, permutation):
    return [sorted(rename_mask(mask, permutation) for mask in row) for row in rows]


def products(rows, n_vars, degree):
    result = []
    for row in rows:
        if not row:
            continue
        gap = degree-max(mask.bit_count() for mask in row)
        for size in range(gap+1):
            for selected in itertools.combinations(range(n_vars), size):
                multiplier = sum(1 << variable for variable in selected)
                terms = set()
                for mask in row:
                    value = mask | multiplier
                    terms.symmetric_difference_update((value,))
                if terms:
                    result.append(sorted(terms))
    return result


def row_space(rows, columns):
    index = {mask: bit for bit, mask in enumerate(columns)}
    require(all(mask in index for row in rows for mask in row),
            'matrix consequence uses a monomial outside the unfiltered F4 support')
    pivots = {}
    for row in rows:
        value = sum(1 << index[mask] for mask in row)
        while value:
            pivot = value.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    # Canonical reduced row space, independent of row order/native kernels.
    for pivot in sorted(pivots):
        for later in sorted(pivots):
            if later > pivot and pivots[later] & (1 << pivot):
                pivots[later] ^= pivots[pivot]
    return [hex(pivots[pivot]) for pivot in sorted(pivots, reverse=True)]


def audit(inputs, exported):
    require(inputs['summands'] == 3 and inputs['matrix_degree'] == 3
            and inputs['fixture']['degree'] == 17 and inputs['fixture']['curve_a'] == 1
            and inputs['basis'] == [1, 2, 4, 8, 16, 32], 'different disclosed encoding control')
    require(exported['schema_version'] == 1 and exported['basis'] == inputs['basis']
            and len(exported['controls']) == len(inputs['controls']), 'native control set changed')
    curve = Curve(inputs['fixture'])
    basis = inputs['basis']
    ell, n_vars = len(basis), 3*len(basis)+curve.n
    order = list(range(2*ell))+list(range(3*ell, n_vars))+list(range(2*ell, 3*ell))
    permutation = [order.index(variable) for variable in range(n_vars)]
    reports = []
    for expected, observed in zip(inputs['controls'], exported['controls']):
        require(observed['trial'] == expected['trial'] and observed['target'] == expected['target']
                and observed['n_vars'] == n_vars and observed['ell'] == ell and observed['m'] == 3,
                'native Boolean control/layout changed')
        rows = equations(curve, basis, expected['target'][0])
        for name in ('direct', 'template', 'reused', 'repeated'):
            require(checked_rows(observed[name], n_vars) == rows,
                    'independent coefficient identity differs from native '+name)
        require(observed['permutation'] == permutation, 'native variable permutation differs')
        interleaved = rename_rows(rows, permutation)
        require(checked_rows(observed['interleaved'], n_vars) == interleaved,
                'renamed Boolean equations differ')
        points = [curve.decode(point) for point in expected['points']]
        target = curve.decode(expected['target'])
        require(curve.mul(curve.g, expected['scalar']) == target,
                'disclosed points/query scalar fail independent group replay')
        assignments = []
        for ordering in itertools.permutations(range(3)):
            first, second, third = [points[index] for index in ordering]
            intermediate = curve.add(first, second)
            require(intermediate is not None and curve.add(intermediate, third) == target,
                    'disclosed full-point chain failed')
            require(all(0 <= point[0] < 1 << ell for point in (first, second, third)),
                    'summand is outside the disclosed standard subspace')
            assignment = first[0] | (second[0] << ell) | (third[0] << (2*ell)) | (intermediate[0] << (3*ell))
            renamed = rename_mask(assignment, permutation)
            require(all(evaluate(row, assignment) == 0 for row in rows)
                    and all(evaluate(row, renamed) == 0 for row in interleaved),
                    'implemented ANF excludes the independently derived chain assignment')
            assignments.append(dict(order=list(ordering), assignment=str(assignment),
                                    renamed_assignment=str(renamed), intermediate=list(intermediate)))
        require([(item['layout'], item['engine']) for item in observed['matrices']]
                == [(layout, engine) for layout in ('original', 'interleaved') for engine in ('f4', 'f5')],
                'missing/reordered native matrix controls')
        matrices = []
        for layout, system in (('original', rows), ('interleaved', interleaved)):
            raw = products(system, n_vars, inputs['matrix_degree'])
            columns = sorted({mask for row in raw for mask in row})
            expected_space = row_space(raw, columns)
            for item in observed['matrices']:
                if item['layout'] != layout:
                    continue
                require(item['status'] == 'REDUCED', 'root matrix control was not reduced')
                reduced = checked_rows(item['rows'], n_vars)
                require(row_space(reduced, columns) == expected_space,
                        'native '+item['engine']+' row space differs from independently constructed unfiltered F4')
                require(type(item['elimination_and_criterion_word_xors']) is int
                        and item['elimination_and_criterion_word_xors'] >= 0, 'missing bounded native work diagnostic')
                for assignment in assignments:
                    mask = int(assignment['assignment' if layout == 'original' else 'renamed_assignment'])
                    require(all(evaluate(row, mask) == 0 for row in reduced),
                            'native matrix consequence excludes the disclosed witness')
                matrices.append(dict(layout=layout, engine=item['engine'], rows=len(reduced),
                    independent_unfiltered_rows=len(raw), columns=len(columns), rank=len(expected_space),
                    row_space_sha256=sha256(expected_space),
                    elimination_and_criterion_word_xors=item['elimination_and_criterion_word_xors']))
        reports.append(dict(trial=expected['trial'], variables=n_vars, equations=len(rows),
            coefficient_identity_verified=True, witness_assignments=assignments, matrices=matrices))
    return dict(schema_version=1, status='AUDITED_DISCLOSED_BOOLEAN_SYSTEM_CONTROL', controls=reports,
        candidate_id=None, run_id=None, measured_costs=None, online_speedup=None,
        complete_dlp=False, natural_yield_estimated=False, promotion_eligible=False,
        production_linear_tail_and_branch_traversal_validated=False,
        scope='two disclosed n17 inputs; exact ANF identity and full-readback root row spaces only')
