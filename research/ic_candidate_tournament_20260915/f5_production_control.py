"""Independent audit of disclosed specialization and production decisive rows.

These guided paths do not estimate natural yield or reproduce search traversal.
The reference enumerates unfiltered Boolean F4 multiples, without native kernels.
"""
import itertools

from f5_boolean_control import (checked_rows, equations, evaluate, products,
                                rename_mask, rename_rows, row_space)
from identity import sha256
from oracle import Curve, require


DEPTHS = [0, 6, 12, 18, 24, 30, 35]


def substitute(rows, variable, value):
    """Square-free GF(2) specialization; cancellations use set parity."""
    bit = 1 << variable
    result = []
    for row in rows:
        terms = set()
        for mask in row:
            if value or not mask & bit:
                terms.symmetric_difference_update((mask & ~bit,))
        result.append(sorted(terms))
    return result


def decisive_reference(rows, n_vars, degree):
    raw = products(rows, n_vars, degree)
    columns = sorted({0, *(1 << v for v in range(n_vars)),
                      *(mask for row in raw for mask in row)})
    index = {mask: bit for bit, mask in enumerate(columns)}
    basis = {int(value, 16).bit_length()-1: int(value, 16)
             for value in row_space(raw, columns)}

    def belongs(row):
        if any(mask not in index for mask in row):
            return False
        value = sum(1 << index[mask] for mask in row)
        while value:
            pivot = value.bit_length()-1
            if pivot not in basis:
                return False
            value ^= basis[pivot]
        return True

    contradiction = belongs([0])
    forced = [sorted([1 << v, *([0] if value else [])])
              for v in range(n_vars) for value in (False, True)
              if belongs(sorted([1 << v, *([0] if value else [])]))]
    return raw, columns, basis, contradiction, forced, belongs


def audit(inputs, exported):
    require(inputs['summands'] == 3 and inputs['matrix_degree'] == 3
            and inputs['basis'] == [1, 2, 4, 8, 16, 32]
            and inputs['fixture']['degree'] == 17 and inputs['fixture']['curve_a'] == 1,
            'different disclosed production control')
    require(exported['schema_version'] == 1 and exported['depths'] == DEPTHS
            and len(exported['controls']) == len(inputs['controls']),
            'native production panel changed')
    curve = Curve(inputs['fixture'])
    n_vars = 35
    order = list(range(12))+list(range(18, 35))+list(range(12, 18))
    permutation = [order.index(v) for v in range(n_vars)]
    reports = []
    for expected, observed in zip(inputs['controls'], exported['controls']):
        require(expected['trial'] == observed['trial'] and observed['n_vars'] == n_vars,
                'different disclosed query/variable width')
        root = equations(curve, inputs['basis'], expected['target'][0])
        require(observed['permutation'] == permutation, 'production permutation changed')
        models = expected['models']
        require(len(models) == 6 and len(observed['paths']) == 12,
                'missing disclosed specialization paths')
        points = [curve.decode(point) for point in expected['points']]
        target = curve.decode(expected['target'])
        require(curve.mul(curve.g, expected['scalar']) == target,
                'disclosed query scalar differs from the group target')
        for model, ordering in zip(models, itertools.permutations(range(3))):
            first, second, third = [points[index] for index in ordering]
            intermediate = curve.add(first, second)
            require(intermediate is not None and curve.add(intermediate, third) == target
                    and all(0 <= point[0] < 64 for point in (first, second, third)),
                    'disclosed full-point chain failed')
            assignment = first[0] | (second[0] << 6) | (third[0] << 12) | (intermediate[0] << 18)
            require(model == dict(order=list(ordering), assignment=str(assignment),
                                  renamed_assignment=str(rename_mask(assignment, permutation))),
                    'frozen guided model differs from independent group chain')
        nodes = observed['nodes']
        used = set()
        path_report = []
        for path, (layout, model) in zip(observed['paths'],
                [(layout, model) for layout in ('original', 'interleaved') for model in models]):
            require(path['layout'] == layout and path['order'] == model['order'],
                    'production specialization paths reordered')
            assignment = int(model['assignment' if layout == 'original' else 'renamed_assignment'])
            system = root if layout == 'original' else rename_rows(root, permutation)
            require(all(evaluate(row, assignment) == 0 for row in system),
                    'disclosed model is not a Boolean root')
            require(len(path['steps']) == n_vars+1 and len(path['node_ids']) == len(DEPTHS),
                    'missing canonical specialization steps')
            for depth, recorded in enumerate(path['steps']):
                require(checked_rows(recorded, n_vars) == system,
                        'production canonical substitution differs from independent specialization')
                require(all(evaluate(row, assignment) == 0 for row in system),
                        'specialization removed a disclosed root')
                if depth in DEPTHS:
                    node_id = path['node_ids'][DEPTHS.index(depth)]
                    require(type(node_id) is int and 0 <= node_id < len(nodes),
                            'invalid production node reference')
                    require(checked_rows(nodes[node_id]['system'], n_vars) == system,
                            'production reduced a different specialized system')
                    used.add(node_id)
                if depth < n_vars:
                    system = substitute(system, depth, bool(assignment & (1 << depth)))
            require(all(not row for row in system), 'fully specialized model is nonzero')
            path_report.append(dict(layout=layout, order=model['order'],
                                    specialization_steps=n_vars, reduced_depths=DEPTHS))
        require(used == set(range(len(nodes))), 'unreferenced or missing production node')
        node_reports = []
        for node in nodes:
            rows = checked_rows(node['system'], n_vars)
            raw, columns, basis, contradiction, forced, belongs = decisive_reference(rows, n_vars, 3)
            require(not contradiction, 'disclosed model path independently contradicts itself')
            require([item['engine'] for item in node['engines']] == ['f4', 'f5'],
                    'missing/reordered production engines')
            for item in node['engines']:
                require(item['status'] == 'REDUCED', 'production matrix failed/oversize')
                reduced = checked_rows(item['rows'], n_vars)
                require(all(belongs(row) for row in reduced),
                        'production row is not an independent F4 consequence')
                require(len(reduced) == len({tuple(row) for row in reduced})
                        and sorted(reduced) == sorted(forced),
                        'production decisive set differs from independent forced-variable set')
            node_reports.append(dict(system_sha256=sha256(rows), independent_rows=len(raw),
                columns=len(columns), rank=len(basis), forced_variables=len(forced),
                independent_row_space_sha256=sha256([hex(basis[p]) for p in sorted(basis, reverse=True)])))
        reports.append(dict(trial=expected['trial'], paths=path_report, nodes=node_reports,
                            unique_systems=len(nodes), production_matrix_calls=2*len(nodes)))
    return dict(schema_version=1, status='AUDITED_DISCLOSED_PRODUCTION_PATH_CONTROL', controls=reports,
        candidate_id=None, run_id=None, measured_costs=None, online_speedup=None,
        complete_dlp=False, natural_yield_estimated=False, actual_solver_traversal_validated=False,
        promotion_eligible=False,
        scope='guided disclosed model paths; canonical substitutions and production F4/F5 decisive rows only')
