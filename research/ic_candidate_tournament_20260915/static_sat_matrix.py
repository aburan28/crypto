"""Auditable prime-subgroup relation rows and exact dense elimination.

The tiny n17 scientific control uses a 29-column matrix. This module keeps
the real relation rows, independent group re-addition and scalar replay; it
does not infer logarithms from a known probe scalar or from a lookup table.
"""
import hashlib
import json

from generic_stages import projected_columns
from oracle import require


class RelationMatrix:
    def __init__(self, curve, base):
        self.curve = curve
        self.base = tuple(base)
        self.columns, self.projected = projected_columns(curve, self.base)
        require(len(self.columns) > 1
                and len(self.columns) == len(set(self.columns))
                and all(location is not None or curve.mul(point, curve.h) is None
                        for point, location in zip(self.base, self.projected)),
                'factor-base orbit projection is incomplete')
        self.pivots = {}
        self.rows = []
        self.keys = set()
        self.duplicate_count = 0
        self.dependent_count = 0

    def row_for(self, scalar, indices):
        c = self.curve
        require(type(scalar) is int and 0 <= scalar < c.r
                and len(indices) == 3
                and all(type(i) is int and 0 <= i < len(self.base)
                        for i in indices),
                'malformed relation scalar or factor-base indices')
        group_sum = None
        entries = [0]*len(self.columns)
        for index in indices:
            group_sum = c.add(group_sum, self.base[index])
            location = self.projected[index]
            if location is not None:
                column, coefficient = location
                entries[column] = (entries[column]+coefficient) % c.r
        require(group_sum == c.mul(c.g, scalar),
                'relation does not re-add to its public group point')
        return entries, c.h*scalar % c.r

    def push(self, scalar, indices):
        row, rhs = self.row_for(scalar, indices)
        key = scalar, tuple(sorted(indices))
        if key in self.keys:
            self.duplicate_count += 1
            return 'duplicate'
        self.keys.add(key)
        self.rows.append((row, rhs, scalar, tuple(indices)))
        reduced = row+[rhs]
        for column in range(len(self.columns)):
            coefficient = reduced[column]
            if coefficient == 0:
                continue
            if column in self.pivots:
                pivot = self.pivots[column]
                reduced = [(a-coefficient*b) % self.curve.r
                           for a,b in zip(reduced, pivot)]
            else:
                inverse = pow(coefficient, -1, self.curve.r)
                self.pivots[column] = [(value*inverse) % self.curve.r
                                       for value in reduced]
                return 'novel_rank'
        require(reduced[-1] == 0, 'verified group relations conflict over subgroup field')
        self.dependent_count += 1
        return 'dependent'

    @property
    def rank(self):
        return len(self.pivots)

    def solve(self):
        size = len(self.columns)
        require(self.rank == size, 'relation matrix is not full column rank')
        logs = [0]*size
        for column in range(size-1, -1, -1):
            pivot = self.pivots[column]
            logs[column] = (pivot[-1]
                            - sum(pivot[later]*logs[later]
                                  for later in range(column+1, size))) % self.curve.r
        require(all(sum(a*b for a,b in zip(row, logs)) % self.curve.r == rhs
                    for row,rhs,_,_ in self.rows),
                'solved logs do not satisfy every retained matrix row')
        require(all(self.curve.mul(self.curve.g, log) == point
                    for point,log in zip(self.columns, logs)),
                'solved logs fail independent scalar replay on an orbit column')
        return logs

    def descent_scalar(self, logs, a, b, indices):
        c = self.curve
        require(len(logs) == len(self.columns)
                and type(a) is int and type(b) is int
                and 0 <= a < c.r and 0 < b < c.r,
                'invalid descent coefficients or factor-base logs')
        projected_log = 0
        for index in indices:
            require(type(index) is int and 0 <= index < len(self.base),
                    'descent factor-base index outside candidate')
            location = self.projected[index]
            if location is not None:
                column, coefficient = location
                projected_log = (projected_log
                                 + coefficient*logs[column]) % c.r
        return (projected_log-c.h*a)*pow(c.h*b % c.r, -1, c.r) % c.r

    def recover_target(self, logs, target, a, b, indices):
        c = self.curve
        require(target is not None and c.mul(target, c.r) is None,
                'invalid supplied subgroup target')
        expected = c.add(c.mul(c.g, a), c.mul(target, b))
        actual = None
        for index in indices:
            require(type(index) is int and 0 <= index < len(self.base),
                    'descent factor-base index outside candidate')
            actual = c.add(actual, self.base[index])
        require(actual == expected,
                'target relation does not re-add to aG+bQ')
        recovered = self.descent_scalar(logs, a, b, indices)
        require(c.mul(c.g, recovered) == target,
                'target logarithm fails independent scalar replay')
        return recovered

    def snapshot(self):
        rows = [dict(entries=[[j, str(value)] for j,value in enumerate(row)
                              if value], rhs=str(rhs), scalar=scalar,
                     indices=list(indices))
                for row,rhs,scalar,indices in self.rows]
        data = json.dumps(rows, sort_keys=True, separators=(',', ':')).encode()
        return dict(modulus=str(self.curve.r),
                    column_points=[list(point) for point in self.columns],
                    columns=len(self.columns), rank=self.rank,
                    accepted_rows=len(self.rows),
                    duplicate_relations=self.duplicate_count,
                    dependent_relations=self.dependent_count,
                    rows_sha256=hashlib.sha256(data).hexdigest(), rows=rows)
