"""Reuse the quadratic dependence of coefficient circuits on an affine b space."""
import block
import run as old


def flatten(circuit):
    return [circuit.constant] + circuit.zLinear + circuit.wLinear + [x for row in circuit.mixed for x in row]


class CachedCircuit:
    def __init__(self, vector, d):
        self.constant = vector[0]
        self.zLinear = vector[1:1+d]
        self.wLinear = vector[1+d:1+2*d]
        self.mixed = [vector[1+2*d+i*d:1+2*d+(i+1)*d] for i in range(d)]


class Search(block.Search):
    def __init__(self, f, curve, v, target, deadline, batch_inverse):
        super().__init__(f, curve, v, target, deadline, True)
        self.batch_inverse = batch_inverse
        self.stats.update(circuit_evaluations=0, affine_dimension=0, cached_field_slots=0)

    def xor(self, a, b):
        return [self.f.add(x, y) for x, y in zip(a, b)]

    def prepare(self):
        f, v = self.f, self.v
        with self.phase('coefficient_table_setup'):
            ordered = []
            for h in v.values:
                self.check()
                ordered += [(b, h) for b in self.c.asRoots(f.add(h, self.r))]
            if not ordered:
                self.ordered, self.cache, self.inverses = [], {}, {}
                return
            b0 = ordered[0][0]
            matrix = block.Linear(f, [f.add(b, b0) for b, _ in ordered])
            basis = [column for _, (column, _) in matrix.ordered]
            k = len(basis)
            if len(ordered) != 2**k or len({b for b, _ in ordered}) != len(ordered):
                raise ArithmeticError('coefficient support is not the claimed affine space')
            self.stats['affine_dimension'] = k

            def evaluate(b):
                self.check()
                h = f.add(f.add(f.sqr(b), b), self.r)
                self.stats['circuit_evaluations'] += 1
                return flatten(block.Circuit(f, self.c, v, self.r, b, h, self.powers))

            constant = evaluate(b0)
            linear = [self.xor(evaluate(f.add(b0, x)), constant) for x in basis]
            mixed = []
            for i, x in enumerate(basis):
                row = []
                for j, y in enumerate(basis[:i]):
                    value = evaluate(f.add(f.add(b0, x), y))
                    row.append(self.xor(self.xor(self.xor(value, constant), linear[i]), linear[j]))
                mixed.append(row)
            vectors, bs = [constant], [b0]
            for i, x in enumerate(basis):
                self.check()
                deltas = [linear[i]]
                for derivative in mixed[i]:
                    deltas += [self.xor(value, derivative) for value in deltas]
                vectors += [self.xor(value, delta) for value, delta in zip(vectors, deltas)]
                bs += [f.add(b, x) for b in bs]
            if set(bs) != {b for b, _ in ordered}:
                raise ArithmeticError('coefficient table misses an admissible b')
            self.cache = dict(zip(bs, vectors))
            self.ordered = [(b, h) for b, h in ordered if b]
            self.stats['cached_field_slots'] = len(vectors) * len(constant)
            self.stats['circuit_evaluations_expected'] = 1 + k + k*(k-1)//2
            if self.stats['circuit_evaluations'] != self.stats['circuit_evaluations_expected']:
                raise ArithmeticError('quadratic interpolation sample count differs')
        self.inverses = None
        if self.batch_inverse:
            with self.phase('coefficient_batch_inverse'):
                self.check()
                self.inverses = dict(zip([b for b, _ in self.ordered],
                                         old.hybrid.batchInverse(f, [b for b, _ in self.ordered])))

    def solve(self, mode):
        self.prepare()
        f, v, seen = self.f, self.v, set()
        for b, h in self.ordered:
            self.check()
            with self.phase('coefficient_lookup'):
                inv_b = self.inverses[b] if self.inverses is not None else f.inv(b)
                circuit = CachedCircuit(self.cache[b], v.d)
                forbidden = {v.indices[x] for x in (0, self.r, h) if x in v.indices}
                self.stats['coefficient_values'] += 1
                self.stats['potential_branches'] += len(v.values) - len(forbidden)
            for zi, wi in self.walk(circuit, b, h, forbidden, v.d, 0, circuit.constant, circuit.wLinear):
                with self.phase('extraction_verification'):
                    result = self.recover(b, inv_b, h, zi, wi)
                if result and result[0] not in seen:
                    seen.add(result[0])
                    yield result
                    if mode == 'first':
                        return
        if self.stats['visited_admissible_leaves'] + self.stats['pruned_admissible_branches'] != self.stats['potential_branches']:
            raise ArithmeticError('cached search branch partition has a gap')
