"""Bilinear support circuit, exact leaves, and certified block rejection."""
from contextlib import contextmanager
import time
import run as old


class Linear:
    def __init__(self, f, columns, lifts=False, stopFull=False):
        self.f, self.rows, self.kernel = f, {}, []
        self.lifts = lifts
        for i, value in enumerate(columns):
            pre = 1 << i if lifts else 0
            for k, (column, lift) in sorted(self.rows.items(), reverse=True):
                if value >> k & 1:
                    value = f.add(value, column)
                    if lifts:
                        pre ^= lift
                        f.word('preimage_xor')
            if value:
                self.rows[value.bit_length() - 1] = value, pre
            elif lifts:
                self.kernel.append(pre)
            if stopFull and len(self.rows) == f.m:
                break
        self.ordered = sorted(self.rows.items(), reverse=True)

    def reduce(self, value):
        pre = 0
        for k, (column, lift) in self.ordered:
            if value >> k & 1:
                value = self.f.add(value, column)
                if self.lifts:
                    pre ^= lift
                    self.f.word('preimage_xor')
        return value, pre

    def separator(self, value):
        residual, _ = self.reduce(value)
        if not residual:
            return None
        mask = 1 << (residual.bit_length() - 1)
        for k, (column, _) in reversed(self.ordered):
            if self.f.parity(column, mask):
                mask ^= 1 << k
                self.f.word('separator_mask_xor')
        if not self.f.parity(value, mask):
            raise ArithmeticError('invalid rejection separator')
        return mask


def direct(f, curve, r, b, z, w):
    """Unexpanded denominator-cleared identity, used by the checker only."""
    b2 = f.sqr(b)
    h2 = f.add(f.add(b2, b), r)
    u, t = f.add(h2, z), f.add(r, z)
    value = f.add(f.sqr(w), f.mul(u, w))
    return f.add(f.add(f.mul(f.sqr(t), f.sqr(value)),
                       f.mul(f.mul(f.mul(b2, r), z), value)),
                 f.add(f.mul(f.mul(f.sqr(r), f.sqr(z)), f.sqr(u)),
                       f.mul(f.sqr(b2), curve.a6)))


class Circuit:
    def __init__(self, f, curve, v, r, b, h2, powers):
        self.f, self.v = f, v
        r2, h22, b2 = f.sqr(r), f.sqr(h2), f.sqr(b)
        rh = f.mul(r2, h22)
        br = f.mul(b2, r)
        self.constant = f.mul(f.sqr(b2), curve.a6)
        self.zLinear, self.wLinear, self.mixed = [], [], []
        for z, z2, z4 in powers:
            self.zLinear.append(f.add(f.mul(rh, z2), f.mul(r2, z4)))
            self.wLinear.append(f.add(f.mul(r2, z4), f.mul(rh, z2)))
            aa = z2
            bb = f.add(f.add(f.mul(f.add(r2, h22), z2), z4), f.mul(br, z))
            cc = f.mul(br, f.add(f.mul(h2, z), z2))
            self.mixed.append([f.add(f.add(f.mul(aa, w4), f.mul(bb, w2)), f.mul(cc, w))
                               for w, w2, w4 in powers])

    def evaluate(self, zIndex, wIndex):
        f = self.f
        out = self.constant
        for i in range(self.v.d):
            if zIndex >> i & 1:
                out = f.add(out, self.zLinear[i])
            if wIndex >> i & 1:
                out = f.add(out, self.wLinear[i])
            if zIndex >> i & 1:
                for j in range(self.v.d):
                    if wIndex >> j & 1:
                        out = f.add(out, self.mixed[i][j])
        return out


class Search:
    def __init__(self, f, curve, v, target, deadline, prune):
        self.f, self.c, self.v, self.target = f, curve, v, target
        self.deadline, self.prune = deadline, prune
        self.times = {}
        self.stats = {'coefficient_values': 0, 'potential_branches': 0,
                      'visited_admissible_leaves': 0, 'pruned_admissible_branches': 0,
                      'block_tests': 0, 'full_rank_blocks': 0, 'rejected_blocks': 0,
                      'rejections_by_free_bits': {}, 'leaf_solutions': 0}
        self.certificates = []
        self.r, s = target
        if not self.r:
            raise ValueError('zero-abscissa target uses the predecessor chart')
        with self.phase('setup'):
            self.powers = [(w, f.sqr(w), f.frob(w, 2)) for w in v.basis]
            self.invR = f.inv(self.r)
            self.dd = f.mul(f.add(self.r, s), self.invR)

    def check(self):
        if time.perf_counter() >= self.deadline:
            raise TimeoutError('block search deadline')

    @contextmanager
    def phase(self, name):
        self.f.phase = name
        start = time.perf_counter()
        try:
            yield
        finally:
            self.times[name] = self.times.get(name, 0.) + time.perf_counter() - start

    def walk(self, circuit, b, h2, forbidden, free, prefix, rhs, columns):
        self.check()
        f, v = self.f, self.v
        if not free:
            if prefix in forbidden:
                return
            self.stats['visited_admissible_leaves'] += 1
            with self.phase('leaf_solve'):
                matrix = Linear(f, columns, True)
                residual, solution = matrix.reduce(rhs)
                if residual:
                    return
                answers = [solution]
                for kernel in matrix.kernel:
                    answers += [x ^ kernel for x in answers]
                    f.word('kernel_enumeration_xor', len(answers) // 2)
            for answer in answers:
                self.check()
                self.stats['leaf_solutions'] += 1
                yield prefix, answer
            return
        if self.prune:
            with self.phase('block_rejection'):
                self.stats['block_tests'] += 1
                relaxed = columns + circuit.zLinear[:free]
                for row in circuit.mixed[:free]:
                    relaxed += row
                matrix = Linear(f, relaxed, False, True)
                if len(matrix.rows) == f.m:
                    self.stats['full_rank_blocks'] += 1
                    mask = None
                else:
                    mask = matrix.separator(rhs)
                if mask is not None:
                    removed = (1 << free) - sum(x >> free == prefix >> free for x in forbidden)
                    self.stats['pruned_admissible_branches'] += removed
                    self.stats['rejected_blocks'] += 1
                    key = str(free)
                    bucket = self.stats['rejections_by_free_bits'].setdefault(key, {'blocks': 0, 'admissible_branches': 0})
                    bucket['blocks'] += 1
                    bucket['admissible_branches'] += removed
                    self.certificates.append({'b': f.toCoords(b), 'h2': f.toCoords(h2),
                                              'prefix': prefix, 'free_bits': free, 'separator': mask})
                    return
        bit = free - 1
        yield from self.walk(circuit, b, h2, forbidden, bit, prefix, rhs, columns)
        with self.phase('partition_updates'):
            nextColumns = [f.add(a, delta) for a, delta in zip(columns, circuit.mixed[bit])]
            nextRhs = f.add(rhs, circuit.zLinear[bit])
        yield from self.walk(circuit, b, h2, forbidden, bit, prefix | (1 << bit), nextRhs, nextColumns)

    def recover(self, b, invB, h2, zIndex, wIndex):
        f, v, c = self.f, self.v, self.c
        z, w = v.values[zIndex], v.values[wIndex]
        u = f.add(h2, z)
        xs = [z, w, f.add(w, u)]
        if len(set(xs)) != 3 or any(x in (0, self.r) for x in xs):
            return None
        value = f.add(f.sqr(w), f.mul(u, w))
        gamma = f.mul(f.mul(f.add(self.r, z), self.invR), invB)
        eta = f.mul(f.mul(z, u), invB)
        delta = f.add(f.add(self.r, f.mul(b, self.dd)), eta)
        a = f.add(f.mul(gamma, value), delta)
        h, constant = old.algebra.residualNorm(f, c, self.target, a, b)
        if h[2] != h2 or old.solvers.scalar.polyEval(f, h, z):
            raise ArithmeticError('bilinear leaf failed norm identity')
        points = [(x, f.mul(f.add(f.add(f.sqr(x), f.mul(a, x)), constant), invB)) for x in xs]
        return old.solvers.checked(f, c, xs, points, self.target)

    def solve(self, mode):
        f, v = self.f, self.v
        seen = set()
        for h2 in v.values:
            self.check()
            with self.phase('coefficient_setup'):
                bs = self.c.asRoots(f.add(h2, self.r))
            for b in bs:
                if not b:
                    continue
                self.check()
                with self.phase('coefficient_setup'):
                    invB = f.inv(b)
                    circuit = Circuit(f, self.c, v, self.r, b, h2, self.powers)
                    forbidden = {v.indices[x] for x in (0, self.r, h2) if x in v.indices}
                    self.stats['coefficient_values'] += 1
                    self.stats['potential_branches'] += len(v.values) - len(forbidden)
                for zIndex, wIndex in self.walk(circuit, b, h2, forbidden, v.d, 0, circuit.constant, circuit.wLinear):
                    with self.phase('extraction_verification'):
                        result = self.recover(b, invB, h2, zIndex, wIndex)
                    if result and result[0] not in seen:
                        seen.add(result[0])
                        yield result
                        if mode == 'first':
                            return
        if self.stats['visited_admissible_leaves'] + self.stats['pruned_admissible_branches'] != self.stats['potential_branches']:
            raise ArithmeticError('complete branch partition has a gap')
