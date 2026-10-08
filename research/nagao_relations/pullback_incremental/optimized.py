"""Normalize the Nagao pullback and optionally cache elimination coordinates."""
import time

import curves
import image_solver
import nagaoannihilator as scalar
import pullback


class NormalizedSearch(pullback.Search):
    def __init__(self, f, curve, d, target, deadline):
        super().__init__(f, curve, d, target, deadline)
        r, s = target
        if r:
            # D=(r+s)/r, J=D^2+D+r. All precomputation is charged.
            self.D = f.mul(f.add(r, s), self.invR)
            self.J = f.add(f.add(f.sqr(self.D), self.D), r)
            self.normalGeometry = []
            for z, _, _ in self.geometry:
                self.checkDeadline()
                c1 = f.mul(z, self.invR)
                self.normalGeometry.append((z, c1, f.add(f.one(), c1)))

    def candidates(self):
        f = self.f
        r, _ = self.target
        if not r:
            yield from super().candidates()
            return
        for bits in range(1 << self.d):
            self.checkDeadline()
            h2 = f.fromCoords(bits)
            for b in image_solver.asRoots(f, self.curve, f.add(h2, r)):
                if not b:
                    continue
                invB = f.inv(b)
                base = f.add(r, f.mul(b, self.D))
                constant = f.mul(f.sqr(b), self.J)
                for z, c1, gamma0 in self.normalGeometry:
                    self.checkDeadline()
                    if z == h2:
                        continue
                    u = f.add(h2, z)
                    eta = f.mul(f.mul(z, u), invB)
                    gamma = f.mul(gamma0, invB)
                    delta = f.add(base, eta)
                    rhs = f.add(f.sqr(eta), constant)
                    for v in self.imagePreimages(self.imageSquares[f.toCoords(u)], f.sqr(gamma), c1, rhs):
                        yield f.add(f.mul(gamma, v), delta), b, invB, z


class CachedSearch(NormalizedSearch):
    def imagePreimages(self, basis, c2, c1, rhs):
        f = self.f
        pivots = []
        kernel = []
        for value, squared in basis:
            self.checkDeadline()
            mapped = f.add(f.mul(c2, squared), f.mul(c1, value))
            bits = f.toCoords(mapped)
            preimage = value
            for pivot, column, columnBits, lift in pivots:
                if bits >> pivot & 1:
                    # Arithmetic is still charged. XOR only updates the cached
                    # representation; it does not replace f.add accounting.
                    mapped = f.add(mapped, column)
                    bits ^= columnBits
                    preimage = f.add(preimage, lift)
            if mapped:
                pivots.append((bits.bit_length()-1, mapped, bits, preimage))
                pivots.sort(reverse=True)
            else:
                kernel.append(preimage)
        if len(kernel) > 1:
            raise ArithmeticError('nonzero quadratic map has kernel dimension at most one')
        bits = f.toCoords(rhs)
        preimage = 0
        for pivot, column, columnBits, lift in pivots:
            if bits >> pivot & 1:
                rhs = f.add(rhs, column)
                bits ^= columnBits
                preimage = f.add(preimage, lift)
        if rhs:
            return []
        return [preimage, f.add(preimage, kernel[0])] if kernel else [preimage]


def cell(n, d, targetCoords, mode, budget, variant='cached-pullback'):
    if variant not in ('normalized-pullback', 'cached-pullback'):
        raise ValueError('unknown candidate variant')
    start = time.perf_counter()
    f = scalar.CountedField(n)
    target = tuple(f.fromCoords(x) for x in targetCoords)
    searchClass = NormalizedSearch if variant == 'normalized-pullback' else CachedSearch
    solutions = set()
    seen = set()
    duplicates = candidates = 0
    first = None
    complete = False
    try:
        f.phase = 'setup'
        search = searchClass(f, curves.Curve(f), d, target, start+budget)
        f.phase = 'search'
        for a, b, invB, z in search.candidates():
            candidates += 1
            if (a, b) in seen:
                duplicates += 1
                continue
            seen.add((a, b))
            f.phase = 'support_extract_verify'
            xs = search.recover(a, b, invB, z)
            f.phase = 'search'
            if xs is not None:
                solutions.add(xs)
                if first is None:
                    first = time.perf_counter()-start
                if mode == 'first':
                    break
        else:
            complete = True
    except TimeoutError:
        pass
    elapsed = time.perf_counter()-start
    status = 'first' if mode == 'first' and solutions else ('complete' if complete else 'timeout')
    return {'n': n, 'd': d, 'target': targetCoords, 'variant': variant,
            'mode': mode, 'status': status, 'solutions': [list(x) for x in sorted(solutions)],
            'verified_unique_relations': len(solutions), 'candidate_functions': candidates,
            'duplicate_functions': duplicates, 'first_verified_seconds': first,
            'all_phase_seconds': elapsed, 'within_budget': elapsed <= budget,
            'field_api_counts': f.report(), 'full_dlp_S': None, 'rho_ratio': None}
