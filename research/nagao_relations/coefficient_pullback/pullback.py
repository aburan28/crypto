"""Intersect a conditioned Nagao coefficient equation with support image space."""
import time
import curves
import nagaoannihilator as scalar
import image_solver


class Search(image_solver.Search):
    def __init__(self, f, curve, d, target, deadline):
        super().__init__(f, curve, d, target, deadline)
        self.invR = f.inv(target[0]) if target[0] else None
        self.imageSquares = {u: [(value, f.sqr(value)) for _, (value, _) in basis]
                             for u, basis in self.images.items()}
        self.checkDeadline()

    def imagePreimages(self, basis, c2, c1, rhs):
        """Solve c2*v^2+c1*v=rhs for v in the supplied support image.

        Both columns and preimage updates use counted field additions. No
        uncharged bit-vector elimination is hidden behind a native XOR table.
        """
        f = self.f
        pivots = {}; kernel = []
        for value, squared in basis:
            self.checkDeadline()
            mapped = f.add(f.mul(c2, squared), f.mul(c1, value))
            preimage = value
            for k, (column, lift) in sorted(pivots.items(), reverse=True):
                if f.toCoords(mapped) >> k & 1:
                    mapped = f.add(mapped, column)
                    preimage = f.add(preimage, lift)
            if mapped:
                pivots[f.toCoords(mapped).bit_length()-1] = (mapped, preimage)
            else:
                kernel.append(preimage)
        if len(kernel) > 1:
            raise ArithmeticError('nonzero quadratic map has kernel dimension at most one')
        preimage = 0
        for k, (column, lift) in sorted(pivots.items(), reverse=True):
            if f.toCoords(rhs) >> k & 1:
                rhs = f.add(rhs, column)
                preimage = f.add(preimage, lift)
        if rhs:
            return []
        return [preimage, f.add(preimage, kernel[0])] if kernel else [preimage]

    def candidates(self):
        f = self.f; r, s = self.target
        if not r:
            # Elimination divides by b*r. Preserve the original exact chart.
            yield from super().candidates()
            return
        for hBits in range(1 << self.d):
            self.checkDeadline(); h2 = f.fromCoords(hBits)
            for b in image_solver.asRoots(f, self.curve, f.add(h2, r)):
                if not b:
                    continue
                h, _ = scalar.residualNorm(f, self.target, 0, b)
                invB = f.inv(b); invBR = f.mul(invB, self.invR)
                for z, _, _ in self.geometry:
                    self.checkDeadline()
                    if z == h2:
                        continue
                    t = f.add(r, z); u = f.add(h2, z)
                    k = scalar.polyEval(f, h, z)
                    c = f.add(f.mul(r, h2), f.mul(z, u))
                    gamma = f.mul(t, invBR)
                    delta = f.mul(f.add(f.mul(t, c), k), invBR)
                    bz = f.mul(b, z)
                    c2 = f.mul(t, f.sqr(gamma))
                    c1 = f.mul(bz, gamma)
                    rhs = f.add(k, f.add(f.mul(t, f.sqr(delta)), f.mul(bz, delta)))
                    for v in self.imagePreimages(self.imageSquares[f.toCoords(u)], c2, c1, rhs):
                        a = f.add(f.mul(gamma, v), delta)
                        yield a, b, invB, z


def cell(n, d, targetCoords, mode, budget):
    start = time.perf_counter(); deadline = start + budget
    f = scalar.CountedField(n); curve = curves.Curve(f)
    target = tuple(f.fromCoords(x) for x in targetCoords)
    solutions = set(); seen = set(); duplicates = candidates = 0
    first = None; complete = False
    try:
        f.phase = 'setup'; search = Search(f, curve, d, target, deadline)
        f.phase = 'search'
        for a, b, invB, z in search.candidates():
            candidates += 1
            if (a, b) in seen:
                duplicates += 1
                continue
            seen.add((a, b)); f.phase = 'support_extract_verify'
            xs = search.recover(a, b, invB, z); f.phase = 'search'
            if xs is not None:
                solutions.add(xs)
                if first is None:
                    first = time.perf_counter() - start
                if mode == 'first':
                    break
        else:
            complete = True
    except TimeoutError:
        pass
    elapsed = time.perf_counter() - start
    status = 'first' if mode == 'first' and solutions else ('complete' if complete else 'timeout')
    return {'n': n, 'd': d, 'target': targetCoords, 'variant': 'coefficient-pullback',
            'mode': mode, 'status': status, 'solutions': [list(x) for x in sorted(solutions)],
            'verified_unique_relations': len(solutions), 'candidate_functions': candidates,
            'duplicate_functions': duplicates, 'first_verified_seconds': first,
            'all_phase_seconds': elapsed, 'within_budget': elapsed <= budget,
            'field_api_counts': f.report(), 'full_dlp_S': None, 'rho_ratio': None}
