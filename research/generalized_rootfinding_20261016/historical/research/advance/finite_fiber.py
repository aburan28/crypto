"""Exact intermediate-fiber recovery for the binary five-source S3 chain.

This research helper is intended to accompany source-projected circuit/SAT
enumeration.  It performs no global source search and makes no speed claim.
For a nonzero binary curve constant, each forward S3 is a nonzero polynomial
of degree at most two, so five fixed sources have at most eight intermediate
chains.  All field arithmetic is performed in the supplied representation.
"""

from rootfinder.binaryfield import BinaryField


class BinaryQuadraticFiber:
    """Deterministic exact quadratic roots with a cached AS inverse map.

    Gaussian elimination proves that t -> t²+t has rank n-1 and kernel
    {0,1}.  The same factorization is reused for every quadratic in the
    supplied field, including even extension degrees.
    """

    def __init__(self, field):
        if not isinstance(field, BinaryField):
            raise TypeError("a verified BinaryField is required")
        self.field = field
        self.pivots = {}
        nullspace = []
        for bit in range(field.degree):
            basis = 1 << bit
            image, coordinates = field.square(basis) ^ basis, basis
            while image:
                pivot = image.bit_length() - 1
                previous = self.pivots.get(pivot)
                if previous is None:
                    self.pivots[pivot] = (image, coordinates)
                    break
                image ^= previous[0]
                coordinates ^= previous[1]
            if not image:
                nullspace.append(coordinates)
        if len(self.pivots) != field.degree - 1 or nullspace != [1]:
            raise ArithmeticError("Artin–Schreier rank/kernel certificate failed")

    def artin_schreier(self, rhs):
        """Return both roots of t²+t=rhs, or the empty tuple."""
        field = self.field
        residual, root = field.validate(rhs), 0
        while residual:
            pivot = residual.bit_length() - 1
            previous = self.pivots.get(pivot)
            if previous is None:
                return ()
            residual ^= previous[0]
            root ^= previous[1]
        roots = tuple(sorted((root, root ^ 1)))
        if any(field.square(value) ^ value != rhs for value in roots):
            raise ArithmeticError("Artin–Schreier solution verification failed")
        return roots

    def roots(self, a, b, c):
        """Return all distinct roots of aX²+bX+c.

        A zero polynomial is explicitly rejected instead of returning a
        truncated field.  It cannot occur in a forward S3 with curve_b != 0.
        """
        field = self.field
        a, b, c = (field.validate(value) for value in (a, b, c))
        if a == 0:
            if b == 0:
                if c == 0:
                    raise ValueError("the zero polynomial has an unrestricted fiber")
                return ()
            roots = (field.mul(c, field.inv(b)),)
        elif b == 0:
            quotient = field.mul(c, field.inv(a))
            roots = (field.pow(quotient, 1 << (field.degree - 1)),)
        else:
            rhs = field.mul(field.mul(a, c), field.inv(field.square(b)))
            scale = field.mul(b, field.inv(a))
            roots = tuple(sorted(field.mul(scale, value)
                                 for value in self.artin_schreier(rhs)))
        for root in roots:
            if field.mul(a, field.square(root)) ^ field.mul(b, root) ^ c:
                raise ArithmeticError("quadratic root verification failed")
        return roots

    def s3_roots(self, x, y, curve_b=1):
        field = self.field
        x, y, curve_b = (field.validate(value) for value in (x, y, curve_b))
        product = field.mul(x, y)
        return self.roots(field.square(x ^ y), product,
                          field.square(product) ^ curve_b)


def evaluate_s3(field, x, y, z, curve_b=1):
    x, y, z, curve_b = (field.validate(value) for value in (x, y, z, curve_b))
    product = field.mul(x, y)
    return (field.mul(field.square(x ^ y), field.square(z))
            ^ field.mul(product, z) ^ field.square(product) ^ curve_b)


def five_source_intermediate_fiber(field, sources, target, *, curve_b=1,
                                   quadratic_solver=None):
    """Recover every (u1,u2,u3) of the original four-link affine chain.

    The local fiber is complete independently of how a source tuple was
    found.  Root discovery, rational-curve lifting and group identities are
    separate contracts; no point or first-relation claim is made here.
    """
    if not isinstance(field, BinaryField):
        raise TypeError("a verified BinaryField is required")
    if len(sources) != 5:
        raise ValueError("five source coordinates are required")
    sources = tuple(field.validate(value) for value in sources)
    target, curve_b = field.validate(target), field.validate(curve_b)
    if curve_b == 0:
        raise ValueError("nonzero curve_b is required for the finite-fiber bound")
    solver = quadratic_solver or BinaryQuadraticFiber(field)
    if not isinstance(solver, BinaryQuadraticFiber) or solver.field is not field:
        raise ValueError("quadratic_solver must use this exact field object")

    x1, x2, x3, x4, x5 = sources
    layers = [()]
    width_by_layer, root_calls = [], 0
    for source_index in range(3):
        next_layer = []
        for chain in layers:
            left, right = ((x1, x2) if not chain
                           else (chain[-1], sources[source_index + 1]))
            root_calls += 1
            next_layer.extend(chain + (root,)
                              for root in solver.s3_roots(left, right, curve_b))
        layers = next_layer
        width_by_layer.append(len(layers))
    complete_chains = [chain for chain in layers
                       if evaluate_s3(field, chain[-1], x5, target, curve_b) == 0]
    for u1, u2, u3 in complete_chains:
        links = ((x1, x2, u1), (u1, x3, u2), (u2, x4, u3), (u3, x5, target))
        if any(evaluate_s3(field, *link, curve_b=curve_b) for link in links):
            raise ArithmeticError("original-chain solution verification failed")
    if len(complete_chains) > 8 or len(set(complete_chains)) != len(complete_chains):
        raise ArithmeticError("finite-fiber count/uniqueness certificate failed")
    return {"status": "complete", "sources": sources,
            "intermediate_chains": sorted(complete_chains),
            "solution_count_exact": len(complete_chains),
            "stats": {"fiber_bound": 8, "quadratic_root_calls": root_calls,
                      "width_by_layer": width_by_layer,
                      "as_rank": len(solver.pivots),
                      "as_kernel": [0, 1]}}
