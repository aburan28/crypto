"""Oriented odd-degree Vélu maps with rational kernel abscissae.

The supported source model is y²+xy=x³+a*x²+b in characteristic two.
The kernel generator can live on the source or a quadratic twist having
the same b.  Only the generator's abscissae enter the map, so extension
arithmetic and arbitrary choices of the image's y sign are unnecessary.

This is a research arithmetic interface. It is not constant time and is
not a point-decomposition solver. See PROTOCOL.md and README.md for scope.
"""
from __future__ import annotations

from pathlib import Path
import sys

DEPENDENCIES = Path(__file__).resolve().parents[1] / "ecc2k130_relations"
if str(DEPENDENCIES) not in sys.path:
    sys.path.insert(0, str(DEPENDENCIES))
from relations import Koblitz


def _element(value, field):
    return type(value) is int and 0 <= value < (1 << field.deg)


def _curve(curve):
    """Copy coefficients so later caller edits do not change the map."""
    if not _element(curve.a, curve.F) or not _element(curve.b, curve.F) or not curve.b:
        raise ValueError("expected nonsingular binary curve with canonical coefficients")
    return Koblitz(curve.F, a=curve.a, b=curve.b)


def _same_field(left, right):
    return (left.deg, left.irr) == (right.deg, right.irr)


def _valid_point(curve, point):
    if point is None:
        return True
    return (isinstance(point, tuple) and len(point) == 2
            and all(_element(coordinate, curve.F) for coordinate in point)
            and curve.on_curve(point))


class BinaryVeluMap:
    """Use ``from_generator`` to construct a checked cyclic quotient.

    ``None`` denotes infinity. Calling the map validates source membership;
    infinity and all rational kernel points return ``None``. Invalid inputs
    raise ValueError before the kernel shortcut. The rational-abscissa
    restriction excludes general kernel polynomials with nonrational roots.
    """

    def __init__(self, *args, **kwargs):
        raise TypeError("use BinaryVeluMap.from_generator")

    @classmethod
    def from_generator(cls, source, kernel_curve, generator, degree):
        """Build from a generator of exact odd order ``degree >= 3``.

        ``kernel_curve`` has the same field and b as ``source``; its a may
        differ. They become isomorphic by y -> y+s*x with
        s²+s=source.a+kernel_curve.a. This transports the cyclic kernel
        without changing its rational x coordinates. Complete enumeration
        checks the generator's exact order, including composite odd degrees.
        """
        source, kernel_curve = _curve(source), _curve(kernel_curve)
        if not _same_field(source.F, kernel_curve.F) or source.b != kernel_curve.b:
            raise ValueError("kernel curve must have the source field and b")
        if type(degree) is not int or degree < 3 or degree % 2 == 0:
            raise ValueError("degree must be an odd integer at least 3")
        if generator is None or not _valid_point(kernel_curve, generator):
            raise ValueError("kernel generator must be a finite point on kernel_curve")
        if kernel_curve.mul(generator, degree) is not None:
            raise ValueError("degree does not annihilate the kernel generator")
        roots, point, total = set(), generator, 0
        for _ in range((degree-1)//2):
            if point is None or point[0] in roots:
                raise ValueError("kernel generator has smaller order than degree")
            roots.add(point[0])
            total ^= point[0]
            previous = point
            point = kernel_curve.add(point, generator)
        if point != kernel_curve.neg(previous):
            raise ValueError("kernel enumeration did not close at the stated odd order")
        b_prime = source.b ^ total ^ source.F.sqr(total)
        if not b_prime:
            raise ValueError("kernel quotient unexpectedly has singular codomain")
        result = object.__new__(cls)
        result.source = source
        result.codomain = Koblitz(source.F, a=source.a, b=b_prime)
        result.degree = degree
        result.kernel_abscissae = frozenset(roots)
        result._roots = tuple(sorted(roots))
        result.half_kernel_sum = total
        return result

    def __call__(self, point):
        if not _valid_point(self.source, point):
            raise ValueError("input must be infinity or a canonical source-curve point")
        if point is None:
            return None
        x, y = point
        if x in self.kernel_abscissae:
            return None
        field = self.source.F
        inverses = field.batch_inv([x ^ u for u in self._roots])
        x_image, y_image = x, y
        x_squared_plus_y = field.sqr(x) ^ y
        for u, inverse in zip(self._roots, inverses):
            c = field.mul(u, inverse)
            c_squared = field.sqr(c)
            x_image ^= c ^ c_squared
            # Paired full Vélu sums, after y -> y+sum(u) normalization:
            # u*(x²+y)/(x+u)² + c³ + c. This is independent of a twist lift.
            y_image ^= (field.mul(field.mul(u, field.sqr(inverse)), x_squared_plus_y)
                        ^ field.mul(c, c_squared) ^ c)
        image = x_image, y_image
        if not self.codomain.on_curve(image):
            raise ArithmeticError("Vélu image failed the codomain equation")
        return image

    def then(self, following):
        """Compose two checked maps when their normalized models match."""
        return ComposedVeluMap(self, following)


class ComposedVeluMap:
    def __init__(self, first, following):
        left, right = first.codomain, following.source
        if not _same_field(left.F, right.F) or (left.a, left.b) != (right.a, right.b):
            raise ValueError("composition requires identical intermediate curve models")
        self.first, self.following = first, following
        self.source, self.codomain = first.source, following.codomain
        self.degree = first.degree * following.degree

    def __call__(self, point):
        return self.following(self.first(point))

    def then(self, following):
        return ComposedVeluMap(self, following)
