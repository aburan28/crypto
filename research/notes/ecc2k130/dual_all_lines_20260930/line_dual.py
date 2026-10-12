"""Build a normalized degree-263 dual from a rational twist-torsion basis.

This research interface accepts one canonical projective line ``(1, t)`` or
``(0, 1)``. Its inputs are full points on the supplied quadratic twist;
there is no implicit torsion discovery or extension-field kernel support.
"""
from __future__ import annotations

from pathlib import Path
import sys

RESEARCH = Path(__file__).resolve().parents[3]
DUAL = RESEARCH / "ecc2k130_dual_transport_20260925"
if str(DUAL) not in sys.path:
    sys.path.insert(0, str(DUAL))
from dual_transport import DualTransport

DEGREE = 263


def _canonical_point(point, curve):
    return (isinstance(point, tuple) and len(point) == 2
            and all(type(z) is int and 0 <= z < (1 << curve.F.deg)
                    for z in point)
            and curve.on_curve(point))


def build_line_dual(source, twist, basis, line, witness):
    """Return the checked dual for one of the 264 rational 263-kernel lines.

    The complement is chosen outside the kernel by projective coordinates:
    ``V`` complements ``U + [t]V`` and ``U`` complements ``V``. The underlying
    builder rejects a dependent basis, wrong order, and ambiguous sign witness.
    """
    if (not isinstance(line, tuple) or len(line) != 2
            or any(type(z) is not int for z in line)
            or not (line[0] == 1 and 0 <= line[1] < DEGREE
                    or line == (0, 1))):
        raise ValueError("line must be canonical (1,t) or (0,1) modulo 263")
    if (not isinstance(basis, tuple) or len(basis) != 2
            or not all(_canonical_point(point, twist) for point in basis)):
        raise ValueError("basis must contain two canonical finite twist points")
    u, v = basis
    if twist.mul(u, DEGREE) is not None or twist.mul(v, DEGREE) is not None:
        raise ValueError("basis points must be 263-torsion")
    if line[0] == 1:
        generator = twist.add(u, twist.mul(v, line[1]))
        complement = v
    else:
        generator, complement = v, u
    if generator is None:
        raise ValueError("line generator is infinity")
    result = DualTransport(source, twist, generator, complement, DEGREE, witness)
    result.line = line
    result.line_generator = generator
    result.line_complement = complement
    return result
