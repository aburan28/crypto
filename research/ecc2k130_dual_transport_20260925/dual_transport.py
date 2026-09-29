"""Research-only normalized dual of a checked odd-degree binary Vélu map.

The caller supplies a rational complementary torsion point on the quadratic
twist. The resulting raw reverse quotient can differ from the dual by
negation, which is fixed using a known non-2-torsion source witness.
"""
from __future__ import annotations

from pathlib import Path
import sys

DEPENDENCIES = Path(__file__).resolve().parents[1] / "ecc2k130_oriented_transport_20260924"
if str(DEPENDENCIES) not in sys.path:
    sys.path.insert(0, str(DEPENDENCIES))
from oriented_velu import BinaryVeluMap


class DualTransport:
    """Forward map, raw reverse quotient and sign-corrected dual.

    This is variable-time research arithmetic, not a general dual builder:
    the kernel and complement must have rational abscissae on a supplied
    twist with the same b as the source. The orientation witness is required
    because the normalized raw reverse map may compose to ``[-degree]``.
    """

    def __init__(self, source, twist, generator, complement, degree, witness):
        if type(degree) is not int or degree < 3 or degree % 2 == 0:
            raise ValueError("degree must be an odd integer at least 3")
        self.forward = BinaryVeluMap.from_generator(source, twist, generator, degree)
        self.forward_twist = BinaryVeluMap.from_generator(twist, twist, generator, degree)
        if (self.forward.codomain.b != self.forward_twist.codomain.b
                or self.forward.codomain.a != source.a
                or self.forward_twist.codomain.a != twist.a):
            raise ArithmeticError("forward curve and twist quotient disagree")
        if (not isinstance(complement, tuple) or len(complement) != 2
                or any(type(z) is not int or z < 0 or z >= (1 << twist.F.deg)
                       for z in complement)
                or not twist.on_curve(complement)):
            raise ValueError("complement must be a canonical finite twist point")
        if twist.mul(complement, degree) is not None:
            raise ValueError("complement must have order dividing degree")
        image_complement = self.forward_twist(complement)
        if image_complement is None:
            raise ValueError("complement lies in the forward kernel")
        if self.forward_twist.codomain.mul(image_complement, degree) is not None:
            raise ArithmeticError("reverse kernel image has wrong order")
        self.complement = complement
        self.reverse_kernel_generator = image_complement
        self.raw_reverse = BinaryVeluMap.from_generator(
            self.forward.codomain, self.forward_twist.codomain, image_complement, degree)
        self.raw_reverse_twist = BinaryVeluMap.from_generator(
            self.forward_twist.codomain, self.forward_twist.codomain, image_complement, degree)
        if ((self.raw_reverse.codomain.a, self.raw_reverse.codomain.b)
                != (self.forward.source.a, self.forward.source.b)):
            raise ArithmeticError("normalized reverse codomain is not the source model")
        if ((self.raw_reverse_twist.codomain.a, self.raw_reverse_twist.codomain.b)
                != (self.forward_twist.source.a, self.forward_twist.source.b)):
            raise ArithmeticError("normalized reverse twist codomain is not the twist model")
        self.source = self.forward.source
        self.codomain = self.forward.codomain
        self.degree = degree
        witness_image = self.forward(witness)  # validates source membership
        expected = self.source.mul(witness, degree)
        if (witness_image is None or expected is None
                or expected == self.source.neg(expected)):
            raise ValueError("orientation witness must separate scalar signs")
        observed = self.raw_reverse(witness_image)
        if observed == expected:
            self.raw_reverse_scalar_sign = 1
        elif observed == self.source.neg(expected):
            self.raw_reverse_scalar_sign = -1
        else:
            raise ArithmeticError("reverse composition is neither scalar sign")

    def dual(self, point):
        """Evaluate the sign-corrected dual from codomain to source."""
        result = self.raw_reverse(point)
        return result if self.raw_reverse_scalar_sign == 1 else self.source.neg(result)

    def compose(self, point):
        """Evaluate dual(forward(point)); expected to equal [degree]point."""
        return self.dual(self.forward(point))
