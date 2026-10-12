from sage.all import GF, PolynomialRing, EllipticCurve, ZZ

CASES = (
    ("Oakley Group 3", 155, 62, 0x7338f, 123, 456,
     3805993847215893016155463826195386266397436443, 12, 4),
    ("Oakley Group 4", 185, 69, 0x1ee9, 24, 13,
     12259964326927110866866776214413170562013096250261263279, 4, 2),
)

for name, degree, middle, b, x, y, r, full_cofactor, generator_cofactor in CASES:
    ring = PolynomialRing(GF(2), "u")
    u = ring.gen()
    modulus = u**degree + u**middle + 1
    assert modulus.is_irreducible()
    field = GF(2**degree, name="u", modulus=modulus)
    def element(bits):
        return sum((field.gen()**i for i in range(bits.bit_length()) if (bits >> i) & 1), field(0))
    curve = EllipticCurve(field, [1, 0, 0, 0, element(b)])
    point = curve(element(x), element(y))
    assert ZZ(r).is_prime(proof=True)
    subgroup_point = generator_cofactor * point
    assert subgroup_point != curve(0)
    assert r * subgroup_point == curve(0)
    assert (generator_cofactor // 2) * r * point != curve(0)
    assert full_cofactor % generator_cofactor == 0
    print(f"{name}: polynomial irreducible; RFC point on curve; "
          f"published point order = {generator_cofactor}*r; "
          f"selected point [{generator_cofactor}]P has prime order r")
