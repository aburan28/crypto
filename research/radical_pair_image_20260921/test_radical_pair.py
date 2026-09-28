#!/usr/bin/env sage -python
"""Small deterministic correctness checks for radical_pair.py."""

from radical_pair import BinaryCell, PrimeCell, radical_pair_basis
from sage.all import GF


def test_pair_basis():
    field = GF(17)
    roots = [field(0), field(3), field(5), field(6), field(7)]
    ring, basis, metadata = radical_pair_basis(field, roots)
    assert len(basis) == 6
    assert {int(g.total_degree()) for g in basis} == {5}
    assert metadata["image_size"] == 15
    leading = ring.ideal([g.lm() for g in basis]).hilbert_series()
    t = leading.parent().gen()
    assert leading == sum((degree + 1) * t**degree for degree in range(5))


def test_prime_cell():
    cell = PrimeCell(17, 2, 2, 3)
    target = list(cell.field)[0]
    certificate = cell.certificate(target)
    assert certificate["presentation_equivalence"]
    assert certificate["root_recovery"]
    assert certificate["nondegenerate_group_equivalence"]
    systems = cell.systems(target)
    assert set(systems) == {"direct", "norm_symmetric", "radical_symmetric"}


def test_binary_cell():
    cell = BinaryCell(2)
    target = cell.roots[0]
    certificate = cell.certificate(target)
    assert certificate["presentation_equivalence"]
    assert certificate["unordered_pair_enumeration"]
    assert certificate["quadratic_root_recovery"]
    assert certificate["polynomial_evaluation"]
    assert certificate["root_recovery"]
    assert certificate["nondegenerate_group_equivalence"]
    systems = cell.systems(target)
    assert set(systems) == {"direct", "radical_symmetric"}


if __name__ == "__main__":
    test_pair_basis()
    test_prime_cell()
    test_binary_cell()
    print("radical pair image tests: pass")
