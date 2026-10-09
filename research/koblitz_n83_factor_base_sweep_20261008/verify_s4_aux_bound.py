#!/usr/bin/env python3
"""Audit the current S4 encoder's unavoidable x-monomial auxiliaries.

This is a source-derived lower bound, not a solver construction or timing.
It deliberately fails if a relied-on source expression changes.
"""

from hashlib import sha256
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DESCENT = ROOT / "src/cryptanalysis/binary_semaev_s4.rs"
ENCODER = ROOT / "src/cryptanalysis/semaev_sat.rs"


def require(source: str, fragments: tuple[str, ...]) -> None:
    for fragment in fragments:
        assert fragment in source, f"source changed: {fragment}"


def direct_x_monomials(l: int) -> tuple[int, int]:
    """Enumerate the distinct block products for small independent checks."""
    blocks = [range(i * l, (i + 1) * l) for i in range(3)]
    quadratics = {
        (a, b)
        for i, j in ((0, 1), (0, 2), (1, 2))
        for a in blocks[i]
        for b in blocks[j]
    }
    cubics = {
        (a, b, c)
        for a in blocks[0]
        for b in blocks[1]
        for c in blocks[2]
    }
    assert quadratics.isdisjoint(cubics)
    return len(quadratics), len(cubics)


def main() -> None:
    descent = DESCENT.read_text()
    encoder = ENCODER.read_text()
    require(
        descent,
        (
            "pub fn n_x_vars(&self) -> u32 {\n        3 * self.l",
            "pub fn n_e_vars(&self) -> u32 {\n        (1..=3).map(|i| e_len(i, self.l) as u32).sum()",
            "i * l as usize - (i - 1)",
            "AnfF2m::from_vars(i as u32 * l, l as usize)",
            "Polynomial multiplication (convolution).  No reduction.",
            "let sigma2 = x0x1.xor(&x0x2).xor(&x1x2);",
            "let sigma3 = x0x1.mul(&xs[2]);",
            "let mut row = sigma.coeffs.clone();",
        ),
    )
    require(
        encoder,
        (
            "break_symmetry: true,",
            "if m.len() >= 2 {\n                aux_of.entry(m.clone()).or_insert(0);",
            "if opts.break_symmetry {\n        next += 2 * l;",
            "for (mono, &z) in aux_of.iter() {",
            "for &v in mono {\n            solver.add_clause(vec![-zi, v as Lit]);",
            "solver.add_clause(big);",
        ),
    )
    for l in range(1, 6):
        assert direct_x_monomials(l) == (3 * l * l, l**3)

    l = 83
    quadratics, cubics = 3 * l * l, l**3
    core = 3 * l + (6 * l - 3) + 2 * l
    variables = core + quadratics + cubics
    and_clauses = 3 * quadratics + 4 * cubics
    assert (core, quadratics, cubics, variables, and_clauses) == (
        910,
        20667,
        571787,
        593364,
        2349149,
    )
    for path in (DESCENT, ENCODER):
        print(f"source {path.relative_to(ROOT)} sha256={sha256(path.read_bytes()).hexdigest()}")
    print("small_block_enumeration l=1..5 PASS")
    print(f"n83 core={core} x_quadratic_aux={quadratics} x_cubic_aux={cubics}")
    print(f"n83 sat_variables_lower_bound={variables}")
    print(f"n83 and_definition_clauses_lower_bound={and_clauses}")
    print("status PASS source-derived bound only; no solver memory or runtime claim")


if __name__ == "__main__":
    main()
