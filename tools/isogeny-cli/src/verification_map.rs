//! Exact rational-map substitution using the replay implementation's polynomial arithmetic.
use super::{
    curve::Model,
    field::Field,
    poly::{self, Poly},
};
pub fn check(f: &Field, source: &Model, target: &Model, n: &Poly, d: &Poly) -> bool {
    let cubic = vec![source.b, source.a, f.zero(), f.one()];
    let derivative = poly::sub(
        f,
        &poly::mul(f, &poly::derivative(f, n), d),
        &poly::mul(f, n, &poly::derivative(f, d)),
    );
    let left = poly::mul(f, &cubic, &poly::mul(f, &derivative, &derivative));
    let d2 = poly::mul(f, d, d);
    let d3 = poly::mul(f, &d2, d);
    let d4 = poly::mul(f, &d2, &d2);
    let n3 = poly::mul(f, &poly::mul(f, n, n), n);
    let right = poly::add(
        f,
        &poly::add(
            f,
            &poly::mul(f, &n3, d),
            &poly::scale(f, &poly::mul(f, n, &d3), &target.a),
        ),
        &poly::scale(f, &d4, &target.b),
    );
    left == right
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;
    #[test]
    fn known_degree_three_map_passes_and_coefficient_mutations_fail() {
        let f = Field::new(&BigUint::from(101u32)).unwrap();
        let source = Model {
            a: f.zero(),
            b: f.from_u64(2),
        };
        let target = Model {
            a: f.zero(),
            b: f.from_u64(47),
        };
        let n = vec![f.from_u64(8), f.zero(), f.zero(), f.one()];
        let d = vec![f.zero(), f.zero(), f.one()];
        assert!(check(&f, &source, &target, &n, &d));
        let mut changed = n.clone();
        changed[0] = f.from_u64(9);
        assert!(!check(&f, &source, &target, &changed, &d));
        let mut changed_target = target;
        changed_target.b = f.from_u64(48);
        assert!(!check(&f, &source, &changed_target, &n, &d));
    }
}
