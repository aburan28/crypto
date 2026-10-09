#[path = "../../src/cryptanalysis/isogeny_walk/curve.rs"]
mod curve;
#[path = "../../src/cryptanalysis/isogeny_walk/field.rs"]
mod field;
#[path = "verification_kernel.rs"]
mod kernel;
#[path = "../../src/cryptanalysis/isogeny_walk/modpoly.rs"]
mod modpoly;
#[path = "../../src/cryptanalysis/isogeny_walk/kernel.rs"]
mod original_kernel;
#[path = "../../src/cryptanalysis/isogeny_walk/poly.rs"]
mod original_poly;
#[path = "verification_poly.rs"]
mod poly;
fn main() {}
#[path = "../../src/cryptanalysis/modular_polynomial.rs"]
pub mod modular_polynomial;
pub mod cryptanalysis {
    pub use crate::modular_polynomial;
}
#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;
    #[test]
    fn selected_recurrence_matches_all_original_indices() {
        let f = field::Field::new(&BigUint::from(1009u32)).unwrap();
        let e = curve::Model {
            a: f.from_u64(7),
            b: f.from_u64(11),
        };
        let mut m: Vec<_> = (0..17).map(|i| f.from_u64((i * i + 3) as u64)).collect();
        m.push(f.one());
        let indices: Vec<_> = (0..=160).collect();
        let full = original_kernel::division_polys_mod(&f, &e, 160, &m);
        let selected = kernel::division_polys_selected(&f, &e, &indices, &m);
        for n in indices {
            assert_eq!(selected[&n], full[n], "index {n}");
        }
    }
    #[test]
    fn blocked_composition_matches_original_horner() {
        let f = field::Field::new(&BigUint::from(1009u32)).unwrap();
        for size in [8, 31, 128] {
            let mut m: Vec<_> = (0..size).map(|i| f.from_u64((i * 19 + 5) as u64)).collect();
            m.push(f.one());
            let a: Vec<_> = (0..size + 5)
                .map(|i| f.from_u64((i * i + 13) as u64))
                .collect();
            let x: Vec<_> = (0..size).map(|i| f.from_u64((i * 7 + 2) as u64)).collect();
            assert_eq!(
                kernel::compose_mod(&f, &a, &x, &m),
                original_poly::eval_mod(&f, &a, &x, &m)
            );
        }
    }
    #[test]
    fn known_kernel_and_mutation_keep_the_original_verdict() {
        let f = field::Field::new(&BigUint::from(101u32)).unwrap();
        let e = curve::Model {
            a: f.zero(),
            b: f.from_u64(2),
        };
        let h = vec![f.zero(), f.one()];
        assert_eq!(
            kernel::verify_kernel(&f, &e, 3, &h).unwrap(),
            original_kernel::verify_kernel(&f, &e, 3, &h).unwrap()
        );
        let changed = vec![f.one(), f.one()];
        assert!(kernel::verify_kernel(&f, &e, 3, &changed).is_err());
        assert!(original_kernel::verify_kernel(&f, &e, 3, &changed).is_err());
    }
}
