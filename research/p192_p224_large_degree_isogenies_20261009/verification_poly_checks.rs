#[path = "../../src/cryptanalysis/isogeny_walk/field.rs"]
mod field;
#[path = "verification_poly.rs"]
mod optimized;
#[path = "../../src/cryptanalysis/isogeny_walk/poly.rs"]
mod reference;
fn main() {}
#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;
    #[test]
    fn agrees_with_unchanged_schoolbook_on_both_source_fields() {
        for p in [
            "6277101735386680763835789423207666416083908700390324961279",
            "26959946667150639794667015087019630673557916260026308143510066298881",
        ] {
            let p: BigUint = p.parse().unwrap();
            let f = field::Field::new(&p).unwrap();
            for (a, b) in [
                (0, 7),
                (1, 200),
                (47, 49),
                (48, 48),
                (63, 64),
                (65, 49),
                (127, 129),
                (256, 31),
                (257, 83),
                (511, 512),
            ] {
                let x: Vec<_> = (0..a)
                    .map(|i| f.from_big(&(&p - BigUint::from((i * i + 7) as u64))))
                    .collect();
                let y: Vec<_> = (0..b)
                    .map(|i| f.from_big(&BigUint::from((i * 71 + 19) as u64)))
                    .collect();
                assert_eq!(
                    optimized::mul(&f, &x, &y),
                    reference::mul(&f, &x, &y),
                    "lengths {a},{b}"
                );
            }
        }
    }

    #[test]
    fn reciprocal_remainder_matches_original_long_division() {
        let p: BigUint = "26959946667150639794667015087019630673557916260026308143510066298881"
            .parse()
            .unwrap();
        let f = field::Field::new(&p).unwrap();
        for size in [128, 129, 191, 256] {
            for extra in [0, 1, size / 2, size - 1] {
                let mut modulus: Vec<_> =
                    (0..size).map(|i| f.from_u64((17 * i + 3) as u64)).collect();
                modulus[size - 1] = f.one();
                let a: Vec<_> = (0..size + extra)
                    .map(|i| f.from_u64((i * i + 31) as u64))
                    .collect();
                assert_eq!(
                    optimized::rem(&f, &a, &modulus),
                    reference::rem(&f, &a, &modulus)
                );
            }
        }
    }
}
