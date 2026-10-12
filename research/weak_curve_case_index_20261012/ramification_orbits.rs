//! Replays Procedure 1's binary ramification-mask screen in IMC 2015, Sec. 3.3.
//! A mask is only a combinatorial candidate; it is not a curve or covering witness.
use std::collections::BTreeSet;

#[derive(Clone, Copy)]
struct Factor {
    d: u8,
    n: u8,
    f: u16,
}

const FACTORS: [Factor; 10] = [
    Factor {
        d: 2,
        n: 2,
        f: 0b101,
    },
    Factor {
        d: 3,
        n: 2,
        f: 0b111,
    },
    Factor {
        d: 3,
        n: 3,
        f: 0b1001,
    },
    Factor {
        d: 5,
        n: 4,
        f: 0b11111,
    },
    Factor {
        d: 5,
        n: 5,
        f: 0b100001,
    },
    Factor {
        d: 7,
        n: 3,
        f: 0b1011,
    },
    Factor {
        d: 7,
        n: 3,
        f: 0b1101,
    },
    Factor {
        d: 7,
        n: 4,
        f: 0b11101,
    },
    Factor {
        d: 7,
        n: 4,
        f: 0b10111,
    },
    Factor {
        d: 7,
        n: 6,
        f: 0b1111111,
    },
];

const FULL_DEGREE_SEVEN: Factor = Factor {
    d: 7,
    n: 7,
    f: 0b10000001,
};

fn degree(x: u16) -> u32 {
    15 - x.leading_zeros()
}

fn rem(mut a: u16, b: u16) -> u16 {
    assert_ne!(b, 0);
    while a != 0 && degree(a) >= degree(b) {
        a ^= b << (degree(a) - degree(b));
    }
    a
}

fn div_exact(mut a: u16, b: u16) -> u16 {
    let mut q = 0;
    while a != 0 && degree(a) >= degree(b) {
        let shift = degree(a) - degree(b);
        q ^= 1 << shift;
        a ^= b << shift;
    }
    assert_eq!(a, 0, "F must divide X^d+1");
    q
}

fn gcd(mut a: u16, mut b: u16) -> u16 {
    while b != 0 {
        (a, b) = (b, rem(a, b));
    }
    a
}

fn mul(mut a: u16, mut b: u16) -> u16 {
    let mut out = 0;
    while b != 0 {
        if b & 1 != 0 {
            out ^= a;
        }
        a <<= 1;
        b >>= 1;
    }
    out
}

fn rotate(mask: u16, d: u8) -> u16 {
    let field_mask = (1u16 << d) - 1;
    ((mask << 1) & field_mask) | (mask >> (d - 1))
}

fn canonical(mask: u16, d: u8) -> u16 {
    let mut smallest = mask;
    let mut current = mask;
    for _ in 1..d {
        current = rotate(current, d);
        smallest = smallest.min(current);
    }
    smallest
}

fn orbits(factor: Factor) -> BTreeSet<u16> {
    assert_eq!(degree(factor.f), u32::from(factor.n));
    let modulus = (1u16 << factor.d) | 1;
    let fhat = div_exact(modulus, factor.f);
    let mut out = BTreeSet::new();
    for a in 1..(1u16 << factor.n) {
        if gcd(a, factor.f) == 1 {
            let phi = mul(a, fhat);
            assert!(phi < 1u16 << factor.d);
            out.insert(canonical(phi, factor.d));
        }
    }
    out
}

fn main() {
    println!("d,n,f_hex,mask,weight");
    for factor in FACTORS.into_iter().chain([FULL_DEGREE_SEVEN]) {
        for mask in orbits(factor) {
            println!(
                "{},{},0x{:x},{},{}",
                factor.d,
                factor.n,
                factor.f,
                mask,
                mask.count_ones()
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_factor_orbit_counts() {
        let expected = [1, 1, 1, 3, 3, 1, 1, 1, 1, 7, 7];
        let found: Vec<usize> = FACTORS
            .into_iter()
            .chain([FULL_DEGREE_SEVEN])
            .map(|f| orbits(f).len())
            .collect();
        assert_eq!(found, expected);
    }

    #[test]
    fn cubic_two_generator_mask_has_two_branch_points() {
        let factor = FACTORS[1];
        assert_eq!(orbits(factor), BTreeSet::from([0b011]));
    }

    #[test]
    fn every_selected_factor_divides_frobenius_cycle_polynomial() {
        for f in FACTORS.into_iter().chain([FULL_DEGREE_SEVEN]) {
            assert_eq!(rem((1u16 << f.d) | 1, f.f), 0);
        }
    }
}
