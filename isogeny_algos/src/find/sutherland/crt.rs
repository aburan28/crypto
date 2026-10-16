//! Streaming explicit CRT (BLS, section 6.2). Integer fixed-point rounding avoids
//! a floating-point precision assumption. Only one residue matrix is needed.
use crate::bigint::Big;
use crate::field::{is_prime, Field, Zp};
use crate::int::Int;
use std::collections::HashSet;

/// Reconstruct signed coefficients with |c| <= 2^bound_bits directly modulo m.
/// The caller must establish this bound for the mathematical coefficients.
/// Primes are consumed once, in the supplied order; m may share factors with M.
pub struct ExplicitCrt {
    primes: Vec<u64>,
    modulus: Int,
    product_mod: Int,
    weights: Vec<Int>,
    inverses: Vec<u64>,
    sums: Vec<Int>,
    fractions: Vec<u128>,
    fractional_bits: u32,
    next: usize,
}

impl ExplicitCrt {
    pub fn new(
        primes: &[u64],
        modulus: &Int,
        count: usize,
        bound_bits: usize,
    ) -> Result<Self, String> {
        if modulus <= &Int::zero() || primes.is_empty() || primes.len() > 1_000_000 {
            return Err("CRT requires m > 0 and between 1 and 1000000 primes".into());
        }
        if count == 0 || bound_bits > 1_000_000 {
            return Err("invalid CRT coefficient count or bound".into());
        }
        let mut seen = HashSet::new();
        let mut product = Big::from_u64(1);
        for &p in primes {
            if p >= (1 << 63) || !is_prime(p) || !seen.insert(p) {
                return Err("CRT primes must be distinct primes below 2^63".into());
            }
            product = product.mul_small(p);
        }
        if product <= Big::from_u64(1).shl(bound_bits + 2) {
            return Err("CRT product must exceed four times the coefficient bound".into());
        }
        let mut weights = Vec::with_capacity(primes.len());
        let mut inverses = Vec::with_capacity(primes.len());
        for &p in primes {
            let mi = product.divrem_small(p).0;
            inverses.push(Zp::new(p).inv(mi.rem_small(p)));
            weights.push(Int::from_big(&mi).modulo(modulus));
        }
        // n/2^k < 1/8. If u_i = a_i*c_i mod p_i, then s=sum u_i/p_i is
        // within 1/4 of the correction integer, because |c| < M/4.
        // Truncation moves s down by <1/8, hence nearest-integer is unchanged.
        let fractional_bits = usize::BITS - primes.len().leading_zeros() + 3;
        Ok(Self {
            primes: primes.to_vec(),
            modulus: modulus.clone(),
            product_mod: Int::from_big(&product).modulo(modulus),
            weights,
            inverses,
            sums: vec![Int::zero(); count],
            fractions: vec![0; count],
            fractional_bits,
            next: 0,
        })
    }

    /// Reject malformed input before updating any accumulator; a rejected push
    /// can be retried with the correct prime and residue vector.
    pub fn push(&mut self, p: u64, residues: &[u64]) -> Result<(), String> {
        if self.primes.get(self.next) != Some(&p)
            || residues.len() != self.sums.len()
            || residues.iter().any(|&r| r >= p)
        {
            return Err("CRT prime/order, coefficient count, or residue is invalid".into());
        }
        let f = Zp::new(p);
        for (k, &c) in residues.iter().enumerate() {
            let u = f.mul(c, self.inverses[self.next]);
            self.sums[k] =
                (&self.sums[k] + &(&self.weights[self.next] * &Int::from(u))).modulo(&self.modulus);
            self.fractions[k] += ((u as u128) << self.fractional_bits) / p as u128;
        }
        self.next += 1;
        Ok(())
    }

    pub fn finish(self) -> Result<Vec<Int>, String> {
        if self.next != self.primes.len() {
            return Err("CRT reconstruction is incomplete".into());
        }
        let half = 1u128 << (self.fractional_bits - 1);
        Ok(self
            .sums
            .into_iter()
            .zip(self.fractions)
            .map(|(sum, fraction)| {
                let correction = (fraction + half) >> self.fractional_bits;
                (sum - &self.product_mod * Int::from(correction as u64)).modulo(&self.modulus)
            })
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn reconstructs_signed_bounds_for_arbitrary_targets() {
        let primes = [101, 103, 107];
        let values: Vec<_> = [-65536i64, -65535, -1, 0, 1, 65535, 65536]
            .into_iter()
            .map(Int::from)
            .collect();
        for modulus in [
            Int::one(),
            Int::from(101u64 * 12),
            Int::one().shl(200) + Int::from(17u64),
        ] {
            let mut crt = ExplicitCrt::new(&primes, &modulus, values.len(), 16).unwrap();
            for p in primes {
                crt.push(p, &values.iter().map(|v| v.mod_u64(p)).collect::<Vec<_>>())
                    .unwrap();
            }
            assert_eq!(
                crt.finish().unwrap(),
                values
                    .iter()
                    .map(|v| v.modulo(&modulus))
                    .collect::<Vec<_>>()
            );
        }
    }
    #[test]
    fn rejects_bad_setup_and_push_is_transactional() {
        for primes in [&[101, 101][..], &[101, 102][..], &[][..]] {
            assert!(ExplicitCrt::new(primes, &Int::from(12u64), 1, 1).is_err());
        }
        assert!(ExplicitCrt::new(&[101], &Int::one(), 1, 10).is_err());
        assert!(ExplicitCrt::new(&[101], &Int::zero(), 1, 1).is_err());
        assert!(ExplicitCrt::new(&[101], &Int::one(), 0, 1).is_err());
        let mut crt = ExplicitCrt::new(&[101, 103], &Int::from(12u64), 2, 5).unwrap();
        assert!(crt.push(103, &[1, 2]).is_err());
        assert!(crt.push(101, &[1]).is_err());
        assert!(crt.push(101, &[1, 101]).is_err());
        crt.push(101, &[100, 1]).unwrap();
        crt.push(103, &[102, 1]).unwrap();
        assert_eq!(crt.finish().unwrap(), vec![Int::from(11u64), Int::one()]);
        assert!(ExplicitCrt::new(&[101], &Int::one(), 1, 1)
            .unwrap()
            .finish()
            .is_err());
    }
}
