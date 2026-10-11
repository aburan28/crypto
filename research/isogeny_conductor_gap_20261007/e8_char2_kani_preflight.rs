//! Exact arithmetic preflight for the odd-torsion C37 Kani example.
//!
//! This checks Frobenius on the surface curve's prime-to-characteristic
//! torsion. It does not construct either isogeny or certify its evaluation.
//!
//! `rustc --edition=2021 --test e8_char2_kani_preflight.rs -o /tmp/e8-preflight-test`
//! `/tmp/e8-preflight-test`
//! `rustc --edition=2021 e8_char2_kani_preflight.rs -o /tmp/e8-preflight`
//! `/tmp/e8-preflight`

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct TauElt {
    a: i128,
    b: i128,
}

impl TauElt {
    fn reduce(self, modulus: i128) -> Self {
        Self {
            a: self.a.rem_euclid(modulus),
            b: self.b.rem_euclid(modulus),
        }
    }

    fn mul(self, rhs: Self, modulus: i128) -> Self {
        // tau^2 + tau + 2 = 0 in the maximal order of Q(sqrt(-7)).
        Self {
            a: self.a * rhs.a - 2 * self.b * rhs.b,
            b: self.a * rhs.b + self.b * rhs.a - self.b * rhs.b,
        }
        .reduce(modulus)
    }

    fn pow(mut self, mut exponent: u32, modulus: i128) -> Self {
        let mut result = Self { a: 1, b: 0 };
        while exponent != 0 {
            if exponent & 1 != 0 {
                result = result.mul(self, modulus);
            }
            self = self.mul(self, modulus);
            exponent >>= 1;
        }
        result
    }

    fn trace(self, modulus: i128) -> i128 {
        (2 * self.a - self.b).rem_euclid(modulus)
    }

    fn norm(self, modulus: i128) -> i128 {
        (self.a * self.a - self.a * self.b + 2 * self.b * self.b).rem_euclid(modulus)
    }
}

fn full_torsion_degree(field_degree: u32, modulus: i128) -> u32 {
    let pi = TauElt { a: 0, b: 1 }.pow(field_degree, modulus);
    let mut power = TauElt { a: 1, b: 0 };
    for degree in 1..=(6 * modulus as u32 * modulus as u32) {
        power = power.mul(pi, modulus);
        if power == (TauElt { a: 1, b: 0 }) {
            return degree;
        }
    }
    panic!("Frobenius order exceeded the explicit search cap");
}

fn pow_mod(mut base: i128, mut exponent: u32, modulus: i128) -> i128 {
    let mut result = 1;
    while exponent != 0 {
        if exponent & 1 != 0 {
            result = result * base % modulus;
        }
        base = base * base % modulus;
        exponent >>= 1;
    }
    result
}

fn main() {
    let n = 37;
    let ell = 73;
    let smooth_m = 81;
    let auxiliary_degree = smooth_m - ell;
    let pi9 = TauElt { a: 0, b: 1 }.pow(n, 9);
    let pi81 = TauElt { a: 0, b: 1 }.pow(n, smooth_m);
    let torsion9 = full_torsion_degree(n, 9);
    let torsion81 = full_torsion_degree(n, smooth_m);
    assert_eq!(auxiliary_degree, 2 * 2 + 2 * 2);
    assert_eq!(pi81.trace(smooth_m), (-534_059_i128).rem_euclid(smooth_m));
    assert_eq!(pi81.norm(smooth_m), pow_mod(2, n, smooth_m));
    println!(
        "{{\"curve\":\"C37 surface K0\",\"n\":{n},\"vertical_prime\":{ell},\
         \"smooth_M\":{smooth_m},\"auxiliary_m\":{auxiliary_degree},\
         \"M1\":9,\"M2\":9,\"pi_mod_9\":[{},{}],\"pi_mod_81\":[{},{}],\
         \"full_9_torsion_extension_degree\":{torsion9},\
         \"full_81_torsion_extension_degree\":{torsion81},\
         \"auxiliary_degree_coprime_to_characteristic\":false,\
         \"constructs_isogeny\":false}}",
        pi9.a, pi9.b, pi81.a, pi81.b
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tau_relation_and_c37_trace() {
        for modulus in [9, 81] {
            let tau = TauElt { a: 0, b: 1 };
            assert_eq!(
                tau.mul(tau, modulus),
                TauElt { a: -2, b: -1 }.reduce(modulus)
            );
            let pi = tau.pow(37, modulus);
            assert_eq!(pi.trace(modulus), (-534_059_i128).rem_euclid(modulus));
            assert_eq!(pi.norm(modulus), pow_mod(2, 37, modulus));
        }
    }

    #[test]
    fn exact_frobenius_orders() {
        assert_eq!(full_torsion_degree(37, 9), 24);
        assert_eq!(full_torsion_degree(37, 81), 216);
        let pi = TauElt { a: 0, b: 1 }.pow(37, 9);
        for divisor in [1, 2, 3, 4, 6, 8, 12] {
            assert_ne!(pi.pow(divisor, 9), TauElt { a: 1, b: 0 });
        }
    }

    #[test]
    fn parity_gate() {
        let vertical_prime = 73;
        for smooth_m in [79, 81, 83] {
            let auxiliary = smooth_m - vertical_prime;
            assert_eq!(auxiliary % 2, 0);
        }
        // If the auxiliary degree is odd, M=N+m has even degree and
        // E[M] cannot have the full etale rank-two 2-primary basis.
        assert_eq!((vertical_prime + 7) % 2, 0);
    }
}
