//! Deterministic algebraic factor-base construction for the frozen order.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::ToPrimitive;
use serde_json::{json, Value};

use super::arithmetic::{integer, D_DEC, P_DEC, T_DEC};
use super::Result;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Entry {
    pub ell: u64,
    /// 0 is split and 1 is ramified, exactly as serialized by the protocol.
    pub kind: u8,
    pub kronecker: i8,
    pub roots: Vec<u64>,
}

pub fn mul_mod(left: u64, right: u64, modulus: u64) -> u64 {
    ((left as u128 * right as u128) % modulus as u128) as u64
}

pub fn pow_mod(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut accumulator = 1 % modulus;
    while exponent != 0 {
        if exponent & 1 == 1 {
            accumulator = mul_mod(accumulator, base, modulus);
        }
        exponent >>= 1;
        if exponent != 0 {
            base = mul_mod(base, base, modulus);
        }
    }
    accumulator
}

pub fn mod_u64(value: &BigInt, modulus: u64) -> u64 {
    value
        .mod_floor(&BigInt::from(modulus))
        .to_u64()
        .expect("reduced nonnegative residue fits u64")
}

pub fn is_characteristic_root(root: u64, ell: u64, trace: &BigInt, field_norm: &BigInt) -> bool {
    if ell < 2 || root >= ell {
        return false;
    }
    let modulus = ell as u128;
    let root = root as u128;
    let trace = mod_u64(trace, ell) as u128;
    let field_norm = mod_u64(field_norm, ell) as u128;
    (root * root + modulus - (trace * root) % modulus + field_norm).is_multiple_of(modulus)
}

/// Tonelli--Shanks for an odd prime, returning the smaller root.
fn sqrt_mod_prime(value: u64, prime: u64) -> Option<u64> {
    debug_assert!(prime > 2 && prime % 2 == 1);
    let value = value % prime;
    if value == 0 {
        return Some(0);
    }
    if pow_mod(value, (prime - 1) / 2, prime) != 1 {
        return None;
    }
    if prime % 4 == 3 {
        let root = pow_mod(value, (prime + 1) / 4, prime);
        return Some(root.min(prime - root));
    }
    let mut odd = prime - 1;
    let mut power = 0u32;
    while odd.is_multiple_of(2) {
        odd /= 2;
        power += 1;
    }
    let non_residue =
        (2..prime).find(|candidate| pow_mod(*candidate, (prime - 1) / 2, prime) == prime - 1)?;
    let mut c = pow_mod(non_residue, odd, prime);
    let mut x = pow_mod(value, odd.div_ceil(2), prime);
    let mut t = pow_mod(value, odd, prime);
    let mut m = power;
    while t != 1 {
        let mut i = 1u32;
        let mut probe = mul_mod(t, t, prime);
        while probe != 1 {
            probe = mul_mod(probe, probe, prime);
            i += 1;
            if i == m {
                return None;
            }
        }
        let b = pow_mod(c, 1u64 << (m - i - 1), prime);
        x = mul_mod(x, b, prime);
        let b2 = mul_mod(b, b, prime);
        t = mul_mod(t, b2, prime);
        c = b2;
        m = i;
    }
    Some(x.min(prime - x))
}

pub fn primes_through(limit: u64) -> Result<Vec<u64>> {
    let length = usize::try_from(limit)
        .map_err(|_| "prime limit does not fit usize".to_owned())?
        .checked_add(1)
        .ok_or_else(|| "prime limit overflow".to_owned())?;
    let mut flags = vec![true; length];
    if !flags.is_empty() {
        flags[0] = false;
    }
    if flags.len() > 1 {
        flags[1] = false;
    }
    let mut prime = 2usize;
    while prime <= limit as usize / prime {
        if flags[prime] {
            let mut multiple = prime * prime;
            while multiple < length {
                flags[multiple] = false;
                multiple += prime;
            }
        }
        prime += 1;
    }
    Ok(flags
        .into_iter()
        .enumerate()
        .filter_map(|(value, flag)| flag.then_some(value as u64))
        .collect())
}

/// Roots of X^2-tX+p modulo an odd prime, in increasing order.
pub fn characteristic_roots(
    ell: u64,
    trace: &BigInt,
    field_norm: &BigInt,
    discriminant: &BigInt,
) -> Vec<u64> {
    if ell == 2 {
        return (0..2)
            .filter(|root| {
                let r = *root as u128;
                let tm = mod_u64(trace, 2) as u128;
                let pm = mod_u64(field_norm, 2) as u128;
                (r * r + 2 - (tm * r) % 2 + pm).is_multiple_of(2)
            })
            .collect();
    }
    let sqrt_d = match sqrt_mod_prime(mod_u64(discriminant, ell), ell) {
        Some(root) => root,
        None => return Vec::new(),
    };
    let tm = mod_u64(trace, ell);
    let inv_two = ell.div_ceil(2);
    let mut roots = vec![
        mul_mod((tm + sqrt_d) % ell, inv_two, ell),
        mul_mod((tm + ell - sqrt_d) % ell, inv_two, ell),
    ];
    roots.sort_unstable();
    roots.dedup();
    roots.retain(|root| is_characteristic_root(*root, ell, trace, field_norm));
    roots
}

pub fn build(bound: u64) -> Result<Vec<Entry>> {
    let p = integer(P_DEC)?;
    let t = integer(T_DEC)?;
    let d = integer(D_DEC)?;
    let mut entries = Vec::new();
    for ell in primes_through(bound)? {
        if mod_u64(&p, ell) == 0 {
            continue;
        }
        let roots = characteristic_roots(ell, &t, &p, &d);
        let kronecker = if ell == 2 {
            match mod_u64(&d, 8) {
                1 | 7 => 1,
                3 | 5 => -1,
                _ => 0,
            }
        } else {
            let residue = mod_u64(&d, ell);
            if residue == 0 {
                0
            } else {
                match pow_mod(residue, (ell - 1) / 2, ell) {
                    1 => 1,
                    value if value == ell - 1 => -1,
                    _ => return Err(format!("invalid Legendre result at ell={ell}")),
                }
            }
        };
        let expected_root_count = match kronecker {
            -1 => 0,
            0 => 1,
            1 => 2,
            _ => unreachable!(),
        };
        if roots.len() != expected_root_count {
            return Err(format!(
                "Kronecker/root-count mismatch at ell={ell}: symbol={kronecker}, roots={roots:?}"
            ));
        }
        if !roots.is_empty() {
            let root_sum = if roots.len() == 1 {
                mul_mod(2, roots[0], ell)
            } else {
                (roots[0] + roots[1]) % ell
            };
            let root_product = if roots.len() == 1 {
                mul_mod(roots[0], roots[0], ell)
            } else {
                mul_mod(roots[0], roots[1], ell)
            };
            if root_sum != mod_u64(&t, ell) || root_product != mod_u64(&p, ell) {
                return Err(format!(
                    "characteristic root sum/product mismatch at ell={ell}"
                ));
            }
        }
        match kronecker {
            -1 => {}
            0 => entries.push(Entry {
                ell,
                kind: 1,
                kronecker: 0,
                roots,
            }),
            1 => entries.push(Entry {
                ell,
                kind: 0,
                kronecker: 1,
                roots,
            }),
            _ => unreachable!(),
        }
    }
    Ok(entries)
}

pub fn manifest(source_commit: &str, algebraic_bound: u64, mappable_bound: u64) -> Result<Value> {
    let entries = build(algebraic_bound)?;
    let records = entries
        .iter()
        .map(|entry| {
            json!({
                "ell": entry.ell,
                "kind": entry.kind,
                "kronecker": entry.kronecker,
                "positive_root_index": 0,
                "roots": entry.roots,
            })
        })
        .collect::<Vec<_>>();
    Ok(json!({
        "D": D_DEC,
        "algebraic_bound": algebraic_bound,
        "curve_uid": super::schema::CURVE_UID,
        "entries": records,
        "entry_count": entries.len() as u64,
        "mappable_bound": mappable_bound,
        "p": P_DEC,
        "schema": "p192-wcm-factor-base-v1",
        "source_commit": source_commit,
        "t": T_DEC,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frozen_mappable_members_and_roots() {
        let entries = build(113).unwrap();
        assert_eq!(
            entries.iter().map(|entry| entry.ell).collect::<Vec<_>>(),
            vec![5, 11, 13, 23, 31, 37, 43, 73, 89, 101, 103, 107, 113]
        );
        assert_eq!(entries[0].kind, 1);
        assert_eq!(entries[0].roots, vec![2]);
        assert_eq!(entries[2].roots, vec![2, 5]);
    }

    #[test]
    fn manifest_first_entry_is_the_frozen_fixture() {
        let manifest = manifest(&"a".repeat(40), 113, 113).unwrap();
        assert_eq!(
            super::super::schema::canonical_json_bytes(&manifest["entries"][0]).unwrap(),
            br#"{"ell":5,"kind":1,"kronecker":0,"positive_root_index":0,"roots":[2]}"#
        );
    }
}
