use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

use super::Result;

pub const P_DEC: &str = "6277101735386680763835789423207666416083908700390324961279";
pub const N_DEC: &str = "6277101735386680763835789423176059013767194773182842284081";
pub const T_DEC: &str = "31607402316713927207482677199";
pub const D_DEC: &str = "-24109379060336110122544161233113975664949272517896865359515";
pub const C_DEC: &str = "14140398275856956083603613626459809774163796198179979683";
pub const LARGE_RESIDUAL_BOUND_SQUARED_DEC: &str = "4611686014132420609";

pub fn integer(text: &str) -> Result<BigInt> {
    BigInt::parse_bytes(text.as_bytes(), 10).ok_or_else(|| format!("invalid integer: {text}"))
}

pub fn unsigned(text: &str) -> Result<BigUint> {
    BigUint::parse_bytes(text.as_bytes(), 10)
        .ok_or_else(|| format!("invalid unsigned integer: {text}"))
}

pub fn p192_parameters() -> Result<(BigInt, BigInt, BigInt)> {
    Ok((integer(P_DEC)?, integer(T_DEC)?, integer(D_DEC)?))
}

pub fn centered_u(x: u64, v: u64, trace: &BigInt) -> Result<BigInt> {
    let numerator = BigInt::from(x) - trace * BigInt::from(v);
    if numerator.is_odd() {
        return Err(format!("x-t*v is odd for (v,x)=({v},{x})"));
    }
    Ok(numerator / 2)
}

pub fn norm(u: &BigInt, v: &BigInt, trace: &BigInt, field_norm: &BigInt) -> BigInt {
    u * u + trace * u * v + field_norm * v * v
}

pub fn primitive(u: &BigInt, v: &BigInt) -> bool {
    u.gcd(v).abs().is_one()
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SmallFactorization {
    pub factors: Vec<(u64, u32)>,
    pub residual: BigUint,
}

pub fn factor_over(value: &BigInt, primes: &[u64]) -> Result<SmallFactorization> {
    let mut residual = value
        .to_biguint()
        .ok_or_else(|| "factorization input must be positive".to_owned())?;
    let mut factors = Vec::new();
    let mut previous = 0u64;
    for &prime in primes {
        if prime <= previous || prime < 2 {
            return Err("factor-base primes must be strictly increasing".to_owned());
        }
        previous = prime;
        let mut exponent = 0u32;
        while (&residual % prime).is_zero() {
            residual /= prime;
            exponent = exponent
                .checked_add(1)
                .ok_or_else(|| "prime exponent overflow".to_owned())?;
        }
        if exponent != 0 {
            factors.push((prime, exponent));
        }
    }
    Ok(SmallFactorization { factors, residual })
}

fn extended_gcd(a: &BigInt, b: &BigInt) -> (BigInt, BigInt, BigInt) {
    if b.is_zero() {
        let sign = if a.is_negative() {
            -BigInt::one()
        } else {
            BigInt::one()
        };
        return (a.abs(), sign, BigInt::zero());
    }
    let (gcd, x1, y1) = extended_gcd(b, &a.mod_floor(b));
    (gcd, y1.clone(), x1 - a.div_floor(b) * y1)
}

/// Canonical column HNF with basis vectors `(a,0)` and `(b,d)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Hnf {
    pub a: BigInt,
    pub b: BigInt,
    pub d: BigInt,
}

impl Hnf {
    pub fn identity() -> Self {
        Self {
            a: BigInt::one(),
            b: BigInt::zero(),
            d: BigInt::one(),
        }
    }

    pub fn prime(ell: u64, root: u64) -> Result<Self> {
        if ell < 2 || root >= ell {
            return Err("invalid prime-ideal HNF parameters".to_owned());
        }
        let a = BigInt::from(ell);
        Ok(Self {
            a: a.clone(),
            b: (-BigInt::from(root)).mod_floor(&a),
            d: BigInt::one(),
        })
    }

    fn vectors(&self) -> [(BigInt, BigInt); 2] {
        [
            (self.a.clone(), BigInt::zero()),
            (self.b.clone(), self.d.clone()),
        ]
    }

    pub fn determinant(&self) -> BigInt {
        &self.a * &self.d
    }
}

pub fn hnf_from_vectors(vectors: &[(BigInt, BigInt)]) -> Result<Hnf> {
    if vectors.len() < 2 {
        return Err("a rank-two HNF needs at least two generators".to_owned());
    }
    let mut index = BigInt::zero();
    for left in 0..vectors.len() {
        for right in (left + 1)..vectors.len() {
            let determinant =
                &vectors[left].0 * &vectors[right].1 - &vectors[right].0 * &vectors[left].1;
            index = index.gcd(&determinant.abs());
        }
    }
    if index.is_zero() {
        return Err("ideal generators have rank below two".to_owned());
    }
    let mut gcd_y = BigInt::zero();
    let mut x_combination = BigInt::zero();
    for (x, y) in vectors {
        let (next_gcd, old_coefficient, new_coefficient) = extended_gcd(&gcd_y, y);
        x_combination = old_coefficient * x_combination + new_coefficient * x;
        gcd_y = next_gcd;
    }
    if gcd_y.is_zero() || !index.mod_floor(&gcd_y).is_zero() {
        return Err("inconsistent HNF y-gcd/index".to_owned());
    }
    let a = &index / &gcd_y;
    let b = x_combination.mod_floor(&a);
    Ok(Hnf { a, b, d: gcd_y })
}

pub fn multiply_hnf(left: &Hnf, right: &Hnf, p: &BigInt, t: &BigInt) -> Result<Hnf> {
    let mut products = Vec::with_capacity(4);
    for (x1, y1) in left.vectors() {
        for (x2, y2) in right.vectors() {
            let x = &x1 * &x2 - p * &y1 * &y2;
            let y = &x1 * &y2 + &x2 * &y1 + t * &y1 * &y2;
            products.push((x, y));
        }
    }
    hnf_from_vectors(&products)
}

pub fn pow_hnf(base: &Hnf, exponent: u32, p: &BigInt, t: &BigInt) -> Result<Hnf> {
    let mut accumulator = Hnf::identity();
    let mut power = base.clone();
    let mut remaining = exponent;
    while remaining != 0 {
        if remaining & 1 == 1 {
            accumulator = multiply_hnf(&accumulator, &power, p, t)?;
        }
        remaining >>= 1;
        if remaining != 0 {
            power = multiply_hnf(&power, &power, p, t)?;
        }
    }
    Ok(accumulator)
}

pub fn principal_hnf(u: &BigInt, v: &BigInt, p: &BigInt, t: &BigInt) -> Result<Hnf> {
    hnf_from_vectors(&[(u.clone(), v.clone()), (-v * p, u + v * t)])
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct D23Control {
    pub norm: BigInt,
    pub prime_ideal_cubed: Hnf,
    pub principal_alpha: Hnf,
}

pub fn d23_control() -> Result<D23Control> {
    let p = BigInt::from(6u8);
    let t = BigInt::one();
    let u = BigInt::one();
    let v = BigInt::one();
    let prime = Hnf::prime(2, 1)?;
    let prime_ideal_cubed = pow_hnf(&prime, 3, &p, &t)?;
    let principal_alpha = principal_hnf(&u, &v, &p, &t)?;
    let alpha_norm = norm(&u, &v, &t, &p);
    if alpha_norm != BigInt::from(8u8)
        || prime_ideal_cubed != principal_alpha
        || prime_ideal_cubed.determinant() != alpha_norm
    {
        return Err("D=-23 positive control failed".to_owned());
    }
    // Independently check the advertised form coefficients and discriminants.
    if BigInt::from(1u8).pow(2) - BigInt::from(4u8) * 2 * 3 != BigInt::from(-23)
        || BigInt::from(1u8).pow(2) - BigInt::from(4u8) * 1 * 6 != BigInt::from(-23)
    {
        return Err("D=-23 form fixture has the wrong discriminant".to_owned());
    }
    Ok(D23Control {
        norm: alpha_norm,
        prime_ideal_cubed,
        principal_alpha,
    })
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct P192SqrtDiscriminantControl {
    pub u: BigInt,
    pub v: BigInt,
    pub norm: BigInt,
    pub factorization: SmallFactorization,
    pub residual_exceeds_bound_squared: bool,
}

pub fn p192_sqrt_discriminant_control() -> Result<P192SqrtDiscriminantControl> {
    let (p, t, d) = p192_parameters()?;
    if &t * &t - &p * 4 != d {
        return Err("P-192 D=t^2-4p identity failed".to_owned());
    }
    let u = -&t;
    let v = BigInt::from(2u8);
    if !primitive(&u, &v) {
        return Err("sqrt(D) fixture must be primitive".to_owned());
    }
    let alpha_norm = norm(&u, &v, &t, &p);
    if alpha_norm != d.abs() {
        return Err("sqrt(D) fixture norm is not |D|".to_owned());
    }
    for (ell, root) in [(5u64, 2u64), (11, 3), (31, 7)] {
        let divisor = BigInt::from(ell);
        if (&u + &v * BigInt::from(root)).mod_floor(&divisor) != BigInt::zero() {
            return Err(format!("sqrt(D) root/orientation mismatch at ell={ell}"));
        }
    }
    let factorization = factor_over(&alpha_norm, &[5, 11, 31])?;
    if factorization.factors != [(5, 1), (11, 1), (31, 1)]
        || factorization.residual != unsigned(C_DEC)?
    {
        return Err("sqrt(D) frozen factorization mismatch".to_owned());
    }
    let residual_exceeds_bound_squared =
        factorization.residual > unsigned(LARGE_RESIDUAL_BOUND_SQUARED_DEC)?;
    if !residual_exceeds_bound_squared {
        return Err("sqrt(D) residual unexpectedly crossed the retention bound".to_owned());
    }
    Ok(P192SqrtDiscriminantControl {
        u,
        v,
        norm: alpha_norm,
        factorization,
        residual_exceeds_bound_squared,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn d23_positive_control_is_principal_by_independent_hnf() {
        let control = d23_control().unwrap();
        assert_eq!(control.norm, BigInt::from(8u8));
        assert_eq!(control.prime_ideal_cubed, control.principal_alpha);
        assert_eq!(control.prime_ideal_cubed.a, BigInt::from(8u8));
        assert_eq!(control.prime_ideal_cubed.b, BigInt::one());
        assert_eq!(control.prime_ideal_cubed.d, BigInt::one());
    }

    #[test]
    fn p192_sqrt_discriminant_stays_a_partial() {
        let control = p192_sqrt_discriminant_control().unwrap();
        assert_eq!(control.u, -integer(T_DEC).unwrap());
        assert_eq!(control.v, BigInt::from(2u8));
        assert_eq!(control.norm, integer(D_DEC).unwrap().abs());
        assert_eq!(control.factorization.residual, unsigned(C_DEC).unwrap());
        assert!(control.residual_exceeds_bound_squared);
    }

    #[test]
    fn malformed_factor_bases_fail_closed() {
        assert!(factor_over(&BigInt::from(30u8), &[2, 2, 3]).is_err());
        assert!(factor_over(&BigInt::from(-1), &[2]).is_err());
    }
}
