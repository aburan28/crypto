//! One entry point for the small-x yield on any SECG curve: picks the limb
//! count for the field and runs the binary or prime backend.
//!
//! Prime fields use [`super::prime_small_x::MontField`] with 3, 4, 6 or 9
//! limbs (`secp112…secp192`, `secp224…secp256`, `secp384r1`, `secp521r1`);
//! binary fields use [`super::binary_small_x::Gf2m`] with 2, 3, 4, 5, 7 or 9
//! limbs (`sect113…sect239` and `sect283`, `sect409`, `sect571`).

use super::binary_small_x::{self, BinCurve};
use super::prime_small_x::{self, Curve, MontField, Options, YieldReport};
use super::secg_standards::{SecgCurve, SecgField};

/// Measure the exact small-x yield of `c`. With `raw_samples > 0` the share
/// of raw small-x pair sums that land in the prime-order subgroup is also
/// sampled (see [`prime_small_x::raw_subgroup_fraction`]).
pub fn measure_secg(
    c: &SecgCurve,
    opts: &Options,
    raw_samples: u64,
) -> Result<YieldReport, String> {
    macro_rules! binary {
        ($l:literal) => {{
            let cv = BinCurve::<$l>::from_secg(c)?;
            let mut r = binary_small_x::measure(&cv, opts)?;
            if raw_samples > 0 {
                r.raw_base = Some(binary_small_x::raw_subgroup_fraction(
                    &cv,
                    opts.size,
                    raw_samples,
                    opts.seed,
                )?);
            }
            Ok(r)
        }};
    }
    macro_rules! prime {
        ($l:literal) => {{
            let cv = Curve::<MontField<$l>>::from_secg(c)?;
            let mut r = prime_small_x::measure(&cv, opts)?;
            if raw_samples > 0 {
                r.raw_base = Some(prime_small_x::raw_subgroup_fraction(
                    &cv,
                    opts.size,
                    raw_samples,
                    opts.seed,
                )?);
            }
            Ok(r)
        }};
    }
    match &c.field {
        SecgField::Binary { degree, .. } => match degree {
            0..=128 => binary!(2),
            129..=192 => binary!(3),
            193..=256 => binary!(4),
            257..=320 => binary!(5),
            321..=448 => binary!(7),
            449..=576 => binary!(9),
            _ => Err(format!("{}: degree {degree} is beyond 576", c.name)),
        },
        SecgField::Prime { p } => match p.bits() {
            0..=192 => prime!(3),
            193..=256 => prime!(4),
            257..=384 => prime!(6),
            385..=576 => prime!(9),
            b => Err(format!("{}: {b} bits is beyond 576", c.name)),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::secg_standards::{secg, secg_names};

    #[test]
    fn every_secg_curve_builds_and_has_its_generator_on_the_curve() {
        for name in secg_names() {
            let c = secg(name).unwrap();
            // A tiny measurement exercises the field, the curve and the base.
            let r = measure_secg(
                c,
                &Options {
                    size: 8,
                    trials: 20,
                    controls: 4,
                    seed: 1,
                },
                0,
            )
            .unwrap_or_else(|e| panic!("{name}: {e}"));
            assert_eq!(r.cofactor, c.cofactor, "{name}");
            assert_eq!(r.controls_recovered, 4, "{name}");
            assert_eq!(r.collisions, 0, "{name}");
        }
    }
}
