//! Workloads: one curve, one prime-order subgroup, one generator, one
//! target.
//!
//! A workload is the unit every comparison is paired on.  It is built
//! from a *construction* (a rule the repository already uses to make a
//! curve) rather than from pasted coordinates, so a run can be replayed
//! from its record on any host, and its identity is computed from what
//! was built — the ICV1 model, the subgroup, the generator and the target
//! point — so two constructions that land on the same curve and point
//! are the same workload.
//!
//! The planted logarithm is derived from the workload seed by SHA-256
//! (`canonical::derive_u64`), never from an RNG stream, and no solver is
//! handed it: the runner verifies the solver's answer afterwards.

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};

use crate::cryptanalysis::ecbench::canonical::{bare_id, derive_u64};
use crate::cryptanalysis::ic_boundary::{
    find_prime_order_curve, koblitz_instance, random_binary_instance, roster_prime_instance,
    BinaryGroup, BinaryInstance, CountedGroup, GroupOps, PrimeCurve, PrimeInstance, PrimePoint,
};

/// How to build the curve.  Each variant is a constructor the repository
/// already has, so the registry's `generator_call` names the same curve.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum CurveSpec {
    /// `ic_boundary::find_prime_order_curve(bits, seed)`: prime order,
    /// cofactor 1, `8 ≤ bits ≤ 32`.
    PrimeSearch { bits: u32, seed: u64 },
    /// `ic_boundary::roster_prime_instance(bits)`: the bench roster.
    PrimeRoster { bits: u32 },
    /// `ic_boundary::koblitz_instance(a, n)`: `K_a / F_{2^n}`, `n ≤ 62`.
    Koblitz { a: u8, n: u32 },
    /// `ic_boundary::random_binary_instance(n, seed, max_cofactor)`.
    BinaryRandom {
        n: u32,
        seed: u64,
        max_cofactor: u64,
    },
    /// A prime curve given outright: `y² = x³ + ax + b` over `F_p`, the
    /// group order, the prime subgroup order `r` and a generator of it.
    /// Checked on build (non-singular, generator on the curve, `[r]G =
    /// O`, `r · cofactor = #E`); the order itself is trusted, as the
    /// registry's recorded orders are.  The runner hands measured
    /// children this form so a child never repeats a curve search.
    PrimeExplicit {
        p: u64,
        a: u64,
        b: u64,
        group_order: u64,
        r: u64,
        gx: u64,
        gy: u64,
    },
}

impl CurveSpec {
    /// Build the instance, or say why it cannot be built.
    pub fn build(&self) -> Result<Instance, String> {
        match *self {
            CurveSpec::PrimeSearch { bits, seed } => {
                if !(8..=32).contains(&bits) {
                    return Err(format!("prime_search takes 8 ≤ bits ≤ 32, not {bits}"));
                }
                Ok(Instance::Prime(find_prime_order_curve(bits, seed)))
            }
            CurveSpec::PrimeRoster { bits } => roster_prime_instance(bits)
                .map(Instance::Prime)
                .ok_or_else(|| format!("the bench roster has no {bits}-bit curve")),
            CurveSpec::Koblitz { a, n } => {
                if a > 1 || !(3..=62).contains(&n) {
                    return Err(format!(
                        "koblitz takes a ∈ {{0,1}} and 3 ≤ n ≤ 62, not a={a} n={n}"
                    ));
                }
                koblitz_instance(a, n).map(|i| Instance::Binary(Box::new(i))).ok_or_else(|| {
                    format!("K_{a} over F_2^{n} has no usable prime-order subgroup in the repository's constructor")
                })
            }
            CurveSpec::BinaryRandom {
                n,
                seed,
                max_cofactor,
            } => random_binary_instance(n, seed, max_cofactor)
                .map(|i| Instance::Binary(Box::new(i)))
                .ok_or_else(|| format!("no binary curve found for n={n} seed={seed}")),
            CurveSpec::PrimeExplicit {
                p,
                a,
                b,
                group_order,
                r,
                gx,
                gy,
            } => {
                if !(5..1 << 62).contains(&p) || a >= p || b >= p {
                    return Err("prime_explicit takes 5 ≤ p < 2^62 and a, b < p".into());
                }
                if r < 2 || group_order % r != 0 {
                    return Err(format!("r = {r} does not divide #E = {group_order}"));
                }
                let curve = PrimeCurve { p, a, b };
                let g = PrimePoint::affine(gx, gy);
                if !curve.is_on_curve(g) {
                    return Err("the generator is not on the curve".into());
                }
                let mut scratch = GroupOps::default();
                if !curve.mul(&mut scratch, g, r).infinity {
                    return Err("[r]G is not the identity".into());
                }
                let mut inst = PrimeInstance {
                    name: String::new(),
                    curve,
                    group_order,
                    r,
                    cofactor: group_order / r,
                    generator: (gx, gy),
                };
                inst.name = inst.curve_id().slug;
                Ok(Instance::Prime(inst))
            }
        }
    }

    /// The same curve as an explicit form, when it has one: what a
    /// measured child rebuilds instead of repeating a search.
    pub fn explicit(inst: &Instance) -> Option<CurveSpec> {
        match inst {
            Instance::Prime(i) => Some(CurveSpec::PrimeExplicit {
                p: i.curve.p,
                a: i.curve.a,
                b: i.curve.b,
                group_order: i.group_order,
                r: i.r,
                gx: i.generator.0,
                gy: i.generator.1,
            }),
            Instance::Binary(_) => None,
        }
    }

    /// The constructor call, in the registry's `generator_call` form.
    pub fn call(&self) -> String {
        match *self {
            CurveSpec::PrimeSearch { bits, seed } => {
                format!("find_prime_order_curve({bits}, {seed})")
            }
            CurveSpec::PrimeRoster { bits } => format!("roster_prime_instance({bits})"),
            CurveSpec::Koblitz { a, n } => format!("KoblitzCurve::new({a}, {n})"),
            CurveSpec::BinaryRandom {
                n,
                seed,
                max_cofactor,
            } => format!("random_binary_instance({n}, {seed}, {max_cofactor})"),
            CurveSpec::PrimeExplicit { p, a, b, .. } => {
                format!("PrimeCurve {{ p: {p}, a: {a}, b: {b} }} (explicit)")
            }
        }
    }
}

/// A built curve instance.
pub enum Instance {
    Prime(PrimeInstance),
    Binary(Box<BinaryInstance>),
}

/// The facts about a built curve every record carries.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CurveFacts {
    /// The ICV1 slug: the curve's name in prose and tables (AGENTS.md §11).
    pub slug: String,
    pub icv1: String,
    /// `prime`, `koblitz` or `binary`.
    pub family: String,
    pub construction: String,
    /// Field degree for binary curves; `None` for prime fields.
    pub field_degree: Option<u32>,
    /// Bits of the field characteristic for prime curves.
    pub field_bits: Option<u32>,
    pub group_order: u64,
    pub r: u64,
    pub cofactor: u64,
    /// Generator coordinates as lower-case hex.
    pub generator: [String; 2],
    /// Automorphisms a generic algorithm may use: 2 (negation), or `2n`
    /// on a Koblitz curve (signed Frobenius).
    pub automorphisms_available: u32,
    /// The slug is in `docs/curves/registry.json`.
    pub registered: bool,
    /// EC1 alias and curve UID of the registered representation whose
    /// subgroup and generator are exactly this one; `None` when the
    /// registry has no such representation (never inferred).
    pub ec1: Option<String>,
    pub curve_uid: Option<String>,
}

const REGISTRY: &str = include_str!("../../../docs/curves/registry.json");

/// The registry entry of `slug`, if any.
fn registry_entry(slug: &str) -> Option<Value> {
    let doc: Value = serde_json::from_str(REGISTRY).ok()?;
    doc.get("curves")?
        .as_array()?
        .iter()
        .find(|c| c.get("slug").and_then(Value::as_str) == Some(slug))
        .cloned()
}

fn parse_hex(s: &str) -> Option<u64> {
    u64::from_str_radix(s.trim_start_matches("0x"), 16).ok()
}

/// `(ec1, curve_uid)` of the representation with this subgroup and
/// generator.
fn matching_representation(entry: &Value, r: u64, gx: u64, gy: u64) -> Option<(String, String)> {
    entry
        .get("representations")?
        .as_array()?
        .iter()
        .find_map(|rep| {
            let c = rep.get("curve")?;
            let order: u64 = c.get("subgroup_order")?.as_str()?.parse().ok()?;
            let g = c.get("generator")?.as_array()?;
            let x = parse_hex(g.first()?.as_str()?)?;
            let y = parse_hex(g.get(1)?.as_str()?)?;
            if order != r || x != gx || y != gy {
                return None;
            }
            Some((
                rep.get("ec1")?.as_str()?.to_string(),
                rep.get("curve_uid")?.as_str()?.to_string(),
            ))
        })
}

impl Instance {
    pub fn r(&self) -> u64 {
        match self {
            Instance::Prime(i) => i.r,
            Instance::Binary(i) => i.r,
        }
    }

    pub fn facts(&self, spec: &CurveSpec) -> CurveFacts {
        let (id, family, degree, bits, order, r, h, gx, gy, aut) = match self {
            Instance::Prime(i) => (
                i.curve_id(),
                "prime",
                None,
                Some(64 - i.curve.p.leading_zeros()),
                i.group_order,
                i.r,
                i.cofactor,
                i.generator.0,
                i.generator.1,
                2,
            ),
            Instance::Binary(i) => (
                i.curve_id(),
                if i.koblitz.is_some() {
                    "koblitz"
                } else {
                    "binary"
                },
                Some(i.n),
                None,
                i.group_order,
                i.r,
                i.cofactor,
                i.generator.x,
                i.generator.y,
                if i.koblitz.is_some() { 2 * i.n } else { 2 },
            ),
        };
        let entry = registry_entry(&id.slug);
        let rep = entry
            .as_ref()
            .and_then(|e| matching_representation(e, r, gx, gy));
        CurveFacts {
            slug: id.slug,
            icv1: id.icv1,
            family: family.into(),
            construction: spec.call(),
            field_degree: degree,
            field_bits: bits,
            group_order: order,
            r,
            cofactor: h,
            generator: [format!("0x{gx:x}"), format!("0x{gy:x}")],
            automorphisms_available: aut,
            registered: entry.is_some(),
            ec1: rep.as_ref().map(|(a, _)| a.clone()),
            curve_uid: rep.map(|(_, u)| u),
        }
    }

    /// `[k]G` as hex coordinates, computed by the instance's own counted
    /// double-and-add on a ledger nobody reads.  Used to build the
    /// target and, separately, to verify an answer.
    pub fn mul_generator_hex(&self, k: u64) -> Option<[String; 2]> {
        let mut scratch = GroupOps::default();
        match self {
            Instance::Prime(i) => {
                let p = i.curve.mul(&mut scratch, i.generator_point(), k);
                (!p.infinity).then(|| [format!("0x{:x}", p.x), format!("0x{:x}", p.y)])
            }
            Instance::Binary(i) => {
                let g = BinaryGroup(&i.fast);
                let p = g.mul(&mut scratch, i.generator, k);
                (!p.infinity).then(|| [format!("0x{:x}", p.x), format!("0x{:x}", p.y)])
            }
        }
    }
}

/// How targets are drawn.  One law today; it is named so a second one
/// (an interval, a structured scalar) cannot be confused with it.
pub const TARGET_LAW: &str = "uniform_scalar_sha256_v1";

/// `k ∈ [1, r−1]` for target `index` under `target_seed`.
pub fn planted_scalar(curve: &CurveFacts, target_seed: u64, index: u64) -> u64 {
    // The curve's model hash and subgroup enter the derivation, so the
    // same seed on two curves gives unrelated scalars.
    // ICV1's last field is the first 12 hex digits of the model hash.
    let model = curve
        .icv1
        .rsplit(':')
        .next()
        .and_then(|h| u64::from_str_radix(h, 16).ok())
        .unwrap_or(0);
    1 + derive_u64(TARGET_LAW, &[model, curve.r, target_seed, index]) % (curve.r - 1)
}

/// One built single-target workload.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Workload {
    /// `W` + 12 hex of the workload record (no `h`, so it composes into
    /// a run id as the repository's convention does).
    pub workload_id: String,
    pub workload_sha256: String,
    pub curve_spec: CurveSpec,
    pub curve: CurveFacts,
    pub target_law: String,
    pub target_seed: u64,
    pub target_index: u64,
    /// `Q = [k]G`, hex.
    pub target: [String; 2],
    /// The planted `k`.  Not in the identity (the target point is) and
    /// never passed to a solver.
    pub planted: u64,
}

impl Workload {
    /// Build workload `index` on `spec` under `target_seed`.  This builds
    /// the curve; for several targets on one curve, build it once and
    /// call [`Workload::on`] for each (a 26-bit prime search counts
    /// points in `O(p)` per candidate curve).
    pub fn build(
        spec: &CurveSpec,
        target_seed: u64,
        index: u64,
    ) -> Result<(Workload, Instance), String> {
        let inst = spec.build()?;
        let w = Self::on(spec, &inst, target_seed, index)?;
        Ok((w, inst))
    }

    /// Workload `index` under `target_seed` on `inst`, already built from
    /// `spec`.
    pub fn on(
        spec: &CurveSpec,
        inst: &Instance,
        target_seed: u64,
        index: u64,
    ) -> Result<Workload, String> {
        let curve = inst.facts(spec);
        if curve.r < 5 {
            return Err(format!(
                "subgroup order {} is too small to measure",
                curve.r
            ));
        }
        let planted = planted_scalar(&curve, target_seed, index);
        let target = inst
            .mul_generator_hex(planted)
            .ok_or("the planted scalar gave the identity; the subgroup order is wrong")?;
        let (workload_id, workload_sha256) = bare_id(
            "W",
            &Self::identity_view(&curve, target_seed, index, &target),
        )?;
        Ok(Workload {
            workload_id,
            workload_sha256,
            curve_spec: spec.clone(),
            curve,
            target_law: TARGET_LAW.into(),
            target_seed,
            target_index: index,
            target,
            planted,
        })
    }

    /// What the workload id is a hash of: the exact curve model, the
    /// subgroup, the generator and the target point, plus the law that
    /// drew it.  One target, cold.
    fn identity_view(
        curve: &CurveFacts,
        target_seed: u64,
        index: u64,
        target: &[String; 2],
    ) -> Value {
        json!({
            "schema": "ecbench.workload/v1",
            "icv1": curve.icv1,
            "subgroup_order": curve.r.to_string(),
            "cofactor": curve.cofactor.to_string(),
            "generator": curve.generator,
            "target": target,
            "target_law": TARGET_LAW,
            "target_seed": target_seed.to_string(),
            "target_index": index.to_string(),
            "target_count": 1,
            "cache_policy": "cold",
        })
    }

    /// Recompute the identity from the fields; an audit compares it.
    pub fn recompute_id(&self) -> Result<String, String> {
        Ok(bare_id(
            "W",
            &Self::identity_view(
                &self.curve,
                self.target_seed,
                self.target_index,
                &self.target,
            ),
        )?
        .0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn registered_prime_curve_resolves_its_slug_and_ec1() {
        let spec = CurveSpec::PrimeSearch {
            bits: 20,
            seed: 59297,
        };
        let (w, _) = Workload::build(&spec, 1, 0).unwrap();
        assert_eq!(w.curve.slug, "icv1-fp20-t727-cd198a38");
        assert!(w.curve.registered);
        assert_eq!(w.curve.ec1.as_deref(), Some("EC1P20Cfphd19e740ea95d"));
        assert!(w.workload_id.starts_with('W') && w.workload_id.len() == 13);
        assert_eq!(w.recompute_id().unwrap(), w.workload_id);
    }

    #[test]
    fn targets_differ_by_index_and_are_stable() {
        let spec = CurveSpec::Koblitz { a: 0, n: 13 };
        let (a, _) = Workload::build(&spec, 7, 0).unwrap();
        let (b, _) = Workload::build(&spec, 7, 1).unwrap();
        let (a2, _) = Workload::build(&spec, 7, 0).unwrap();
        assert_ne!(a.workload_id, b.workload_id);
        assert_eq!(a, a2);
        assert_eq!(a.curve.family, "koblitz");
        assert_eq!(a.curve.automorphisms_available, 26);
    }

    #[test]
    fn explicit_form_rebuilds_the_same_workload() {
        let spec = CurveSpec::PrimeSearch {
            bits: 16,
            seed: 59297,
        };
        let (w, inst) = Workload::build(&spec, 3, 2).unwrap();
        let ex = CurveSpec::explicit(&inst).unwrap();
        let (w2, _) = Workload::build(&ex, 3, 2).unwrap();
        assert_eq!(w.workload_id, w2.workload_id);
        assert_eq!(w.planted, w2.planted);
        let CurveSpec::PrimeExplicit {
            p,
            a,
            b,
            group_order,
            r,
            gx,
            ..
        } = ex
        else {
            unreachable!()
        };
        let bad = CurveSpec::PrimeExplicit {
            p,
            a,
            b,
            group_order,
            r,
            gx,
            gy: 1,
        };
        assert!(bad.build().is_err());
    }

    #[test]
    fn planted_scalars_are_uniform() {
        // Kolmogorov–Smirnov against U(0, 1) for k / r over 20 000 indices,
        // and the mean distance |k − r/2| / r against its 0.25.  The draw
        // is deterministic, so this is a fixed check of the derivation,
        // not a flaky statistical test.
        let (w, _) = Workload::build(&CurveSpec::Koblitz { a: 0, n: 13 }, 0, 0).unwrap();
        let n = 20_000u64;
        let mut xs: Vec<f64> = (0..n)
            .map(|i| planted_scalar(&w.curve, 20261002, i) as f64 / w.curve.r as f64)
            .collect();
        xs.sort_by(f64::total_cmp);
        let d = xs
            .iter()
            .enumerate()
            .map(|(i, x)| ((i + 1) as f64 / n as f64 - x).max(x - i as f64 / n as f64))
            .fold(0.0, f64::max);
        // 1.63 / √n is the 1 % critical value.
        assert!(d < 1.63 / (n as f64).sqrt(), "KS D = {d}");
        let gap = xs.iter().map(|x| (x - 0.5).abs()).sum::<f64>() / n as f64;
        assert!((gap - 0.25).abs() < 0.005, "mean gap {gap}");
    }

    #[test]
    fn unusable_koblitz_curve_is_refused() {
        assert!(CurveSpec::Koblitz { a: 0, n: 17 }.build().is_err());
    }
}
