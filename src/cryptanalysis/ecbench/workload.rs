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

use crate::binary_ecc::IrreduciblePoly;
use crate::cryptanalysis::ecbench::canonical::{bare_id, compat_u128, derive_u128, derive_u64};
use crate::cryptanalysis::ecbench::conductor_gap::{self, GapLevel};
use crate::cryptanalysis::ecbench::fp4::Fp4Instance;
use crate::cryptanalysis::ecbench_large_prime::binary_explicit_instance;
use crate::cryptanalysis::ic_boundary::{
    find_prime_order_curve, koblitz_instance, random_binary_instance, roster_prime_instance,
    BinaryGroup, BinaryInstance, CountedGroup, GroupOps, PrimeCurve, PrimeInstance, PrimePoint,
};
use crate::cryptanalysis::koblitz_strong_rho::{RawPointG, WideStrongRho};
use crate::cryptanalysis::koblitz_wide::WideKoblitz;

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
    /// `ic_boundary::koblitz_instance(a, n)`: `K_a / F_{2^n}`, `n ≤ 62`;
    /// for `64 ≤ n ≤ 127`, `koblitz_wide::WideKoblitz` over the registry's
    /// modulus for that `K_a` (the m = 83 gate curve of AGENTS.md §8a).
    Koblitz { a: u8, n: u32 },
    /// `ic_boundary::random_binary_instance(n, seed, max_cofactor)`.
    BinaryRandom {
        n: u32,
        seed: u64,
        max_cofactor: u64,
    },
    /// `conductor_gap::{crater, floor}(params(ell, t), seed)`: a curve of
    /// the CM pair over `GF((t² + 19ℓ²)/4)` with prime order `q + 1 − t`,
    /// at the crater (endomorphism conductor 1) or the floor (conductor
    /// `ℓ`), certified from the Frobenius discriminant `−19ℓ²` alone.
    /// `seed` drives the floor search and is ignored at the crater.
    ConductorGap {
        ell: u64,
        t: u64,
        level: GapLevel,
        seed: u64,
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
    /// An externally supplied public point, with no logarithm in the fixture.
    /// Restricted to cofactor-one prime-order regression cohorts. Primality
    /// and point-count certificates accompany the input catalog.
    PrimePublicFixture {
        p: u64,
        a: u64,
        b: u64,
        group_order: u64,
        r: u64,
        gx: u64,
        gy: u64,
        qx: u64,
        qy: u64,
    },
    /// A binary curve given outright in polynomial basis. `modulus` is
    /// the complete integer bitmask, including its degree-`n` leading bit.
    /// The constructor checks irreducibility, subgroup divisibility,
    /// generator membership, and `[r]G = O`.
    BinaryExplicit {
        n: u32,
        modulus: u64,
        a: u64,
        b: u64,
        group_order: u64,
        r: u64,
        gx: u64,
        gy: u64,
    },
    /// `ecbench::fp4::Fp4Instance::build(p, seed)`: a curve over
    /// `F_{p⁴}` with rational 2-torsion `T = (x₀, 0)`, `f'(x₀) = c²`
    /// (`quartic_folds::generate_translation_instance(p, seed, 16)`),
    /// `p ≡ 1 (mod 4)`, `p < 2^15`.  Unregistered: ICV1 has no
    /// extension-field form.
    Fp4Translation { p: u64, seed: u64 },
}

impl CurveSpec {
    /// Preserve an external public point when making a measured-child input.
    pub fn for_child(&self, inst: &Instance) -> CurveSpec {
        if matches!(self, Self::PrimePublicFixture { .. }) {
            self.clone()
        } else {
            Self::explicit(inst).unwrap_or_else(|| self.clone())
        }
    }
    /// Build the instance, or say why it cannot be built.
    pub fn build(&self) -> Result<Instance, String> {
        match *self {
            CurveSpec::PrimePublicFixture {
                p,
                a,
                b,
                group_order,
                r,
                gx,
                gy,
                qx,
                qy,
            } => {
                if r != group_order {
                    return Err("public fixture requires cofactor one".into());
                }
                if p < 5
                    || !crate::prime_hyperelliptic::fp2::is_prime_u64(p)
                    || !crate::prime_hyperelliptic::fp2::is_prime_u64(r)
                {
                    return Err(
                        "public fixture requires a prime field and prime group order".into(),
                    );
                }
                let trace = p as i128 + 1 - r as i128;
                if trace.unsigned_abs().pow(2) > 4 * p as u128
                    || (r as u128).pow(2) <= 16 * p as u128
                {
                    return Err(
                        "public fixture order needs a Hasse-unique prime-order witness".into(),
                    );
                }
                let a2 = (a as u128 * a as u128) % p as u128;
                let a3 = (a2 * a as u128) % p as u128;
                let b2 = (b as u128 * b as u128) % p as u128;
                if (4 * a3 + 27 * b2).is_multiple_of(p as u128) || gx >= p || gy >= p {
                    return Err(
                        "public fixture requires a nonsingular curve and canonical generator"
                            .into(),
                    );
                }
                let inst = CurveSpec::PrimeExplicit {
                    p,
                    a,
                    b,
                    group_order,
                    r,
                    gx,
                    gy,
                }
                .build()?;
                let Instance::Prime(ref i) = inst else {
                    unreachable!()
                };
                let q = PrimePoint::affine(qx, qy);
                if qx >= p
                    || qy >= p
                    || !i.curve.is_on_curve(q)
                    || !i.curve.mul(&mut GroupOps::default(), q, r).infinity
                {
                    return Err("fixture target is not in the prime subgroup".into());
                }
                Ok(inst)
            }
            CurveSpec::PrimeSearch { bits, seed } => {
                if !(8..=32).contains(&bits) {
                    return Err(format!("prime_search takes 8 ≤ bits ≤ 32, not {bits}"));
                }
                Ok(Instance::Prime(find_prime_order_curve(bits, seed)))
            }
            CurveSpec::PrimeRoster { bits } => roster_prime_instance(bits)
                .map(Instance::Prime)
                .ok_or_else(|| format!("the bench roster has no {bits}-bit curve")),
            CurveSpec::Koblitz { a, n } if (64..=127).contains(&n) && a <= 1 => {
                let irr = registry_koblitz_modulus(a, n).ok_or_else(|| {
                    format!("K_{a} over F_2^{n} is not in docs/curves/registry.json; register it first (AGENTS.md §11)")
                })?;
                let kc = WideKoblitz::new(a, &irr)?;
                let rho = kc.strong_rho();
                Ok(Instance::Wide(Box::new(WideInstance { kc, rho })))
            }
            CurveSpec::Koblitz { a, n } => {
                if a > 1 || !(3..=62).contains(&n) {
                    return Err(format!(
                        "koblitz takes a ∈ {{0,1}} and 3 ≤ n ≤ 62, or 64 ≤ n ≤ 127 (wide), not a={a} n={n}"
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
            CurveSpec::ConductorGap {
                ell,
                t,
                level,
                seed,
            } => {
                let p = conductor_gap::params(ell, t)?;
                match level {
                    GapLevel::Crater => conductor_gap::crater(p),
                    GapLevel::Floor => conductor_gap::floor(p, seed),
                }
                .map(Instance::Prime)
            }
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
            CurveSpec::BinaryExplicit {
                n,
                modulus,
                a,
                b,
                group_order,
                r,
                gx,
                gy,
            } => binary_explicit_instance(n, modulus, a, b, group_order, r, gx, gy)
                .map(|i| Instance::Binary(Box::new(i))),
            CurveSpec::Fp4Translation { p, seed } => {
                Fp4Instance::build(p, seed).map(|i| Instance::Fp4(Box::new(i)))
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
            Instance::Binary(i) if i.koblitz.is_none() => {
                let modulus = i
                    .irreducible
                    .low_terms
                    .iter()
                    .fold(1u64 << i.n, |v, &term| v | (1u64 << term));
                Some(CurveSpec::BinaryExplicit {
                    n: i.n,
                    modulus,
                    a: i.a,
                    b: i.b,
                    group_order: i.group_order,
                    r: i.r,
                    gx: i.generator.x,
                    gy: i.generator.y,
                })
            }
            // Preserve the original Koblitz constructor: it carries the
            // certified Frobenius action and endomorphism identity that a
            // generic explicit binary instance intentionally does not claim.
            Instance::Binary(_) | Instance::Wide(_) | Instance::Fp4(_) => None,
        }
    }

    /// The constructor call, in the registry's `generator_call` form.
    pub fn call(&self) -> String {
        match *self {
            CurveSpec::PrimeSearch { bits, seed } => {
                format!("find_prime_order_curve({bits}, {seed})")
            }
            CurveSpec::PrimeRoster { bits } => format!("roster_prime_instance({bits})"),
            CurveSpec::Koblitz { a, n } if n >= 64 => {
                format!("WideKoblitz::new({a}, registry modulus for n = {n})")
            }
            CurveSpec::Koblitz { a, n } => format!("KoblitzCurve::new({a}, {n})"),
            CurveSpec::BinaryRandom {
                n,
                seed,
                max_cofactor,
            } => format!("random_binary_instance({n}, {seed}, {max_cofactor})"),
            CurveSpec::ConductorGap {
                ell,
                t,
                level,
                seed,
            } => match level {
                GapLevel::Crater => format!("conductor_gap::crater(params({ell}, {t}))"),
                GapLevel::Floor => format!("conductor_gap::floor(params({ell}, {t}), {seed})"),
            },
            CurveSpec::PrimeExplicit { p, a, b, .. }
            | CurveSpec::PrimePublicFixture { p, a, b, .. } => {
                format!("PrimeCurve {{ p: {p}, a: {a}, b: {b} }} (explicit)")
            }
            CurveSpec::BinaryExplicit {
                n, modulus, a, b, ..
            } => format!("binary_explicit_instance({n}, 0x{modulus:x}, 0x{a:x}, 0x{b:x}, ...)"),
            CurveSpec::Fp4Translation { p, seed } => {
                format!("quartic_folds::generate_translation_instance({p}, {seed}, 16)")
            }
        }
    }
}

/// A built curve instance.
pub enum Instance {
    Prime(PrimeInstance),
    Binary(Box<BinaryInstance>),
    /// A Koblitz curve past one word (`64 ≤ n ≤ 127`).
    Wide(Box<WideInstance>),
    /// A curve over `F_{p⁴}` (`ecbench::fp4`); a point is one 64-bit word
    /// per coordinate, four 16-bit limbs.
    Fp4(Box<Fp4Instance>),
}

/// A wide Koblitz curve and the strong reference's arithmetic on it, which
/// every wide method and the workload's own `[k]G` use.
pub struct WideInstance {
    pub kc: WideKoblitz,
    pub rho: WideStrongRho,
}

/// The registry's modulus for `K_a` over `F_{2^n}`, if it lists one.
pub fn registry_koblitz_modulus(a: u8, n: u32) -> Option<IrreduciblePoly> {
    let doc: Value = serde_json::from_str(REGISTRY).ok()?;
    let entry = doc.get("curves")?.as_array()?.iter().find(|c| {
        let p = c.get("params");
        c.get("family").and_then(Value::as_str) == Some("koblitz")
            && p.and_then(|p| p.get("a")).and_then(Value::as_u64) == Some(u64::from(a))
            && p.and_then(|p| p.get("n")).and_then(Value::as_u64) == Some(u64::from(n))
    })?;
    let hex = entry.get("params")?.get("modulus")?.as_str()?;
    let bits = u128::from_str_radix(hex.trim_start_matches("0x"), 16).ok()?;
    if 127 - bits.leading_zeros() != n {
        return None;
    }
    let low_terms = (0..n).filter(|&i| (bits >> i) & 1 == 1).collect();
    Some(IrreduciblePoly {
        degree: n,
        low_terms,
    })
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
    /// Orders as `u128` for wide curves; a value that fits a `u64` is
    /// written as the number it always was (`canonical::compat_u128`).
    #[serde(with = "compat_u128")]
    pub group_order: u128,
    #[serde(with = "compat_u128")]
    pub r: u128,
    #[serde(with = "compat_u128")]
    pub cofactor: u128,
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

fn parse_hex(s: &str) -> Option<u128> {
    u128::from_str_radix(s.trim_start_matches("0x"), 16).ok()
}

/// `(ec1, curve_uid)` of the representation with this subgroup and
/// generator.
fn matching_representation(entry: &Value, r: u128, gx: u128, gy: u128) -> Option<(String, String)> {
    entry
        .get("representations")?
        .as_array()?
        .iter()
        .find_map(|rep| {
            let c = rep.get("curve")?;
            let order: u128 = c.get("subgroup_order")?.as_str()?.parse().ok()?;
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
    pub fn r(&self) -> u128 {
        match self {
            Instance::Prime(i) => u128::from(i.r),
            Instance::Binary(i) => u128::from(i.r),
            Instance::Wide(i) => i.kc.r,
            Instance::Fp4(i) => u128::from(i.r),
        }
    }

    pub fn facts(&self, spec: &CurveSpec) -> CurveFacts {
        if let Instance::Fp4(i) = self {
            let (icv1, slug) = i.names();
            let g = Fp4Instance::encode(&i.g).unwrap_or_default();
            return CurveFacts {
                slug,
                icv1,
                family: "fp4".into(),
                construction: spec.call(),
                field_degree: None,
                field_bits: Some(64 - i.p.leading_zeros()),
                group_order: u128::from(i.group_order),
                r: u128::from(i.r),
                cofactor: u128::from(i.cofactor),
                generator: g,
                automorphisms_available: 2,
                registered: false,
                ec1: None,
                curve_uid: None,
            };
        }
        let (id, family, degree, bits, order, r, h, gx, gy, aut) = match self {
            Instance::Prime(i) => (
                i.curve_id(),
                "prime",
                None,
                Some(64 - i.curve.p.leading_zeros()),
                u128::from(i.group_order),
                u128::from(i.r),
                u128::from(i.cofactor),
                u128::from(i.generator.0),
                u128::from(i.generator.1),
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
                u128::from(i.group_order),
                u128::from(i.r),
                u128::from(i.cofactor),
                u128::from(i.generator.x),
                u128::from(i.generator.y),
                if i.koblitz.is_some() { 2 * i.n } else { 2 },
            ),
            Instance::Fp4(_) => unreachable!("returned above"),
            Instance::Wide(i) => (
                i.kc.curve_id(),
                "koblitz",
                Some(i.kc.n),
                None,
                i.kc.group_order,
                i.kc.r,
                i.kc.cofactor,
                i.kc.generator.0,
                i.kc.generator.1,
                2 * i.kc.n,
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
    pub fn mul_generator_hex(&self, k: u128) -> Option<[String; 2]> {
        let mut scratch = GroupOps::default();
        if let Instance::Wide(i) = self {
            return match i.rho.scalar_mul(i.rho.generator(), k % i.kc.r) {
                RawPointG::Affine { x, y } => Some([format!("0x{x:x}"), format!("0x{y:x}")]),
                RawPointG::Infinity => None,
            };
        }
        // A one-word instance's subgroup order fits a word, so a reduced
        // scalar does too.
        let k = u64::try_from(k % self.r()).ok()?;
        if let Instance::Fp4(i) = self {
            return i.mul_generator_hex(k);
        }
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
            Instance::Wide(_) | Instance::Fp4(_) => unreachable!("returned above"),
        }
    }
}

impl Instance {
    /// A public target: for `counter = 0, 1, …`, an abscissa derived from
    /// `(model, r, seed, index, counter)` by SHA-256, lifted to a point
    /// (the ordinate's choice also hashed), times the cofactor; the first
    /// that is not the identity.  When `r ∤ h` (so `r² ∤ #E`, which
    /// `Workload::on` checks) the r-torsion is cyclic and that point lies
    /// in `⟨G⟩`.
    pub fn hash_to_subgroup(
        &self,
        curve: &CurveFacts,
        seed: u64,
        index: u64,
    ) -> Option<([String; 2], u64)> {
        let model = model_hash(curve);
        let mut scratch = GroupOps::default();
        if let Instance::Fp4(_) = self {
            // Planted targets only on E(F_{p⁴}): no public-target law yet.
            return None;
        }
        if let Instance::Wide(i) = self {
            // PUBLIC_TARGET_LAW_WIDE: the same inputs as decimal strings, a
            // 128-bit draw, its top bit choosing the ordinate.
            for counter in 0..10_000u64 {
                let parts = [
                    model.to_string(),
                    curve.r.to_string(),
                    seed.to_string(),
                    index.to_string(),
                    counter.to_string(),
                ];
                let h = derive_u128(PUBLIC_TARGET_LAW_WIDE, &parts);
                let flip = h >> 127 == 1;
                let pts = i.kc.points_with_x(h & ((1u128 << i.kc.n) - 1));
                let Some(&(x, y)) = pts.get(usize::from(flip && pts.len() > 1)) else {
                    continue;
                };
                if let RawPointG::Affine { x, y } =
                    i.rho.scalar_mul(RawPointG::Affine { x, y }, i.kc.cofactor)
                {
                    return Some(([format!("0x{x:x}"), format!("0x{y:x}")], counter));
                }
            }
            return None;
        }
        for counter in 0..10_000u64 {
            // One-word curves only: their order fits the law's u64 parts.
            let r = u64::try_from(curve.r).ok()?;
            let h = derive_u64(PUBLIC_TARGET_LAW, &[model, r, seed, index, counter]);
            let flip = (h >> 63) & 1 == 1;
            match self {
                Instance::Prime(i) => {
                    let c = &i.curve;
                    let x = h % c.p;
                    let Some(y) = c.sqrt(c.rhs(x)) else { continue };
                    let y = if flip && y != 0 { c.p - y } else { y };
                    let p = PrimePoint::affine(x, y);
                    if !c.is_on_curve(p) {
                        continue;
                    }
                    let q = c.mul(&mut scratch, p, i.cofactor);
                    if !q.infinity {
                        return Some(([format!("0x{:x}", q.x), format!("0x{:x}", q.y)], counter));
                    }
                }
                Instance::Wide(_) | Instance::Fp4(_) => unreachable!("handled above"),
                Instance::Binary(i) => {
                    let mask = if i.n >= 64 {
                        u64::MAX
                    } else {
                        (1u64 << i.n) - 1
                    };
                    let pts = i.points_with_x(h & mask);
                    let Some(&p) = pts.get(usize::from(flip && pts.len() > 1)) else {
                        continue;
                    };
                    let g = BinaryGroup(&i.fast);
                    let q = g.mul(&mut scratch, p, i.cofactor);
                    if !q.infinity {
                        return Some(([format!("0x{:x}", q.x), format!("0x{:x}", q.y)], counter));
                    }
                }
            }
        }
        None
    }
}

/// ICV1's last field: the first 12 hex digits of the model hash.
fn model_hash(curve: &CurveFacts) -> u64 {
    curve
        .icv1
        .rsplit(':')
        .next()
        .and_then(|h| u64::from_str_radix(h, 16).ok())
        .unwrap_or(0)
}

/// How planted targets are drawn: a scalar derived by SHA-256, `Q = [k]G`.
pub const TARGET_LAW: &str = "uniform_scalar_sha256_v1";

/// How public targets are drawn: an abscissa derived by SHA-256, lifted to
/// a point, the cofactor cleared.  Nobody ever knows the logarithm, which
/// is what the IC claim rules ask of a target ("previously unseen public
/// target"; `identity.py` refuses a planted workload).
pub const PUBLIC_TARGET_LAW: &str = "hash_to_subgroup_v1";

/// The public-target law for a subgroup order past `2^64`: the same
/// inputs, hashed as decimal strings, with a 128-bit draw.
pub const PUBLIC_TARGET_LAW_WIDE: &str = "hash_to_subgroup_v2";

/// Whether a workload's target was built from a known scalar.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TargetKind {
    /// `Q = [k]G` for a derived `k`, kept to check the answer against.
    #[default]
    Planted,
    /// The same `k` on every curve of one order: the derivation hashes
    /// `[r, target_seed, index]` only, so a paired comparison can
    /// divide the same planted scalar's cost on two isogenous curves
    /// (their `r` is equal) and resample targets jointly.  A paired
    /// quotient is a different statistic from the independent-target
    /// law, not a replacement for it.
    PlantedShared,
    /// `Q` hashed to the subgroup; the answer is checked by `[d]G = Q`.
    Public,
}

impl TargetKind {
    pub fn is_planted(&self) -> bool {
        matches!(self, TargetKind::Planted | TargetKind::PlantedShared)
    }
    /// Only the default kind stays unserialised.
    pub fn is_plain_planted(&self) -> bool {
        *self == TargetKind::Planted
    }
}

/// The planted-scalar law for a subgroup order past `2^64`: the same
/// inputs, hashed as decimal strings, with a 128-bit draw.  A one-word
/// order keeps [`TARGET_LAW`], so every existing workload is unchanged.
pub const TARGET_LAW_WIDE: &str = "uniform_scalar_sha256_v2";

/// The paired planted-scalar law: `k` derives from `[r, target_seed,
/// index]` alone — no curve model — so every curve of the same subgroup
/// order draws the same scalar list.  One-word orders only (the paired
/// comparisons this serves run on one-word instances).
pub const SHARED_TARGET_LAW: &str = "uniform_scalar_sha256_shared_v1";

/// `k ∈ [1, r−1]` for target `index` under `target_seed`, the same on
/// every curve of subgroup order `r`.
pub fn shared_planted_scalar(r: u64, target_seed: u64, index: u64) -> u128 {
    u128::from(1 + derive_u64(SHARED_TARGET_LAW, &[r, target_seed, index]) % (r - 1))
}

/// `k ∈ [1, r−1]` for target `index` under `target_seed`.
pub fn planted_scalar(curve: &CurveFacts, target_seed: u64, index: u64) -> u128 {
    // The curve's model hash and subgroup enter the derivation, so the
    // same seed on two curves gives unrelated scalars.
    let model = model_hash(curve);
    match u64::try_from(curve.r) {
        Ok(r) => u128::from(1 + derive_u64(TARGET_LAW, &[model, r, target_seed, index]) % (r - 1)),
        Err(_) => {
            let parts = [model, target_seed, index].map(|v| v.to_string());
            let parts = [
                parts[0].clone(),
                curve.r.to_string(),
                parts[1].clone(),
                parts[2].clone(),
            ];
            // The modulo bias of a 128-bit draw over r < 2^127 is below 2^-1.
            1 + derive_u128(TARGET_LAW_WIDE, &parts) % (curve.r - 1)
        }
    }
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
    /// `Q`, hex.
    pub target: [String; 2],
    /// The planted `k` with `Q = [k]G`, for a planted target.  Not in the
    /// identity (the target point is) and never passed to a solver.
    /// `None` for a public target, whose logarithm nobody knows.
    #[serde(default, with = "compat_u128::option")]
    pub planted: Option<u128>,
    /// The hash counter that landed a public target in the subgroup.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_counter: Option<u64>,
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
        kind: TargetKind,
    ) -> Result<(Workload, Instance), String> {
        let inst = spec.build()?;
        let w = Self::on(spec, &inst, target_seed, index, kind)?;
        Ok((w, inst))
    }

    /// Workload `index` under `target_seed` on `inst`, already built from
    /// `spec`.
    pub fn on(
        spec: &CurveSpec,
        inst: &Instance,
        target_seed: u64,
        index: u64,
        kind: TargetKind,
    ) -> Result<Workload, String> {
        let curve = inst.facts(spec);
        if curve.r < 5 {
            return Err(format!(
                "subgroup order {} is too small to measure",
                curve.r
            ));
        }
        let (law, target, planted, counter) =
            if let CurveSpec::PrimePublicFixture { qx, qy, .. } = *spec {
                if kind != TargetKind::Public || index != 0 {
                    return Err("external fixture is exactly one public target per curve".into());
                }
                (
                    "ecbench.public-fixture/v1",
                    [format!("0x{qx:x}"), format!("0x{qy:x}")],
                    None,
                    None,
                )
            } else {
                match kind {
                    TargetKind::Planted => {
                        let k = planted_scalar(&curve, target_seed, index);
                        let q = inst.mul_generator_hex(k).ok_or(
                            "the planted scalar gave the identity; the subgroup order is wrong",
                        )?;
                        let law = if u64::try_from(curve.r).is_ok() {
                            TARGET_LAW
                        } else {
                            TARGET_LAW_WIDE
                        };
                        (law, q, Some(k), None)
                    }
                    TargetKind::PlantedShared => {
                        let r = u64::try_from(curve.r)
                            .map_err(|_| "planted_shared targets need a one-word subgroup order")?;
                        let k = shared_planted_scalar(r, target_seed, index);
                        let q = inst.mul_generator_hex(k).ok_or(
                            "the planted scalar gave the identity; the subgroup order is wrong",
                        )?;
                        (SHARED_TARGET_LAW, q, Some(k), None)
                    }
                    TargetKind::Public => {
                        // With r | h the r-torsion can be Z/r × Z/r, and a cleared
                        // point need not lie in <G>.
                        if curve.cofactor.is_multiple_of(curve.r) {
                            return Err(format!(
                        "r = {} divides the cofactor {}: a hashed target need not lie in <G>",
                        curve.r, curve.cofactor
                    ));
                        }
                        let (q, c) = inst
                            .hash_to_subgroup(&curve, target_seed, index)
                            .ok_or("no public target found in 10 000 hash attempts")?;
                        let law = if u64::try_from(curve.r).is_ok() {
                            PUBLIC_TARGET_LAW
                        } else {
                            PUBLIC_TARGET_LAW_WIDE
                        };
                        (law, q, None, Some(c))
                    }
                }
            };
        let (workload_id, workload_sha256) = bare_id(
            "W",
            &Self::identity_view(&curve, law, target_seed, index, &target),
        )?;
        Ok(Workload {
            workload_id,
            workload_sha256,
            curve_spec: spec.clone(),
            curve,
            target_law: law.into(),
            target_seed,
            target_index: index,
            target,
            planted,
            target_counter: counter,
        })
    }

    /// The target's kind, from its law.
    pub fn kind(&self) -> TargetKind {
        if self.target_law == PUBLIC_TARGET_LAW
            || self.target_law == PUBLIC_TARGET_LAW_WIDE
            || self.target_law == "ecbench.public-fixture/v1"
        {
            TargetKind::Public
        } else if self.target_law == SHARED_TARGET_LAW {
            TargetKind::PlantedShared
        } else {
            TargetKind::Planted
        }
    }

    /// What the workload id is a hash of: the exact curve model, the
    /// subgroup, the generator and the target point, plus the law that
    /// drew it.  One target, cold.
    fn identity_view(
        curve: &CurveFacts,
        law: &str,
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
            "target_law": law,
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
                &self.target_law,
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

    fn public_fixture() -> CurveSpec {
        CurveSpec::PrimePublicFixture {
            p: 343073,
            a: 36086,
            b: 72172,
            group_order: 342833,
            r: 342833,
            gx: 0,
            gy: 2882,
            qx: 314327,
            qy: 50562,
        }
    }

    #[test]
    fn public_fixture_survives_measured_child_normalization() {
        let spec = public_fixture();
        let (parent, inst) = Workload::build(&spec, 7, 0, TargetKind::Public).unwrap();
        let child = spec.for_child(&inst);
        let rebuilt = Workload::build(&child, 7, 0, TargetKind::Public).unwrap().0;
        assert_eq!(parent.workload_id, rebuilt.workload_id);
        assert_eq!(parent.target, rebuilt.target);
        assert_eq!(rebuilt.kind(), TargetKind::Public);
        assert!(rebuilt.planted.is_none());
    }

    #[test]
    fn fixture_rejects_invalid_points_and_reuse_as_planted_or_multiple_targets() {
        let spec = public_fixture();
        assert!(Workload::build(&spec, 7, 0, TargetKind::Planted).is_err());
        assert!(Workload::build(&spec, 7, 1, TargetKind::Public).is_err());
        let mut value = serde_json::to_value(spec).unwrap();
        value["qx"] = json!(343073);
        assert!(serde_json::from_value::<CurveSpec>(value)
            .unwrap()
            .build()
            .is_err());
    }

    #[test]
    fn public_fixture_rejects_composite_or_impossible_order_claims() {
        for (key, value) in [("p", 343075u64), ("r", 342835), ("r", 17)] {
            let mut v = serde_json::to_value(public_fixture()).unwrap();
            v[key] = json!(value);
            if key == "r" {
                v["group_order"] = json!(value);
            }
            assert!(serde_json::from_value::<CurveSpec>(v)
                .unwrap()
                .build()
                .is_err());
        }
    }

    #[test]
    fn public_fixture_rejects_a_singular_additive_group_despite_hasse_witness() {
        let mut v = serde_json::to_value(public_fixture()).unwrap();
        v["a"] = json!(0);
        v["b"] = json!(0);
        v["r"] = v["p"].clone();
        v["group_order"] = v["p"].clone();
        for key in ["gx", "gy", "qx", "qy"] {
            v[key] = json!(1);
        }
        assert!(serde_json::from_value::<CurveSpec>(v)
            .unwrap()
            .build()
            .is_err());
    }

    #[test]
    fn registered_prime_curve_resolves_its_slug_and_ec1() {
        let spec = CurveSpec::PrimeSearch {
            bits: 20,
            seed: 59297,
        };
        let (w, _) = Workload::build(&spec, 1, 0, TargetKind::Planted).unwrap();
        assert_eq!(w.curve.slug, "icv1-fp20-t727-cd198a38");
        assert!(w.curve.registered);
        assert_eq!(w.curve.ec1.as_deref(), Some("EC1P20Cfphd19e740ea95d"));
        assert!(w.workload_id.starts_with('W') && w.workload_id.len() == 13);
        assert_eq!(w.recompute_id().unwrap(), w.workload_id);
    }

    #[test]
    fn targets_differ_by_index_and_are_stable() {
        let spec = CurveSpec::Koblitz { a: 0, n: 13 };
        let (a, _) = Workload::build(&spec, 7, 0, TargetKind::Planted).unwrap();
        let (b, _) = Workload::build(&spec, 7, 1, TargetKind::Planted).unwrap();
        let (a2, _) = Workload::build(&spec, 7, 0, TargetKind::Planted).unwrap();
        assert_ne!(a.workload_id, b.workload_id);
        assert_eq!(a, a2);
        assert_eq!(a.curve.family, "koblitz");
        assert_eq!(a.curve.automorphisms_available, 26);
    }

    #[test]
    fn public_targets_lie_in_the_subgroup_and_have_no_planted_answer() {
        for spec in [
            CurveSpec::Koblitz { a: 0, n: 13 },
            CurveSpec::PrimeSearch {
                bits: 16,
                seed: 59297,
            },
        ] {
            let (w, inst) = Workload::build(&spec, 7, 0, TargetKind::Public).unwrap();
            let (w1, _) = Workload::build(&spec, 7, 1, TargetKind::Public).unwrap();
            let (p, _) = Workload::build(&spec, 7, 0, TargetKind::Planted).unwrap();
            assert_eq!(w.target_law, PUBLIC_TARGET_LAW);
            assert_eq!(w.kind(), TargetKind::Public);
            assert!(w.planted.is_none() && w.target_counter.is_some());
            assert_ne!(w.workload_id, w1.workload_id);
            assert_ne!(w.workload_id, p.workload_id);
            assert_eq!(w.recompute_id().unwrap(), w.workload_id);
            assert_eq!(
                Workload::build(&spec, 7, 0, TargetKind::Public).unwrap().0,
                w
            );
            // In <G>, by exhaustion on a subgroup this small.
            assert!((1..w.curve.r).any(|k| inst.mul_generator_hex(k).as_ref() == Some(&w.target)));
        }
    }

    #[test]
    fn explicit_form_rebuilds_the_same_workload() {
        let spec = CurveSpec::PrimeSearch {
            bits: 16,
            seed: 59297,
        };
        let (w, inst) = Workload::build(&spec, 3, 2, TargetKind::Planted).unwrap();
        let ex = CurveSpec::explicit(&inst).unwrap();
        let (w2, _) = Workload::build(&ex, 3, 2, TargetKind::Planted).unwrap();
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
        let (w, _) = Workload::build(
            &CurveSpec::Koblitz { a: 0, n: 13 },
            0,
            0,
            TargetKind::Planted,
        )
        .unwrap();
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
    fn conductor_gap_pair_shares_order_and_differs_in_j() {
        let crater = CurveSpec::ConductorGap {
            ell: 263,
            t: 241,
            level: GapLevel::Crater,
            seed: 0,
        };
        let floor = CurveSpec::ConductorGap {
            ell: 263,
            t: 241,
            level: GapLevel::Floor,
            seed: 3,
        };
        let (Instance::Prime(c), Instance::Prime(f)) =
            (crater.build().unwrap(), floor.build().unwrap())
        else {
            panic!("conductor_gap builds prime instances");
        };
        assert_eq!(
            (c.curve.p, c.group_order, c.cofactor),
            (343_073, 342_833, 1)
        );
        assert_eq!(
            (f.curve.p, f.group_order, f.cofactor),
            (343_073, 342_833, 1)
        );
        assert_ne!(
            conductor_gap::j_invariant(&c.curve),
            conductor_gap::j_invariant(&f.curve)
        );
        assert_eq!(
            conductor_gap::j_invariant(&c.curve),
            343_073 - 884_736 % 343_073
        );
        // The explicit form a measured child receives rebuilds the same curve.
        let ex = CurveSpec::explicit(&Instance::Prime(f.clone())).unwrap();
        let Instance::Prime(f2) = ex.build().unwrap() else {
            panic!()
        };
        assert_eq!(
            (f2.curve.a, f2.curve.b, f2.generator),
            (f.curve.a, f.curve.b, f.generator)
        );
        assert_eq!(crater.call(), "conductor_gap::crater(params(263, 241))");
        assert!(CurveSpec::ConductorGap {
            ell: 263,
            t: 3,
            level: GapLevel::Crater,
            seed: 0
        }
        .build()
        .is_err());
    }

    #[test]
    fn unusable_koblitz_curve_is_refused() {
        assert!(CurveSpec::Koblitz { a: 0, n: 17 }.build().is_err());
    }
}
