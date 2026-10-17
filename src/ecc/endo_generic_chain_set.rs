//! Curve-generic loader for the frozen chain sets of the native GLV-chain
//! benchmark (`research/endo_native_bench_20261016`): the
//! `endosweep-chainsweep/1` JSON the autoresearcher's chain sweep writes
//! (`harness/endosweep/chainsweep.py`, generalised per curve), for any short
//! Weierstrass curve over a prime field, the curve's own `p`, `a`, `b`, `n`,
//! cofactor, generator, `λ`, reduced basis and test vectors included.
//!
//! The file is frozen input with Python provenance; nothing here runs Python.
//! [`ChainSet::parse`] checks the integers ([`ChainSet::check_integers`]: `p`
//! and `n` probable primes, the cofactor consistent with Hasse's bound,
//! `ω`'s root a root of its minimal polynomial mod `n`, every norm
//! `N(a + bω)` and the product of its steps, every eigenvalue equal to
//! `a' + b'·ω_root` for the matched element `±(a + bω)`), and
//! [`ChainSet::prepare`] builds the curve and every chain in Montgomery form
//! at a limb count `L`, which checks every field element is canonical, the
//! steps compose from `E` back to `E`, the polynomials are monic of the right
//! degrees, and each GLV basis has determinant `±n` and lies in `λ`'s
//! lattice ([`PreparedChain::from_steps`], [`GlvContext::from_parts`]).
//! The generator is checked on the curve and of order `n` by the crate's
//! textbook `BigUint` arithmetic.  Test vectors are checked by the tests and
//! by the benchmark before any timing.

use super::endo_generic_point::{
    textbook_mul, Affine, Curve, GlvContext, PreparedChain, StepInput,
};
use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use serde_json::Value;
use std::path::Path;

/// The format tag the chain sweep writes.
pub const FORMAT: &str = "endosweep-chainsweep/1";

/// A test vector: `φ(P)`, and `k = k1 + k2·λ` with `k·P`.
#[derive(Clone, Debug)]
pub struct VectorSpec {
    pub p: (BigUint, BigUint),
    pub phi_p: (BigUint, BigUint),
    pub k: BigUint,
    pub k1: BigInt,
    pub k2: BigInt,
    pub k_p: (BigUint, BigUint),
}

/// The chain sweep's own count for one evaluator.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ModelOps {
    pub mul: u64,
    pub sqr: u64,
}

/// One (element, step order) chain.
#[derive(Clone, Debug)]
pub struct ChainSpec {
    pub id: String,
    pub element: (i64, i64),
    pub matched_element: (i64, i64),
    pub norm: u64,
    pub order: Vec<u32>,
    pub eigenvalue: BigUint,
    pub glv_basis: [[BigInt; 2]; 2],
    pub babai_bound_bits: u32,
    pub model_generic: ModelOps,
    pub model_optimised: ModelOps,
    pub steps: Vec<StepInput>,
    pub iso_u: BigUint,
    pub vectors: Vec<VectorSpec>,
}

/// The whole frozen set.
#[derive(Clone, Debug)]
pub struct ChainSet {
    pub target: String,
    pub p: BigUint,
    pub a: BigUint,
    pub b: BigUint,
    pub n: BigUint,
    pub cofactor: BigUint,
    pub d_k: i64,
    pub class_number: Option<u64>,
    pub omega_root: BigUint,
    /// The generator the file names (`curve.G`), if any.
    pub generator: Option<(BigUint, BigUint)>,
    pub generator_source: String,
    pub chains: Vec<ChainSpec>,
}

/// A chain set built at limb count `L`: the curve, its generator, and one GLV
/// context per chain (same order as [`ChainSet::chains`]).
pub struct PreparedSet<const L: usize> {
    pub curve: Curve<L>,
    pub generator: Affine<L>,
    pub contexts: Vec<GlvContext<L>>,
}

fn field<'a>(v: &'a Value, key: &str, ctx: &str) -> Result<&'a Value, String> {
    v.get(key)
        .ok_or_else(|| format!("{ctx}: missing field {key:?}"))
}

fn as_str<'a>(v: &'a Value, ctx: &str) -> Result<&'a str, String> {
    v.as_str()
        .ok_or_else(|| format!("{ctx}: expected a string"))
}

fn as_u64(v: &Value, ctx: &str) -> Result<u64, String> {
    v.as_u64()
        .ok_or_else(|| format!("{ctx}: expected a non-negative integer"))
}

fn as_i64(v: &Value, ctx: &str) -> Result<i64, String> {
    v.as_i64()
        .ok_or_else(|| format!("{ctx}: expected an integer"))
}

fn as_array<'a>(v: &'a Value, ctx: &str) -> Result<&'a Vec<Value>, String> {
    v.as_array()
        .ok_or_else(|| format!("{ctx}: expected an array"))
}

fn hex(s: &str, ctx: &str) -> Result<BigUint, String> {
    let d = s
        .strip_prefix("0x")
        .ok_or_else(|| format!("{ctx}: {s:?} is not 0x-prefixed hex"))?;
    BigUint::parse_bytes(d.as_bytes(), 16).ok_or_else(|| format!("{ctx}: {s:?} is not hex"))
}

fn hex_below(v: &Value, bound: &BigUint, ctx: &str) -> Result<BigUint, String> {
    let x = hex(as_str(v, ctx)?, ctx)?;
    if &x >= bound {
        return Err(format!("{ctx}: value is not below the modulus"));
    }
    Ok(x)
}

/// A non-negative integer given as 0x-hex, decimal string or JSON number.
fn uint_any(v: &Value, ctx: &str) -> Result<BigUint, String> {
    if let Some(u) = v.as_u64() {
        return Ok(BigUint::from(u));
    }
    if let Value::Number(num) = v {
        // the crate builds serde_json with `arbitrary_precision`, so a large
        // JSON integer keeps every digit in its text
        let t = num.to_string();
        return BigUint::parse_bytes(t.as_bytes(), 10)
            .ok_or_else(|| format!("{ctx}: {t} is not a non-negative integer"));
    }
    let s = as_str(v, ctx)?;
    if s.starts_with("0x") {
        hex(s, ctx)
    } else {
        BigUint::parse_bytes(s.as_bytes(), 10)
            .ok_or_else(|| format!("{ctx}: {s:?} is not an integer"))
    }
}

fn decimal(v: &Value, ctx: &str) -> Result<BigInt, String> {
    let s = as_str(v, ctx)?;
    BigInt::parse_bytes(s.as_bytes(), 10)
        .ok_or_else(|| format!("{ctx}: {s:?} is not a signed decimal"))
}

fn pair(v: &Value, ctx: &str) -> Result<(i64, i64), String> {
    let a = as_array(v, ctx)?;
    if a.len() != 2 {
        return Err(format!("{ctx}: expected two integers"));
    }
    Ok((as_i64(&a[0], ctx)?, as_i64(&a[1], ctx)?))
}

fn point(v: &Value, p: &BigUint, ctx: &str) -> Result<(BigUint, BigUint), String> {
    let a = as_array(v, ctx)?;
    if a.len() != 2 {
        return Err(format!("{ctx}: expected [x, y]"));
    }
    Ok((hex_below(&a[0], p, ctx)?, hex_below(&a[1], p, ctx)?))
}

fn poly(v: &Value, p: &BigUint, ctx: &str) -> Result<Vec<BigUint>, String> {
    as_array(v, ctx)?
        .iter()
        .map(|c| hex_below(c, p, ctx))
        .collect()
}

fn model(v: &Value, ctx: &str) -> Result<ModelOps, String> {
    Ok(ModelOps {
        mul: as_u64(field(v, "M", ctx)?, ctx)?,
        sqr: as_u64(field(v, "S", ctx)?, ctx)?,
    })
}

/// Miller–Rabin with the first 40 prime bases: a probable-prime test (no
/// certificate), enough to refuse a mistyped modulus.
pub fn is_probable_prime(n: &BigUint) -> bool {
    let small = [
        2u32, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83,
        89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    ];
    if n < &BigUint::from(2u8) {
        return false;
    }
    for &q in &small {
        let q = BigUint::from(q);
        if n == &q {
            return true;
        }
        if (n % &q).is_zero() {
            return false;
        }
    }
    let nm1 = n - 1u8;
    let s = nm1.trailing_zeros().unwrap_or(0);
    let d = &nm1 >> s;
    'base: for &q in &small {
        let mut x = BigUint::from(q).modpow(&d, n);
        if x.is_one() || x == nm1 {
            continue;
        }
        for _ in 1..s {
            x = (&x * &x) % n;
            if x == nm1 {
                continue 'base;
            }
        }
        return false;
    }
    true
}

/// The trace `t` and norm `N` of `ω` for the fundamental discriminant `D`:
/// `ω = (1 + √D)/2` (`t = 1`, `N = (1 − D)/4`) when `D ≡ 1 (mod 4)`, and
/// `ω = √(D/4)` (`t = 0`, `N = −D/4`) otherwise.
pub fn omega_trace_norm(d_k: i64) -> Result<(i64, i64), String> {
    if d_k >= 0 {
        return Err("D_K must be negative".into());
    }
    match d_k.rem_euclid(4) {
        1 => Ok((1, (1 - d_k) / 4)),
        0 => Ok((0, -d_k / 4)),
        _ => Err(format!("{d_k} is not a discriminant")),
    }
}

impl ChainSet {
    /// Parse and validate the JSON text of a chain set (the integer checks;
    /// see [`Self::prepare`] for the field-level ones).
    pub fn parse(text: &str) -> Result<Self, String> {
        let root: Value = serde_json::from_str(text).map_err(|e| format!("not JSON: {e}"))?;
        let format = as_str(field(&root, "format", "chain set")?, "format")?;
        if format != FORMAT {
            return Err(format!("format {format:?}, expected {FORMAT:?}"));
        }
        let gethex = |k: &str| -> Result<BigUint, String> {
            hex(as_str(field(&root, k, "chain set")?, k)?, k)
        };
        let p = gethex("p")?;
        let a = gethex("a")?;
        let b = gethex("b")?;
        let n = gethex("n")?;
        if a >= p || b >= p {
            return Err("a and b must be canonical (below p)".into());
        }
        let cofactor = uint_any(field(&root, "cofactor", "chain set")?, "cofactor")?;
        let d_k = as_i64(field(&root, "D_K", "chain set")?, "D_K")?;
        let omega_root = gethex("omega_root")?;
        if omega_root >= n {
            return Err("omega_root must be below n".into());
        }
        let (generator, generator_source) = match root.get("curve").and_then(|c| c.get("G")) {
            Some(g) => (
                Some(point(g, &p, "curve.G")?),
                root["curve"]
                    .get("G_source")
                    .and_then(Value::as_str)
                    .unwrap_or("curve.G")
                    .to_string(),
            ),
            None => (None, String::new()),
        };
        let mut chains = Vec::new();
        for (idx, c) in as_array(field(&root, "chains", "chain set")?, "chains")?
            .iter()
            .enumerate()
        {
            let id = as_str(field(c, "id", "chain")?, "id")?.to_string();
            let ctx = format!("chain {idx} ({id})");
            let mut steps = Vec::new();
            for (sidx, st) in as_array(field(c, "steps", &ctx)?, &ctx)?.iter().enumerate() {
                let sctx = format!("{ctx} step {sidx}");
                let ell = u32::try_from(as_u64(field(st, "ell", &sctx)?, &sctx)?)
                    .map_err(|_| format!("{sctx}: ell too large"))?;
                let dom = field(st, "domain", &sctx)?;
                let cod = field(st, "codomain", &sctx)?;
                let velu2 = match st.get("velu2") {
                    Some(v) => Some((
                        hex_below(field(v, "x0", &sctx)?, &p, &sctx)?,
                        hex_below(field(v, "v", &sctx)?, &p, &sctx)?,
                    )),
                    None => None,
                };
                steps.push(StepInput {
                    ell,
                    domain_a: hex_below(field(dom, "a", &sctx)?, &p, &sctx)?,
                    domain_b: hex_below(field(dom, "b", &sctx)?, &p, &sctx)?,
                    codomain_a: hex_below(field(cod, "a", &sctx)?, &p, &sctx)?,
                    codomain_b: hex_below(field(cod, "b", &sctx)?, &p, &sctx)?,
                    psi: poly(field(st, "psi", &sctx)?, &p, &sctx)?,
                    n: poly(field(st, "N", &sctx)?, &p, &sctx)?,
                    m: poly(field(st, "M", &sctx)?, &p, &sctx)?,
                    velu2,
                });
            }
            let order: Vec<u32> = as_array(field(c, "order", &ctx)?, &ctx)?
                .iter()
                .map(|v| as_u64(v, &ctx).map(|x| x as u32))
                .collect::<Result<_, _>>()?;
            if order != steps.iter().map(|s| s.ell).collect::<Vec<_>>() {
                return Err(format!("{ctx}: order does not match the steps"));
            }
            let basis = as_array(field(c, "glv_basis", &ctx)?, &ctx)?;
            if basis.len() != 2 {
                return Err(format!("{ctx}: glv_basis must have two rows"));
            }
            let mut glv_basis: [[BigInt; 2]; 2] = Default::default();
            for (r, row) in basis.iter().enumerate() {
                let row = as_array(row, &ctx)?;
                if row.len() != 2 {
                    return Err(format!("{ctx}: glv_basis rows have two entries"));
                }
                for (j, e) in row.iter().enumerate() {
                    glv_basis[r][j] = decimal(e, &ctx)?;
                }
            }
            let mops = field(c, "model_ops", &ctx)?;
            let mut vectors = Vec::new();
            for tv in as_array(field(c, "test_vectors", &ctx)?, &ctx)? {
                vectors.push(VectorSpec {
                    p: point(field(tv, "P", &ctx)?, &p, &ctx)?,
                    phi_p: point(field(tv, "phiP", &ctx)?, &p, &ctx)?,
                    k: hex(as_str(field(tv, "k", &ctx)?, &ctx)?, &ctx)?,
                    k1: decimal(field(tv, "k1", &ctx)?, &ctx)?,
                    k2: decimal(field(tv, "k2", &ctx)?, &ctx)?,
                    k_p: point(field(tv, "kP", &ctx)?, &p, &ctx)?,
                });
            }
            chains.push(ChainSpec {
                id,
                element: pair(field(c, "element", &ctx)?, &ctx)?,
                matched_element: pair(field(c, "matched_element", &ctx)?, &ctx)?,
                norm: as_u64(field(c, "norm", &ctx)?, &ctx)?,
                order,
                eigenvalue: hex_below(field(c, "eigenvalue", &ctx)?, &n, &ctx)?,
                glv_basis,
                babai_bound_bits: as_u64(field(c, "babai_bound_bits", &ctx)?, &ctx)? as u32,
                model_generic: model(field(mops, "generic", &ctx)?, &ctx)?,
                model_optimised: model(field(mops, "optimised", &ctx)?, &ctx)?,
                steps,
                iso_u: hex_below(field(c, "isomorphism_u", &ctx)?, &p, &ctx)?,
                vectors,
            });
        }
        if chains.is_empty() {
            return Err("the chain set is empty".into());
        }
        let mut ids: Vec<&str> = chains.iter().map(|c| c.id.as_str()).collect();
        ids.sort_unstable();
        if ids.windows(2).any(|w| w[0] == w[1]) {
            return Err("duplicate chain id".into());
        }
        let set = ChainSet {
            target: as_str(field(&root, "target", "chain set")?, "target")?.to_string(),
            p,
            a,
            b,
            n,
            cofactor,
            d_k,
            class_number: root.get("class_number").and_then(Value::as_u64),
            omega_root,
            generator,
            generator_source,
            chains,
        };
        set.check_integers()?;
        Ok(set)
    }

    /// Read and parse a chain-set file.
    pub fn load(path: &Path) -> Result<Self, String> {
        let text =
            std::fs::read_to_string(path).map_err(|e| format!("read {}: {e}", path.display()))?;
        Self::parse(&text)
    }

    /// The integer facts the set rests on (see the module docs).
    pub fn check_integers(&self) -> Result<(), String> {
        if !is_probable_prime(&self.p) {
            return Err("p is not a probable prime".into());
        }
        if !is_probable_prime(&self.n) || self.n == BigUint::from(2u8) {
            return Err("n is not an odd probable prime".into());
        }
        if self.cofactor.is_zero() {
            return Err("cofactor must be positive".into());
        }
        // Hasse: |p + 1 − h·n| ≤ 2√p, i.e. (p + 1 − h·n)² ≤ 4p
        let hn = BigInt::from(&self.cofactor * &self.n);
        let t = BigInt::from(&self.p + 1u8) - hn;
        if &t * &t > BigInt::from(&self.p * 4u8) {
            return Err("h·n is outside the Hasse interval of p".into());
        }
        let (tw, nw) = omega_trace_norm(self.d_k)?;
        let nb = BigInt::from(self.n.clone());
        let w = BigInt::from(self.omega_root.clone());
        if !(&w * &w - &w * tw + nw).mod_floor(&nb).is_zero() {
            return Err("omega_root is not a root of omega's minimal polynomial mod n".into());
        }
        for c in &self.chains {
            let (a, b) = c.element;
            let norm = a * a + tw * a * b + nw * b * b;
            if norm < 0 || norm as u64 != c.norm {
                return Err(format!("{}: norm is not N(a + b·omega)", c.id));
            }
            let product: u64 = c.order.iter().map(|&l| u64::from(l)).product();
            if product != c.norm {
                return Err(format!(
                    "{}: the steps' degrees do not multiply to the norm",
                    c.id
                ));
            }
            let (ma, mb) = c.matched_element;
            if (ma, mb) != (a, b) && (ma, mb) != (-a, -b) {
                return Err(format!("{}: matched element is not ±element", c.id));
            }
            let lam = (BigInt::from(ma) + BigInt::from(mb) * &w).mod_floor(&nb);
            if lam != BigInt::from(c.eigenvalue.clone()) {
                return Err(format!(
                    "{}: eigenvalue is not a' + b'·omega_root mod n",
                    c.id
                ));
            }
            if c.vectors.is_empty() {
                return Err(format!("{}: no test vectors", c.id));
            }
            for v in &c.vectors {
                let diff =
                    &v.k1 + &v.k2 * BigInt::from(c.eigenvalue.clone()) - BigInt::from(v.k.clone());
                if !diff.mod_floor(&nb).is_zero() {
                    return Err(format!("{}: a test vector's k1 + k2·λ is not k", c.id));
                }
            }
        }
        Ok(())
    }

    /// The least number of 64-bit limbs that holds `p`.
    pub fn limbs(&self) -> usize {
        (self.p.bits() as usize).div_ceil(64)
    }

    pub fn get(&self, id: &str) -> Option<&ChainSpec> {
        self.chains.iter().find(|c| c.id == id)
    }

    /// Build the curve, the generator and every chain at `L` limbs.  The
    /// generator is `curve.G` from the file, or else `fallback` (e.g. the
    /// crate's own generator of a named curve), or else the first test
    /// vector's `P`; it must lie on the curve and have order `n` under the
    /// textbook arithmetic.
    pub fn prepare<const L: usize>(
        &self,
        fallback: Option<(BigUint, BigUint)>,
    ) -> Result<PreparedSet<L>, String> {
        let curve = Curve::<L>::new(&self.p, &self.a, &self.b)?;
        let (gx, gy) = self
            .generator
            .clone()
            .or(fallback)
            .unwrap_or_else(|| self.chains[0].vectors[0].p.clone());
        let generator = curve.affine(&gx, &gy)?;
        if !curve.is_on_curve(&generator) {
            return Err("the generator is not on the curve".into());
        }
        if textbook_mul(&curve, &generator, &self.n).is_some() {
            return Err("the generator does not have order n".into());
        }
        let mut contexts = Vec::with_capacity(self.chains.len());
        for c in &self.chains {
            let chain = PreparedChain::from_steps(&curve, &c.id, c.norm, &c.steps, &c.iso_u)
                .map_err(|e| format!("{}: {e}", c.id))?;
            if chain.step_degrees() != c.order {
                return Err(format!("{}: prepared degrees differ", c.id));
            }
            let ctx = GlvContext::from_parts(
                chain,
                &c.eigenvalue,
                &self.n,
                c.glv_basis.clone(),
                c.babai_bound_bits,
            )
            .map_err(|e| format!("{}: {e}", c.id))?;
            contexts.push(ctx);
        }
        Ok(PreparedSet {
            curve,
            generator,
            contexts,
        })
    }
}

impl ChainSpec {
    /// `7·5·5`.
    pub fn shape(&self) -> String {
        self.order
            .iter()
            .map(u32::to_string)
            .collect::<Vec<_>>()
            .join("·")
    }

    /// `4 + 1ω`.
    pub fn element_label(&self) -> String {
        let (a, b) = self.element;
        if b >= 0 {
            format!("{a} + {b}ω")
        } else {
            format!("{a} - {}ω", -b)
        }
    }

    /// The largest `|k1|`, `|k2|` of the test vectors, in bits.
    pub fn vector_bits(&self) -> u64 {
        self.vectors
            .iter()
            .map(|v| v.k1.abs().bits().max(v.k2.abs().bits()))
            .max()
            .unwrap_or(0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ecc::endo_generic_point::{
        jacobian_to_biguint, scalar_mul_glv, scalar_mul_wnaf, ChainEval, TableMode,
    };
    use crate::ecc::CurveParams;
    use num_bigint::RandBigInt;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    const CPB: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/research/cryptopro_b_chain_sweep_20261006/frozen_inputs/chains.constants.json"
    );
    const NATIVE_DIR: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/research/endo_native_bench_20261016/frozen_inputs"
    );

    fn cpb() -> ChainSet {
        ChainSet::load(Path::new(CPB)).expect("the CryptoPro-B chain set parses")
    }

    fn cpb_generator() -> Option<(BigUint, BigUint)> {
        let c = CurveParams::gost_cryptopro_b();
        Some((c.gx, c.gy))
    }

    /// Every vector of every chain, both evaluators, both table modes, and
    /// `φ = [λ]` on random points of the group; generic over `L`.
    fn check_set<const L: usize>(
        set: &ChainSet,
        fallback: Option<(BigUint, BigUint)>,
        extra: usize,
    ) {
        let ps = set.prepare::<L>(fallback).expect("prepares");
        let c = &ps.curve;
        let mut rng = StdRng::seed_from_u64(0x5ee9);
        for (spec, ctx) in set.chains.iter().zip(ps.contexts.iter()) {
            for v in &spec.vectors {
                let p = c.affine(&v.p.0, &v.p.1).unwrap();
                assert!(c.is_on_curve(&p), "{}", spec.id);
                for eval in [ChainEval::Generic, ChainEval::Optimised] {
                    let phi = ctx.chain.apply_with(c, &p, eval);
                    assert_eq!(
                        jacobian_to_biguint(c, &phi),
                        Some(v.phi_p.clone()),
                        "{} {eval:?}",
                        spec.id
                    );
                }
                assert_eq!(
                    ctx.decompose(&v.k),
                    (v.k1.clone(), v.k2.clone()),
                    "{}",
                    spec.id
                );
                for t in [TableMode::Affine, TableMode::Jacobian] {
                    let got = scalar_mul_glv(c, ctx, &p, &v.k, 4, t, ChainEval::Optimised);
                    assert_eq!(
                        jacobian_to_biguint(c, &got),
                        Some(v.k_p.clone()),
                        "{}",
                        spec.id
                    );
                }
            }
            for _ in 0..extra {
                let r = rng.gen_biguint_below(&set.n);
                let p = c
                    .to_affine(&scalar_mul_wnaf(c, &ps.generator, &r, 5, TableMode::Affine))
                    .unwrap();
                let want = jacobian_to_biguint(
                    c,
                    &scalar_mul_wnaf(c, &p, &ctx.lambda, 5, TableMode::Jacobian),
                );
                for eval in [ChainEval::Generic, ChainEval::Optimised] {
                    assert_eq!(
                        jacobian_to_biguint(c, &ctx.chain.apply_with(c, &p, eval)),
                        want,
                        "{}",
                        spec.id
                    );
                }
            }
        }
    }

    #[test]
    fn cryptopro_b_set_loads_and_reproduces_its_vectors_on_the_generic_path() {
        let set = cpb();
        assert_eq!(set.limbs(), 4);
        assert_eq!(set.chains.len(), 22);
        check_set::<4>(&set, cpb_generator(), 2);
    }

    #[test]
    fn native_sets_present_load_and_reproduce_their_vectors() {
        // The frozen chain sets of research/endo_native_bench_20261016, when
        // committed; each at its own limb count.
        let Ok(dir) = std::fs::read_dir(NATIVE_DIR) else {
            return;
        };
        for e in dir.flatten() {
            let path = e.path();
            if !path.to_string_lossy().ends_with(".constants.json") {
                continue;
            }
            let set = ChainSet::load(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
            match set.limbs() {
                4 => check_set::<4>(&set, None, 1),
                6 => check_set::<6>(&set, None, 1),
                13 => check_set::<13>(&set, None, 1),
                l => panic!("{}: {l} limbs not built", path.display()),
            }
        }
    }

    #[test]
    fn tampered_sets_are_refused() {
        let text = std::fs::read_to_string(CPB).unwrap();
        let v: Value = serde_json::from_str(&text).unwrap();
        let refuse = |v: &Value, what: &str| {
            let r = ChainSet::parse(&v.to_string());
            assert!(r.is_err(), "{what} was accepted");
        };
        let mut w = v.clone();
        w["format"] = Value::String("endosweep-chainsweep/0".into());
        refuse(&w, "format");
        let mut w = v.clone();
        w["p"] = Value::String(
            "0x8000000000000000000000000000000000000000000000000000000000000c9b".into(),
        );
        refuse(&w, "a composite p");
        let mut w = v.clone();
        w["cofactor"] = Value::from(2u8);
        refuse(&w, "a cofactor outside Hasse");
        let mut w = v.clone();
        w["chains"][0]["norm"] = Value::from(156u32);
        refuse(&w, "a wrong norm");
        let mut w = v.clone();
        w["chains"][0]["eigenvalue"] = Value::String("0x5".into());
        refuse(&w, "a wrong eigenvalue");
        let mut w = v.clone();
        let first = w["chains"][0].clone();
        w["chains"].as_array_mut().unwrap().push(first);
        refuse(&w, "a duplicate id");
        let mut w = v.clone();
        w["chains"][0]["steps"][0]["N"][0] = v["p"].clone();
        refuse(&w, "a non-canonical coefficient");
        // field-level refusals happen at prepare()
        let mut w = v.clone();
        let a0 = w["chains"][0]["steps"][0]["domain"]["a"].clone();
        if w["chains"][0]["steps"].as_array().unwrap().len() > 1 {
            w["chains"][0]["steps"][1]["domain"]["a"] = a0;
            let set = ChainSet::parse(&w.to_string()).unwrap();
            assert!(set.prepare::<4>(cpb_generator()).is_err());
        }
        let mut w = v.clone();
        w["chains"][0]["isomorphism_u"] = Value::String("0x2".into());
        let set = ChainSet::parse(&w.to_string()).unwrap();
        let err = set.prepare::<4>(cpb_generator()).err().unwrap();
        assert!(err.contains("iso_u"), "{err}");
    }

    #[test]
    fn omega_trace_and_norm() {
        assert_eq!(omega_trace_norm(-619).unwrap(), (1, 155));
        assert_eq!(omega_trace_norm(-8).unwrap(), (0, 2));
        assert_eq!(omega_trace_norm(-339).unwrap(), (1, 85));
        assert!(omega_trace_norm(-2).is_err());
        assert!(is_probable_prime(&BigUint::from(1_000_000_007u64)));
        assert!(!is_probable_prime(&BigUint::from(1_000_000_007u64 * 3)));
    }

    /// The chain sweep's model and the generic harness's counted operations
    /// agree for every CryptoPro-B chain and both evaluators — the same
    /// assertion the specialised harness's counting build makes.
    #[cfg(feature = "endo-chain-opcount")]
    #[test]
    fn cryptopro_b_phi_counts_match_the_model_on_the_generic_path() {
        use crate::ecc::endo_generic_field::opcount;
        let set = cpb();
        let ps = set.prepare::<4>(cpb_generator()).unwrap();
        for (spec, ctx) in set.chains.iter().zip(ps.contexts.iter()) {
            for (eval, want) in [
                (ChainEval::Generic, spec.model_generic),
                (ChainEval::Optimised, spec.model_optimised),
            ] {
                opcount::reset();
                let _ = ctx.chain.apply_with(&ps.curve, &ps.generator, eval);
                let got = opcount::read();
                assert_eq!(
                    (got.mul, got.sqr, got.inv),
                    (want.mul, want.sqr, 0),
                    "{}",
                    spec.id
                );
            }
        }
    }
}
