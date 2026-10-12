//! The frozen chain set of the CryptoPro-B chain sweep
//! (`research/cryptopro_b_chain_sweep_20261006`): every isogeny-chain
//! realisation, within a stated cost factor of the cheapest, of the cheap
//! endomorphisms of GOST CryptoPro-B, loaded at run time from the JSON
//! file the autoresearcher's chain sweep writes (format
//! [`FORMAT`]).
//!
//! The file is frozen input with Python provenance (the autoresearcher's
//! `harness/endosweep/chainsweep.py`); nothing here runs Python.  Loading
//! validates it against this crate: the curve constants must equal those
//! of [`crate::ecc::cryptopro_b_chain_consts`], every field element must be
//! a canonical value below `p`, and each chain goes through the same checks
//! as the compiled-in chains ([`PreparedChain::from_steps`],
//! [`GlvContext::from_parts`]): degrees, monic polynomials, steps that
//! compose from `E` back to `E`, a non-singular basis inside `λ`'s lattice.
//! The test vectors are checked separately by the tests and the benchmark.

use super::cryptopro_b_chain_consts::{A, B, N, P};
use super::cryptopro_b_point::{CryptoProBAffine, GlvContext, PreparedChain, StepConstants};
use crate::ct_bignum::{Uint, U256};
use num_bigint::{BigInt, BigUint};
use serde_json::Value;
use std::path::Path;

/// The format tag the chain sweep writes.
pub const FORMAT: &str = "endosweep-chainsweep/1";

/// One step: an `ell`-isogeny in Kohel form, `x' = N(x)/ψ(x)²`,
/// `y' = y·M(x)/ψ(x)³`, between the declared curves.
#[derive(Clone, Debug)]
pub struct StepSpec {
    pub ell: u32,
    pub domain_a: [u64; 4],
    pub domain_b: [u64; 4],
    pub codomain_a: [u64; 4],
    pub codomain_b: [u64; 4],
    pub psi: Vec<[u64; 4]>,
    pub n: Vec<[u64; 4]>,
    pub m: Vec<[u64; 4]>,
}

/// A test vector: `φ(P)`, and `k = k1 + k2·λ` with `k·P`.
#[derive(Clone, Debug)]
pub struct VectorSpec {
    pub p: [[u64; 4]; 2],
    pub phi_p: [[u64; 4]; 2],
    pub k: BigUint,
    pub k1: BigInt,
    pub k2: BigInt,
    pub k_p: [[u64; 4]; 2],
}

/// The chain sweep's own operation count for one evaluator (field
/// multiplications and squarings), carried so that the benchmark can set
/// its counted values beside the model's.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ModelOps {
    pub mul: u64,
    pub sqr: u64,
}

/// One (element, step order) chain.
#[derive(Clone, Debug)]
pub struct ChainSpec {
    /// `"<a>+<b>w/<l1>.<l2>…"`, e.g. `4+1w/7.5.5`.
    pub id: String,
    /// The element `a + b·ω` the chain realises.
    pub element: (i64, i64),
    /// The element it was matched to on points (the same up to sign).
    pub matched_element: (i64, i64),
    pub norm: u64,
    /// The step degrees in evaluation order.
    pub order: Vec<u32>,
    /// `λ` with `φ = [λ]` on the prime-order group.
    pub eigenvalue: [u64; 4],
    pub glv_basis: [[String; 2]; 2],
    pub babai_bound_bits: u32,
    pub model_generic: ModelOps,
    pub model_optimised: ModelOps,
    pub steps: Vec<StepSpec>,
    pub iso_u: [u64; 4],
    pub vectors: Vec<VectorSpec>,
}

/// The whole frozen set.
#[derive(Clone, Debug)]
pub struct ChainSet {
    pub target: String,
    pub d_k: i64,
    pub class_number: Option<u64>,
    pub omega_root: BigUint,
    pub nmax: u64,
    pub lmax: u64,
    pub export_factor: f64,
    pub chains: Vec<ChainSpec>,
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

fn hex_biguint(s: &str, ctx: &str) -> Result<BigUint, String> {
    let digits = s
        .strip_prefix("0x")
        .ok_or_else(|| format!("{ctx}: {s:?} is not 0x-prefixed hex"))?;
    BigUint::parse_bytes(digits.as_bytes(), 16).ok_or_else(|| format!("{ctx}: {s:?} is not hex"))
}

/// A canonical value below `bound` as little-endian limbs.
fn hex_limbs_below(s: &str, bound: &BigUint, ctx: &str) -> Result<[u64; 4], String> {
    let v = hex_biguint(s, ctx)?;
    if &v >= bound {
        return Err(format!("{ctx}: {s} is not below the modulus"));
    }
    Ok(U256::from_biguint(&v).0)
}

fn pair(v: &Value, ctx: &str) -> Result<(i64, i64), String> {
    let a = as_array(v, ctx)?;
    if a.len() != 2 {
        return Err(format!("{ctx}: expected two integers"));
    }
    Ok((as_i64(&a[0], ctx)?, as_i64(&a[1], ctx)?))
}

fn point(v: &Value, p: &BigUint, ctx: &str) -> Result<[[u64; 4]; 2], String> {
    let a = as_array(v, ctx)?;
    if a.len() != 2 {
        return Err(format!("{ctx}: expected [x, y]"));
    }
    Ok([
        hex_limbs_below(as_str(&a[0], ctx)?, p, ctx)?,
        hex_limbs_below(as_str(&a[1], ctx)?, p, ctx)?,
    ])
}

fn poly(v: &Value, p: &BigUint, ctx: &str) -> Result<Vec<[u64; 4]>, String> {
    as_array(v, ctx)?
        .iter()
        .map(|c| hex_limbs_below(as_str(c, ctx)?, p, ctx))
        .collect()
}

fn decimal(s: &str, ctx: &str) -> Result<BigInt, String> {
    BigInt::parse_bytes(s.as_bytes(), 10)
        .ok_or_else(|| format!("{ctx}: {s:?} is not a signed decimal"))
}

fn model(v: &Value, ctx: &str) -> Result<ModelOps, String> {
    Ok(ModelOps {
        mul: as_u64(field(v, "M", ctx)?, ctx)?,
        sqr: as_u64(field(v, "S", ctx)?, ctx)?,
    })
}

impl ChainSet {
    /// Parse and validate the JSON text of a chain set.
    pub fn parse(text: &str) -> Result<Self, String> {
        let root: Value = serde_json::from_str(text).map_err(|e| format!("not JSON: {e}"))?;
        let format = as_str(field(&root, "format", "chain set")?, "format")?;
        if format != FORMAT {
            return Err(format!("format {format:?}, expected {FORMAT:?}"));
        }
        let p = Uint(P).to_biguint();
        let n = Uint(N).to_biguint();
        for (key, want) in [("p", P), ("a", A), ("b", B), ("n", N)] {
            let got = hex_biguint(as_str(field(&root, key, "chain set")?, key)?, key)?;
            if got != Uint(want).to_biguint() {
                return Err(format!(
                    "curve constant {key} differs from GOST CryptoPro-B"
                ));
            }
        }
        if as_u64(field(&root, "cofactor", "chain set")?, "cofactor")? != 1 {
            return Err("cofactor must be 1".into());
        }
        let bounds = field(&root, "bounds", "chain set")?;
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
                let ell = as_u64(field(st, "ell", &sctx)?, &sctx)?;
                let ell = u32::try_from(ell).map_err(|_| format!("{sctx}: ell too large"))?;
                let dom = field(st, "domain", &sctx)?;
                let cod = field(st, "codomain", &sctx)?;
                steps.push(StepSpec {
                    ell,
                    domain_a: hex_limbs_below(as_str(field(dom, "a", &sctx)?, &sctx)?, &p, &sctx)?,
                    domain_b: hex_limbs_below(as_str(field(dom, "b", &sctx)?, &sctx)?, &p, &sctx)?,
                    codomain_a: hex_limbs_below(
                        as_str(field(cod, "a", &sctx)?, &sctx)?,
                        &p,
                        &sctx,
                    )?,
                    codomain_b: hex_limbs_below(
                        as_str(field(cod, "b", &sctx)?, &sctx)?,
                        &p,
                        &sctx,
                    )?,
                    psi: poly(field(st, "psi", &sctx)?, &p, &sctx)?,
                    n: poly(field(st, "N", &sctx)?, &p, &sctx)?,
                    m: poly(field(st, "M", &sctx)?, &p, &sctx)?,
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
            let mut glv_basis: [[String; 2]; 2] = Default::default();
            for (r, row) in basis.iter().enumerate() {
                let row = as_array(row, &ctx)?;
                if row.len() != 2 {
                    return Err(format!("{ctx}: glv_basis rows have two entries"));
                }
                for (j, e) in row.iter().enumerate() {
                    let s = as_str(e, &ctx)?;
                    decimal(s, &ctx)?;
                    glv_basis[r][j] = s.to_string();
                }
            }
            let mops = field(c, "model_ops", &ctx)?;
            let mut vectors = Vec::new();
            for tv in as_array(field(c, "test_vectors", &ctx)?, &ctx)? {
                vectors.push(VectorSpec {
                    p: point(field(tv, "P", &ctx)?, &p, &ctx)?,
                    phi_p: point(field(tv, "phiP", &ctx)?, &p, &ctx)?,
                    k: hex_biguint(as_str(field(tv, "k", &ctx)?, &ctx)?, &ctx)?,
                    k1: decimal(as_str(field(tv, "k1", &ctx)?, &ctx)?, &ctx)?,
                    k2: decimal(as_str(field(tv, "k2", &ctx)?, &ctx)?, &ctx)?,
                    k_p: point(field(tv, "kP", &ctx)?, &p, &ctx)?,
                });
            }
            chains.push(ChainSpec {
                id,
                element: pair(field(c, "element", &ctx)?, &ctx)?,
                matched_element: pair(field(c, "matched_element", &ctx)?, &ctx)?,
                norm: as_u64(field(c, "norm", &ctx)?, &ctx)?,
                order,
                eigenvalue: hex_limbs_below(
                    as_str(field(c, "eigenvalue", &ctx)?, &ctx)?,
                    &n,
                    &ctx,
                )?,
                glv_basis,
                babai_bound_bits: as_u64(field(c, "babai_bound_bits", &ctx)?, &ctx)? as u32,
                model_generic: model(field(mops, "generic", &ctx)?, &ctx)?,
                model_optimised: model(field(mops, "optimised", &ctx)?, &ctx)?,
                steps,
                iso_u: hex_limbs_below(as_str(field(c, "isomorphism_u", &ctx)?, &ctx)?, &p, &ctx)?,
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
        Ok(ChainSet {
            target: as_str(field(&root, "target", "chain set")?, "target")?.to_string(),
            d_k: as_i64(field(&root, "D_K", "chain set")?, "D_K")?,
            class_number: root.get("class_number").and_then(Value::as_u64),
            omega_root: hex_biguint(
                as_str(field(&root, "omega_root", "chain set")?, "omega_root")?,
                "omega_root",
            )?,
            nmax: as_u64(field(bounds, "nmax", "bounds")?, "nmax")?,
            lmax: as_u64(field(bounds, "lmax", "bounds")?, "lmax")?,
            export_factor: field(bounds, "export_factor", "bounds")?
                .as_f64()
                .ok_or("bounds: export_factor is not a number")?,
            chains,
        })
    }

    /// Read and parse a chain-set file.
    pub fn load(path: &Path) -> Result<Self, String> {
        let text =
            std::fs::read_to_string(path).map_err(|e| format!("read {}: {e}", path.display()))?;
        Self::parse(&text)
    }

    pub fn get(&self, id: &str) -> Option<&ChainSpec> {
        self.chains.iter().find(|c| c.id == id)
    }
}

impl ChainSpec {
    /// `7·5·5` for the step degrees in order.
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

    /// The chain with its constants in Montgomery form, validated.
    pub fn prepared(&self) -> PreparedChain {
        let steps: Vec<StepConstants<'_>> = self
            .steps
            .iter()
            .map(|s| StepConstants {
                ell: s.ell,
                domain_a: &s.domain_a,
                domain_b: &s.domain_b,
                codomain_a: &s.codomain_a,
                codomain_b: &s.codomain_b,
                psi: &s.psi,
                n: &s.n,
                m: &s.m,
            })
            .collect();
        PreparedChain::from_steps(&self.id, self.norm, &steps, &self.iso_u)
    }

    /// The GLV context of this chain (its own `λ` and reduced basis).
    pub fn glv_context(&self) -> GlvContext {
        let basis = [
            [self.glv_basis[0][0].as_str(), self.glv_basis[0][1].as_str()],
            [self.glv_basis[1][0].as_str(), self.glv_basis[1][1].as_str()],
        ];
        GlvContext::from_parts(
            self.prepared(),
            &self.eigenvalue,
            basis,
            self.babai_bound_bits,
        )
    }

    /// The test vector's affine points.
    pub fn vector_points(v: &VectorSpec) -> (CryptoProBAffine, CryptoProBAffine, CryptoProBAffine) {
        (
            CryptoProBAffine::from_limbs(&v.p[0], &v.p[1]),
            CryptoProBAffine::from_limbs(&v.phi_p[0], &v.phi_p[1]),
            CryptoProBAffine::from_limbs(&v.k_p[0], &v.k_p[1]),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ecc::cryptopro_b_point::{
        scalar_mul_glv_with, scalar_mul_wnaf, ChainEval, CryptoProBJacobian, TableMode,
    };
    use crate::ecc::curve::CurveParams;
    use num_bigint::RandBigInt;
    use num_traits::Zero;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    const FROZEN: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/research/cryptopro_b_chain_sweep_20261006/frozen_inputs/chains.constants.json"
    );

    fn frozen_text() -> String {
        std::fs::read_to_string(FROZEN).expect("the frozen chain set is committed")
    }

    fn frozen() -> ChainSet {
        ChainSet::parse(&frozen_text()).expect("the frozen chain set parses")
    }

    #[test]
    fn frozen_set_loads_and_every_chain_prepares() {
        let set = frozen();
        assert_eq!(set.d_k, -619);
        assert_eq!(set.class_number, Some(5));
        assert!(!set.chains.is_empty());
        for c in &set.chains {
            let product: u64 = c.order.iter().map(|&l| u64::from(l)).product();
            assert_eq!(product, c.norm, "{}", c.id);
            let (a, b) = c.element;
            // N(a + bω) = a² + ab + 155b² for D = -619
            assert_eq!((a * a + a * b + 155 * b * b) as u64, c.norm, "{}", c.id);
            assert!(!c.vectors.is_empty(), "{}", c.id);
            let ctx = c.glv_context();
            assert_eq!(ctx.chain.step_degrees(), c.order, "{}", c.id);
        }
    }

    #[test]
    fn frozen_set_reproduces_its_test_vectors() {
        let set = frozen();
        for c in &set.chains {
            let ctx = c.glv_context();
            for v in &c.vectors {
                let (p, phi_p, k_p) = ChainSpec::vector_points(v);
                let phi_want = CryptoProBJacobian::from_affine(&phi_p);
                for eval in [ChainEval::Generic, ChainEval::Optimised] {
                    assert!(
                        ctx.chain.apply_with(&p, eval).eq_point(&phi_want),
                        "{} {eval:?}",
                        c.id
                    );
                }
                let (k1, k2) = ctx.decompose(&v.k);
                assert_eq!((k1, k2), (v.k1.clone(), v.k2.clone()), "{}", c.id);
                for table in [TableMode::Affine, TableMode::Jacobian] {
                    let got = scalar_mul_glv_with(&ctx, &p, &v.k, 4, table, ChainEval::Optimised);
                    assert!(
                        got.eq_point(&CryptoProBJacobian::from_affine(&k_p)),
                        "{}",
                        c.id
                    );
                }
            }
        }
    }

    #[test]
    fn frozen_chains_are_multiplication_by_lambda_and_orderings_agree() {
        let set = frozen();
        let curve = CurveParams::gost_cryptopro_b();
        let g = CryptoProBAffine::from_textbook(&curve.generator()).expect("generator");
        let n = Uint(N).to_biguint();
        let mut rng = StdRng::seed_from_u64(0x5ee9_c4a1);
        let mut by_element: std::collections::BTreeMap<(i64, i64), BigUint> = Default::default();
        for c in &set.chains {
            let ctx = c.glv_context();
            // every ordering of one element is the same endomorphism up to sign
            let lam = ctx.lambda.clone();
            let canon = std::cmp::min(lam.clone(), &n - &lam);
            let prev = by_element.entry(c.element).or_insert_with(|| canon.clone());
            assert_eq!(*prev, canon, "{}: orderings disagree", c.id);
            for _ in 0..3 {
                let r = loop {
                    let r = rng.gen_biguint_below(&n);
                    if !r.is_zero() {
                        break r;
                    }
                };
                let p = scalar_mul_wnaf(&g, &r, 5)
                    .to_affine()
                    .expect("not the identity");
                let want = scalar_mul_wnaf(&p, &ctx.lambda, 5);
                for eval in [ChainEval::Generic, ChainEval::Optimised] {
                    assert!(ctx.chain.apply_with(&p, eval).eq_point(&want), "{}", c.id);
                }
            }
        }
    }

    #[test]
    fn tampered_sets_are_refused() {
        let text = frozen_text();
        let wrong_format = text.replacen(FORMAT, "endosweep-chainsweep/0", 1);
        assert!(ChainSet::parse(&wrong_format)
            .unwrap_err()
            .contains("format"));
        let p_hex = "0x8000000000000000000000000000000000000000000000000000000000000c99";
        assert!(text.contains(p_hex));
        let wrong_p = text.replacen(
            p_hex,
            "0x8000000000000000000000000000000000000000000000000000000000000c9b",
            1,
        );
        assert!(ChainSet::parse(&wrong_p).is_err());
        // a coefficient at or above p is not canonical
        let mut v: Value = serde_json::from_str(&text).expect("JSON");
        v["chains"][0]["steps"][0]["N"][0] = Value::String(p_hex.to_string());
        assert!(ChainSet::parse(&v.to_string())
            .unwrap_err()
            .contains("not below the modulus"));
        // duplicate ids
        let mut v: Value = serde_json::from_str(&text).expect("JSON");
        let first = v["chains"][0].clone();
        v["chains"].as_array_mut().expect("chains").push(first);
        assert!(ChainSet::parse(&v.to_string())
            .unwrap_err()
            .contains("duplicate"));
        // an order that is not the steps' order
        let mut v: Value = serde_json::from_str(&text).expect("JSON");
        v["chains"][0]["order"] = serde_json::json!([3]);
        assert!(ChainSet::parse(&v.to_string()).is_err());
    }

    #[test]
    #[should_panic(expected = "does not compose")]
    fn a_step_that_does_not_compose_is_refused() {
        let mut set = frozen();
        let c = set
            .chains
            .iter_mut()
            .find(|c| c.steps.len() > 1)
            .expect("a multi-step chain");
        c.steps[1].domain_a = c.steps[0].domain_a;
        let _ = c.prepared();
    }

    /// The chain sweep's Python model and this crate's counted operations
    /// agree exactly, for every chain in the set and both evaluators.
    #[cfg(feature = "cryptopro-b-opcount")]
    #[test]
    fn frozen_set_counts_match_the_chain_sweep_model() {
        use crate::ecc::cryptopro_b_field::opcount;
        let set = frozen();
        let curve = CurveParams::gost_cryptopro_b();
        let g = CryptoProBAffine::from_textbook(&curve.generator()).expect("generator");
        for c in &set.chains {
            let ctx = c.glv_context();
            for (eval, want) in [
                (ChainEval::Generic, c.model_generic),
                (ChainEval::Optimised, c.model_optimised),
            ] {
                opcount::reset();
                let _ = ctx.chain.apply_with(&g, eval);
                let got = opcount::read();
                assert_eq!(
                    (got.mul, got.sqr, got.inv),
                    (want.mul, want.sqr, 0),
                    "{} {eval:?}",
                    c.id
                );
            }
        }
    }
}
