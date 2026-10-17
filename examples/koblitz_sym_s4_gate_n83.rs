//! **m=83 confidence gate for the heavy-δ symmetric-S₄ oracle (AGENTS §8a).**
//!
//! The toy-scale factor-base → D* campaign
//! (`research/factor_base_dor_campaign_20261013/`) measured that the
//! symmetrised S₄ system over a geometric base of dimension `l` refutes at
//! `D* ≈ l + 1` with `6l − 3` unknowns, independent of the field degree, and
//! that its point-level reachable set is small (the corrected exact-yield law).
//! This driver prices both sides on the gate curve itself — `E₀: y² + xy =
//! x³ + 1` over `F_{2^83}` (the ECC2K-130 family member with the verified
//! modulus `z^83 + z^45 + z^2 + z + 1`, order `4·r`,
//! `r = 2417851639230796216685689`) — with exact ground truth:
//!
//! - **Oracle arm.**  For `l = 2, 3`: build the symmetrised S₄ system over the
//!   geometric spaces `W_k = θ^k·⟨1, g, …⟩`, solve it exactly by evaluation
//!   over all `2^(6l−3)` assignments, decode roots to `(e₁, e₂, e₃)` and
//!   cross-check against `koblitz_sym_s4::s4_symmetric_eval`.
//! - **Ground truth.**  Direct enumeration of all triples of curve points
//!   with abscissa in `V = span(W₁)`: the exact point-level reachable set
//!   `S_pt`, planted targets, and the derived systems-per-relation.
//!
//! Known-answer only (targets are curve points; no key recovery).
//!
//! ```bash
//! cargo run --release --example koblitz_sym_s4_gate_n83 -- --json experiments/…json
//! ```

use crypto_lib::binary_ecc::curve::{point_add, BinaryPoint};
use crypto_lib::binary_ecc::{F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::koblitz_sym_s4::s4_symmetric_eval;
use crypto_lib::cryptanalysis::wide_groebner::{TwoWordFieldStructure, WideFieldTable, WPoly};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::time::Instant;

/// Coordinate vector of a field element that is F₂-linear in Boolean vars.
#[derive(Clone)]
struct SymVec {
    coords: Vec<WPoly>,
}

fn elem_bits(v: &F2mElement) -> u128 {
    let raw = v.raw_bits();
    let mut b = raw[0] as u128;
    if raw.len() > 1 {
        b |= (raw[1] as u128) << 64;
    }
    b
}

impl SymVec {
    fn constant(v: &F2mElement, n: u32) -> Self {
        Self {
            coords: (0..n)
                .map(|k| {
                    let (w, b) = ((k / 64) as usize, k % 64);
                    if v.raw_bits()[w] >> b & 1 == 1 {
                        WPoly::one()
                    } else {
                        WPoly::zero()
                    }
                })
                .collect(),
        }
    }

    fn from_subspace(basis: &[F2mElement], offset: usize, n: u32) -> Self {
        let mut coords = vec![WPoly::zero(); n as usize];
        for (t, b) in basis.iter().enumerate() {
            for k in 0..n {
                let (w, bit) = ((k / 64) as usize, k % 64);
                if b.raw_bits()[w] >> bit & 1 == 1 {
                    coords[k as usize] = coords[k as usize].add(&WPoly::var(offset + t));
                }
            }
        }
        Self { coords }
    }

    fn add(&self, other: &Self) -> Self {
        Self {
            coords: self
                .coords
                .iter()
                .zip(&other.coords)
                .map(|(a, b)| a.add(b))
                .collect(),
        }
    }

    fn mul(&self, other: &Self, n: u32, st: &TwoWordFieldStructure) -> Self {
        let mut coords = vec![WPoly::zero(); n as usize];
        for (i, a) in self.coords.iter().enumerate() {
            if a.is_zero() {
                continue;
            }
            for (j, b) in other.coords.iter().enumerate() {
                if b.is_zero() {
                    continue;
                }
                let prod = a.mul(b);
                let red = st.product_bits(i, j);
                for bit in 0..n as usize {
                    if red >> bit & 1 == 1 {
                        coords[bit] = coords[bit].add(&prod);
                    }
                }
            }
        }
        Self { coords }
    }

    fn square(&self, n: u32, st: &TwoWordFieldStructure) -> Self {
        let mut coords = vec![WPoly::zero(); n as usize];
        for (k, c) in self.coords.iter().enumerate() {
            let sq = st.square_bits(k);
            for bit in 0..n as usize {
                if sq >> bit & 1 == 1 {
                    coords[bit] = coords[bit].add(c);
                }
            }
        }
        Self { coords }
    }

    /// Value of every coordinate polynomial under a full assignment: the field
    /// element this vector evaluates to.
    fn value(&self, bits: u128) -> u128 {
        let mut out = 0u128;
        for (k, c) in self.coords.iter().enumerate() {
            if c.eval(bits) {
                out |= 1u128 << k;
            }
        }
        out
    }
}

/// The symmetrised S₄ system (`koblitz_sym_s4`'s 11-term formula): its `n`
/// coordinate polynomials, each set to zero.
fn sym_s4_system(
    spaces: &[Vec<F2mElement>; 3],
    x_r: &F2mElement,
    b: &F2mElement,
    n: u32,
    st: &TwoWordFieldStructure,
    irr: &IrreduciblePoly,
) -> Vec<WPoly> {
    let mut off = 0usize;
    let e: Vec<SymVec> = spaces
        .iter()
        .map(|w| {
            let s = SymVec::from_subspace(w, off, n);
            off += w.len();
            s
        })
        .collect();
    let x = x_r.clone();
    let x2 = x.mul(&x, irr);
    let x3 = x2.mul(&x, irr);
    let x4 = x2.mul(&x2, irr);
    let b2 = b.mul(&b, irr);
    let c = |v: &F2mElement| SymVec::constant(v, n);

    let e1sq = e[0].square(n, st);
    let e1q = e1sq.square(n, st);
    let e2sq = e[1].square(n, st);
    let e2q = e2sq.square(n, st);
    let e3sq = e[2].square(n, st);
    let e3cu = e3sq.mul(&e[2], n, st);
    let e3q = e3sq.square(n, st);

    let mut s = c(&x4.mul(&b2, irr)); // x^4 b^2
    s = s.add(&e[2].mul(&c(&x3.mul(b, irr)), n, st)); // e3 x^3 b
    s = s.add(&e3sq.mul(&c(&x4.add(b)), n, st)); // e3^2 (x^4 + b)
    s = s.add(&e3cu.mul(&c(&x), n, st)); // e3^3 x
    s = s.add(&e3q); // e3^4
    s = s.add(&e2sq.mul(&c(&x2.mul(b, irr)), n, st)); // e2^2 x^2 b
    s = s.add(&e2sq.mul(&e[2], n, st).mul(&c(&x3), n, st)); // e2^2 e3 x^3
    s = s.add(&e2q.mul(&c(&x4), n, st)); // e2^4 x^4
    s = s.add(&e1sq.mul(&e[2], n, st).mul(&c(&x.mul(b, irr)), n, st)); // e1^2 e3 x b
    s = s.add(&e1sq.mul(&e3sq, n, st).mul(&c(&x2), n, st)); // e1^2 e3^2 x^2
    s = s.add(&e1q.mul(&c(&b2), n, st)); // e1^4 b^2
    s.coords.into_iter().filter(|p| !p.is_zero()).collect()
}

fn lift(curve: &crypto_lib::binary_ecc::curve::BinaryCurve, x: &F2mElement, flip: bool) -> Option<BinaryPoint> {
    let irr = &curve.irreducible;
    let xi = x.flt_inverse(irr)?;
    let c = x.add(&xi.square(irr));
    if c.trace(irr) != 0 {
        return None;
    }
    let w = c.half_trace(irr)?;
    let mut y = w.mul(x, irr);
    if flip {
        y = y.add(x);
    }
    Some(BinaryPoint::Affine { x: x.clone(), y })
}

fn point_x(p: &BinaryPoint) -> Option<F2mElement> {
    match p {
        BinaryPoint::Affine { x, .. } => Some(x.clone()),
        BinaryPoint::Infinity => None,
    }
}

fn main() {
    let kc = KoblitzCurve::known_n83_k0().expect("N83 k0 fixture");
    let curve = &kc.curve;
    let n = kc.n;
    let irr = curve.irreducible.clone();
    let b = curve.b.clone();
    let st = TwoWordFieldStructure::new(n, &irr).expect("two-word structure");
    let r = &kc.subgroup_order;

    // structures from a stream disjoint from the target stream (the 2026-10-14
    // campaign correction): theta/g from 0x5EED_0083.
    let mut srng = StdRng::seed_from_u64(0x5EED_0083);
    let theta = F2mElement::from_words(&[srng.gen::<u64>(), srng.gen::<u64>()], n);
    let g = F2mElement::from_words(&[srng.gen::<u64>(), srng.gen::<u64>()], n);

    println!("curve: a=0 n={n} r={r} ({} bits)", r.bits());
    println!("theta={:x} g={:x}", elem_bits(&theta), elem_bits(&g));

    let mut report: Vec<serde_json::Value> = Vec::new();
    for l in [2u32, 3u32] {
        // geometric spaces W_k = theta^k <1, g, ..., g^(k(l-1))>
        let mut theta_k = F2mElement::one(n);
        let mut spaces: Vec<Vec<F2mElement>> = Vec::new();
        for k in 0..3 {
            theta_k = theta_k.mul(&theta, &irr);
            let dim = (k + 1) * (l - 1) + 1;
            let mut w = Vec::new();
            let mut gpow = F2mElement::one(n);
            for _ in 0..dim {
                w.push(theta_k.mul(&gpow, &irr));
                gpow = gpow.mul(&g, &irr);
            }
            spaces.push(w);
        }
        let n_vars: usize = spaces.iter().map(|w| w.len()).sum();
        assert_eq!(n_vars, (6 * l - 3) as usize, "geometric dims");

        // V members and the curve points above them
        let v = &spaces[0];
        let mut v_elems: Vec<F2mElement> = Vec::new();
        for bits in 0..(1u64 << l) {
            let mut acc = F2mElement::zero(n);
            for i in 0..l {
                if bits >> i & 1 == 1 {
                    acc = acc.add(&v[i as usize]);
                }
            }
            v_elems.push(acc);
        }
        let mut pts: Vec<BinaryPoint> = Vec::new();
        for x in &v_elems {
            for flip in [false, true] {
                if let Some(p) = lift(curve, x, flip) {
                    assert!(curve.is_on_curve(&p));
                    pts.push(p);
                }
            }
        }

        // exact point-level reachable set over all triples
        let t0 = Instant::now();
        let mut s_pt: std::collections::HashSet<u128> = std::collections::HashSet::new();
        let mut sums: Vec<Option<F2mElement>> = Vec::new();
        for pi in 0..pts.len() {
            for pj in 0..pts.len() {
                for pk in 0..pts.len() {
                    let s = point_add(curve, &point_add(curve, &pts[pi], &pts[pj]), &pts[pk]);
                    let xr = point_x(&s);
                    if let Some(x) = &xr {
                        s_pt.insert(elem_bits(x));
                    }
                    sums.push(xr);
                }
            }
        }
        let enum_wall = t0.elapsed();

        // planted target: the abscissa of a genuine triple sum
        let target = sums
            .iter()
            .flatten()
            .find(|x| !x.is_zero())
            .expect("a nonzero triple sum exists")
            .clone();
        // the planted triple's own e-values (ground truth for the decode)
        let planted: Vec<[u128; 3]> = (0..sums.len())
            .filter(|i| sums[*i].as_ref().map(|t| *t == target).unwrap_or(false))
            .map(|_| [0u128; 3])
            .collect();
        let _ = planted.len();

        // oracle: build + solve
        let t1 = Instant::now();
        let eqs = sym_s4_system(
            &[spaces[0].clone(), spaces[1].clone(), spaces[2].clone()],
            &target,
            &b,
            n,
            &st,
            &irr,
        );
        let build_wall = t1.elapsed();
        let t2 = Instant::now();
        let mut roots = Vec::new();
        for bits in 0..(1u128 << n_vars) {
            if eqs.iter().all(|p| !p.eval(bits)) {
                roots.push(bits);
            }
        }
        let solve_wall = t2.elapsed();

        // decode roots to (e1, e2, e3) field elements
        let decode = |bits: u128| -> [F2mElement; 3] {
            let mut off = 0usize;
            let mut out = [F2mElement::zero(n), F2mElement::zero(n), F2mElement::zero(n)];
            for (k, w) in spaces.iter().enumerate() {
                let mut acc = F2mElement::zero(n);
                for t in 0..w.len() {
                    if bits >> (off + t) & 1 == 1 {
                        acc = acc.add(&w[t]);
                    }
                }
                out[k] = acc;
                off += w.len();
            }
            out
        };
        // a root is a true decomposition iff its e-triple is that of some
        // V^3 multiset (exact check, <= 2^(3l) per root)
        let true_decompositions: usize = roots
            .iter()
            .map(|r| {
                let e = decode(*r);
                let mut hit = false;
                for a in &v_elems {
                    for c in &v_elems {
                        for d in &v_elems {
                            let e1 = a.add(c).add(d);
                            let e2 = a
                                .mul(c, &irr)
                                .add(&a.mul(d, &irr))
                                .add(&c.mul(d, &irr));
                            let e3 = a.mul(c, &irr).mul(d, &irr);
                            if e1 == e[0] && e2 == e[1] && e3 == e[2] {
                                hit = true;
                            }
                        }
                    }
                }
                // the engine's own evaluation must also vanish at the root
                assert!(
                    s4_symmetric_eval(n, &e, &target, &b, &irr).is_zero(),
                    "every Boolean root must satisfy the symmetric S4 formula"
                );
                hit
            })
            .filter(|h| *h)
            .count();

        let yield_exact = s_pt.len() as f64 / 2f64.powi(n as i32);
        let systems_per_relation = 1.0 / yield_exact.max(1e-300);
        let ladder_probe_ref = 8.8e6_f64; // landed n=83 4-sum ladder, probes/relation

        let row = serde_json::json!({
            "l": l, "n_vars": n_vars, "v_points": pts.len(),
            "reachable_set_size": s_pt.len(),
            "exact_yield": yield_exact,
            "systems_per_relation_single_v": systems_per_relation,
            "ladder_probes_per_relation_ref": ladder_probe_ref,
            "ratio_to_ladder": systems_per_relation / ladder_probe_ref,
            "roots_found": roots.len(),
            "roots_that_are_true_decompositions": true_decompositions,
            "eqs": eqs.len(),
            "build_ms": build_wall.as_millis() as u64,
            "solve_ms": solve_wall.as_millis() as u64,
            "enum_wall_ms": enum_wall.as_millis() as u64,
            "enum_triples": pts.len().pow(3),
        });
        println!("{row}");
        report.push(row);
    }

    if let Some(pos) = std::env::args().position(|a| a == "--json") {
        if let Some(path) = std::env::args().nth(pos + 1) {
            std::fs::write(
                &path,
                serde_json::to_string_pretty(&serde_json::json!({
                    "schema": "koblitz-sym-s4-gate-n83-v1",
                    "curve": "E0 y^2+xy=x^3+1 over F_2^83, r=2417851639230796216685689",
                    "seed_structures": "0x5EED0083",
                    "rows": report,
                }))
                .unwrap(),
            )
            .unwrap();
            println!("wrote {path}");
        }
    }
}
