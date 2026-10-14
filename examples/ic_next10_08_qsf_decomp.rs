//! IC next-10 item 8, arm D: relation yield and refutation degree of a
//! quasi-subfield (Frobenius-stable) subspace factor base against random
//! subspaces of the same dimension, in the `m = 2` binary decomposition system
//! `S3(x1, x2, x_R) = 0` on `y^2 + xy = x^3 + 1` over `F_{2^n}`.
//!
//! A base is either `g=<hex>` (the subspace `ker g(sigma)`, `sigma` = squaring,
//! `g | t^n - 1`; `g = t^k + 1` is the subfield `F_{2^k}`) or `random=<seed>`.
//!
//! Yield: for each uniformly random target abscissa `x_R` of a curve point,
//! count ordered pairs `(x1, x2) in V^2`, both curve abscissae, with
//! `S3(x1, x2, x_R) = 0`, exactly: loop over `x1 in V`, solve the quadratic in
//! `x2`, test `x2 in V`.
//!
//! Refutation degree: for targets with no algebraic solution in `V^2` (verified
//! by the same exact sweep, curve test dropped), the least `d` at which the
//! degree-`d` matrix-F4 closure of the Weil-descended system contains `1`.
//!
//! ```sh
//! cargo run --release --example ic_next10_08_qsf_decomp -- \
//!     <n> <k> <base> <yield_targets> <refute_targets> <d_max> <seed>
//! ```
//! Prints one JSON line.
use crypto_lib::binary_ecc::F2mElement;
use crypto_lib::cryptanalysis::koblitz_groebner::{
    build_decomposition_system, matrix_f4_f2_counted, solving_profile, FieldStructure,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::find_irreducible_sparse;
use crypto_lib::cryptanalysis::pq_groebner_f2::F2BoolPoly;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::time::Instant;

/// Echelon form of an F_2 subspace of u64 bitmasks, with a membership test
/// and a solver for the linear map whose images were inserted.
struct Echelon {
    rows: Vec<(u32, u64, u64)>, // (pivot bit, image, preimage combination)
}
impl Echelon {
    fn new() -> Self {
        Echelon { rows: Vec::new() }
    }
    fn reduce(&self, mut v: u64) -> (u64, u64) {
        let mut comb = 0u64;
        for &(p, r, c) in &self.rows {
            if (v >> p) & 1 == 1 {
                v ^= r;
                comb ^= c;
            }
        }
        (v, comb)
    }
    /// insert image `v` with preimage `pre`; returns false if dependent
    fn insert(&mut self, v: u64, pre: u64) -> bool {
        let (r, c) = self.reduce(v);
        if r == 0 {
            return false;
        }
        let p = 63 - r.leading_zeros();
        let c = c ^ pre;
        // keep fully reduced: clear bit p from existing rows
        for row in self.rows.iter_mut() {
            if (row.1 >> p) & 1 == 1 {
                row.1 ^= r;
                row.2 ^= c;
            }
        }
        self.rows.push((p, r, c));
        true
    }
    fn solve(&self, v: u64) -> Option<u64> {
        let (r, c) = self.reduce(v);
        if r == 0 {
            Some(c)
        } else {
            None
        }
    }
}

fn elem(bits: u64, n: u32) -> F2mElement {
    F2mElement::from_bit_positions(&(0..n).filter(|i| (bits >> i) & 1 == 1).collect::<Vec<_>>(), n)
}

fn trace(gf: &Gf2, a: u64, n: u32) -> u64 {
    let mut acc = 0u64;
    let mut t = a;
    for _ in 0..n {
        acc ^= t;
        t = gf.sqr(t);
    }
    acc
}

struct Ctx {
    gf: Gf2,
    n: u32,
    qsolve: Echelon, // solves z^2 + z = c
}

impl Ctx {
    /// roots x of  A x^2 + B x + C = 0 where S3(u, x, w) with
    /// A = (u+w)^2, B = u w, C = (u w)^2 + b, b = 1.
    fn s3_roots(&self, u: u64, w: u64, out: &mut Vec<u64>) {
        out.clear();
        let gf = &self.gf;
        let a = gf.sqr(u ^ w);
        let bb = gf.mul(u, w);
        let c = gf.sqr(bb) ^ 1;
        if a == 0 {
            if bb != 0 {
                out.push(gf.mul(c, gf.inv(bb)));
            }
            return;
        }
        if bb == 0 {
            // x^2 = C/A, unique square root
            let r = gf.mul(c, gf.inv(a));
            out.push(gf.sqr_k(r, self.n - 1));
            return;
        }
        // x = (B/A) z, z^2 + z = C A / B^2
        let t = gf.mul(bb, gf.inv(a));
        let rhs = gf.mul(gf.mul(c, a), gf.inv(gf.sqr(bb)));
        if let Some(z) = self.qsolve.solve(rhs) {
            out.push(gf.mul(t, z));
            out.push(gf.mul(t, z ^ 1));
        }
    }
    fn on_curve_x(&self, x: u64) -> bool {
        // y^2 + xy = x^3 + 1 has a solution iff x = 0 or Tr(x + 1/x^2) = 0
        if x == 0 {
            return true;
        }
        let gf = &self.gf;
        trace(gf, x ^ gf.inv(gf.sqr(x)), self.n) == 0
    }
}

fn kernel_basis(cols: &[u64], n: u32) -> Vec<u64> {
    // kernel of the map e_k -> cols[k]
    let mut ech = Echelon::new();
    let mut ker = Vec::new();
    for k in 0..n {
        let (r, c) = ech.reduce(cols[k as usize]);
        if r == 0 {
            ker.push(c ^ (1u64 << k));
        } else {
            ech.insert(cols[k as usize], 1u64 << k);
        }
    }
    ker
}

fn signature(polys: &[F2BoolPoly]) -> Vec<Vec<u64>> {
    let mut s: Vec<Vec<u64>> = polys
        .iter()
        .map(|p| {
            let mut m: Vec<u64> = p.terms.iter().map(|t| t.mask).collect();
            m.sort_unstable();
            m
        })
        .filter(|m: &Vec<u64>| !m.is_empty())
        .collect();
    s.sort();
    s
}

fn is_one(p: &F2BoolPoly) -> bool {
    p.terms.len() == 1 && p.terms[0].mask == 0
}

/// degree-d closure; returns (refuted, converged, xors, size_cap)
fn closure(polys: &[F2BoolPoly], n_vars: usize, d: u32, max_rounds: usize) -> (bool, bool, u64, bool) {
    let mut gens = polys.to_vec();
    let mut prev = signature(&gens);
    let mut xors = 0u64;
    for _ in 0..max_rounds {
        let (reduced, x) = match matrix_f4_f2_counted(&gens, n_vars, d) {
            Some(r) => r,
            None => return (false, false, xors, true),
        };
        xors += x;
        let next: Vec<_> = reduced.into_iter().filter(|p| !p.terms.is_empty()).collect();
        if next.iter().any(is_one) {
            return (true, true, xors, false);
        }
        let sig = signature(&next);
        if sig == prev {
            return (false, true, xors, false);
        }
        prev = sig;
        gens = next;
    }
    (false, false, xors, false)
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let n: u32 = a[0].parse().unwrap();
    let k: u32 = a[1].parse().unwrap();
    let base = a[2].clone();
    let n_yield: usize = a[3].parse().unwrap();
    let n_ref: usize = a[4].parse().unwrap();
    let d_max: u32 = a[5].parse().unwrap();
    let seed: u64 = a[6].parse().unwrap();
    let onepass = std::env::var("QSF_ONEPASS").map(|v| v == "1").unwrap_or(false);
    let irr = find_irreducible_sparse(n).unwrap();
    let gf = Gf2::new(&irr);
    let fmask = (1u64 << n) - 1;
    let mut rng = StdRng::seed_from_u64(seed);

    // quadratic solver z^2 + z
    let mut qsolve = Echelon::new();
    for i in 0..n {
        let z = 1u64 << i;
        qsolve.insert(gf.sqr(z) ^ z, z);
    }
    let ctx = Ctx { gf, n, qsolve };
    let gf = &ctx.gf;

    // basis
    let basis: Vec<u64> = if let Some(h) = base.strip_prefix("g=") {
        let g = u128::from_str_radix(h.trim_start_matches("0x"), 16).unwrap();
        let cols: Vec<u64> = (0..n)
            .map(|kk| {
                let x = 1u64 << kk;
                let mut acc = 0u64;
                for i in 0..=127u32 {
                    if (g >> i) & 1 == 1 {
                        acc ^= gf.sqr_k(x, i % n);
                    }
                }
                acc
            })
            .collect();
        kernel_basis(&cols, n)
    } else {
        let s: u64 = base.strip_prefix("random=").unwrap().parse().unwrap();
        let mut r2 = StdRng::seed_from_u64(s);
        let mut e = Echelon::new();
        let mut b = Vec::new();
        while b.len() < k as usize {
            let v = r2.gen::<u64>() & fmask;
            if v != 0 && e.insert(v, 0) {
                b.push(v);
            }
        }
        b
    };
    assert_eq!(basis.len(), k as usize, "basis dimension {} != k {}", basis.len(), k);
    let mut vech = Echelon::new();
    for &b in &basis {
        vech.insert(b, 0);
    }
    let in_v = |x: u64| vech.reduce(x).0 == 0;
    let v: Vec<u64> = (0..(1u64 << k))
        .map(|m| {
            let mut acc = 0;
            for (i, &b) in basis.iter().enumerate() {
                if (m >> i) & 1 == 1 {
                    acc ^= b;
                }
            }
            acc
        })
        .collect();
    let v_on_curve: Vec<u64> = v.iter().copied().filter(|&x| ctx.on_curve_x(x)).collect();

    // exact per-target sweep: (#ordered curve decompositions, #algebraic solutions)
    let mut roots = Vec::new();
    let mut sweep = |xr: u64| -> (u64, u64) {
        let (mut curve, mut alg) = (0u64, 0u64);
        for &x1 in &v {
            ctx.s3_roots(x1, xr, &mut roots);
            let x1c = ctx.on_curve_x(x1);
            for &x2 in roots.iter() {
                if in_v(x2) {
                    alg += 1;
                    if x1c && ctx.on_curve_x(x2) {
                        curve += 1;
                    }
                }
            }
        }
        (curve, alg)
    };

    // ---- yield ----
    let t0 = Instant::now();
    let mut decomposed = 0usize;
    let mut total_dec = 0u64;
    let mut done = 0usize;
    while done < n_yield {
        let xr = rng.gen::<u64>() & fmask;
        if xr == 0 || !ctx.on_curve_x(xr) {
            continue;
        }
        let (c, _) = sweep(xr);
        if c > 0 {
            decomposed += 1;
        }
        total_dec += c;
        done += 1;
    }
    let yield_secs = t0.elapsed().as_secs_f64();
    eprintln!(
        "{}",
        serde_json::json!({"partial": "yield", "n": n, "k": k, "base": base, "decomposed": decomposed,
            "yield_targets": n_yield, "mean_ordered_decomps": total_dec as f64 / n_yield as f64})
    );

    // ---- refutation degree ----
    let st = FieldStructure::new(n, &irr);
    let basis_e: Vec<F2mElement> = basis.iter().map(|&b| elem(b, n)).collect();
    let b_e = F2mElement::one(n);
    let mut ref_rows = Vec::new();
    let mut tries = 0;
    while ref_rows.len() < n_ref && tries < 100000 {
        tries += 1;
        let xr = rng.gen::<u64>() & fmask;
        if xr == 0 || !ctx.on_curve_x(xr) {
            continue;
        }
        let (_, alg) = sweep(xr);
        if alg > 0 {
            continue; // need a verified non-decomposable target
        }
        let sys = build_decomposition_system(&basis_e, &elem(xr, n), &b_e, 2, &st).unwrap();
        let t1 = Instant::now();
        let mut first = None;
        let mut xors_total = 0u64;
        let mut notes = Vec::new();
        for d in 2..=d_max {
            // QSF_ONEPASS=1: one degree-d Macaulay matrix (no closure rounds),
            // which fits the memory cap at n = 31, d = 4.
            let (refuted, converged, xors, cap) = if onepass {
                match solving_profile(&sys.equations, sys.n_vars, d) {
                    Some(p) => (p.refuted, true, 0, false),
                    None => (false, false, 0, true),
                }
            } else {
                closure(&sys.equations, sys.n_vars, d, 12)
            };
            xors_total += xors;
            notes.push(format!(
                "d{d}:{}",
                if refuted { "R" } else if cap { "CAP" } else if converged { "noR" } else { "nc" }
            ));
            if refuted {
                first = Some(d);
                break;
            }
            if cap {
                break;
            }
        }
        ref_rows.push(serde_json::json!({
            "x_r": format!("{xr:x}"), "refute_degree": first, "steps": notes.join(" "),
            "word_xors": xors_total, "secs": t1.elapsed().as_secs_f64(), "vars": sys.n_vars,
        }));
    }

    println!(
        "{}",
        serde_json::json!({
            "n": n, "k": k, "base": base, "seed": seed,
            "basis_hex": basis.iter().map(|b| format!("{b:x}")).collect::<Vec<_>>(),
            "v_size": v.len(), "v_on_curve": v_on_curve.len(),
            "yield_targets": n_yield, "decomposed": decomposed,
            "yield": decomposed as f64 / n_yield as f64,
            "mean_ordered_decomps": total_dec as f64 / n_yield as f64,
            "yield_secs": yield_secs,
            "refute": ref_rows, "refute_tries": tries, "d_max": d_max,
            "refute_mode": if onepass { "one-pass Macaulay" } else { "F4 closure" },
        })
    );
}
