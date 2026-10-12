//! Independent cross-check of `examples/hky_last_fall_degree.rs` on the
//! committed `dreg_ladder` draws (`research/fall_degree_bounds_20261008`).
//!
//! Written separately on purpose. It shares only the system construction
//! (`koblitz_groebner::build_decomposition_system`), and computes with a
//! different algorithm:
//!
//! - **Mutant closure the way `lfd_fast.sage` does it.** In rounds: take the
//!   whole current basis, add every variable multiple of every basis element
//!   of degree `< c`, and recompute a fully reduced row echelon form with
//!   columns in descending monomial order (so each pivot is the row's leading
//!   monomial and `V ∩ R_{<c}` is spanned by the rows of degree `< c`).
//!   Repeat until the dimension stops growing. The engine instead inserts
//!   one row at a time from a work queue into a non-reduced echelon basis.
//!   (A first version pivoted on the lowest monomial; that breaks the graded
//!   property, under-generates the closure, and disagreed with the engine.)
//! - **Plain refutation degree** (no mutants): the least `c` at which the
//!   degree-`c` Macaulay space contains `1`. It should reproduce the ladder's
//!   recorded `D*` when the ladder's outcome is a plain refutation.
//!
//! Monomials are ordered by descending `(degree, mask)`, unlike the engine's
//! deglex. Any degree-compatible order gives the same closure.
//!
//! Usage: `hky_crosscheck <ladder.jsonl> --m M --cmax C [--plain-max P]`
use crypto_lib::binary_ecc::F2mElement;
use crypto_lib::cryptanalysis::koblitz_groebner::{build_decomposition_system, FieldStructure};
use crypto_lib::cryptanalysis::koblitz_index_calculus::find_irreducible_sparse;
use serde_json::{json, Value};
use std::collections::HashMap;

struct Space {
    monos: Vec<u64>,
    index: HashMap<u64, usize>,
    w: usize,
}

impl Space {
    fn new(nv: u32, c: u32) -> Space {
        let mut monos: Vec<u64> = (0u64..1 << nv).filter(|m| m.count_ones() <= c).collect();
        // descending (degree, mask): column 0 is the largest monomial, so a
        // lowest-column pivot is the leading monomial, as in lfd_fast.sage
        // (which sorts its columns descending before echelon_form)
        monos.sort_by_key(|&m| std::cmp::Reverse((m.count_ones(), m)));
        let index = monos.iter().enumerate().map(|(i, &m)| (m, i)).collect();
        let w = monos.len().div_ceil(64);
        Space { monos, index, w }
    }
    fn vec(&self, poly: &[u64]) -> Vec<u64> {
        let mut v = vec![0u64; self.w];
        for m in poly {
            let i = self.index[m];
            v[i / 64] ^= 1 << (i % 64);
        }
        v
    }
    fn terms(&self, v: &[u64]) -> Vec<u64> {
        let mut out = vec![];
        for (k, &word) in v.iter().enumerate() {
            let mut b = word;
            while b != 0 {
                out.push(self.monos[k * 64 + b.trailing_zeros() as usize]);
                b &= b - 1;
            }
        }
        out
    }
}

/// Fully reduced row echelon form, pivoting on the lowest set column.
fn rref(mut rows: Vec<Vec<u64>>) -> Vec<Vec<u64>> {
    let mut out: Vec<Vec<u64>> = vec![];
    let mut pivots: Vec<usize> = vec![];
    for mut r in rows.drain(..) {
        for (p, row) in pivots.iter().zip(&out) {
            if r[p / 64] >> (p % 64) & 1 == 1 {
                for (x, y) in r.iter_mut().zip(row) {
                    *x ^= y;
                }
            }
        }
        if let Some(k) = r.iter().position(|&x| x != 0) {
            let p = k * 64 + r[k].trailing_zeros() as usize;
            // keep the basis fully reduced: clear p from earlier rows
            for row in out.iter_mut() {
                if row[p / 64] >> (p % 64) & 1 == 1 {
                    for (x, y) in row.iter_mut().zip(&r) {
                        *x ^= y;
                    }
                }
            }
            out.push(r);
            pivots.push(p);
        }
    }
    out
}

fn mul_mono(poly: &[u64], m: u64) -> Vec<u64> {
    let mut acc: HashMap<u64, bool> = HashMap::new();
    for &t in poly {
        let e = acc.entry(t | m).or_insert(false);
        *e = !*e;
    }
    acc.into_iter().filter(|&(_, b)| b).map(|(t, _)| t).collect()
}

fn deg(poly: &[u64]) -> u32 {
    poly.iter().map(|m| m.count_ones()).max().unwrap_or(0)
}

/// Macaulay rows `m·f`, `deg m ≤ c − deg f`.
fn macaulay(eqs: &[Vec<u64>], sp: &Space, c: u32) -> Vec<Vec<u64>> {
    let mut rows = vec![];
    for f in eqs {
        let d = deg(f);
        if d > c {
            continue;
        }
        for &m in sp.monos.iter().filter(|m| m.count_ones() <= c - d) {
            let p = mul_mono(f, m);
            if deg(&p) <= c {
                rows.push(sp.vec(&p));
            }
        }
    }
    rows
}

fn contains_one(basis: &[Vec<u64>], sp: &Space) -> bool {
    // the constant monomial is the last column; in RREF `1` is a row with only that bit
    let last = sp.monos.len() - 1;
    basis.iter().any(|r| {
        r.iter().enumerate().all(|(k, &x)| {
            if k == last / 64 { x == 1u64 << (last % 64) } else { x == 0 }
        })
    })
}

fn closure_round_based(eqs: &[Vec<u64>], nv: u32, c: u32) -> (usize, bool) {
    let sp = Space::new(nv, c);
    let mut basis = rref(macaulay(eqs, &sp, c));
    loop {
        if contains_one(&basis, &sp) {
            return (sp.monos.len(), true);
        }
        let mut next = basis.clone();
        for r in &basis {
            let t = sp.terms(r);
            if deg(&t) < c {
                for i in 0..nv {
                    next.push(sp.vec(&mul_mono(&t, 1 << i)));
                }
            }
        }
        let nb = rref(next);
        if nb.len() == basis.len() {
            return (nb.len(), contains_one(&nb, &sp));
        }
        basis = nb;
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |k: &str, d: u32| {
        args.iter()
            .position(|a| a == k)
            .map_or(d, |i| args[i + 1].parse().unwrap())
    };
    let (m, cmax, plain_max) = (get("--m", 3) as usize, get("--cmax", 5), get("--plain-max", 0));
    let text = std::fs::read_to_string(&args[1]).expect("read ladder records");
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let rec: Value = serde_json::from_str(line).unwrap();
        if rec.get("draw").is_none() || rec["solutions"].as_u64() != Some(0) {
            continue;
        }
        let n = rec["n"].as_u64().unwrap() as u32;
        let irr = find_irreducible_sparse(n).unwrap();
        let st = FieldStructure::new(n, &irr);
        let el = |w: u64| F2mElement::from_biguint(&num_bigint::BigUint::from(w), n);
        let basis: Vec<F2mElement> = rec["v_basis"]
            .as_array()
            .unwrap()
            .iter()
            .map(|w| el(w.as_u64().unwrap()))
            .collect();
        let sys = build_decomposition_system(&basis, &el(rec["x_r"].as_u64().unwrap()), &F2mElement::one(n), m, &st)
            .unwrap();
        let eqs: Vec<Vec<u64>> = sys
            .equations
            .iter()
            .map(|p| p.terms.iter().map(|t| t.mask).collect())
            .collect();
        let nv = sys.n_vars as u32;
        let mut d_f = None;
        for c in 2..=cmax {
            if closure_round_based(&eqs, nv, c).1 {
                d_f = Some(c);
                break;
            }
        }
        let mut plain = None;
        for c in 2..=plain_max {
            let sp = Space::new(nv, c);
            if contains_one(&rref(macaulay(&eqs, &sp, c)), &sp) {
                plain = Some(c);
                break;
            }
        }
        println!(
            "{}",
            json!({"cell": rec["cell"], "draw": rec["draw"], "n_vars": nv, "d_star_recorded": rec["outcome"],
                   "d_f_crosscheck": d_f, "plain_refutation_degree": plain, "plain_max": plain_max, "cmax": cmax})
        );
    }
}
