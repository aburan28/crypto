//! How expensive is the `Parameterized` (Gröbner) cache layer?
//!
//! `algebra_cache` has three layers. `Preprocessing` is the one on the
//! relation-search path and it is cheap -- see `preprocessing_cost.rs`, which
//! measures it. `Parameterized` is the other extreme: `parameter_basis_cached`
//! runs Buchberger over F2 on `n_vars + n` Boolean variables, is capped at 16,
//! and is marked offline-experiment-only. This probe puts a number on it,
//! because "is any cached artifact expensive enough to deserve a durable tier"
//! cannot be answered without one.
//!
//! **Run every point under an external `timeout`.** The Buchberger engine has
//! no interrupt hook, so a point past the feasible edge does not come back:
//!
//!     for pt in "3 2 2" "5 2 2" "5 3 2"; do
//!       timeout 45 cargo run --release --example gb_probe -- $pt
//!     done
//!
//! Measured 2026-09-21, 4-core x86_64 container, rustc 1.94.1, release:
//!
//!     n=3 ell=2 m=2   7 vars      49 ms
//!     n=5 ell=2 m=2   9 vars    5340 ms
//!     n=5 ell=3 m=2  11 vars   did not finish in 45 s
//!     n=7 ell=2 m=2, n=5 ell=2 m=3, n=7 ell=3 m=2, n=9 ell=2 m=2
//!                                did not finish in 45 s
//!
//! Two more variables cost about 100x, which is what Buchberger over F2 should
//! do and why the 16-variable cap is where it is. Scoped to these points on
//! this machine; nothing here is a claim about any attack's cost.

use crypto_lib::binary_ecc::F2mElement;
use crypto_lib::cryptanalysis::algebra_cache::AlgebraCache;
use crypto_lib::cryptanalysis::koblitz_groebner::FieldStructure;
use crypto_lib::cryptanalysis::koblitz_index_calculus::find_irreducible;
use crypto_lib::cryptanalysis::polynomial_reuse::{parameter_basis_cached, DecompositionTemplate};
use num_bigint::BigUint;
use std::time::Instant;

fn fe(x: u64, n: u32) -> F2mElement {
    F2mElement::from_biguint(&BigUint::from(x), n)
}

fn main() {
    // One point per process, so a point that does not terminate bounds only
    // itself. With no arguments, the smallest known-feasible point.
    let a: Vec<u32> = std::env::args()
        .skip(1)
        .filter_map(|x| x.parse().ok())
        .collect();
    let (n, ell, m) = if a.len() == 3 {
        (a[0], a[1] as usize, a[2] as usize)
    } else {
        println!("usage: gb_probe <n> <ell> <m>   (run under `timeout`; defaulting to 3 2 2)");
        (3, 2, 2)
    };

    println!(
        "{:>4} {:>4} {:>3} {:>7} {:>7}  {:>14}",
        "n", "ell", "m", "n_vars", "total", "groebner"
    );

    let Some(irr) = find_irreducible(n) else {
        println!("{n:>4} {ell:>4} {m:>3}   no irreducible of that degree");
        return;
    };
    let st = FieldStructure::new(n, &irr);
    let basis: Vec<_> = (0..ell).map(|k| fe(1u64 << k, n)).collect();
    let Some(t) = DecompositionTemplate::build(&basis, &fe(1, n), m, &st) else {
        println!("{n:>4} {ell:>4} {m:>3}   layout refused (MAX_VARS, or n > 64)");
        return;
    };
    let total = t.n_vars + n as usize;
    if total > 16 {
        println!(
            "{n:>4} {ell:>4} {m:>3} {:>7} {total:>7}   over the 16-var cap",
            t.n_vars
        );
        return;
    }

    let mut cache = AlgebraCache::local(64 * 1024 * 1024);
    let start = Instant::now();
    let got = parameter_basis_cached(&t, &mut cache);
    let ns = start.elapsed().as_nanos();
    match got {
        Some(g) => println!(
            "{n:>4} {ell:>4} {m:>3} {:>7} {total:>7}  {:>11.1} ms   ({} generators)",
            t.n_vars,
            ns as f64 / 1e6,
            g.len(),
        ),
        None => println!("{n:>4} {ell:>4} {m:>3} {:>7} {total:>7}   refused", t.n_vars),
    }
}
