//! Analysis of an `fp_backend_bench` run: per-arm medians and the
//! pre-registered claims of `research/fp_backend_20261015/PROTOCOL.md`.
//!
//! ```text
//! cargo run --release --example fp_backend_report -- RUN_DIR
//! ```
//!
//! Reads `RUN_DIR/raw.csv`; writes `RUN_DIR/claims.csv` and
//! `RUN_DIR/tables.md`. A claim compares a baseline arm `B` with a
//! candidate arm `C`: for each repetition `r` the paired ratio is
//! `B_r / C_r` (ns per op, so > 1 means `C` is faster); the verdict uses
//! the median paired ratio against the frozen threshold, and the minimum
//! and maximum paired ratios are reported as the spread.

use std::collections::BTreeMap;
use std::fmt::Write as _;

/// (id, description, baseline arm, candidate arm, threshold): PASS when
/// the median paired ratio baseline/candidate is at least the threshold.
const CLAIMS: &[(&str, &str, &str, &str, f64)] = &[
    // C1–C4: the f4_fp_tower u32 path.
    (
        "C1a",
        "Fp32 dense row kernel, p31, 1e4: no regression vs f4_fp_tower",
        "f4_fp_tower/p31/axpy_dense_1e4",
        "Fp32/p31/axpy_dense_1e4",
        0.95,
    ),
    (
        "C1b",
        "Fp32 dense row kernel, p31, 1e5: no regression vs f4_fp_tower",
        "f4_fp_tower/p31/axpy_dense_1e5",
        "Fp32/p31/axpy_dense_1e5",
        0.95,
    ),
    (
        "C2a",
        "Fp32 sparse scatter, p31: no regression vs f4_fp_tower",
        "f4_fp_tower/p31/scatter_sparse",
        "Fp32/p31/scatter_sparse",
        0.95,
    ),
    (
        "C2b",
        "Fp32 sparse scatter, p32: no regression vs f4_fp_tower",
        "f4_fp_tower/p32/scatter_sparse",
        "Fp32/p32/scatter_sparse",
        0.95,
    ),
    (
        "C3a",
        "Fp32 dense row kernel, p32, 1e4: >= 2x f4_fp_tower (eager above 2^31)",
        "f4_fp_tower/p32/axpy_dense_1e4",
        "Fp32/p32/axpy_dense_1e4",
        2.0,
    ),
    (
        "C3b",
        "Fp32 dense row kernel, p32, 1e5: >= 2x f4_fp_tower",
        "f4_fp_tower/p32/axpy_dense_1e5",
        "Fp32/p32/axpy_dense_1e5",
        2.0,
    ),
    (
        "C4",
        "Fp32 mul throughput, p31: no regression vs f4_fp_tower Fp::mul",
        "f4_fp_tower/p31/mul_tp",
        "Fp32/p31/mul_tp",
        0.95,
    ),
    // C5: prime_fast::Fp64.
    (
        "C5a",
        "Fp64 mul throughput, p61: no regression vs prime_fast::Fp64",
        "prime_fast::Fp64/p61/mul_tp",
        "Fp64/p61/mul_tp",
        0.95,
    ),
    (
        "C5b",
        "Fp64 mul latency, p61: no regression vs prime_fast::Fp64",
        "prime_fast::Fp64/p61/mul_lat",
        "Fp64/p61/mul_lat",
        0.95,
    ),
    (
        "C5c",
        "Fp64 inverse, p61: >= 1.4x prime_fast::Fp64",
        "prime_fast::Fp64/p61/inv",
        "Fp64/p61/inv",
        1.4,
    ),
    // C6: ct_bignum MontgomeryContext.
    (
        "C6a",
        "FpP192 mul throughput >= 1.5x ct_bignum<3>, P-192",
        "ct_bignum/P-192/mul_tp",
        "FpP192/P-192/mul_tp",
        1.5,
    ),
    (
        "C6b",
        "FpP192 mul latency >= 2x ct_bignum<3>, P-192",
        "ct_bignum/P-192/mul_lat",
        "FpP192/P-192/mul_lat",
        2.0,
    ),
    (
        "C6c",
        "FpWide<3> mul throughput >= 1.0x ct_bignum<3>, P-192",
        "ct_bignum/P-192/mul_tp",
        "FpWide<3>/P-192/mul_tp",
        1.0,
    ),
    (
        "C6d",
        "FpWide<3> mul latency >= 1.25x ct_bignum<3>, P-192",
        "ct_bignum/P-192/mul_lat",
        "FpWide<3>/P-192/mul_lat",
        1.25,
    ),
    (
        "C6e",
        "FpWide<4+spare> mul throughput >= 1.3x ct_bignum<4>, P-224",
        "ct_bignum/P-224/mul_tp",
        "FpWide<4+spare>/P-224/mul_tp",
        1.3,
    ),
    (
        "C6f",
        "FpWide<4+spare> mul latency >= 1.8x ct_bignum<4>, P-224",
        "ct_bignum/P-224/mul_lat",
        "FpWide<4+spare>/P-224/mul_lat",
        1.8,
    ),
    (
        "C6g",
        "FpWide<4> mul throughput >= 1.1x ct_bignum<4>, P-256",
        "ct_bignum/P-256/mul_tp",
        "FpWide<4>/P-256/mul_tp",
        1.1,
    ),
    (
        "C6h",
        "FpWide<4> mul latency >= 1.3x ct_bignum<4>, P-256",
        "ct_bignum/P-256/mul_lat",
        "FpWide<4>/P-256/mul_lat",
        1.3,
    ),
    // C7: BigUint.
    (
        "C7a",
        "FpWide<3> mul >= 5x BigUint (a*b) % p, P-192",
        "BigUint/P-192/mul_tp",
        "FpWide<3>/P-192/mul_tp",
        5.0,
    ),
    (
        "C7b",
        "FpWide<4+spare> mul >= 5x BigUint, P-224",
        "BigUint/P-224/mul_tp",
        "FpWide<4+spare>/P-224/mul_tp",
        5.0,
    ),
    (
        "C7c",
        "FpWide<4> mul >= 4x BigUint, P-256",
        "BigUint/P-256/mul_tp",
        "FpWide<4>/P-256/mul_tp",
        4.0,
    ),
    (
        "C7d",
        "FpWide<3> inverse >= 5x BigUint modpow(p-2), P-192",
        "BigUint/P-192/inv",
        "FpWide<3>/P-192/inv",
        5.0,
    ),
    (
        "C7e",
        "FpWide<4+spare> inverse >= 5x BigUint modpow(p-2), P-224",
        "BigUint/P-224/inv",
        "FpWide<4+spare>/P-224/inv",
        5.0,
    ),
    // C8–C9: lazy against eager reduction.
    (
        "C8a",
        "Fp32 lazy row update >= 5x eager, p31, 1e4",
        "Fp32/p31/axpy_eager_1e4",
        "Fp32/p31/axpy_lazy_1e4",
        5.0,
    ),
    (
        "C8b",
        "Fp64 lazy row update >= 2x eager, p61, 1e4",
        "Fp64/p61/axpy_eager_1e4",
        "Fp64/p61/axpy_lazy_1e4",
        2.0,
    ),
    (
        "C8c",
        "FpWide<3> lazy row update >= 1.5x eager, P-192, 1e4",
        "FpWide<3>/P-192/axpy_eager_1e4",
        "FpWide<3>/P-192/axpy_lazy_1e4",
        1.5,
    ),
    (
        "C8d",
        "FpWide<3> lazy row update >= 1.5x eager, P-192, 1e5",
        "FpWide<3>/P-192/axpy_eager_1e5",
        "FpWide<3>/P-192/axpy_lazy_1e5",
        1.5,
    ),
    (
        "C8e",
        "FpWide<4+spare> lazy row update >= 1.5x eager, P-224, 1e4",
        "FpWide<4+spare>/P-224/axpy_eager_1e4",
        "FpWide<4+spare>/P-224/axpy_lazy_1e4",
        1.5,
    ),
    (
        "C8f",
        "FpWide<4+spare> lazy row update >= 1.5x eager, P-224, 1e5",
        "FpWide<4+spare>/P-224/axpy_eager_1e5",
        "FpWide<4+spare>/P-224/axpy_lazy_1e5",
        1.5,
    ),
    (
        "C8g",
        "FpWide<4> lazy row update >= 1.5x eager, P-256, 1e4",
        "FpWide<4>/P-256/axpy_eager_1e4",
        "FpWide<4>/P-256/axpy_lazy_1e4",
        1.5,
    ),
    (
        "C9a",
        "FpWide<3> lazy dot >= 1.5x eager, P-192",
        "FpWide<3>/P-192/dot_eager_1e4",
        "FpWide<3>/P-192/dot_lazy_1e4",
        1.5,
    ),
    (
        "C9b",
        "FpWide<4+spare> lazy dot >= 1.5x eager, P-224",
        "FpWide<4+spare>/P-224/dot_eager_1e4",
        "FpWide<4+spare>/P-224/dot_lazy_1e4",
        1.5,
    ),
    // C10–C12: inversion, Legendre, square root choices.
    (
        "C10a",
        "binary-GCD (Pornin) inverse >= 1.3x plain binary, FpWide<4>, P-224",
        "FpWide<4>/P-224/inv_plain_binary",
        "FpWide<4>/P-224/inv",
        1.3,
    ),
    (
        "C10b",
        "binary-GCD (Pornin) inverse >= 2.5x Fermat, FpWide<4>, P-224",
        "FpWide<4>/P-224/inv_fermat",
        "FpWide<4>/P-224/inv",
        2.5,
    ),
    (
        "C11a",
        "binary Jacobi >= 3x Euler, FpWide<3>, P-192",
        "FpWide<3>/P-192/legendre_euler",
        "FpWide<3>/P-192/legendre",
        3.0,
    ),
    (
        "C11b",
        "binary Jacobi >= 3x Euler, FpWide<4>, P-224",
        "FpWide<4>/P-224/legendre_euler",
        "FpWide<4>/P-224/legendre",
        3.0,
    ),
    (
        "C12",
        "windowed Tonelli-Shanks >= 3x plain, FpWide<4>, P-224 (s = 96)",
        "FpWide<4>/P-224/sqrt_ts_plain",
        "FpWide<4>/P-224/sqrt",
        3.0,
    ),
    // C13–C14: the specialisations earn their place.
    (
        "C13",
        "Solinas FpP192 mul >= 1.3x Montgomery FpWide<3>, P-192",
        "FpWide<3>/P-192/mul_tp",
        "FpP192/P-192/mul_tp",
        1.3,
    ),
    (
        "C14",
        "spare-bit FpWide<4+spare> mul >= 1.15x FpWide<4>, P-224",
        "FpWide<4>/P-224/mul_tp",
        "FpWide<4+spare>/P-224/mul_tp",
        1.15,
    ),
];

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

fn main() {
    let dir = std::env::args().nth(1).expect("RUN_DIR");
    let raw = std::fs::read_to_string(format!("{dir}/raw.csv")).expect("raw.csv");
    // arm -> rep -> ns
    let mut by_arm: BTreeMap<String, BTreeMap<usize, f64>> = BTreeMap::new();
    let mut meta: BTreeMap<String, (String, String, String)> = BTreeMap::new();
    for line in raw.lines().skip(1) {
        let f: Vec<&str> = line.split(',').collect();
        assert_eq!(f.len(), 9, "bad line {line}");
        let rep: usize = f[0].parse().unwrap();
        let ns: f64 = f[8].parse().unwrap();
        by_arm.entry(f[2].to_string()).or_default().insert(rep, ns);
        meta.insert(f[2].to_string(), (f[3].into(), f[4].into(), f[5].into()));
    }
    let med = |arm: &str| -> Option<f64> {
        by_arm
            .get(arm)
            .map(|m| median(&mut m.values().copied().collect::<Vec<_>>()))
    };

    let mut claims = String::from(
        "id,description,baseline,candidate,baseline_median_ns,candidate_median_ns,ratio_median,ratio_min,ratio_max,n_pairs,threshold,verdict\n",
    );
    let mut md = String::from("## Pre-registered claims\n\n| id | claim | baseline ns | candidate ns | ratio (median of paired) | paired min–max | threshold | verdict |\n|:--|:--|--:|--:|--:|:--|--:|:--|\n");
    let (mut pass, mut fail, mut missing) = (0, 0, 0);
    for &(id, desc, base, cand, thr) in CLAIMS {
        let (Some(b), Some(c)) = (by_arm.get(base), by_arm.get(cand)) else {
            missing += 1;
            writeln!(claims, "{id},\"{desc}\",{base},{cand},,,,,,0,{thr},MISSING").unwrap();
            writeln!(md, "| {id} | {desc} | – | – | – | – | {thr} | MISSING |").unwrap();
            continue;
        };
        let mut ratios: Vec<f64> = b
            .iter()
            .filter_map(|(r, &bv)| c.get(r).map(|&cv| bv / cv))
            .collect();
        let n = ratios.len();
        let rmed = median(&mut ratios);
        let (rmin, rmax) = (ratios[0], ratios[n - 1]);
        let verdict = if rmed >= thr { "PASS" } else { "FAIL" };
        if rmed >= thr {
            pass += 1;
        } else {
            fail += 1;
        }
        let (bm, cm) = (med(base).unwrap(), med(cand).unwrap());
        writeln!(
            claims,
            "{id},\"{desc}\",{base},{cand},{bm:.4},{cm:.4},{rmed:.4},{rmin:.4},{rmax:.4},{n},{thr},{verdict}"
        )
        .unwrap();
        writeln!(
            md,
            "| {id} | {desc} | {bm:.2} | {cm:.2} | {rmed:.3} | {rmin:.3}–{rmax:.3} | {thr} | {verdict} |"
        )
        .unwrap();
    }
    writeln!(
        md,
        "\n{pass} PASS, {fail} FAIL, {missing} MISSING of {}.\n",
        CLAIMS.len()
    )
    .unwrap();

    // Per-backend tables: rows = backend/prime, columns = ops.
    let ops = [
        "mul_tp",
        "mul_lat",
        "sqr_lat",
        "add_tp",
        "inv",
        "batch_inv",
        "sqrt",
        "legendre",
        "axpy_lazy_1e4",
        "axpy_eager_1e4",
        "axpy_lazy_1e5",
        "axpy_eager_1e5",
        "dot_lazy_1e4",
        "dot_eager_1e4",
    ];
    let mut rows: Vec<(String, String)> = Vec::new();
    for (b, p, _) in meta.values() {
        if !rows.contains(&(b.clone(), p.clone())) {
            rows.push((b.clone(), p.clone()));
        }
    }
    md.push_str("## Median ns per operation (all arms)\n\n| backend | prime |");
    for o in ops {
        write!(md, " {o} |").unwrap();
    }
    md.push_str("\n|:--|:--|");
    for _ in ops {
        md.push_str("--:|");
    }
    md.push('\n');
    for (b, p) in &rows {
        if !ops
            .iter()
            .any(|o| by_arm.contains_key(&format!("{b}/{p}/{o}")))
        {
            continue;
        }
        write!(md, "| {b} | {p} |").unwrap();
        for o in ops {
            match med(&format!("{b}/{p}/{o}")) {
                Some(v) if v >= 100.0 => write!(md, " {v:.0} |").unwrap(),
                Some(v) => write!(md, " {v:.2} |").unwrap(),
                None => md.push_str(" – |"),
            }
        }
        md.push('\n');
    }
    md.push_str("\n## Other arms (median ns per op, min–max over repetitions)\n\n| arm | median | min | max |\n|:--|--:|--:|--:|\n");
    for (arm, m) in &by_arm {
        let op = &meta[arm].2;
        if ops.contains(&op.as_str()) {
            continue;
        }
        let mut v: Vec<f64> = m.values().copied().collect();
        let md_ = median(&mut v);
        writeln!(
            md,
            "| `{arm}` | {md_:.2} | {:.2} | {:.2} |",
            v[0],
            v[v.len() - 1]
        )
        .unwrap();
    }
    std::fs::write(format!("{dir}/claims.csv"), claims).unwrap();
    std::fs::write(format!("{dir}/tables.md"), &md).unwrap();
    println!("{md}");
}
