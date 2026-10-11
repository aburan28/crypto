//! Export the Gröbner-basis workloads this repository's engines are timed
//! on, in the input languages of msolve, Singular and Magma, so the same
//! systems can be run through those systems and their bases checked
//! against ours.
//!
//! ```text
//! cargo run --release --example gb_cas_export -- fp <outdir> [repeats]
//! cargo run --release --example gb_cas_export -- fp-xl <outdir> [repeats]
//! cargo run --release --example gb_cas_export -- bool <node_dump.jsonl> <outdir>
//! ```
//!
//! `fp`: the random dense systems of `examples/f4_fp_bench.rs` (same
//! generator, same seeds; the basis fingerprints printed here equal the
//! bench's).  For each case it times `f4_fp::f4` (grevlex, the bench's
//! degree bound; median of `repeats`) and writes
//!
//! - `<case>.ms`: msolve input;
//! - `<case>.sing`: a Singular script that times `std` (reduced, `dp`)
//!   and checks that our basis is a reduced Gröbner basis of the same
//!   ideal (`ours ⊆ std(I)`, `I ⊆ ours` by normal forms, equal size and
//!   leading terms);
//! - `<case>.magma`: a Magma script that times `GroebnerBasis` (F4, the
//!   default) and compares with our basis;
//!
//! and prints one JSON line with our timing and basis shape.
//!
//! `bool`: the root systems of a decomposition search, read from a node
//! dump (`KIC_F4_NODE_DUMP=<file>`, one line per reduced node; a root is
//! a line on which every variable occurs).  Each is written with the
//! field equations `x² + x` (msolve, Singular) or in Magma's Boolean
//! polynomial ring; the decision to compare is whether the ideal is the
//! unit ideal (no decomposition of the target).
//!
//! The exporter only writes inputs and our own measured timings; the
//! external tools are run by `research/gb_cas_comparison_20261011/run.sh`.

use crypto_lib::cryptanalysis::f4_fp::{f4, F4Options, Ordering, Poly};
use crypto_lib::cryptanalysis::pq_groebner_f2::F2BoolPoly;
use serde_json::json;
use std::collections::hash_map::DefaultHasher;
use std::fmt::Write as _;
use std::hash::{Hash, Hasher};
use std::path::Path;

// ── the f4_fp_bench systems (copied verbatim: same generator, same seeds) ──

fn monomials(n: usize, d: u32) -> Vec<Vec<u32>> {
    fn rec(n: usize, d: u32, cur: &mut Vec<u32>, out: &mut Vec<Vec<u32>>) {
        if cur.len() == n {
            out.push(cur.clone());
            return;
        }
        let used: u32 = cur.iter().sum();
        for k in 0..=d - used {
            cur.push(k);
            rec(n, d, cur, out);
            cur.pop();
        }
    }
    let mut out = Vec::new();
    rec(n, d, &mut Vec::new(), &mut out);
    out
}

fn system(n: usize, m: usize, d: u32, p: u64, seed: u64) -> Vec<Poly> {
    let mut s = seed | 1;
    let mut rnd = move || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        s
    };
    let monos = monomials(n, d);
    (0..m)
        .map(|_| monos.iter().map(|e| (e.clone(), rnd() % p)).collect())
        .collect()
}

fn fingerprint<T: Hash>(x: &T) -> String {
    let mut h = DefaultHasher::new();
    x.hash(&mut h);
    format!("{:016x}", h.finish())
}

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.total_cmp(b));
    v[v.len() / 2]
}

/// (name, n, m, deg, p, max_degree): `f4_fp_bench`'s cases, standard then
/// large.
const FP_CASES: &[(&str, usize, usize, u32, u64, u32)] = &[
    ("quad_n4_p31", 4, 4, 2, 31, 12),
    ("quad_n5_p31", 5, 5, 2, 31, 12),
    ("quad_n6_p65521", 6, 6, 2, 65521, 12),
    ("quad_n7_p65521", 7, 7, 2, 65521, 12),
    ("cubic_n3_p29", 3, 3, 3, 29, 12),
    ("cubic_n4_p31", 4, 4, 3, 31, 14),
    ("quart_n3_p31", 3, 3, 4, 31, 16),
    ("overdet_n6_m9_p101", 6, 9, 2, 101, 10),
    ("quad_n8_p65521", 8, 8, 2, 65521, 12),
    ("quad_n9_p65521", 9, 9, 2, 65521, 14),
    ("cubic_n5_p31", 5, 5, 3, 31, 16),
    ("overdet_n10_m15_p65521", 10, 15, 2, 65521, 10),
    ("overdet_n12_m24_p65521", 12, 24, 2, 65521, 10),
];

/// Larger cases, exported with `fp-xl`: where the linear algebra of a
/// full-size F4 dominates (not in `f4_fp_bench`'s fixed suite).
const FP_XL_CASES: &[(&str, usize, usize, u32, u64, u32)] = &[
    ("quad_n10_p65521", 10, 10, 2, 65521, 16),
    ("quad_n11_p65521", 11, 11, 2, 65521, 18),
    ("cubic_n6_p31", 6, 6, 3, 31, 18),
];

// ── writers ────────────────────────────────────────────────────────────

/// A term `c·x1^e1·…` in the shared syntax of msolve, Singular and Magma
/// (`*` and `^`); variables are `x1 … xn`.
fn term(c: u64, e: &[u32], out: &mut String) {
    let mut first = true;
    if c != 1 || e.iter().all(|&k| k == 0) {
        write!(out, "{c}").unwrap();
        first = false;
    }
    for (i, &k) in e.iter().enumerate() {
        if k == 0 {
            continue;
        }
        if !first {
            out.push('*');
        }
        first = false;
        write!(out, "x{}", i + 1).unwrap();
        if k > 1 {
            write!(out, "^{k}").unwrap();
        }
    }
}

fn poly_fp(p: &Poly) -> String {
    let mut s = String::new();
    let mut any = false;
    for (e, c) in p {
        if *c == 0 {
            continue;
        }
        if any {
            s.push('+');
        }
        term(*c, e, &mut s);
        any = true;
    }
    if !any {
        s.push('0');
    }
    s
}

fn poly_bool(p: &F2BoolPoly, n: usize) -> String {
    let mut s = String::new();
    for (k, t) in p.terms.iter().enumerate() {
        if k > 0 {
            s.push('+');
        }
        let e: Vec<u32> = (0..n).map(|i| ((t.mask >> i) & 1) as u32).collect();
        term(1, &e, &mut s);
    }
    if p.terms.is_empty() {
        s.push('0');
    }
    s
}

fn vars(n: usize) -> String {
    (1..=n)
        .map(|i| format!("x{i}"))
        .collect::<Vec<_>>()
        .join(",")
}

fn write_msolve(path: &Path, n: usize, p: u64, polys: &[String]) {
    let mut s = format!("{}\n{p}\n", vars(n));
    s.push_str(&polys.join(",\n"));
    s.push('\n');
    std::fs::write(path, s).expect("write msolve input");
}

/// Singular: time `std` of the input (reduced, `dp`), then check `ours`
/// against it.  Prints `std_ms`, the basis size and `same_ideal` /
/// `same_reduced_basis` (1 or 0).
fn write_singular(path: &Path, n: usize, p: u64, input: &[String], ours: &[String]) {
    let s = format!(
        "ring r = {p},({vars}),dp;\noption(redSB);\nsystem(\"--ticks-per-sec\",1000);\n\
         ideal I = {input};\n\
         int t0 = rtimer;\nideal G = std(I);\nint t = rtimer - t0;\n\
         ideal O = {ours};\n\
         int sub1 = (size(reduce(O, G)) == 0);\n\
         int sub2 = (size(reduce(I, std(O))) == 0);\n\
         int same = 0;\nif (size(G) == size(O)) {{ same = (size(reduce(O, G)) == 0) && (size(reduce(G, std(O))) == 0) && (size(simplify(lead(G),1)) == size(simplify(lead(O),1))); }}\n\
         print(\"std_ms \" + string(t));\nprint(\"std_size \" + string(size(G)));\n\
         print(\"ours_size \" + string(size(O)));\n\
         print(\"same_ideal \" + string(sub1 && sub2));\nprint(\"same_size \" + string(same));\nquit;\n",
        vars = vars(n),
        input = input.join(",\n  "),
        ours = ours.join(",\n  "),
    );
    std::fs::write(path, s).expect("write Singular script");
}

/// Magma: time `GroebnerBasis` (F4 is Magma's default for these
/// coefficient rings) and compare with our basis.  Not run here: Magma
/// is not installed in this container.
fn write_magma(path: &Path, n: usize, p: u64, input: &[String], ours: &[String]) {
    let names = (1..=n)
        .map(|i| format!("x{i}"))
        .collect::<Vec<_>>()
        .join(",");
    let s = format!(
        "// {name}: run as `magma -b {name}`; prints wall and CPU seconds of GroebnerBasis.\n\
         K := GF({p});\nP<{names}> := PolynomialRing(K, {n}, \"grevlex\");\n\
         I := ideal<P | {input}>;\n\
         t := Cputime(); w := Realtime();\nG := GroebnerBasis(I);\n\
         printf \"magma_gb_cpu_s %o\\nmagma_gb_wall_s %o\\n\", Cputime(t), Realtime(w);\n\
         O := [P | {ours}];\n\
         printf \"magma_basis_size %o\\nours_size %o\\n\", #G, #O;\n\
         printf \"same_ideal %o\\n\", ideal<P | O> eq I;\n\
         printf \"same_reduced_basis %o\\n\", Set(G) eq Set(O);\n\
         quit;\n",
        name = path.file_name().unwrap().to_string_lossy(),
        input = input.join(",\n  "),
        ours = ours.join(",\n  "),
    );
    std::fs::write(path, s).expect("write Magma script");
}

fn run_fp(outdir: &Path, repeats: usize, cases: &[(&str, usize, usize, u32, u64, u32)]) {
    for &(name, n, m, d, p, dmax) in cases {
        let sys = system(
            n,
            m,
            d,
            p,
            0x9e37_79b9_7f4a_7c15 ^ (n as u64 * 131 + d as u64),
        );
        let opts = F4Options::new(Ordering::Grevlex, dmax);
        let mut times = Vec::new();
        let mut report = None;
        for _ in 0..repeats {
            let r = f4(&sys, n, p, &opts);
            times.push(r.ms);
            report = Some(r);
        }
        let r = report.expect("repeats ≥ 1");
        let input: Vec<String> = sys.iter().map(poly_fp).collect();
        let ours: Vec<String> = r.basis.iter().map(poly_fp).collect();
        write_msolve(&outdir.join(format!("{name}.ms")), n, p, &input);
        write_singular(&outdir.join(format!("{name}.sing")), n, p, &input, &ours);
        write_magma(&outdir.join(format!("{name}.magma")), n, p, &input, &ours);
        println!(
            "{}",
            json!({
                "case": name, "n": n, "m": m, "deg": d, "p": p, "max_degree": dmax,
                "ours_f4_ms_median": median(times.clone()), "ours_f4_ms_all": times,
                "repeats": repeats, "basis_len": r.basis.len(), "basis_fp": fingerprint(&r.basis),
                "steps": r.steps, "max_rows": r.max_rows, "max_cols": r.max_cols,
                "inconsistent": r.inconsistent,
            })
        );
    }
}

fn run_bool(dump: &Path, outdir: &Path) {
    let text = std::fs::read_to_string(dump).expect("read node dump");
    let mut k = 0;
    for line in text.lines() {
        let v: serde_json::Value = serde_json::from_str(line).expect("dump line");
        let n = v["n_vars"].as_u64().expect("n_vars") as usize;
        let system: Vec<F2BoolPoly> = serde_json::from_value(v["system"].clone()).expect("system");
        let occurring = system
            .iter()
            .flat_map(|p| p.terms.iter())
            .fold(0u64, |a, t| a | t.mask);
        if occurring.count_ones() as usize != n {
            continue; // not a root: some variable already assigned
        }
        k += 1;
        let name = format!("root{k}");
        let input: Vec<String> = system.iter().map(|p| poly_bool(p, n)).collect();
        let mut with_field = input.clone();
        with_field.extend((1..=n).map(|i| format!("x{i}^2+x{i}")));
        write_msolve(&outdir.join(format!("{name}.ms")), n, 2, &with_field);
        let sing = format!(
            "ring r = 2,({vars}),dp;\noption(redSB);\nsystem(\"--ticks-per-sec\",1000);\n\
             ideal I = {input};\n\
             int t0 = rtimer;\nideal G = std(I);\nint t = rtimer - t0;\n\
             print(\"std_ms \" + string(t));\nprint(\"std_size \" + string(size(G)));\n\
             print(\"unit_ideal \" + string(G[1] == 1));\nquit;\n",
            vars = vars(n),
            input = with_field.join(",\n  "),
        );
        std::fs::write(outdir.join(format!("{name}.sing")), sing).expect("write");
        let names = vars(n);
        let magma = format!(
            "// {name}.magma: root decomposition system, Boolean ring (field equations implicit).\n\
             B<{names}> := BooleanPolynomialRing({n}, \"grevlex\");\n\
             I := ideal<B | {input}>;\n\
             t := Cputime(); w := Realtime();\nG := GroebnerBasis(I);\n\
             printf \"magma_gb_cpu_s %o\\nmagma_gb_wall_s %o\\n\", Cputime(t), Realtime(w);\n\
             printf \"unit_ideal %o\\n\", G eq [B!1];\nquit;\n",
            input = input.join(",\n  "),
        );
        std::fs::write(outdir.join(format!("{name}.magma")), magma).expect("write");
        println!(
            "{}",
            json!({"root": name, "n_vars": n, "generators": system.len(),
                   "terms": system.iter().map(|p| p.terms.len()).sum::<usize>(),
                   "max_degree": system.iter().filter_map(|p| p.terms.first()).map(|t| t.mask.count_ones()).max()})
        );
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some(tier @ ("fp" | "fp-xl")) => {
            let outdir = Path::new(args.get(2).expect("outdir"));
            std::fs::create_dir_all(outdir).expect("outdir");
            let repeats = args.get(3).map_or(3, |s| s.parse().expect("repeats"));
            run_fp(
                outdir,
                repeats,
                if tier == "fp" { FP_CASES } else { FP_XL_CASES },
            );
        }
        Some("bool") => {
            let dump = Path::new(args.get(2).expect("node dump"));
            let outdir = Path::new(args.get(3).expect("outdir"));
            std::fs::create_dir_all(outdir).expect("outdir");
            run_bool(dump, outdir);
        }
        _ => eprintln!("usage: gb_cas_export fp <outdir> [repeats] | bool <dump.jsonl> <outdir>"),
    }
}
