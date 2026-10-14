//! Linear-algebra engine benchmark for the index-calculus logarithm solve.
//!
//! Generates a deterministic synthetic relation system shaped like an
//! index-calculus relation matrix (each row is `m` factor-base summands with
//! coefficients ±1 or uniform residues, columns uniform, rows drawn until
//! every column is covered and the system has full column rank, plus a small
//! excess), plants a logarithm vector, and solves it with one engine:
//!
//! | arm | engine |
//! |:--|:--|
//! | `igauss` | `IncrementalGauss` (dense RREF, rows fed one at a time; ecbench default) |
//! | `sgauss` | `StructuredGauss` (sparse RREF, Markowitz pivot, rows fed one at a time) |
//! | `full-bw` | `solve_sparse_system`: filter (dup/singleton/clique/merge) + block Wiedemann + reconstruct |
//! | `full-bl` | the same with block Lanczos |
//! | `core-bw` / `core-bl` | the engine alone on the *core* (`export-core` below), filter disabled |
//!
//! `export-core` writes the singleton- and clique-filtered, *unmerged*,
//! square core (coefficients stay ±1/±2) in the CADO-NFS `bwc` binary format
//! (`matrix.bin`, `rhs.txt`), so CADO's block Wiedemann runs on exactly the
//! matrix `core-bw` sees.  `check-cado` verifies a CADO solution file
//! against that core and against the planted logarithms.
//!
//! Every arm checks its answer against the planted vector; a run that does
//! not reproduce it is reported `verified: false`.  One JSON line per run.
//!
//! ```text
//! cargo run --release --example la_engine_bench -- ours --arm full-bw --cols 16384 --weight 3
//! cargo run --release --example la_engine_bench -- export-core --cols 16384 --weight 3 --dir /tmp/c
//! cargo run --release --example la_engine_bench -- check-cado --cols 16384 --weight 3 --dir /tmp/c --sol /tmp/c/wdir/K.sols0-1.0
//! ```

use crypto_lib::cryptanalysis::ic_framework::linalg::matrix_by_name;
use crypto_lib::cryptanalysis::ic_boundary::RowStatus;
use crypto_lib::cryptanalysis::{
    filter_relations, solve_sparse_system, FilterOptions, SparseCoreSolver, SparseRow,
    SparseSolveOptions, SparseSolveOutcome,
};
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::collections::HashMap;
use std::io::{Read, Write};
use std::time::Instant;

/// 2^62 − 57, prime; fits the `u64` word path and CADO's one-limb `p1` layer.
const DEFAULT_P: u64 = 4_611_686_018_427_387_847;

fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}
fn addmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 + b as u128) % m as u128) as u64
}

#[derive(Clone, Debug)]
struct Params {
    cols: usize,
    weight: usize,
    excess: usize,
    seed: u64,
    p: u64,
    wide: bool,
    /// All coefficients +1 (a repeated column gives +2): the shape CADO's
    /// relation format can carry (nonnegative valuations).  Draws the same
    /// random stream as `unit`, so the column pattern is identical.
    positive: bool,
    block: usize,
    /// `--merge-max-weight` / `--max-row-weight`: override the filter's merge
    /// thresholds (defaults 8 and 32) for the `full-*` arms.
    merge_max_weight: Option<usize>,
    max_row_weight: Option<usize>,
}

struct System {
    rows: Vec<SparseRow<u64>>,
    planted: Vec<u64>,
    /// Rows drawn before every column was covered (coupon-collector point).
    rows_to_cover: usize,
}

/// Draw rows until every column occurs, then `excess` more.  Full column
/// rank is not forced here; every arm reports it if the rows fall short.
fn generate(p: &Params) -> System {
    let mut rng = StdRng::seed_from_u64(p.seed);
    let planted: Vec<u64> = (0..p.cols).map(|_| rng.gen_range(1..p.p)).collect();
    let mut covered = vec![false; p.cols];
    let mut uncovered = p.cols;
    let mut rows = Vec::new();
    let mut rows_to_cover = 0;
    let mut extra = 0usize;
    while extra < p.excess {
        let mut entries = Vec::with_capacity(p.weight);
        for _ in 0..p.weight {
            let c = rng.gen_range(0..p.cols) as u32;
            let v = if p.wide {
                rng.gen_range(1..p.p)
            } else if rng.gen::<bool>() || p.positive {
                1
            } else {
                p.p - 1
            };
            entries.push((c, v));
        }
        let row = SparseRow::new(entries, 0, p.p);
        if row.entries.is_empty() {
            continue;
        }
        let mut rhs = 0u64;
        for &(c, v) in &row.entries {
            rhs = addmod(rhs, mulmod(v, planted[c as usize], p.p), p.p);
            if !covered[c as usize] {
                covered[c as usize] = true;
                uncovered -= 1;
            }
        }
        rows.push(SparseRow {
            entries: row.entries,
            rhs,
        });
        if uncovered == 0 {
            if rows_to_cover == 0 {
                rows_to_cover = rows.len();
            }
            extra += 1;
        }
    }
    System {
        rows,
        planted,
        rows_to_cover,
    }
}

/// Singleton + clique filtering only (no merge), trimmed to square, so the
/// coefficients stay small integers that CADO's format can carry.
fn unmerged_core(p: &Params, sys: &System) -> (Vec<SparseRow<u64>>, Vec<u32>) {
    let opts = FilterOptions {
        remove_duplicates: true,
        remove_singletons: true,
        target_excess: 0,
        merge_max_weight: 0,
        max_row_weight: usize::MAX,
    };
    let f = filter_relations(sys.rows.clone(), p.cols, p.p, &opts).expect("consistent");
    let index: HashMap<u32, u32> = f
        .core_columns
        .iter()
        .enumerate()
        .map(|(i, &c)| (c, i as u32))
        .collect();
    let core = f
        .core_rows
        .iter()
        .map(|r| SparseRow {
            entries: r.entries.iter().map(|&(c, v)| (index[&c], v)).collect(),
            rhs: r.rhs,
        })
        .collect();
    (core, f.core_columns.clone())
}

fn sparse_opts(p: &Params, solver: SparseCoreSolver, filter: bool) -> SparseSolveOptions {
    let mut o = SparseSolveOptions {
        solver,
        ..SparseSolveOptions::default()
    };
    o.wiedemann.block_m = p.block;
    o.wiedemann.block_n = p.block;
    o.lanczos.block_size = p.block;
    if let Some(w) = p.merge_max_weight {
        o.filter.merge_max_weight = w;
    }
    if let Some(w) = p.max_row_weight {
        o.filter.max_row_weight = w;
    }
    if !filter {
        o.filter = FilterOptions {
            remove_duplicates: false,
            remove_singletons: false,
            target_excess: usize::MAX,
            merge_max_weight: 0,
            max_row_weight: usize::MAX,
        };
    }
    o
}

fn threads() -> usize {
    rayon::current_num_threads()
}

fn run_ours(p: &Params, arm: &str) -> serde_json::Value {
    let t0 = Instant::now();
    let sys = generate(p);
    let gen_s = t0.elapsed().as_secs_f64();
    let nnz: usize = sys.rows.iter().map(|r| r.entries.len()).sum();
    let mut out = serde_json::json!({
        "arm": arm, "cols": p.cols, "weight": p.weight, "excess": p.excess,
        "seed": p.seed, "modulus": p.p.to_string(), "coeffs": coeffs_name(p),
        "block": p.block, "threads": threads(),
        "merge_max_weight": p.merge_max_weight, "max_row_weight": p.max_row_weight,
        "rows": sys.rows.len(), "rows_to_cover": sys.rows_to_cover, "nonzeros": nnz,
        "generate_s": gen_s,
    });
    let o = out.as_object_mut().unwrap();
    match arm {
        "igauss" | "sgauss" => {
            let name = if arm == "igauss" { "incremental-gauss" } else { "structured-gauss" };
            let mut m = matrix_by_name(name, p.cols, p.p).unwrap();
            let t = Instant::now();
            let mut fed = 0usize;
            for row in &sys.rows {
                let mut dense = vec![0u64; p.cols];
                for &(c, v) in &row.entries {
                    dense[c as usize] = v;
                }
                fed += 1;
                if m.add_row(dense, row.rhs) == RowStatus::Inconsistent {
                    break;
                }
                if m.rank() == p.cols {
                    break;
                }
            }
            let solve_s = t.elapsed().as_secs_f64();
            let solved: Vec<Option<u64>> = (0..p.cols).map(|c| m.pinned(c)).collect();
            let ok = m.rank() == p.cols
                && solved.iter().zip(&sys.planted).all(|(s, &x)| *s == Some(x));
            let (work, unit) = m.work();
            o.insert("rows_fed".into(), fed.into());
            o.insert("rank".into(), m.rank().into());
            o.insert("work".into(), work.into());
            o.insert("work_unit".into(), unit.into());
            o.insert("solve_s".into(), solve_s.into());
            o.insert("verified".into(), ok.into());
        }
        "full-bw" | "full-bl" => {
            let solver = if arm == "full-bw" {
                SparseCoreSolver::BlockWiedemann
            } else {
                SparseCoreSolver::BlockLanczos
            };
            let opts = sparse_opts(p, solver, true);
            let tf = Instant::now();
            let f = filter_relations(sys.rows.clone(), p.cols, p.p, &opts.filter);
            let filter_s = tf.elapsed().as_secs_f64();
            drop(f);
            let t = Instant::now();
            let (outcome, report) = solve_sparse_system(&sys.rows, p.cols, p.p, &opts);
            let solve_s = t.elapsed().as_secs_f64();
            let ok = matches!(&outcome, SparseSolveOutcome::Solved(x) if *x == sys.planted);
            o.insert("filter_s".into(), filter_s.into());
            o.insert("solve_s".into(), solve_s.into());
            o.insert("core_s_derived".into(), (solve_s - filter_s).max(0.0).into());
            o.insert("report".into(), serde_json::to_value(&report).unwrap());
            o.insert("verified".into(), ok.into());
        }
        "core-bw" | "core-bl" => {
            let solver = if arm == "core-bw" {
                SparseCoreSolver::BlockWiedemann
            } else {
                SparseCoreSolver::BlockLanczos
            };
            let tf = Instant::now();
            let (core, core_columns) = unmerged_core(p, &sys);
            let filter_s = tf.elapsed().as_secs_f64();
            let n = core_columns.len();
            let core_nnz: usize = core.iter().map(|r| r.entries.len()).sum();
            let opts = sparse_opts(p, solver, false);
            let t = Instant::now();
            let (outcome, report) = solve_sparse_system(&core, n, p.p, &opts);
            let solve_s = t.elapsed().as_secs_f64();
            let ok = matches!(&outcome, SparseSolveOutcome::Solved(x)
                if x.iter().zip(&core_columns).all(|(&v, &c)| v == sys.planted[c as usize]));
            o.insert("core_rows".into(), core.len().into());
            o.insert("core_cols".into(), n.into());
            o.insert("core_nonzeros".into(), core_nnz.into());
            o.insert("filter_s".into(), filter_s.into());
            o.insert("solve_s".into(), solve_s.into());
            o.insert("report".into(), serde_json::to_value(&report).unwrap());
            o.insert("verified".into(), ok.into());
        }
        other => panic!("unknown arm {other}"),
    }
    out
}

fn u32le(w: &mut impl Write, v: u32) {
    w.write_all(&v.to_le_bytes()).unwrap();
}

/// CADO `bwc` GF(p) input: per row `len`, then `len` × (`col`, `i32 coeff`),
/// all 32-bit little endian; rhs as ASCII (below).
fn export_core(p: &Params, dir: &str) -> serde_json::Value {
    std::fs::create_dir_all(dir).unwrap();
    let sys = generate(p);
    let (core, cols) = unmerged_core(p, &sys);
    assert_eq!(core.len(), cols.len(), "core must be square");
    let mut w = std::io::BufWriter::new(std::fs::File::create(format!("{dir}/matrix.bin")).unwrap());
    let half = p.p / 2;
    for r in &core {
        u32le(&mut w, r.entries.len() as u32);
        for &(c, v) in &r.entries {
            let s: i64 = if v > half { -((p.p - v) as i64) } else { v as i64 };
            assert!(s.abs() < (1 << 31), "coefficient too large for CADO int32");
            u32le(&mut w, c);
            w.write_all(&(s as i32).to_le_bytes()).unwrap();
        }
    }
    drop(w);
    // `bwc.pl` takes the right-hand side as ASCII: `rows nrhs p`, then one
    // decimal residue per row.
    let mut w = std::io::BufWriter::new(std::fs::File::create(format!("{dir}/rhs.txt")).unwrap());
    writeln!(w, "{} 1 {}", core.len(), p.p).unwrap();
    for r in &core {
        writeln!(w, "{}", r.rhs).unwrap();
    }
    drop(w);
    let nnz: usize = core.iter().map(|r| r.entries.len()).sum();
    serde_json::json!({"cols": p.cols, "weight": p.weight, "seed": p.seed,
        "core_dim": cols.len(), "core_nonzeros": nnz, "modulus": p.p.to_string(),
        "dir": dir})
}

/// Read a CADO solution: `core_dim` + `nrhs` residues of two LE words each
/// (`K.sols*` holds `x || y` with `(M || RHS)(x; y) = 0`, so `x / (−y)`
/// solves `M x = RHS`).  Some builds pad `x` to a multiple of the block
/// size; the trailing `y` words are always last.
fn check_cado(p: &Params, dir: &str, sol: &str) -> serde_json::Value {
    let sys = generate(p);
    let (core, cols) = unmerged_core(p, &sys);
    let n = cols.len();
    let mut bytes = Vec::new();
    std::fs::File::open(sol).unwrap().read_to_end(&mut bytes).unwrap();
    let words: Vec<u64> = bytes
        .chunks_exact(8)
        .map(|c| u64::from_le_bytes(c.try_into().unwrap()))
        .collect();
    let y = *words.last().unwrap() % p.p;
    let mut result = serde_json::json!({"dir": dir, "sol": sol, "core_dim": n,
        "words": words.len()});
    let o = result.as_object_mut().unwrap();
    if y == 0 || words.len() < n + 1 {
        o.insert("verified".into(), false.into());
        o.insert("reason".into(), "zero rhs coordinate or short file".into());
        return result;
    }
    // x / (−y)
    let inv = modinv(p.p - y, p.p);
    let x: Vec<u64> = words[..n].iter().map(|&v| mulmod(v % p.p, inv, p.p)).collect();
    let rows_ok = core.iter().all(|r| {
        let mut s = 0u64;
        for &(c, v) in &r.entries {
            s = addmod(s, mulmod(v, x[c as usize], p.p), p.p);
        }
        s == r.rhs
    });
    let planted_ok = x
        .iter()
        .zip(&cols)
        .all(|(&v, &c)| v == sys.planted[c as usize]);
    o.insert("rows_ok".into(), rows_ok.into());
    o.insert("verified".into(), (rows_ok && planted_ok).into());
    result
}

fn coeffs_name(p: &Params) -> &'static str {
    if p.wide {
        "wide"
    } else if p.positive {
        "positive"
    } else {
        "unit"
    }
}

// ── CADO's own filter: purge → merge-dl → replay-dl → bwc ─────────────

/// Relations in CADO's indexed format, `a,b:i1,…,ik` (hex), an index
/// repeated once per unit of its (nonnegative) exponent.  `a` is the row
/// number, so the purged file maps back to the original rows.
fn cado_export_rels(p: &Params, dir: &str) -> serde_json::Value {
    assert!(p.positive, "CADO relations need --coeffs positive");
    std::fs::create_dir_all(dir).unwrap();
    let sys = generate(p);
    let mut w = std::io::BufWriter::new(std::fs::File::create(format!("{dir}/rels")).unwrap());
    for (i, r) in sys.rows.iter().enumerate() {
        let mut ids = Vec::new();
        for &(c, v) in &r.entries {
            assert!(v <= 64, "positive coefficients stay small");
            for _ in 0..v {
                ids.push(format!("{c:x}"));
            }
        }
        writeln!(w, "{i:x},1:{}", ids.join(",")).unwrap();
    }
    serde_json::json!({"rows": sys.rows.len(), "cols": p.cols, "dir": dir})
}

/// Original row of each line of a purged file (its `a`, in hex).
fn purged_rows(path: &str) -> Vec<usize> {
    let text = std::fs::read_to_string(path).unwrap();
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .map(|l| usize::from_str_radix(l.split(',').next().unwrap(), 16).unwrap())
        .collect()
}

/// Right-hand side of each row of the replayed matrix: a relation set
/// `n id:e …` (ids index the purged file, exponents signed) sums
/// `e · rhs(row)`.  Written as `bwc.pl`'s ASCII rhs.
fn cado_rhs(p: &Params, dir: &str) -> serde_json::Value {
    let t = Instant::now();
    let sys = generate(p);
    let gen_s = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let purged = purged_rows(&format!("{dir}/purged"));
    let index = std::fs::read_to_string(format!("{dir}/small.index")).unwrap();
    let mut lines = index.lines();
    let n: usize = lines.next().unwrap().trim().parse().unwrap();
    let mut w = std::io::BufWriter::new(std::fs::File::create(format!("{dir}/rhs.txt")).unwrap());
    writeln!(w, "{n} 1 {}", p.p).unwrap();
    for line in lines.take(n) {
        let mut acc = 0u64;
        for item in line.split_whitespace().skip(1) {
            let (id, e) = item.split_once(':').unwrap_or((item, "1"));
            let id = usize::from_str_radix(id, 16).unwrap();
            let e: i64 = e.parse().unwrap();
            let e = if e < 0 { p.p - ((-e) as u64 % p.p) } else { e as u64 % p.p };
            acc = addmod(acc, mulmod(e, sys.rows[purged[id]].rhs, p.p), p.p);
        }
        writeln!(w, "{acc}").unwrap();
    }
    serde_json::json!({"small_rows": n, "rhs_s": t.elapsed().as_secs_f64(), "generate_s": gen_s})
}

/// Logs of the replayed matrix's columns from a `K.sols` file (column `i`
/// is ideal `h` per `small.ideal`), then iterative propagation through
/// every original row with exactly one unknown — what CADO's
/// `reconstructlog-dl` does, minus the NFS-specific parts.
fn cado_reconstruct(p: &Params, dir: &str, sol: &str) -> serde_json::Value {
    let sys = generate(p);
    let t = Instant::now();
    let ideals = std::fs::read_to_string(format!("{dir}/small.ideal")).unwrap();
    let col_of: Vec<usize> = ideals
        .lines()
        // A one-token line is the header (column count), when present.
        .filter(|l| !l.starts_with('#') && l.split_whitespace().count() == 2)
        .map(|l| {
            let mut it = l.split_whitespace();
            let _i: usize = it.next().unwrap().parse().unwrap();
            usize::from_str_radix(it.next().unwrap(), 16).unwrap()
        })
        .collect();
    let mut bytes = Vec::new();
    std::fs::File::open(sol).unwrap().read_to_end(&mut bytes).unwrap();
    let words: Vec<u64> = bytes
        .chunks_exact(8)
        .map(|c| u64::from_le_bytes(c.try_into().unwrap()))
        .collect();
    let k = col_of.len();
    let y = *words.last().unwrap() % p.p;
    let mut value: Vec<Option<u64>> = vec![None; p.cols];
    if y != 0 && words.len() > k {
        let inv = modinv(p.p - y, p.p);
        for (i, &c) in col_of.iter().enumerate() {
            value[c] = Some(mulmod(words[i] % p.p, inv, p.p));
        }
    }
    // Propagate: a row with one unknown determines it.
    let mut col_rows: Vec<Vec<u32>> = vec![Vec::new(); p.cols];
    for (i, r) in sys.rows.iter().enumerate() {
        for &(c, _) in &r.entries {
            col_rows[c as usize].push(i as u32);
        }
    }
    let mut queue: Vec<u32> = (0..sys.rows.len() as u32).collect();
    while let Some(i) = queue.pop() {
        let r = &sys.rows[i as usize];
        let unknown: Vec<&(u32, u64)> =
            r.entries.iter().filter(|(c, _)| value[*c as usize].is_none()).collect();
        if unknown.len() != 1 {
            continue;
        }
        let (c, v) = *unknown[0];
        let mut s = r.rhs;
        for &(c2, v2) in &r.entries {
            if let Some(x) = value[c2 as usize] {
                s = (s + p.p - mulmod(v2, x, p.p)) % p.p;
            }
        }
        value[c as usize] = Some(mulmod(s, modinv(v, p.p), p.p));
        queue.extend(col_rows[c as usize].iter().copied());
    }
    let reconstruct_s = t.elapsed().as_secs_f64();
    let known = value.iter().filter(|v| v.is_some()).count();
    let ok = value
        .iter()
        .zip(&sys.planted)
        .all(|(v, &x)| *v == Some(x));
    serde_json::json!({"small_cols": k, "determined": known, "cols": p.cols,
        "reconstruct_s": reconstruct_s, "verified": ok})
}

fn modinv(a: u64, m: u64) -> u64 {
    let (mut r0, mut r1) = (m as i128, a as i128);
    let (mut t0, mut t1) = (0i128, 1i128);
    while r1 != 0 {
        let q = r0 / r1;
        (r0, r1) = (r1, r0 - q * r1);
        (t0, t1) = (t1, t0 - q * t1);
    }
    ((t0 % m as i128 + m as i128) % m as i128) as u64
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let cmd = args.first().cloned().unwrap_or_default();
    let get = |k: &str| {
        args.iter()
            .position(|a| a == k)
            .and_then(|i| args.get(i + 1).cloned())
    };
    let num = |k: &str, d: u64| get(k).map(|v| v.parse::<u64>().unwrap()).unwrap_or(d);
    let p = Params {
        cols: num("--cols", 4096) as usize,
        weight: num("--weight", 3) as usize,
        excess: num("--excess", 32) as usize,
        seed: num("--seed", 1),
        p: num("--modulus", DEFAULT_P),
        wide: get("--coeffs").as_deref() == Some("wide"),
        positive: get("--coeffs").as_deref() == Some("positive"),
        block: num("--block", 4) as usize,
        merge_max_weight: get("--merge-max-weight").map(|v| v.parse().unwrap()),
        max_row_weight: get("--max-row-weight").map(|v| v.parse().unwrap()),
    };
    let v = match cmd.as_str() {
        "ours" => run_ours(&p, &get("--arm").expect("--arm")),
        "export-core" => export_core(&p, &get("--dir").expect("--dir")),
        "check-cado" => check_cado(
            &p,
            &get("--dir").expect("--dir"),
            &get("--sol").expect("--sol"),
        ),
        "cado-export-rels" => cado_export_rels(&p, &get("--dir").expect("--dir")),
        "cado-rhs" => cado_rhs(&p, &get("--dir").expect("--dir")),
        "cado-reconstruct" => cado_reconstruct(
            &p,
            &get("--dir").expect("--dir"),
            &get("--sol").expect("--sol"),
        ),
        _ => {
            eprintln!("usage: la_engine_bench ours|export-core|check-cado|cado-export-rels|cado-rhs|cado-reconstruct [options]");
            std::process::exit(2);
        }
    };
    println!("{}", serde_json::to_string(&v).unwrap());
}
