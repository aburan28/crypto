//! Dense, bit-packed Boolean Macaulay elimination: the refutation degree of an affine
//! Boolean system at one degree `D`, from the Singular ideal files that
//! `rr_degree_ladder --dump-dir` writes.
//!
//!     macaulay_dense --file cell-dN-arm.sing --degree D [--threads 4] [--mem-gb 12]
//!
//! The degree-`D` Macaulay matrix of the system is the set of products `m·f` reduced
//! modulo `x_i² = x_i`, over every generator `f` and every multilinear monomial `m` with
//! `deg m ≤ D − deg f`; its columns are the multilinear monomials of degree `≤ D`.  The
//! system is refuted at `D` exactly when the constant `1` lies in the row space.  Rows are
//! generated in batches and reduced against a growing echelon basis whose pivot is the
//! lowest set column; the constant monomial has the highest column index, so a row that
//! reduces to that single bit is the polynomial `1`.  Batches are reduced in parallel
//! against the basis that exists when the batch starts, then serially among themselves.
//!
//! Output: one JSON line `{"file","degree","n_vars","n_cols","n_rows","rank","refuted",
//! "secs"}`.  Exit 0 whether or not refuted; the memory cap (`--mem-gb`, default 12) is a
//! refusal before any allocation, never a kill.

use rayon::prelude::*;
use std::time::Instant;

fn parse_sing(path: &str) -> (usize, Vec<Vec<u64>>) {
    let s = std::fs::read_to_string(path).expect("read .sing");
    let n_vars: usize = s
        .split("x(1..")
        .nth(1)
        .and_then(|t| t.split(')').next())
        .expect("ring header")
        .parse()
        .expect("n_vars");
    let body = s
        .split("ideal I =")
        .nth(1)
        .expect("ideal")
        .trim()
        .trim_end_matches(';');
    let polys = body
        .split(',')
        .map(|p| {
            let mut terms: Vec<u64> = p
                .split('+')
                .map(|mono| {
                    let mono = mono.trim();
                    if mono == "1" {
                        0u64
                    } else {
                        mono.split('*').fold(0u64, |acc, v| {
                            let i: usize = v
                                .trim()
                                .trim_start_matches("x(")
                                .trim_end_matches(')')
                                .parse()
                                .expect("var");
                            acc | 1u64 << (i - 1)
                        })
                    }
                })
                .collect();
            terms.sort_unstable();
            terms.dedup(); // the exporter never repeats a monomial; keep the invariant explicit
            terms
        })
        .filter(|t| !t.is_empty())
        .collect();
    (n_vars, polys)
}

/// Multilinear monomials of degree ≤ d over `v` variables, as masks, grouped by degree.
fn monomials_upto(v: usize, d: usize) -> Vec<u64> {
    let mut out = Vec::new();
    for mask in 0u64..(1u64 << v) {
        if (mask.count_ones() as usize) <= d {
            out.push(mask);
        }
    }
    out
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |name: &str| {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1))
            .cloned()
    };
    let file = get("--file").expect("--file");
    let degree: usize = get("--degree").expect("--degree").parse().expect("degree");
    let threads: usize = get("--threads").map_or(4, |t| t.parse().expect("threads"));
    let mem_gb: f64 = get("--mem-gb").map_or(12.0, |t| t.parse().expect("mem-gb"));
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .expect("pool");
    let started = Instant::now();

    let (v, polys) = parse_sing(&file);
    assert!(v <= 24, "the column table is indexed by mask; v ≤ 24");
    // columns: monomials of degree ≤ D, the constant one last
    let mut cols = monomials_upto(v, degree);
    cols.retain(|&m| m != 0);
    cols.push(0);
    let n_cols = cols.len();
    let words = n_cols.div_ceil(64);
    let mut col_index = vec![u32::MAX; 1usize << v];
    for (i, &m) in cols.iter().enumerate() {
        col_index[m as usize] = i as u32;
    }
    let const_col = n_cols - 1;
    // rows: (f, m) with deg m ≤ D − deg f
    let degs: Vec<usize> = polys
        .iter()
        .map(|p| p.iter().map(|t| t.count_ones() as usize).max().unwrap_or(0))
        .collect();
    let mut jobs: Vec<(usize, u64)> = Vec::new();
    for (fi, &df) in degs.iter().enumerate() {
        if df > degree {
            continue;
        }
        for &m in monomials_upto(v, degree - df).iter() {
            jobs.push((fi, m));
        }
    }
    let n_rows = jobs.len();
    let basis_bytes = (n_cols as f64) * (words as f64) * 8.0;
    if basis_bytes > mem_gb * 1e9 {
        println!(
            r#"{{"file":"{file}","degree":{degree},"n_vars":{v},"n_cols":{n_cols},"n_rows":{n_rows},"refused":"basis would need {:.1} GB"}}"#,
            basis_bytes / 1e9
        );
        return;
    }
    eprintln!(
        "{file}: D={degree} v={v} cols={n_cols} rows={n_rows} words={words} basis≤{:.1} GB",
        basis_bytes / 1e9
    );

    // basis[pivot] = row with lowest set bit at `pivot`
    let mut basis: Vec<Option<Box<[u64]>>> = (0..n_cols).map(|_| None).collect();
    let mut rank = 0usize;
    let mut refuted = false;
    let lowest = |row: &[u64]| -> Option<usize> {
        row.iter()
            .enumerate()
            .find(|(_, w)| **w != 0)
            .map(|(i, w)| i * 64 + w.trailing_zeros() as usize)
    };
    let batch = 2048usize;
    let mut pos = 0usize;
    while pos < n_rows && !refuted {
        let end = (pos + batch).min(n_rows);
        // build and reduce the batch in parallel against the current basis
        let basis_ref = &basis;
        let col_index_ref = &col_index;
        let polys_ref = &polys;
        let mut reduced: Vec<Vec<u64>> = jobs[pos..end]
            .par_iter()
            .map(|&(fi, m)| {
                let mut row = vec![0u64; words];
                for &t in &polys_ref[fi] {
                    let c = col_index_ref[(t | m) as usize] as usize;
                    row[c / 64] ^= 1u64 << (c % 64);
                }
                loop {
                    match lowest(&row) {
                        None => break,
                        Some(p) => match &basis_ref[p] {
                            Some(b) => {
                                for (r, bw) in row.iter_mut().zip(b.iter()) {
                                    *r ^= *bw;
                                }
                            }
                            None => break,
                        },
                    }
                }
                row
            })
            .filter(|row| row.iter().any(|w| *w != 0))
            .collect();
        // serial pass: reduce against pivots inserted by earlier rows of this batch
        for row in reduced.iter_mut() {
            loop {
                match lowest(row) {
                    None => break,
                    Some(p) => {
                        if p == const_col {
                            refuted = true;
                            break;
                        }
                        match &basis[p] {
                            Some(b) => {
                                for (r, bw) in row.iter_mut().zip(b.iter()) {
                                    *r ^= *bw;
                                }
                            }
                            None => {
                                basis[p] = Some(row.clone().into_boxed_slice());
                                rank += 1;
                                break;
                            }
                        }
                    }
                }
            }
            if refuted {
                break;
            }
        }
        pos = end;
        if pos.is_multiple_of(batch * 32) {
            eprintln!(
                "  rows {pos}/{n_rows}, rank {rank}, {:.0} s",
                started.elapsed().as_secs_f64()
            );
        }
    }
    println!(
        r#"{{"file":"{file}","degree":{degree},"n_vars":{v},"n_cols":{n_cols},"n_rows":{n_rows},"rank":{rank},"refuted":{refuted},"secs":{:.3}}}"#,
        started.elapsed().as_secs_f64()
    );
}
