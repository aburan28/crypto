//! GF(2) rank and span membership for packed bit matrices.
//!
//! The linear-algebra kernel of `research/ic_tree_split_cutout_20260930/cutout.py`.
//! Reads from stdin three little-endian `u64` values `s k w`, then `s` rows and
//! `k` sample rows of `w` words each.  Brings the `s` rows to row echelon form
//! and counts the samples that lie in their row space.
//!
//! Output: one JSON line `{"rank": r, "in_span": c}`.
//!
//!     cutout.py ... | gf2_span        (spawned by cutout.py, one call per arm and degree)
use std::io::{Read, Write};

fn main() {
    let mut buf = Vec::new();
    std::io::stdin().read_to_end(&mut buf).expect("read stdin");
    let words: Vec<u64> = buf
        .chunks_exact(8)
        .map(|c| u64::from_le_bytes(c.try_into().unwrap()))
        .collect();
    let (s, k, w) = (words[0] as usize, words[1] as usize, words[2] as usize);
    assert_eq!(words.len(), 3 + (s + k) * w, "length does not match header");
    let mut rows: Vec<Vec<u64>> = (0..s).map(|i| words[3 + i * w..3 + (i + 1) * w].to_vec()).collect();
    let samples = &words[3 + s * w..];

    // Row echelon form: row r has its first set bit at pivots[r], and every
    // later row is zero there.
    let mut pivots = Vec::new();
    let mut r = 0;
    'cols: for col in 0..w * 64 {
        if r == s {
            break;
        }
        let (cw, bit) = (col / 64, 1u64 << (col % 64));
        let p = match (r..s).find(|&i| rows[i][cw] & bit != 0) {
            Some(p) => p,
            None => continue 'cols,
        };
        rows.swap(r, p);
        let (head, tail) = rows.split_at_mut(r + 1);
        let pivot = &head[r];
        for row in tail.iter_mut() {
            if row[cw] & bit != 0 {
                for j in cw..w {
                    row[j] ^= pivot[j];
                }
            }
        }
        pivots.push(col);
        r += 1;
    }

    let mut in_span = 0usize;
    for i in 0..k {
        let mut v = samples[i * w..(i + 1) * w].to_vec();
        for (idx, &col) in pivots.iter().enumerate() {
            let (cw, bit) = (col / 64, 1u64 << (col % 64));
            if v[cw] & bit != 0 {
                for j in cw..w {
                    v[j] ^= rows[idx][j];
                }
            }
        }
        in_span += v.iter().all(|&x| x == 0) as usize;
    }
    let mut out = std::io::stdout();
    writeln!(out, "{{\"rank\": {}, \"in_span\": {}}}", pivots.len(), in_span).unwrap();
}
