//! IC next-10 item 8, arm A2: exhaustive search for linearised quasi-subfield
//! polynomials with coefficients in `F_{2^n}` (not only `F_2`) at prime `n`.
//!
//! `L(X) = X^{2^{n0}} + c_j X^{2^j} + ... + c_0 X` splits completely in
//! `F_{2^n}` iff the `F_2`-linear map `x -> L(x)` has kernel of dimension `n0`.
//! Scaling `X -> gamma X` multiplies `c_j` by `gamma^{2^j - 2^{n0}}`, which runs
//! over all of `F_{2^n}^*` when `gcd(2^{n0 - j} - 1, 2^n - 1) = 1`, i.e. when
//! `n` does not divide `n0 - j` (always, for prime `n` and `0 <= j < n0 < n`).
//! So `c_j = 1` loses nothing, and the search runs over `c_0 .. c_{j-1}` only:
//! `2^{j n}` candidates per `(n, n0, j)`.
//!
//! ```sh
//! cargo run --release --example ic_next10_08_qsf_general -- <j_max> <n> [<n> ...]
//! ```
//! Prints one CSV row per `(n, n0, j)`.
use crypto_lib::cryptanalysis::koblitz_index_calculus::find_irreducible;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;

fn kernel_dim(cols: &mut [u64], n: u32) -> u32 {
    // rank of the n x n matrix whose columns are bitmasks
    let mut rank = 0u32;
    let len = cols.len();
    for bit in 0..n {
        let mask = 1u64 << bit;
        let mut piv = None;
        for (c, item) in cols.iter().enumerate().skip(rank as usize) {
            if item & mask != 0 {
                piv = Some(c);
                break;
            }
        }
        if let Some(p) = piv {
            cols.swap(rank as usize, p);
            let pv = cols[rank as usize];
            for c in (rank as usize + 1)..len {
                if cols[c] & mask != 0 {
                    cols[c] ^= pv;
                }
            }
            rank += 1;
        }
    }
    n - rank
}

fn main() {
    let args: Vec<u32> = std::env::args().skip(1).map(|a| a.parse().unwrap()).collect();
    let j_max = args[0];
    println!("n,n0,j,examined,found,beta_if_found,witness_lambda_hex");
    for &n in &args[1..] {
        let irr = find_irreducible(n).expect("irreducible");
        let gf = Gf2::new(&irr);
        for n0 in 2..n {
            // precompute (z^k)^{2^i} for i <= n0
            let pw: Vec<Vec<u64>> = (0..n)
                .map(|k| (0..=n0).map(|i| gf.sqr_k(1u64 << k, i)).collect())
                .collect();
            for j in 0..=j_max.min(n0 - 1) {
                let examined: u64 = 1u64 << (j * n);
                let mut found = 0u64;
                let mut witness = String::new();
                let mask = (1u64 << n) - 1;
                let mut cols = vec![0u64; n as usize];
                for code in 0..examined {
                    for k in 0..n as usize {
                        let mut acc = pw[k][n0 as usize] ^ pw[k][j as usize]; // c_j = 1
                        for i in 0..j as usize {
                            let c = (code >> (i as u32 * n)) & mask;
                            if c != 0 {
                                acc ^= gf.mul(c, pw[k][i]);
                            }
                        }
                        cols[k] = acc;
                    }
                    if kernel_dim(&mut cols, n) == n0 {
                        found += 1;
                        if witness.is_empty() {
                            let cs: Vec<String> = (0..j)
                                .map(|i| format!("{:x}", (code >> (i * n)) & mask))
                                .collect();
                            witness = format!("c0..c{}=[{}];c{}=1", j, cs.join(" "), j);
                        }
                    }
                }
                let beta = (j as f64) * (n as f64) / ((n0 * n0) as f64);
                println!(
                    "{n},{n0},{j},{examined},{found},{},{}",
                    if found > 0 { format!("{beta:.4}") } else { String::new() },
                    witness
                );
            }
        }
    }
}
