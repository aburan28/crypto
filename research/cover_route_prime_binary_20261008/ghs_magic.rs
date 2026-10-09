//! GHS magic number for the two RFC 2409 EC2N curves (standalone; `rustc -O`).
//! `m = dim_F2 Span{(1, sigma^i(gamma))}` with gamma = sqrt(b), sigma the q-Frobenius.
//! Field elements are bit vectors over F_2 modulo the RFC's trinomial.

#[derive(Clone)]
struct F2n { n: usize, red: Vec<u64> } // red: bits of the reduction polynomial (degree n)

type El = Vec<u64>;

impl F2n {
    fn new(n: usize, mid: usize) -> Self { F2n::from_bits(n, &[0usize, mid, n]) }
    fn from_bits(n: usize, bits: &[usize]) -> Self {
        let w = (n + 64) / 64 + 1;
        let mut red = vec![0u64; w];
        for &k in bits { red[k / 64] |= 1u64 << (k % 64); }
        F2n { n, red }
    }
    fn from_hex(n: usize, hex: &str) -> (Self, Vec<usize>) {
        let h = hex.trim_start_matches("0x");
        let mut bits = Vec::new();
        for (i, ch) in h.chars().rev().enumerate() {
            let v = ch.to_digit(16).unwrap();
            for j in 0..4 { if (v >> j) & 1 == 1 { bits.push(4 * i + j); } }
        }
        assert_eq!(*bits.iter().max().unwrap(), n, "modulus degree");
        (F2n::from_bits(n, &bits), bits)
    }
    fn el_from_hex(&self, hex: &str) -> El {
        let h = hex.trim_start_matches("0x");
        let mut e = self.zero();
        for (i, ch) in h.chars().rev().enumerate() {
            let v = ch.to_digit(16).unwrap() as u64;
            for j in 0..4 { if (v >> j) & 1 == 1 { let k = 4 * i + j; e[k / 64] |= 1u64 << (k % 64); } }
        }
        e
    }
    fn words(&self) -> usize { (self.n + 63) / 64 }
    fn zero(&self) -> El { vec![0; self.words()] }
    fn from_u64(&self, v: u64) -> El { let mut e = self.zero(); e[0] = v; e }
    fn bit(e: &El, i: usize) -> bool { i / 64 < e.len() && (e[i / 64] >> (i % 64)) & 1 == 1 }
    fn mul(&self, a: &El, b: &El) -> El {
        let w = self.words();
        let mut t = vec![0u64; 2 * w];
        for i in 0..self.n {
            if F2n::bit(a, i) {
                // t ^= b << i
                let (ws, bs) = (i / 64, i % 64);
                for j in 0..w {
                    t[j + ws] ^= b[j] << bs;
                    if bs > 0 && j + ws + 1 < 2 * w { t[j + ws + 1] ^= b[j] >> (64 - bs); }
                }
            }
        }
        // reduce degrees 2n-2 .. n
        for d in (self.n..2 * self.n - 1).rev() {
            if (t[d / 64] >> (d % 64)) & 1 == 1 {
                let s = d - self.n;
                for k in 0..self.red.len() * 64 {
                    if (self.red[k / 64] >> (k % 64)) & 1 == 1 {
                        let p = k + s;
                        t[p / 64] ^= 1u64 << (p % 64);
                    }
                }
            }
        }
        let mut r = t[..w].to_vec();
        // mask high bits
        let hb = self.n % 64;
        if hb > 0 { r[w - 1] &= (1u64 << hb) - 1; }
        r
    }
    fn sq(&self, a: &El) -> El { self.mul(a, a) }
    fn pow2k(&self, a: &El, k: usize) -> El { let mut r = a.clone(); for _ in 0..k { r = self.sq(&r); } r }
    fn sqrt(&self, a: &El) -> El { self.pow2k(a, self.n - 1) }
}

fn rank_f2(rows: &mut Vec<Vec<u64>>, nbits: usize) -> usize {
    let mut rank = 0;
    for bit in 0..nbits {
        let (wi, bi) = (bit / 64, bit % 64);
        let piv = (rank..rows.len()).find(|&r| (rows[r][wi] >> bi) & 1 == 1);
        if let Some(p) = piv {
            rows.swap(rank, p);
            let pr = rows[rank].clone();
            for r in 0..rows.len() {
                if r != rank && (rows[r][wi] >> bi) & 1 == 1 {
                    for j in 0..pr.len() { rows[r][j] ^= pr[j]; }
                }
            }
            rank += 1;
        }
    }
    rank
}

fn main() {
    let mut curves: Vec<(String, F2n, El)> = Vec::new();
    for (name, n, mid, b) in [("RFC 2409 Oakley group 3", 155usize, 62usize, 0x7338fu64), ("RFC 2409 Oakley group 4", 185, 69, 0x1ee9)] {
        let f = F2n::new(n, mid); let bb = f.from_u64(b); curves.push((name.to_string(), f, bb));
    }
    if let Some(path) = std::env::args().nth(1) {
        for line in std::fs::read_to_string(path).unwrap().lines() {
            let t: Vec<&str> = line.split_whitespace().collect();
            if t.len() < 4 { continue; }
            let n: usize = t[1].parse().unwrap();
            let (f, _) = F2n::from_hex(n, t[2]);
            let bb = f.el_from_hex(t[3]);
            curves.push((t[0].to_string(), f, bb));
        }
    }
    for (name, f, bb) in curves {
        let n = f.n;
        let b = &bb;
        let gamma = f.sqrt(&bb);
        assert_eq!(f.sq(&gamma), bb, "sqrt check");
        // smallest proper subfield containing b
        let mut sub = "none";
        let mut subdeg = n;
        for l in 1..n { if n % l == 0 && f.pow2k(&bb, l) == bb { subdeg = l; sub = "yes"; break; } }
        println!("{name}: n = {n}, b = {}..., b in proper subfield: {sub} (degree {subdeg})", b[0] & 0xffff);
        for l in std::iter::once(1usize).chain((2..n).filter(|l| n % l == 0)) {
            let np = n / l;
            // rows (1, sigma^i gamma), i < np, as bit vectors of length n+1
            let mut rows: Vec<Vec<u64>> = Vec::new();
            let mut cur = gamma.clone();
            for _ in 0..np {
                let mut v = vec![0u64; (n + 1 + 63) / 64 + 1];
                v[0] |= 1; // the "1" coordinate at bit 0
                for i in 0..n { if F2n::bit(&cur, i) { v[(i + 1) / 64] |= 1u64 << ((i + 1) % 64); } }
                rows.push(v);
                cur = f.pow2k(&cur, l);
            }
            let m = rank_f2(&mut rows, n + 1);
            if m >= 1 && m - 1 < 64 {
                let genus_hi = 1u128 << (m - 1);
                println!("  base F_2^{l:<3} (n' = {np:<3}): magic number m = {m:<3} genus {genus_hi} or {}", genus_hi - 1);
            } else {
                println!("  base F_2^{l:<3} (n' = {np:<3}): magic number m = {m:<3} genus 2^{} or 2^{} - 1", m - 1, m - 1);
            }
        }
    }
}
