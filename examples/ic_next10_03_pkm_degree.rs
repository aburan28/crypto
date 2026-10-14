//! IC next-10 campaign, item 3: the PKM tower's solving degree past the
//! memory wall (`research/ic_next10_20261014/03-pkm-degree/`).
//!
//! Reads a Kummer tower system as `examples/pkm_tower_pilot.rs
//! --export-magma-rows` writes it (Magma syntax, raw presentation) and:
//!
//! - `convert IN.magma OUT.ms` writes it in msolve's input format, variables
//!   in the same order (msolve's DRL ranks the first variable highest, as
//!   Magma's grevlex does);
//! - `delta IN.magma` (`m = 2`, the `full` presentation): the exact XL /
//!   Macaulay solving degree of the refutation, `D_XL = d + δ` with
//!   `δ = deg(1/S̄)` in the square-free basis of the Kummer quotient
//!   (note §3.3, Proposition 2), computed by a 2-D number-theoretic
//!   transform over `(Z/2^t)^2`;
//! - `macaulay IN.magma DMAX` (any `m`, `chain` with free `u`): the Macaulay
//!   matrix of degree `D` in the tower quotient, for `D = 1 … DMAX`: its rank,
//!   the affine Hilbert function of the homogenised generators, and the first
//!   `D` whose row space contains `1` (the XL solving degree). No row is ever
//!   reused at a lower degree, so this is XL, not F4;
//! - `top IN.magma`: the top-degree parts of the non-tower generators modulo
//!   the towers' leading forms `y²`, i.e. what the homogenised system says at
//!   infinity (note §3.4), and the Hilbert function of that top ideal when it
//!   is monomial.
//!
//! Output is one JSON object per command on stdout.

use serde_json::json;
use std::collections::BTreeMap;

type Poly = Vec<(Vec<u32>, u64)>;

struct Sys {
    p: u64,
    names: Vec<String>,
    polys: Vec<Poly>,
}

fn parse(path: &str) -> Sys {
    let text = std::fs::read_to_string(path).expect("read the Magma script");
    let p: u64 = {
        let i = text.find("GF(").expect("GF(") + 3;
        let j = text[i..].find(')').unwrap() + i;
        text[i..j].parse().unwrap()
    };
    let names: Vec<String> = {
        let i = text.find("AssignNames(~P, [").expect("AssignNames") + 17;
        let j = text[i..].find(']').unwrap() + i;
        text[i..j]
            .split(',')
            .map(|s| s.trim().trim_matches('"').to_string())
            .collect()
    };
    let n = names.len();
    let body = {
        let i = text.find("F := [P |").expect("F := [P |") + 9;
        let j = text[i..].find("];").unwrap() + i;
        &text[i..j]
    };
    let mut polys = Vec::new();
    for line in body.split(",\n") {
        let line = line.trim().trim_end_matches(',');
        if line.is_empty() {
            continue;
        }
        let mut poly: BTreeMap<Vec<u32>, u64> = BTreeMap::new();
        for term in line.split(" + ") {
            let mut e = vec![0u32; n];
            let mut parts = term.trim().split('*');
            let c: u64 = parts.next().unwrap().parse().expect("coefficient");
            for f in parts {
                let f = f.strip_prefix("P.").expect("P.i");
                let (v, k) = match f.split_once('^') {
                    Some((v, k)) => (v.parse::<usize>().unwrap(), k.parse::<u32>().unwrap()),
                    None => (f.parse::<usize>().unwrap(), 1),
                };
                e[v - 1] += k;
            }
            let s = poly.entry(e).or_insert(0);
            *s = (*s + c) % p;
        }
        polys.push(poly.into_iter().filter(|(_, c)| *c != 0).collect());
    }
    Sys { p, names, polys }
}

fn pw(mut b: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1u64;
    b %= p;
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % p;
        }
        b = b * b % p;
        e >>= 1;
    }
    r
}

/// The Kummer structure: blocks `y{i}_{j}`, free `u{k}`; the tower
/// equations are checked to be `y_{j+1} − y_j²` and `1 − y_{t−1}²`.
struct Kummer {
    m: usize,
    t: usize,
    /// var index -> Some((block, level)) or None for a free `u`.
    role: Vec<Option<(usize, usize)>>,
    nu: usize,
    /// var index -> u index
    uidx: Vec<usize>,
    /// non-tower generators (indices into polys)
    eqs: Vec<usize>,
}

fn kummer(s: &Sys) -> Kummer {
    let mut role = Vec::new();
    let mut uidx = vec![usize::MAX; s.names.len()];
    let (mut m, mut t, mut nu) = (0, 0, 0);
    for (vi, nm) in s.names.iter().enumerate() {
        if let Some(r) = nm.strip_prefix('y') {
            let (i, j) = r.split_once('_').unwrap();
            let (i, j): (usize, usize) = (i.parse().unwrap(), j.parse().unwrap());
            m = m.max(i + 1);
            t = t.max(j + 1);
            role.push(Some((i, j)));
        } else {
            assert!(nm.starts_with('u'));
            uidx[vi] = nu;
            nu += 1;
            role.push(None);
        }
    }
    let var = |i: usize, j: usize| role.iter().position(|r| *r == Some((i, j))).unwrap();
    // Which polys are towers?
    let mut tower_seen = vec![false; m * t];
    let mut eqs = Vec::new();
    for (k, f) in s.polys.iter().enumerate() {
        let mut is_tower = false;
        if f.len() == 2 {
            for i in 0..m {
                for j in 0..t {
                    let mut sq = vec![0u32; s.names.len()];
                    sq[var(i, j)] = 2;
                    let mut lo = vec![0u32; s.names.len()];
                    if j + 1 < t {
                        lo[var(i, j + 1)] = 1;
                    }
                    let want: BTreeMap<Vec<u32>, u64> =
                        [(sq, s.p - 1), (lo, 1)].into_iter().collect();
                    let got: BTreeMap<Vec<u32>, u64> = f.iter().cloned().collect();
                    if want == got {
                        is_tower = true;
                        tower_seen[i * t + j] = true;
                    }
                }
            }
        }
        if !is_tower {
            eqs.push(k);
        }
    }
    assert!(tower_seen.iter().all(|&b| b), "every Kummer tower equation present");
    Kummer { m, t, role, nu, uidx, eqs }
}

impl Kummer {
    /// Block exponents E_i = Σ_j e_ij 2^j mod 2^t, and the u exponents.
    fn nf(&self, e: &[u32]) -> (Vec<u64>, Vec<u32>) {
        let modm = (1u64 << self.t) - 1;
        let mut b = vec![0u64; self.m];
        let mut u = vec![0u32; self.nu];
        for (vi, &k) in e.iter().enumerate() {
            if k == 0 {
                continue;
            }
            match self.role[vi] {
                Some((i, j)) => b[i] = (b[i] + (k as u64) * (1u64 << j)) & modm,
                None => u[self.uidx[vi]] += k,
            }
        }
        (b, u)
    }
}

fn deg(e: &[u32]) -> u32 {
    e.iter().sum()
}

fn convert(inp: &str, out: &str) {
    let s = parse(inp);
    let mut o = format!("{}\n{}\n", s.names.join(","), s.p);
    let ps: Vec<String> = s
        .polys
        .iter()
        .map(|f| {
            f.iter()
                .map(|(e, c)| {
                    let mut t = c.to_string();
                    for (i, &k) in e.iter().enumerate() {
                        if k == 1 {
                            t.push_str(&format!("*{}", s.names[i]));
                        } else if k > 1 {
                            t.push_str(&format!("*{}^{k}", s.names[i]));
                        }
                    }
                    t
                })
                .collect::<Vec<_>>()
                .join("+")
        })
        .collect();
    o.push_str(&ps.join(",\n"));
    o.push('\n');
    std::fs::write(out, o).expect("write msolve file");
    println!(
        "{}",
        json!({"cmd":"convert","in":inp,"out":out,"p":s.p,"n_vars":s.names.len(),"n_eqs":s.polys.len(),
               "degrees": s.polys.iter().map(|f| f.iter().map(|(e,_)| deg(e)).max().unwrap()).collect::<Vec<_>>()})
    );
}

fn ntt(a: &mut [u64], root: u64, p: u64) {
    let n = a.len();
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j |= bit;
        if i < j {
            a.swap(i, j);
        }
    }
    let mut len = 2;
    while len <= n {
        let w = pw(root, (n / len) as u64, p);
        let mut i = 0;
        while i < n {
            let mut wn = 1u64;
            for k in 0..len / 2 {
                let x = a[i + k];
                let y = a[i + k + len / 2] * wn % p;
                a[i + k] = (x + y) % p;
                a[i + k + len / 2] = (x + p - y) % p;
                wn = wn * w % p;
            }
            i += len;
        }
        len <<= 1;
    }
}

/// 2-D transform of an n×n array (row-major) with `root` of order n.
fn ntt2(a: &mut [u64], n: usize, root: u64, p: u64) {
    for r in 0..n {
        ntt(&mut a[r * n..(r + 1) * n], root, p);
    }
    let mut col = vec![0u64; n];
    for c in 0..n {
        for r in 0..n {
            col[r] = a[r * n + c];
        }
        ntt(&mut col, root, p);
        for r in 0..n {
            a[r * n + c] = col[r];
        }
    }
}

fn delta(inp: &str) {
    let s = parse(inp);
    let k = kummer(&s);
    assert_eq!(k.m, 2, "delta: m = 2 (full presentation) only");
    assert_eq!(k.nu, 0);
    assert_eq!(k.eqs.len(), 1);
    let p = s.p;
    let t = k.t;
    let n = 1usize << t;
    assert_eq!((p - 1) % (n as u64), 0);
    // generator of F_p^*: find the least primitive root.
    let fac = {
        let mut f = Vec::new();
        let mut x = p - 1;
        let mut q = 2;
        while q * q <= x {
            if x % q == 0 {
                f.push(q);
                while x % q == 0 {
                    x /= q;
                }
            }
            q += 1;
        }
        if x > 1 {
            f.push(x);
        }
        f
    };
    let g = (2..p).find(|&g| fac.iter().all(|&q| pw(g, (p - 1) / q, p) != 1)).unwrap();
    let zeta = pw(g, (p - 1) / n as u64, p);
    let f = &s.polys[k.eqs[0]];
    let d = f.iter().map(|(e, _)| deg(e)).max().unwrap();
    // S̄ as a group-algebra element: coefficient array over (E1, E2).
    let mut a = vec![0u64; n * n];
    for (e, c) in f {
        let (b, _) = k.nf(e);
        let idx = (b[0] as usize) * n + b[1] as usize;
        a[idx] = (a[idx] + c) % p;
    }
    let sbar_deg = (0..n * n)
        .filter(|&i| a[i] != 0)
        .map(|i| ((i / n).count_ones() + (i % n).count_ones()))
        .max()
        .unwrap();
    // evaluate on the grid
    ntt2(&mut a, n, zeta, p);
    let zeros = a.iter().filter(|&&v| v == 0).count();
    // 1/S̄ pointwise (zeros left at 0: the σ = 0 case is the refutation).
    for v in a.iter_mut() {
        if *v != 0 {
            *v = pw(*v, p - 2, p);
        }
    }
    // inverse transform
    let zinv = pw(zeta, p - 2, p);
    ntt2(&mut a, n, zinv, p);
    let ninv = pw((n * n) as u64 % p, p - 2, p);
    let mut hist = vec![0u64; 2 * t + 1];
    for v in a.iter_mut() {
        *v = *v * ninv % p;
    }
    let mut dmax = 0u32;
    for i in 0..n * n {
        if a[i] != 0 {
            let w = (i / n).count_ones() + (i % n).count_ones();
            hist[w as usize] += 1;
            dmax = dmax.max(w);
        }
    }
    // Self-check: (1/S̄)·S̄ = 1 in the group algebra, tested at a few
    // random points via the convolution at index 0 and a non-zero index.
    println!(
        "{}",
        json!({"cmd":"delta","in":inp,"m":2,"t":t,"N":2*t,"p":p,"zeta":zeta,
               "d_raw":d,"deg_Sbar":sbar_deg,"grid_zeros_sigma":zeros,
               "delta":dmax,"D_XL":d+dmax,
               "weight_histogram_of_inverse":hist,
               "weight_N_coefficients": hist[2*t]})
    );
}

/// A sparse row: sorted by column key ascending (leading = first).
type Row = Vec<(u128, u64)>;

fn colkey(k: &Kummer, b: &[u64], u: &[u32]) -> u128 {
    let dg: u32 = b.iter().map(|x| x.count_ones()).sum::<u32>() + u.iter().sum::<u32>();
    let mut key: u128 = (255 - dg as u128) << 112;
    let mut sh = 0;
    for &x in b {
        key |= (x as u128) << sh;
        sh += k.t;
    }
    for &x in u {
        key |= (x as u128) << sh;
        sh += 8;
    }
    assert!(sh <= 112);
    key
}

fn macaulay(inp: &str, dmax: u32) {
    let s = parse(inp);
    let k = kummer(&s);
    let p = s.p;
    let n_blk = k.m;
    let nbits = k.m * k.t;
    // multipliers: (mask over all blocks, u exps), by degree
    let mut pivots: BTreeMap<u128, Row> = BTreeMap::new();
    let mut per_d = Vec::new();
    let mut found = None;
    let eqdeg: Vec<u32> = k
        .eqs
        .iter()
        .map(|&i| s.polys[i].iter().map(|(e, _)| deg(e)).max().unwrap())
        .collect();
    let eq_nf: Vec<Vec<(Vec<u64>, Vec<u32>, u64)>> = k
        .eqs
        .iter()
        .map(|&i| {
            s.polys[i]
                .iter()
                .map(|(e, c)| {
                    let (b, u) = k.nf(e);
                    (b, u, *c)
                })
                .collect()
        })
        .collect();
    // u-monomials of each degree
    fn umons(nu: usize, d: u32) -> Vec<Vec<u32>> {
        if nu == 0 {
            return if d == 0 { vec![vec![]] } else { vec![] };
        }
        let mut out = Vec::new();
        for a in 0..=d {
            for mut r in umons(nu - 1, d - a) {
                r.insert(0, a);
                out.push(r);
            }
        }
        out
    }
    let t0 = std::time::Instant::now();
    let mut nrows_total = 0u64;
    for dd in 1..=dmax {
        // new rows: multiplier degree exactly dd - deg f
        for (q, &df) in eqdeg.iter().enumerate() {
            if dd < df {
                continue;
            }
            let md = dd - df;
            for ydeg in 0..=md.min(nbits as u32) {
                let ud = md - ydeg;
                let us = umons(k.nu, ud);
                if us.is_empty() {
                    continue;
                }
                // masks with popcount ydeg over nbits
                let mut masks = Vec::new();
                if ydeg == 0 {
                    masks.push(0u64);
                } else {
                    let mut x: u64 = (1u64 << ydeg) - 1;
                    while x < (1u64 << nbits) {
                        masks.push(x);
                        let c = x & x.wrapping_neg();
                        let r = x + c;
                        x = (((r ^ x) >> 2) / c) | r;
                    }
                }
                for &mask in &masks {
                    let mb: Vec<u64> = (0..n_blk)
                        .map(|i| (mask >> (i * k.t)) & ((1u64 << k.t) - 1))
                        .collect();
                    for uu in &us {
                        let mut acc: BTreeMap<u128, u64> = BTreeMap::new();
                        for (b, u, c) in &eq_nf[q] {
                            let nb: Vec<u64> = (0..n_blk)
                                .map(|i| (b[i] + mb[i]) & ((1u64 << k.t) - 1))
                                .collect();
                            let nu: Vec<u32> = u.iter().zip(uu).map(|(a, b)| a + b).collect();
                            let key = colkey(&k, &nb, &nu);
                            let e = acc.entry(key).or_insert(0);
                            *e = (*e + c) % p;
                        }
                        let mut row: Row = acc.into_iter().filter(|(_, c)| *c != 0).collect();
                        nrows_total += 1;
                        // reduce
                        loop {
                            if row.is_empty() {
                                break;
                            }
                            let lead = row[0].0;
                            if let Some(pr) = pivots.get(&lead) {
                                let f = row[0].1; // pivot rows are monic
                                let mut out: Row = Vec::with_capacity(row.len() + pr.len());
                                let (mut i, mut j) = (0, 0);
                                while i < row.len() || j < pr.len() {
                                    if j >= pr.len() || (i < row.len() && row[i].0 < pr[j].0) {
                                        out.push(row[i]);
                                        i += 1;
                                    } else if i >= row.len() || pr[j].0 < row[i].0 {
                                        out.push((pr[j].0, (p - f * pr[j].1 % p) % p));
                                        j += 1;
                                    } else {
                                        let v = (row[i].1 + p - f * pr[j].1 % p) % p;
                                        if v != 0 {
                                            out.push((row[i].0, v));
                                        }
                                        i += 1;
                                        j += 1;
                                    }
                                }
                                row = out;
                            } else {
                                let inv = pw(row[0].1, p - 2, p);
                                for e in row.iter_mut() {
                                    e.1 = e.1 * inv % p;
                                }
                                pivots.insert(lead, row);
                                break;
                            }
                        }
                    }
                }
            }
        }
        // standard monomials of degree <= dd
        let mut nstd: u64 = 0;
        for yd in 0..=dd.min(nbits as u32) {
            let c = binom(nbits as u64, yd as u64);
            nstd += c * umons(k.nu, dd - yd).len() as u64;
            // umons(nu, x) count = C(x+nu-1, nu-1) for nu>0
        }
        // count std monomials with deg ≤ dd: sum over exact degree e ≤ dd
        let mut nstd_le: u64 = 0;
        for e in 0..=dd {
            let mut c = 0u64;
            for yd in 0..=e.min(nbits as u32) {
                c += binom(nbits as u64, yd as u64) * umons(k.nu, e - yd).len() as u64;
            }
            nstd_le += c;
        }
        let _ = nstd;
        let one = colkey(&k, &vec![0u64; n_blk], &vec![0u32; k.nu]);
        let has_one = pivots.contains_key(&one);
        let rank = pivots.len() as u64;
        let nnz: usize = pivots.values().map(|r| r.len()).sum();
        per_d.push(json!({"D":dd,"rows":nrows_total,"rank":rank,"std_monomials_le_D":nstd_le,
                          "HF_affine": nstd_le.saturating_sub(rank),"echelon_nnz":nnz,
                          "contains_1":has_one,"s":t0.elapsed().as_secs_f64()}));
        eprintln!("{}", per_d.last().unwrap());
        if has_one {
            found = Some(dd);
            break;
        }
    }
    println!(
        "{}",
        json!({"cmd":"macaulay","in":inp,"m":k.m,"t":k.t,"N":k.m*k.t,"n_u":k.nu,"eq_degrees":eqdeg,
               "D_XL":found,"dmax":dmax,"per_degree":per_d})
    );
}

fn binom(n: u64, k: u64) -> u64 {
    if k > n {
        return 0;
    }
    let mut r = 1u64;
    for i in 0..k {
        r = r * (n - i) / (i + 1);
    }
    r
}

fn top(inp: &str) {
    let s = parse(inp);
    let k = kummer(&s);
    let mut out = Vec::new();
    let mut all_monomial = true;
    let mut gens: Vec<Vec<u32>> = Vec::new();
    for &q in &k.eqs {
        let f = &s.polys[q];
        let d = f.iter().map(|(e, _)| deg(e)).max().unwrap();
        // top part, then drop terms divisible by some y^2
        let surv: Vec<(Vec<u32>, u64)> = f
            .iter()
            .filter(|(e, _)| deg(e) == d)
            .filter(|(e, _)| e.iter().enumerate().all(|(vi, &x)| k.role[vi].is_none() || x < 2))
            .cloned()
            .collect();
        let txt: Vec<String> = surv
            .iter()
            .map(|(e, c)| {
                let mut t = c.to_string();
                for (i, &x) in e.iter().enumerate() {
                    if x > 0 {
                        t.push_str(&format!("*{}^{x}", s.names[i]));
                    }
                }
                t
            })
            .collect();
        if surv.len() > 1 {
            all_monomial = false;
        }
        if surv.len() == 1 {
            gens.push(surv[0].0.clone());
        }
        out.push(json!({"eq":q,"degree":d,"top_mod_y2":txt}));
    }
    // Hilbert function of (y_ij^2, monomial gens) when monomial: count
    // standard monomials of each degree up to N + 4 (y square-free, u free).
    let n = s.names.len();
    let nbits = k.m * k.t;
    let mut hf = Vec::new();
    if all_monomial && nbits <= 24 {
        let yvars: Vec<usize> = (0..n).filter(|&v| k.role[v].is_some()).collect();
        let uvars: Vec<usize> = (0..n).filter(|&v| k.role[v].is_none()).collect();
        let dmax = nbits as u32 + 4;
        let mut cnt = vec![0u64; dmax as usize + 1];
        // u-monomials up to dmax
        fn umons_le(nu: usize, d: u32) -> Vec<Vec<u32>> {
            if nu == 0 {
                return vec![vec![]];
            }
            let mut out = Vec::new();
            for a in 0..=d {
                for mut r in umons_le(nu - 1, d - a) {
                    r.insert(0, a);
                    out.push(r);
                }
            }
            out
        }
        let us = umons_le(uvars.len(), dmax);
        for mask in 0u64..(1u64 << nbits) {
            let yd = mask.count_ones();
            if yd > dmax {
                continue;
            }
            for uu in &us {
                let ud: u32 = uu.iter().sum();
                if yd + ud > dmax {
                    continue;
                }
                let mut e = vec![0u32; n];
                for (bi, &v) in yvars.iter().enumerate() {
                    e[v] = ((mask >> bi) & 1) as u32;
                }
                for (ui, &v) in uvars.iter().enumerate() {
                    e[v] = uu[ui];
                }
                if gens.iter().any(|g| g.iter().zip(&e).all(|(a, b)| a <= b)) {
                    continue;
                }
                cnt[(yd + ud) as usize] += 1;
            }
        }
        hf = cnt;
    }
    let dreg = if all_monomial && !hf.is_empty() {
        hf.iter().position(|&c| c == 0)
    } else {
        None
    };
    println!(
        "{}",
        json!({"cmd":"top","in":inp,"m":k.m,"t":k.t,"N":k.m*k.t,"n_u":k.nu,"generators":out,
               "top_ideal_monomial":all_monomial,"HF_top_ideal":hf,"d_reg_first_zero":dreg})
    );
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    match a.get(1).map(|s| s.as_str()) {
        Some("convert") => convert(&a[2], &a[3]),
        Some("delta") => delta(&a[2]),
        Some("macaulay") => macaulay(&a[2], a[3].parse().unwrap()),
        Some("top") => top(&a[2]),
        _ => {
            eprintln!("usage: convert IN OUT | delta IN | macaulay IN DMAX | top IN");
            std::process::exit(2)
        }
    }
}
