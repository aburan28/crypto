// 09-cm-norm: norm-form / CM-ideal factor bases. Exact checks on pair 46 + toy coverage measurement.
use std::fmt::Write as _;

fn mulm(a: u64, b: u64, m: u64) -> u64 { ((a as u128 * b as u128) % m as u128) as u64 }
fn powm(mut a: u64, mut e: u64, m: u64) -> u64 { let mut r = 1 % m; a %= m; while e > 0 { if e & 1 == 1 { r = mulm(r, a, m); } a = mulm(a, a, m); e >>= 1; } r }
fn is_prime(n: u64) -> bool {
    if n < 2 { return false; }
    for &q in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] { if n % q == 0 { return n == q; } }
    let (mut d, mut s) = (n - 1, 0); while d % 2 == 0 { d /= 2; s += 1; }
    'w: for &a in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powm(a, d, n); if x == 1 || x == n - 1 { continue; }
        for _ in 1..s { x = mulm(x, x, n); if x == n - 1 { continue 'w; } }
        return false;
    } true
}
fn inv(a: u64, m: u64) -> u64 { powm(a, m - 2, m) }
fn sqrt_mod(a: u64, p: u64) -> Option<u64> { // Tonelli-Shanks
    let a = a % p; if a == 0 { return Some(0); }
    if powm(a, (p - 1) / 2, p) != 1 { return None; }
    let (mut q, mut s) = (p - 1, 0u32); while q % 2 == 0 { q /= 2; s += 1; }
    let mut z = 2; while powm(z, (p - 1) / 2, p) != p - 1 { z += 1; }
    let (mut m, mut c, mut t, mut r) = (s, powm(z, q, p), powm(a, q, p), powm(a, (q + 1) / 2, p));
    while t != 1 { let mut i = 0; let mut tt = t; while tt != 1 { tt = mulm(tt, tt, p); i += 1; }
        let b = powm(c, 1u64 << (m - i - 1), p); m = i; c = mulm(b, b, p); t = mulm(t, c, p); r = mulm(r, b, p); }
    Some(r)
}
fn isqrt(n: u64) -> u64 { let mut x = (n as f64).sqrt() as u64; while x * x > n { x -= 1; } while (x + 1) * (x + 1) <= n { x += 1; } x }

#[derive(Clone, Copy, PartialEq)] enum Pt { Inf, A(u64, u64) }
struct Curve { p: u64, a: u64, b: u64 }
impl Curve {
    fn add(&self, P: Pt, Q: Pt) -> Pt { let p = self.p; match (P, Q) {
        (Pt::Inf, _) => Q, (_, Pt::Inf) => P,
        (Pt::A(x1, y1), Pt::A(x2, y2)) => {
            let l = if x1 == x2 { if (y1 + y2) % p == 0 { return Pt::Inf; }
                mulm((3 * mulm(x1, x1, p) + self.a) % p, inv(2 * y1 % p, p), p)
            } else { mulm((y2 + p - y1) % p, inv((x2 + p - x1) % p, p), p) };
            let x3 = (mulm(l, l, p) + 2 * p - x1 - x2) % p; let y3 = (mulm(l, (x1 + p - x3) % p, p) + p - y1) % p; Pt::A(x3, y3) } } }
    fn mul(&self, mut k: u64, P: Pt) -> Pt { let mut r = Pt::Inf; let mut q = P; while k > 0 { if k & 1 == 1 { r = self.add(r, q); } q = self.add(q, q); k >>= 1; } r }
    fn lift(&self, x: u64) -> Option<Pt> { let p = self.p; let r = (mulm(mulm(x, x, p), x, p) + mulm(self.a, x, p) + self.b) % p; sqrt_mod(r, p).map(|y| Pt::A(x, y)) }
}

struct Toy { p: u64, t: i64, u: u64, n: u64, a: u64, b: u64 }
fn find_toys(start: u64, count: usize) -> Vec<Toy> {
    let mut out = vec![]; let mut p = start;
    while out.len() < count { p += 1; if !is_prime(p) { continue; }
        let mut found = None; let mut u = 1u64;
        while 19 * u * u <= 4 * p { let r = 4 * p - 19 * u * u; let s = isqrt(r); if s * s == r { found = Some((s, u)); break; } u += 1; }
        let Some((t, u)) = found else { continue };
        for &tt in &[t as i64, -(t as i64)] { let n = (p as i64 + 1 - tt) as u64; if !is_prime(n) { continue; }
            // j = -96^3: k = j/(1728-j), a = 3k, b = 2k
            let jm = (p - (884736 % p)) % p; let k = mulm(jm, inv((1728 + 884736) % p, p), p);
            let (mut a, mut b) = (3 * k % p, 2 * k % p);
            let c = Curve { p, a, b }; let mut x = 0; let pt = loop { if let Some(q) = c.lift(x) { break q; } x += 1; };
            if c.mul(n, pt) != Pt::Inf { let mut d = 2; while powm(d, (p - 1) / 2, p) == 1 { d += 1; } a = mulm(a, mulm(d, d, p), p); b = mulm(b, powm(d, 3, p), p); }
            let c2 = Curve { p, a, b }; let mut x = 0; let pt = loop { if let Some(q) = c2.lift(x) { break q; } x += 1; };
            if c2.mul(n, pt) != Pt::Inf { continue; }
            out.push(Toy { p, t: tt, u, n, a, b }); break; }
    } out
}

struct Rng(u64); impl Rng { fn next(&mut self) -> u64 { self.0 ^= self.0 << 13; self.0 ^= self.0 >> 7; self.0 ^= self.0 << 17; self.0 } }

fn coverage(s: &[u64], n: u64, m: usize) -> f64 { // fraction of Z/N that are sums of <= m elements of s (s symmetric)
    let mut bits = vec![0u64; (n as usize + 63) / 64]; let mut set = |x: u64| bits[(x / 64) as usize] |= 1 << (x % 64);
    let mut e: Vec<u64> = vec![0]; e.extend_from_slice(s); let l = e.len();
    match m { 1 => for &x in &e { set(x) },
        2 => for i in 0..l { for j in i..l { set((e[i] + e[j]) % n) } },
        _ => for i in 0..l { for j in i..l { let y = (e[i] + e[j]) % n; for k in j..l { set((y + e[k]) % n) } } } }
    let c: u64 = bits.iter().map(|w| w.count_ones() as u64).sum(); c as f64 / n as f64
}
fn sym(v: &[u64], n: u64) -> Vec<u64> { let mut s = vec![]; for &x in v { s.push(x % n); s.push((n - x % n) % n); } s }

fn ideal_set(lam: u64, n: u64, bcount: usize) -> (Vec<u64>, i64) { // B smallest-norm a+bw up to sign
    let mut el: Vec<(i64, i64, i64)> = vec![]; let r = 200i64;
    for b in 0..=r { for a in -r..=r { if b == 0 && a <= 0 { continue; } let nm = a * a + a * b + 5 * b * b; el.push((nm, a, b)); } }
    el.sort(); el.truncate(bcount); let maxn = el.last().unwrap().0;
    let v: Vec<u64> = el.iter().map(|&(_, a, b)| { let a = ((a % n as i64) + n as i64) as u64 % n; (a + mulm(b as u64, lam, n)) % n }).collect();
    (v, maxn)
}

fn main() {
    let mut js = String::from("{\n");
    // ---- C1: exact pair-46 checks ----
    let p: u64 = 5947698030651343; let n: u64 = 5947697996161643; let t: u64 = 34489701; let v: u64 = 34489703; let lam: u64 = 2973849084305079;
    let c_n = p + 1 - t == n; let c_p = 5 * t * t + 19 * t + 19 == p; let c_disc = (4 * p as u128) == (t as u128 * t as u128 + 19 * v as u128 * v as u128);
    let c_lam = (mulm(lam, lam, n) + n - lam + 5) % n == 0; let c_vlam = mulm(v, lam, n) == 2; let c_vlam_t = mulm(v, lam, n) == t; let c_vconj_2 = mulm(v, (n + 1 - lam) % n, n) == 2;
    let leg = powm(n_mod(v), (v - 1) / 2, v); let kron: i64 = if leg == 1 { 1 } else { -1 };
    fn n_mod(v: u64) -> u64 { v - 19 % v }
    let h_floor = (v as i64 - kron) as u64; let minnorm_nonint = 19u128 * (v as u128) * (v as u128) / 4;
    let other_root = (n + 1 - lam) % n;
    let c_prime = (is_prime(p), is_prime(n), is_prime(v));
    // crater curve check: lambda acts on a point? we can only check that lam is root; also check E(F_p) order N on crater/floor via random point
    let crater = Curve { p, a: 5757116383569871, b: 3560173328417442 }; let floor = Curve { p, a: 2087370704375809, b: 1402887937802154 };
    let mut x = 1; let pc = loop { if let Some(q) = crater.lift(x) { break q; } x += 1; }; let mut x = 1; let pf = loop { if let Some(q) = floor.lift(x) { break q; } x += 1; };
    let c_ordc = crater.mul(n, pc) == Pt::Inf; let c_ordf = floor.mul(n, pf) == Pt::Inf;
    println!("pair46: N=p+1-t {c_n}, p=5t^2+19t+19 {c_p}, 4p=t^2+19v^2 {c_disc}, lam^2-lam+5=0 {c_lam}, v*lam=2 {c_vlam}, v*lam=t {c_vlam_t}, v*(1-lam)=2 {c_vconj_2}, primes {:?}, [N]P=O crater {c_ordc} floor {c_ordf}", c_prime);
    println!("pair46: (-19/v)={kron}, h(Z[pi])={h_floor}, min norm of non-integer in Z[pi] >= {minnorm_nonint}, other root {other_root}");
    write!(js, "\"pair46_exact\": {{\"N_eq_p1t\": {c_n}, \"p_eq_5t2_19t_19\": {c_p}, \"4p_eq_t2_19v2\": {c_disc}, \"lambda_root_x2_x_5\": {c_lam}, \"v_lambda_eq_2_modN\": {c_vlam}, \"v_lambda_eq_t_modN\": {c_vlam_t}, \"v_conjlambda_eq_2_modN\": {c_vconj_2}, \"p_prime\": {}, \"N_prime\": {}, \"v_prime\": {}, \"NP_inf_crater\": {c_ordc}, \"NP_inf_floor\": {c_ordf}, \"kronecker_m19_v\": {kron}, \"h_floor_Zpi\": {h_floor}, \"min_norm_noninteger_Zpi_ge\": \"{minnorm_nonint}\", \"lambda_conjugate\": {other_root}}},\n", c_prime.0, c_prime.1, c_prime.2).unwrap();

    // ---- C3: toy coverage ----
    let mut toys = find_toys(1 << 20, 6); toys.extend(find_toys(1 << 22, 2));
    let bs = [32usize, 64, 128, 256, 512, 1024];
    js.push_str("\"toys\": [\n");
    for (ci, ty) in toys.iter().enumerate() {
        let c = Curve { p: ty.p, a: ty.a, b: ty.b }; let nn = ty.n;
        let mut x = 0; let g = loop { if let Some(q) = c.lift(x) { break q; } x += 1; };
        // lambda roots mod N
        let disc = (nn - 19 % nn) % nn; let s = sqrt_mod(disc, nn).unwrap(); let i2 = inv(2, nn);
        let l1 = mulm((1 + s) % nn, i2, nn); let l2 = mulm((1 + nn - s) % nn, i2, nn);
        // enumerate multiples to get logs of small-x points
        let xmax = 6000u64; let mut xlog = vec![0u64; xmax as usize]; let mut q = g;
        for k in 1..nn { if let Pt::A(xx, _) = q { if xx < xmax && xlog[xx as usize] == 0 { xlog[xx as usize] = k; } } q = c.add(q, g); }
        let smallx: Vec<u64> = (0..xmax).filter(|&x| xlog[x as usize] != 0).map(|x| xlog[x as usize]).collect();
        println!("curve {ci}: p={} t={} u={} N={} a={} b={} G={:?} lambdas {l1} {l2}", ty.p, ty.t, ty.u, nn, ty.a, ty.b, match g { Pt::A(x, y) => (x, y), _ => (0, 0) });
        write!(js, "{{\"idx\": {ci}, \"p\": {}, \"t\": {}, \"u\": {}, \"N\": {}, \"a\": {}, \"b\": {}, \"lambda1\": {l1}, \"lambda2\": {l2}, \"rows\": [\n", ty.p, ty.t, ty.u, nn, ty.a, ty.b).unwrap();
        let mut rows = vec![];
        for &bc in &bs { for m in 1..=3usize { if m == 3 && bc > 128 { continue; }
            let (i1, maxn) = ideal_set(l1, nn, bc); let (i2s, _) = ideal_set(l2, nn, bc);
            let ap: Vec<u64> = (1..=bc as u64).collect();
            let cov_i1 = coverage(&sym(&i1, nn), nn, m); let cov_i2 = coverage(&sym(&i2s, nn), nn, m);
            let cov_ap = coverage(&sym(&ap, nn), nn, m); let cov_sx = coverage(&sym(&smallx[..bc], nn), nn, m);
            let mut nulls = vec![]; for rep in 0..20u64 { let mut r = Rng(0x9E3779B97F4A7C15 ^ (1000 * ci as u64 + rep + 1) * 0x2545F4914F6CDD1D);
                let mut vv = vec![]; let mut seen = std::collections::HashSet::new();
                while vv.len() < bc { let x = 1 + r.next() % (nn - 1); let key = x.min(nn - x); if seen.insert(key) { vv.push(x); } }
                nulls.push(coverage(&sym(&vv, nn), nn, m)); }
            let ns: Vec<String> = nulls.iter().map(|x| format!("{x:.8}")).collect();
            println!("  B={bc} m={m} maxnorm={maxn} ideal1={cov_i1:.6} ideal2={cov_i2:.6} ap={cov_ap:.6} smallx={cov_sx:.6} null_mean={:.6}", nulls.iter().sum::<f64>() / 20.0);
            rows.push(format!("{{\"B\": {bc}, \"m\": {m}, \"max_norm\": {maxn}, \"ideal_l1\": {cov_i1:.8}, \"ideal_l2\": {cov_i2:.8}, \"ap_pos_and_floor\": {cov_ap:.8}, \"smallx\": {cov_sx:.8}, \"null\": [{}]}}", ns.join(", ")));
        } }
        js.push_str(&rows.join(",\n")); js.push_str(if ci + 1 < toys.len() { "]},\n" } else { "]}\n" });
    }
    js.push_str("]\n}\n");
    std::fs::write("raw.json", js).unwrap();
}
