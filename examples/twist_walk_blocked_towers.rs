//! Can a short 3-isogeny walk leave a GHS tower that MMT Lemma 6 condition (8)
//! blocks on the quadratic twist?
//!
//! ```text
//! cargo run --release --example twist_walk_blocked_towers -- control
//! cargo run --release --example twist_walk_blocked_towers -- catalog \
//!     --registry docs/curves/registry.json --out research/twist_walk_blocked_towers_20261014/catalog_towers.csv
//! cargo run --release --example twist_walk_blocked_towers -- null \
//!     --registry docs/curves/registry.json --samples 300 --seed 20261014 \
//!     --out research/twist_walk_blocked_towers_20261014/null_towers.csv
//! ```
//!
//! A quadratic twist keeps `b` and flips `Tr_{K/F_2}(a)`.  The twist of a tower
//! `(N, l, n = N/l)` is blocked iff `n` is even and `Tr_{K/k}(b) = 0`.  The walk
//! is over `j`-invariants with edges the roots of `Φ_3(X, j) mod 2`; Frobenius
//! `j ↦ j²` is free (it preserves `m` and whether `Tr_{K/k}(b) = 0`), so nodes
//! are Frobenius orbits.  Degree-2 steps are Frobenius in characteristic 2 and
//! are not walked.  Protocol: `research/twist_walk_blocked_towers_20261014/PROTOCOL.md`.

use crypto_lib::binary_ecc::{F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::binary_isogeny::{
    find_roots_in_f2m, phi_l_mod2_in_x, phi_l_mod2_monomials,
};
use crypto_lib::cryptanalysis::ec_trapdoor::{
    absolute_trace, magic_number, mmt_condition_8, relative_trace, FieldTower,
};
use num_bigint::BigUint;
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::collections::{HashMap, VecDeque};
use std::{fs, io::Write, time::Instant};

/// Walk depth, in 3-isogeny steps.
const DEPTH: usize = 4;

/// The 13 blocked twist towers of the pilot (`catalog_towers.csv`): registry
/// slug and base degree `l`.
const FROZEN_TOWERS: [(&str, u32); 13] = [
    ("icv1-f2m18-t393-d2c04152", 1),
    ("icv1-f2m18-t393-d2c04152", 3),
    ("icv1-f2m24-t1545-54539071", 1),
    ("icv1-f2m24-tm3047-2d0c0e8c", 1),
    ("icv1-f2m24-tm3047-2d0c0e8c", 2),
    ("icv1-f2m208-t32421385432645202936390897948121-225f5036", 1),
    ("icv1-f2m208-t32421385432645202936390897948121-225f5036", 2),
    ("icv1-f2m208-t32421385432645202936390897948121-225f5036", 13),
    ("icv1-f2m208-t32421385432645202936390897948121-225f5036", 26),
    (
        "icv1-f2m368-tm15466245949934416937841361495935610942842793802409489935-43ee62b3",
        1,
    ),
    (
        "icv1-f2m368-tm15466245949934416937841361495935610942842793802409489935-43ee62b3",
        2,
    ),
    (
        "icv1-f2m368-tm15466245949934416937841361495935610942842793802409489935-43ee62b3",
        23,
    ),
    (
        "icv1-f2m368-tm15466245949934416937841361495935610942842793802409489935-43ee62b3",
        46,
    ),
];

struct Field {
    m: u32,
    irr: IrreduciblePoly,
}

fn parse_hex(s: &str) -> BigUint {
    BigUint::parse_bytes(s.trim_start_matches("0x").as_bytes(), 16).expect("hex")
}

fn field_from_modulus(modulus: &BigUint) -> Field {
    let m = (modulus.bits() - 1) as u32;
    let low_terms = (0..m).filter(|&i| modulus.bit(u64::from(i))).collect();
    Field {
        m,
        irr: IrreduciblePoly {
            degree: m,
            low_terms,
        },
    }
}

/// Smallest raw-bit representative of the Frobenius orbit `{j^(2^i)}`.
fn canon(j: &F2mElement, f: &Field) -> Vec<u64> {
    let mut cur = j.clone();
    let mut best = j.raw_bits().to_vec();
    for _ in 1..f.m {
        cur = cur.square(&f.irr);
        if cur.raw_bits() < best.as_slice() {
            best = cur.raw_bits().to_vec();
        }
    }
    best
}

/// Breadth-first search of the `Φ_3` graph on Frobenius orbits of `j`, to
/// `DEPTH`, with edges supplied by `neighbours`.  Returns `(j, depth)` for each
/// node, the start first.
fn walk_with<F: Fn(&F2mElement) -> Vec<F2mElement>>(
    j0: &F2mElement,
    f: &Field,
    neighbours: F,
) -> Vec<(F2mElement, usize)> {
    let mut seen: HashMap<Vec<u64>, usize> = HashMap::new();
    seen.insert(canon(j0, f), 0);
    let mut nodes = vec![(j0.clone(), 0usize)];
    let mut queue = VecDeque::from([(j0.clone(), 0usize)]);
    while let Some((j, d)) = queue.pop_front() {
        if d == DEPTH {
            continue;
        }
        let mut roots = neighbours(&j);
        roots.sort_by_key(|r| r.to_biguint());
        for r in roots {
            if r.is_zero() {
                continue;
            }
            let key = canon(&r, f);
            if seen.contains_key(&key) {
                continue;
            }
            seen.insert(key, d + 1);
            nodes.push((r.clone(), d + 1));
            queue.push_back((r, d + 1));
        }
    }
    nodes
}

/// The walk the experiment uses: edges are the roots of `Φ_3(X, j) mod 2`.
fn walk(j0: &F2mElement, f: &Field) -> Vec<(F2mElement, usize)> {
    walk_with(j0, f, |j| {
        let poly = phi_l_mod2_in_x(3, j, f.m, &f.irr);
        find_roots_in_f2m(&poly, f.m, &f.irr)
    })
}

/// All `x` with `Φ_3(x, j) = 0` by scanning every element of the field and
/// evaluating `Φ_3` from its tabled monomials.  Independent of the root finder
/// and of the polynomial packing in `phi_l_mod2_in_x`.
fn exhaustive_roots(j: &F2mElement, f: &Field, monos: &[(u32, u32)]) -> Vec<F2mElement> {
    let mut jp = vec![F2mElement::one(f.m)];
    for _ in 0..5 {
        jp.push(jp.last().unwrap().mul(j, &f.irr));
    }
    let mut out = Vec::new();
    for v in 0u32..(1u32 << f.m) {
        let pos: Vec<u32> = (0..f.m).filter(|i| (v >> i) & 1 == 1).collect();
        let x = F2mElement::from_bit_positions(&pos, f.m);
        let mut xp = vec![F2mElement::one(f.m)];
        for _ in 0..5 {
            xp.push(xp.last().unwrap().mul(&x, &f.irr));
        }
        let mut acc = F2mElement::zero(f.m);
        for &(i, k) in monos {
            acc = acc.add(&xp[i as usize].mul(&jp[k as usize], &f.irr));
        }
        if acc.is_zero() {
            out.push(x);
        }
    }
    out
}

/// `z^i` with absolute trace 1, the `a` of the quadratic twist class.
fn twist_a(f: &Field) -> F2mElement {
    let tower = FieldTower::new(f.m, f.m, 1, f.irr.clone());
    for i in 0..f.m {
        let x = F2mElement::from_bit_positions(&[i], f.m);
        if absolute_trace(&tower, &x) {
            return x;
        }
    }
    unreachable!("the absolute trace is onto F_2")
}

fn random_element(rng: &mut StdRng, m: u32) -> F2mElement {
    let pos: Vec<u32> = (0..m).filter(|_| rng.gen_bool(0.5)).collect();
    F2mElement::from_bit_positions(&pos, m)
}

/// A uniform-ish nonzero `b` with `Tr_{K/k}(b) = 0`: a random `b`, minus a
/// `w` with the same relative trace (`w = c w0`, `c = Tr(b) / Tr(w0) ∈ k`).
fn random_blocked_b(rng: &mut StdRng, tower: &FieldTower) -> F2mElement {
    let f_irr = &tower.irr;
    loop {
        let b = random_element(rng, tower.big_n);
        let tb = relative_trace(tower, &b);
        let candidate = if tb.is_zero() {
            b
        } else {
            let w0 = loop {
                let w = random_element(rng, tower.big_n);
                if !relative_trace(tower, &w).is_zero() {
                    break w;
                }
            };
            let theta = relative_trace(tower, &w0);
            let c = tb.mul(&theta.flt_inverse(f_irr).expect("theta != 0"), f_irr);
            b.add(&c.mul(&w0, f_irr))
        };
        if !candidate.is_zero() && relative_trace(tower, &candidate).is_zero() {
            return candidate;
        }
    }
}

struct TowerRow {
    nodes_by_depth: [usize; DEPTH + 1],
    unblocked_by_depth: [usize; DEPTH + 1],
    min_m_unblocked: [Option<u32>; DEPTH + 1],
    start_m: u32,
    start_blocked: bool,
    first_unblocked_depth: Option<usize>,
}

fn eval_tower(nodes: &[(F2mElement, usize)], tower: &FieldTower, a_tw: &F2mElement) -> TowerRow {
    let zero_a = F2mElement::zero(tower.big_n);
    let mut row = TowerRow {
        nodes_by_depth: [0; DEPTH + 1],
        unblocked_by_depth: [0; DEPTH + 1],
        min_m_unblocked: [None; DEPTH + 1],
        start_m: 0,
        start_blocked: false,
        first_unblocked_depth: None,
    };
    for (j, d) in nodes {
        let b = j.flt_inverse(&tower.irr).expect("j != 0");
        let m_prime = magic_number(tower, &zero_a, &b);
        let ok = mmt_condition_8(tower, a_tw, &b);
        row.nodes_by_depth[*d] += 1;
        if *d == 0 {
            row.start_m = m_prime;
            row.start_blocked = !ok;
        }
        if ok {
            row.unblocked_by_depth[*d] += 1;
            for dd in *d..=DEPTH {
                row.min_m_unblocked[dd] =
                    Some(row.min_m_unblocked[dd].map_or(m_prime, |x| x.min(m_prime)));
            }
            if row.first_unblocked_depth.is_none_or(|x| *d < x) {
                row.first_unblocked_depth = Some(*d);
            }
        }
    }
    row
}

const HEADER: &str = "population,sample,N,l,n,start_m,start_blocked,nodes_d0,nodes_d1,nodes_d2,nodes_d3,nodes_d4,unblocked_d0,unblocked_d1,unblocked_d2,unblocked_d3,unblocked_d4,min_m_unblocked_le1,min_m_unblocked_le2,min_m_unblocked_le3,min_m_unblocked_le4,first_unblocked_depth";

fn csv_row(population: &str, sample: &str, n_total: u32, l: u32, r: &TowerRow) -> String {
    let opt = |x: Option<u32>| x.map_or(String::new(), |v| v.to_string());
    format!(
        "{population},{sample},{n_total},{l},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
        n_total / l,
        r.start_m,
        r.start_blocked,
        r.nodes_by_depth[0],
        r.nodes_by_depth[1],
        r.nodes_by_depth[2],
        r.nodes_by_depth[3],
        r.nodes_by_depth[4],
        r.unblocked_by_depth[0],
        r.unblocked_by_depth[1],
        r.unblocked_by_depth[2],
        r.unblocked_by_depth[3],
        r.unblocked_by_depth[4],
        opt(r.min_m_unblocked[1]),
        opt(r.min_m_unblocked[2]),
        opt(r.min_m_unblocked[3]),
        opt(r.min_m_unblocked[4]),
        r.first_unblocked_depth
            .map_or("none".to_string(), |d| d.to_string()),
    )
}

fn load_registry(path: &str) -> serde_json::Value {
    serde_json::from_str(&fs::read_to_string(path).expect("registry")).expect("json")
}

fn find_curve<'a>(reg: &'a serde_json::Value, slug: &str) -> &'a serde_json::Value {
    reg["curves"]
        .as_array()
        .expect("curves")
        .iter()
        .find(|c| c["slug"] == slug)
        .unwrap_or_else(|| panic!("slug {slug} not in registry"))
}

fn arg(args: &[String], k: &str) -> Option<String> {
    args.iter()
        .position(|a| a == k)
        .map(|i| args[i + 1].clone())
}

fn write_out(path: &str, rows: &[String]) {
    let mut f = fs::File::create(path).expect("out");
    writeln!(f, "{HEADER}").unwrap();
    for r in rows {
        writeln!(f, "{r}").unwrap();
    }
}

/// Positive control A: at `N = 20 > 16`, where `find_roots_in_f2m` takes its
/// Cantor-Zassenhaus path, the roots of `Φ_3(X, j)` equal an exhaustive search.
/// Draws continue until 8 values of `j` have at least two roots (splitting is
/// what the root finder must get right), capped at 150 draws.
///
/// Positive control B: at `N = 18`, the whole walk plus tower evaluation
/// (`walk` + `eval_tower`) equals the same BFS fed by exhaustive neighbour
/// scans, on 10 random blocked starts that have at least one neighbour.
fn control(seed: u64) -> bool {
    let monos = phi_l_mod2_monomials(3).expect("Φ_3 tabled");
    let mut rng = StdRng::seed_from_u64(seed);
    let mut all_ok = true;

    let f20 = Field {
        m: 20,
        irr: IrreduciblePoly {
            degree: 20,
            low_terms: vec![0, 3],
        },
    };
    let (mut multi, mut trial, mut total_roots) = (0usize, 0usize, 0usize);
    while multi < 8 && trial < 150 {
        let j = loop {
            let x = random_element(&mut rng, f20.m);
            if !x.is_zero() {
                break x;
            }
        };
        let poly = phi_l_mod2_in_x(3, &j, f20.m, &f20.irr);
        let mut got: Vec<BigUint> = find_roots_in_f2m(&poly, f20.m, &f20.irr)
            .iter()
            .map(|r| r.to_biguint())
            .collect();
        got.sort();
        // cheap filter: only j with >= 2 roots found, plus every 10th draw
        if got.len() < 2 && trial % 10 != 0 {
            trial += 1;
            continue;
        }
        let mut want: Vec<BigUint> = exhaustive_roots(&j, &f20, &monos)
            .iter()
            .map(|r| r.to_biguint())
            .collect();
        want.sort();
        let ok = got == want;
        total_roots += want.len();
        if want.len() >= 2 {
            multi += 1;
        }
        println!(
            "control A draw {trial}: roots found {}, exhaustive {} -> {}",
            got.len(),
            want.len(),
            if ok { "PASS" } else { "FAIL" }
        );
        all_ok &= ok;
        trial += 1;
    }
    all_ok &= multi >= 8;
    println!(
        "control A: {} ({total_roots} roots checked; {multi} draws with >= 2 roots)",
        if all_ok { "PASS" } else { "FAIL" }
    );

    let f18 = Field {
        m: 18,
        irr: IrreduciblePoly {
            degree: 18,
            low_terms: vec![0, 7],
        },
    };
    let tower = FieldTower::new(18, 6, 3, f18.irr.clone());
    let a_tw = twist_a(&f18);
    let (mut done, mut draws, mut b_ok) = (0usize, 0usize, true);
    while done < 10 && draws < 600 {
        draws += 1;
        let b = random_blocked_b(&mut rng, &tower);
        let j0 = b.flt_inverse(&f18.irr).expect("b != 0");
        let reference = walk_with(&j0, &f18, |j| exhaustive_roots(j, &f18, &monos));
        if reference.len() < 2 {
            continue;
        }
        let a = eval_tower(&walk(&j0, &f18), &tower, &a_tw);
        let r = eval_tower(&reference, &tower, &a_tw);
        let same = csv_row("c", "0", 18, 3, &a) == csv_row("c", "0", 18, 3, &r);
        println!(
            "control B start {draws}: {} nodes, first unblocked depth {:?} -> {}",
            reference.len(),
            r.first_unblocked_depth,
            if same { "PASS" } else { "FAIL" }
        );
        b_ok &= same;
        done += 1;
    }
    b_ok &= done >= 10;
    println!(
        "control B: {} ({done} starts with neighbours, {draws} draws)",
        if b_ok { "PASS" } else { "FAIL" }
    );
    all_ok & b_ok
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mode = args.get(1).map(String::as_str).unwrap_or("");
    let t0 = Instant::now();
    match mode {
        "control" => {
            let seed: u64 = arg(&args, "--seed").map_or(20261014, |s| s.parse().unwrap());
            if !control(seed) {
                std::process::exit(1);
            }
        }
        "catalog" => {
            let reg = load_registry(&arg(&args, "--registry").expect("--registry"));
            let out = arg(&args, "--out").expect("--out");
            let mut rows = Vec::new();
            let mut walks: HashMap<&str, Vec<(F2mElement, usize)>> = HashMap::new();
            for &(slug, l) in &FROZEN_TOWERS {
                let c = find_curve(&reg, slug);
                let f = field_from_modulus(&parse_hex(c["params"]["modulus"].as_str().unwrap()));
                let b =
                    F2mElement::from_biguint(&parse_hex(c["params"]["b"].as_str().unwrap()), f.m);
                let j0 = b.flt_inverse(&f.irr).expect("b != 0");
                let nodes = walks.entry(slug).or_insert_with(|| walk(&j0, &f));
                let tower = FieldTower::new(f.m, f.m / l, l, f.irr.clone());
                let r = eval_tower(nodes, &tower, &twist_a(&f));
                assert!(
                    r.start_blocked,
                    "{slug} l={l}: the pilot says the twist is blocked, the native check disagrees"
                );
                println!(
                    "catalog {} N={} l={l}: start m={}, nodes/depth {:?}, first unblocked depth {:?}",
                    &slug[..slug.len().min(24)],
                    f.m,
                    r.start_m,
                    r.nodes_by_depth,
                    r.first_unblocked_depth
                );
                rows.push(csv_row("catalog", slug, f.m, l, &r));
            }
            write_out(&out, &rows);
        }
        "null" => {
            let reg = load_registry(&arg(&args, "--registry").expect("--registry"));
            let out = arg(&args, "--out").expect("--out");
            let samples: usize = arg(&args, "--samples").map_or(300, |s| s.parse().unwrap());
            let seed: u64 = arg(&args, "--seed").map_or(20261014, |s| s.parse().unwrap());
            let mut rng = StdRng::seed_from_u64(seed);
            // One null population per distinct (N, l) among the frozen towers.
            let mut configs: Vec<(&str, u32)> = Vec::new();
            for &(slug, l) in &FROZEN_TOWERS {
                let n_total = field_degree(&reg, slug);
                if !configs
                    .iter()
                    .any(|&(s, ll)| ll == l && field_degree(&reg, s) == n_total)
                {
                    configs.push((slug, l));
                }
            }
            let mut rows = Vec::new();
            for (slug, l) in configs {
                let c = find_curve(&reg, slug);
                let f = field_from_modulus(&parse_hex(c["params"]["modulus"].as_str().unwrap()));
                let tower = FieldTower::new(f.m, f.m / l, l, f.irr.clone());
                let a_tw = twist_a(&f);
                let t1 = Instant::now();
                for s in 0..samples {
                    let b = random_blocked_b(&mut rng, &tower);
                    let j0 = b.flt_inverse(&f.irr).expect("b != 0");
                    let nodes = walk(&j0, &f);
                    let r = eval_tower(&nodes, &tower, &a_tw);
                    assert!(r.start_blocked, "constructed start must be blocked");
                    rows.push(csv_row("null", &s.to_string(), f.m, l, &r));
                }
                println!(
                    "null N={} l={l}: {samples} starts in {:.1}s",
                    f.m,
                    t1.elapsed().as_secs_f64()
                );
            }
            write_out(&out, &rows);
        }
        "classes" => {
            let lo: u32 = arg(&args, "--min-n").map_or(8, |s| s.parse().unwrap());
            let hi: u32 = arg(&args, "--max-n").map_or(20, |s| s.parse().unwrap());
            let out = arg(&args, "--out").expect("--out");
            class_ceiling(lo, hi, &out);
        }
        "summary" => {
            let catalog = arg(&args, "--catalog").expect("--catalog CSV");
            let null = arg(&args, "--null").expect("--null CSV");
            let out = arg(&args, "--out").expect("--out");
            summary(&catalog, &null, &out);
        }
        _ => {
            eprintln!("usage: twist_walk_blocked_towers control|catalog|null|summary|classes [--registry F] [--out F] [--samples N] [--seed S] [--catalog F --null F]");
            std::process::exit(2);
        }
    }
    eprintln!("elapsed {:.1}s", t0.elapsed().as_secs_f64());
}

fn field_degree(reg: &serde_json::Value, slug: &str) -> u64 {
    find_curve(reg, slug)["params"]["m"].as_u64().expect("m")
}

// ── Summary: counts and Clopper-Pearson intervals ───────────────────────

fn ln_gamma(x: f64) -> f64 {
    // Lanczos approximation (g = 7, n = 9), accurate to ~1e-13 for x > 0.
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = C[0];
    let t = x + 7.5;
    for (i, c) in C.iter().enumerate().skip(1) {
        a += c / (x + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

/// `P(X <= k)` for `X ~ Binomial(n, p)`, summed in log space.
fn binom_cdf(k: u64, n: u64, p: f64) -> f64 {
    if p <= 0.0 {
        return 1.0;
    }
    if p >= 1.0 {
        return if k >= n { 1.0 } else { 0.0 };
    }
    let (lp, lq) = (p.ln(), (1.0 - p).ln());
    let lnn = ln_gamma(n as f64 + 1.0);
    let mut total = 0.0;
    for i in 0..=k.min(n) {
        let ln_pmf = lnn - ln_gamma(i as f64 + 1.0) - ln_gamma((n - i) as f64 + 1.0)
            + i as f64 * lp
            + (n - i) as f64 * lq;
        total += ln_pmf.exp();
    }
    total.min(1.0)
}

/// Two-sided 95% Clopper-Pearson interval for `k` successes in `n` trials.
fn clopper_pearson(k: u64, n: u64) -> (f64, f64) {
    let alpha = 0.05;
    let bisect = |f: &dyn Fn(f64) -> f64, target: f64, increasing: bool| {
        // find p in [0, 1] with f(p) = target; f monotone
        let (mut lo, mut hi) = (0.0f64, 1.0f64);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if (f(mid) < target) == increasing {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    };
    // lower: P(X >= k | p) = alpha/2, increasing in p
    let lower = if k == 0 {
        0.0
    } else {
        bisect(&|p| 1.0 - binom_cdf(k - 1, n, p), alpha / 2.0, true)
    };
    // upper: P(X <= k | p) = alpha/2, decreasing in p
    let upper = if k == n {
        1.0
    } else {
        bisect(&|p| binom_cdf(k, n, p), alpha / 2.0, false)
    };
    (lower, upper)
}

struct Cfg {
    n_total: u32,
    l: u32,
    null_n: u64,
    has_neighbour: u64,
    escape: [u64; DEPTH + 1],
    escape_generic_m: u64,
}

fn summary(catalog_csv: &str, null_csv: &str, out: &str) {
    let read = |path: &str| -> Vec<Vec<String>> {
        let text = fs::read_to_string(path).expect("csv");
        text.lines()
            .skip(1)
            .map(|l| l.split(',').map(str::to_string).collect())
            .collect()
    };
    let col = |name: &str| HEADER.split(',').position(|c| c == name).expect(name);
    let (cn, cl, cnn) = (col("N"), col("l"), col("n"));
    let (cnodes1, cfirst) = (col("nodes_d1"), col("first_unblocked_depth"));
    let cmin4 = col("min_m_unblocked_le4");
    let mut cfgs: Vec<Cfg> = Vec::new();
    for r in read(null_csv) {
        let (n_total, l): (u32, u32) = (r[cn].parse().unwrap(), r[cl].parse().unwrap());
        let i = match cfgs.iter().position(|c| c.n_total == n_total && c.l == l) {
            Some(i) => i,
            None => {
                cfgs.push(Cfg {
                    n_total,
                    l,
                    null_n: 0,
                    has_neighbour: 0,
                    escape: [0; DEPTH + 1],
                    escape_generic_m: 0,
                });
                cfgs.len() - 1
            }
        };
        let c = &mut cfgs[i];
        c.null_n += 1;
        if r[cnodes1].parse::<u64>().unwrap() > 0 {
            c.has_neighbour += 1;
        }
        if let Ok(first) = r[cfirst].parse::<usize>() {
            for d in first..=DEPTH {
                c.escape[d] += 1;
            }
            let n_rel: u32 = r[cnn].parse().unwrap();
            if r[cmin4].parse::<u32>().ok() == Some(n_rel) {
                c.escape_generic_m += 1;
            }
        }
    }
    let mut lines = vec![
        "kind,N,l,n,sample,start_m,has_neighbour,escape_le1,escape_le2,escape_le3,escape_le4,escape_le4_ci_lo,escape_le4_ci_hi,has_neighbour_ci_lo,has_neighbour_ci_hi,escape_with_generic_m,first_unblocked_depth,min_m_unblocked_le4".to_string(),
    ];
    for c in &cfgs {
        let (elo, ehi) = clopper_pearson(c.escape[DEPTH], c.null_n);
        let (nlo, nhi) = clopper_pearson(c.has_neighbour, c.null_n);
        lines.push(format!(
            "null,{},{},{},{},,{}/{},{}/{},{}/{},{}/{},{}/{},{elo:.4},{ehi:.4},{nlo:.4},{nhi:.4},{}/{},,",
            c.n_total, c.l, c.n_total / c.l, c.null_n, c.has_neighbour, c.null_n,
            c.escape[1], c.null_n, c.escape[2], c.null_n, c.escape[3], c.null_n, c.escape[4], c.null_n,
            c.escape_generic_m, c.escape[DEPTH].max(1),
        ));
    }
    let cn_idx = (col("sample"), col("start_m"));
    for r in read(catalog_csv) {
        lines.push(format!(
            "catalog,{},{},{},{},{},{},,,,,,,,,,{},{}",
            r[cn],
            r[cl],
            r[cnn],
            r[cn_idx.0],
            r[cn_idx.1],
            r[cnodes1].parse::<u64>().unwrap() > 0,
            r[cfirst],
            r[cmin4],
        ));
    }
    fs::write(out, lines.join("\n") + "\n").expect("summary out");
    println!("wrote {out}");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clopper_pearson_matches_known_values() {
        // 0 of 2000: upper = 1 - (alpha/2)^(1/n)
        let (lo, hi) = clopper_pearson(0, 2000);
        assert_eq!(lo, 0.0);
        assert!(
            (hi - (1.0 - 0.025f64.powf(1.0 / 2000.0))).abs() < 1e-9,
            "hi = {hi}"
        );
        // 10 of 20: the standard exact interval (0.2719, 0.7281)
        let (lo, hi) = clopper_pearson(10, 20);
        assert!(
            (lo - 0.2719).abs() < 5e-4 && (hi - 0.7281).abs() < 5e-4,
            "({lo}, {hi})"
        );
        // 20 of 20: lower = (alpha/2)^(1/n)
        let (lo, hi) = clopper_pearson(20, 20);
        assert_eq!(hi, 1.0);
        assert!((lo - 0.025f64.powf(1.0 / 20.0)).abs() < 1e-9, "lo = {lo}");
    }
}

/// Exact class-level ceiling.  The twist-side isogeny class of `b` is
/// `{b' : t(b') = t(b)}` (the `a = 0` family has trace `t(b)`, its twist
/// `-t(b)`, and `b` is unchanged by the twist).  For every `N` in `lo..=hi` and
/// every `l | N` with `n = N/l` even, count the trace classes with a blocked
/// twist tower (`Tr_{K/k}(b) = 0`), those with no unblocked curve at all, and
/// the mixed ones, where some isogenous curve is unblocked.  A walk can leave a
/// blocked tower only inside a mixed class.
fn class_ceiling(lo: u32, hi: u32, out: &str) {
    use crypto_lib::cryptanalysis::ghs_class_coverage::{traces_a0, Gf2Log};
    let mut lines = vec!["N,l,n,curves,classes,classes_with_blocked,classes_all_blocked,classes_mixed,curves_blocked,curves_blocked_in_mixed_classes".to_string()];
    for n_total in lo..=hi {
        let f = Gf2Log::new(n_total);
        let t = traces_a0(&f);
        for l in (1..=n_total).filter(|l| n_total % l == 0 && (n_total / l) % 2 == 0) {
            // per trace value: (blocked count, unblocked count)
            let mut per: HashMap<i32, (u64, u64)> = HashMap::new();
            for b in 1..f.size() as u32 {
                let e = per.entry(t[b as usize]).or_insert((0, 0));
                if f.rel_trace(b, l) == 0 {
                    e.0 += 1;
                } else {
                    e.1 += 1;
                }
            }
            let classes = per.len();
            let with_blocked = per.values().filter(|v| v.0 > 0).count();
            let all_blocked = per.values().filter(|v| v.0 > 0 && v.1 == 0).count();
            let mixed = per.values().filter(|v| v.0 > 0 && v.1 > 0).count();
            let blocked: u64 = per.values().map(|v| v.0).sum();
            let blocked_mixed: u64 = per.values().filter(|v| v.1 > 0).map(|v| v.0).sum();
            lines.push(format!(
                "{n_total},{l},{},{},{classes},{with_blocked},{all_blocked},{mixed},{blocked},{blocked_mixed}",
                n_total / l,
                f.size() - 1
            ));
        }
        eprintln!("N = {n_total} done");
    }
    fs::write(out, lines.join("\n") + "\n").expect("classes out");
    println!("wrote {out}");
}
