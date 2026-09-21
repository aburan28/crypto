//! FFD measurement on the subspace-restricted, Weil-descended `S₄`
//! system at `l = 8` (`n = 3l = 24`) — the missing "FFD required"
//! datum for the binary decomposition beat in
//! `docs/ic/boundary_targets.json` (next_target: sub-2^(2l) oracle at
//! l=8, "ffd min/max/mean over >=16 draws present").
//!
//! # What is measured
//!
//! The decomposition problem: given a target `x_R ∈ F_{2^24}`, do
//! `X₁, X₂, X₃` in the factor-base subspace `V = ⟨1, z, …, z^{l−1}⟩`
//! exist with `S₄(X₁, X₂, X₃, x_R) = 0` (Koblitz `b = 1`)?
//!
//! `binary_semaev_s4::weil_descend_s4` builds the descended system in
//! two spaces: `n = 24` quadratic equations in the `e`-variables
//! (elementary symmetric functions, `6l − 3 = 45` of them) plus the
//! correspondence constraints `e_{i,d} = σ_{i,d}(X)` (linear /
//! quadratic / cubic in the `3l = 24` factor-base bits).  This program
//! **eliminates the `e`-variables exactly** (substituting the
//! correspondence ANFs — they are independent of `x_R`, so all
//! `e`-monomial images are computed once and cached), yielding `24`
//! equations of degree ≤ 6 in `24` Boolean variables, and measures the
//! first fall degree of that system with the degree-general Macaulay
//! harness (`ffd_harness::measure_monomial_system`).
//!
//! Two variants per draw:
//!   * **full** — all `24` x-bits free (`X₁` swept by the oracle);
//!   * **x1-fixed** — `X₁` frozen to a random subspace element, `16`
//!     unknowns (`X₂, X₃` bits).  This is the per-`X₁` system a
//!     Gröbner-based oracle would solve inside its `2^l` sweep, so its
//!     fall degree is the one that decides whether a sub-2^(2l) oracle
//!     of that shape can exist.
//!
//! # Correctness gates (fail-closed)
//!
//! 1. **Sanity**: the degree-general Macaulay builder must reproduce
//!    the battle-tested quadratic path (`run_sweep`, seed `0xFFDDEAD`)
//!    row-for-row, rank-for-rank on the descended `S₃` systems at
//!    `n ∈ {5,6,7}`, including the documented FFD = 3.
//! 2. **Elimination exactness**: on 64 random points of `V³` per draw,
//!    the eliminated system is satisfied **iff**
//!    `symmetrised_s4_eval(σ(X), x_R) = 0` — zero mismatches tolerated.
//! 3. **Witness**: whenever `semaev_decomp::decompose` finds a
//!    decomposition of the draw's target, the eliminated system (and
//!    its `X₁`-fixed fold at the witness's `X₁`) must vanish on it.
//!    Extra draws are taken until at least one decomposable target has
//!    been witness-checked.
//!
//! # Scope and honesty
//!
//! This is a MEASUREMENT run: no oracle is claimed and no cost win is
//! asserted.  Falls at degree `D ≥ 2·δ_min` may be Koszul syzygies
//! (`f_i f_j − f_j f_i` with `δ_i = δ_j = δ_min`) rather than
//! structure; each fall is flagged `structural` only when
//! `D < 2·δ_min`.  A run with no structural fall up to the censoring
//! degrees is evidence *against* the low-degree Gröbner route at this
//! rung — recorded as such, closing only the tested scope.
//!
//! ```bash
//! cargo run --release --example ffd_s4_subspace -- \
//!     --draws 16 --seed 0x54FFD518 --dmax-full 7 --dmax-fix 9 \
//!     --mem-cap-gb 1.2 --out <path>.json
//! ```

use crypto_lib::binary_ecc::{F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::binary_semaev_s4::{
    e_len, elementary_symmetric_3, symmetrised_s4_eval, weil_descend_s4, S4System,
};
use crypto_lib::cryptanalysis::ffd_harness::{
    choose_irreducible, macaulay_memory_estimate_monomial, measure_monomial_system, measure_one,
    quad_monomial_index, random_nonzero_f2m, run_sweep, weil_descend_s3, F2BoolPoly,
};
use crypto_lib::cryptanalysis::semaev_decomp::{decompose, Gf2};
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde_json::{json, Value};
use std::collections::HashMap;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

// ── Fixed parameters (match the recorded l=8 baseline arm) ─────────

const N: u32 = 24;
const L: u32 = 8;
/// Irreducible used by the recorded baseline arm
/// (`decomp_bench_to_l8.stdout.txt`: `n = 24  l = 8   [0, 1, 3, 4]`),
/// i.e. `z²⁴ + z⁴ + z³ + z + 1`.
const LOW_TERMS: [u32; 4] = [0, 1, 3, 4];
const D_MIN: u32 = 4;
/// Seed of the documented quadratic FFD sweep (`ffd_sweep_demo`).
const S3_SANITY_SEED: u64 = 0x_FFD_DEAD;

// ── Deterministic RNG (splitmix64) ─────────────────────────────────

struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

// ── Monomial packing (≤6 slots × 5 bits; var id + 1 per slot) ──────

fn pack_mono(m: &[u32]) -> u64 {
    let mut p = 0u64;
    for (k, &v) in m.iter().enumerate() {
        debug_assert!(v < 31 && k < 8);
        p |= ((v as u64) + 1) << (5 * k);
    }
    p
}

fn unpack_mono(p: u64) -> Vec<u32> {
    let mut out = Vec::new();
    for k in 0..8 {
        let slot = ((p >> (5 * k)) & 31) as u32;
        if slot == 0 {
            break;
        }
        out.push(slot - 1);
    }
    out
}

fn merge_squarefree(a: &[u32], b: &[u32]) -> Vec<u32> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => {
                out.push(a[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                out.push(b[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                out.push(a[i]);
                i += 1;
                j += 1;
            }
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
}

/// XOR-product of two squarefree-monomial sets, returned packed.
fn fast_mul_packed(a: &[Vec<u32>], b: &[Vec<u32>]) -> Vec<u64> {
    let mut acc: HashMap<u64, ()> = HashMap::new();
    for ma in a {
        for mb in b {
            let p = pack_mono(&merge_squarefree(ma, mb));
            if acc.remove(&p).is_none() {
                acc.insert(p, ());
            }
        }
    }
    let mut out: Vec<u64> = acc.into_keys().collect();
    out.sort_unstable();
    out
}

fn toggle(acc: &mut HashMap<u64, ()>, p: u64) {
    if acc.remove(&p).is_none() {
        acc.insert(p, ());
    }
}

// ── Sanity: degree-general builder vs the quadratic path ───────────

fn boolpolys_to_monomials(eqs: &[F2BoolPoly], v: u32) -> Vec<Vec<Vec<u32>>> {
    eqs.iter()
        .map(|p| {
            let mut monos = Vec::new();
            if p.coeffs[0] {
                monos.push(vec![]);
            }
            for i in 0..v {
                if p.coeffs[1 + i as usize] {
                    monos.push(vec![i]);
                }
            }
            for i in 0..v {
                for j in (i + 1)..v {
                    if p.coeffs[quad_monomial_index(i, j, v)] {
                        monos.push(vec![i, j]);
                    }
                }
            }
            monos
        })
        .collect()
}

fn sanity_s3() -> Value {
    println!("── sanity: degree-general Macaulay vs quadratic path (S₃, seed 0x{:X}) ──", S3_SANITY_SEED);
    let reference = run_sweep(3..=7, 4, S3_SANITY_SEED);
    let mut rng = StdRng::seed_from_u64(S3_SANITY_SEED);
    let mut checked = Vec::new();
    for n in 3..=7u32 {
        let irr = choose_irreducible(n);
        let b = random_nonzero_f2m(&mut rng, n);
        let x3 = random_nonzero_f2m(&mut rng, n);
        if n < 5 {
            continue; // keep the RNG stream in lockstep with run_sweep
        }
        // Independent reference for this exact input.
        let one = measure_one(n, &irr, &b, &x3, 4);
        let eqs = weil_descend_s3(n, &irr, &b, &x3);
        let mono = boolpolys_to_monomials(&eqs, 2 * n);
        let (per_deg, fall) = measure_monomial_system(&mono, 2 * n, 2, 4);
        assert_eq!(fall, one.fall_degree, "fall mismatch at n={n}");
        assert_eq!(per_deg.len(), one.per_degree.len(), "degree count mismatch at n={n}");
        for (g, r) in per_deg.iter().zip(&one.per_degree) {
            assert_eq!(g.degree, r.degree);
            assert_eq!(
                g.rows_constructed, r.rows_constructed,
                "rows mismatch n={n} D={}", r.degree
            );
            assert_eq!(g.cols, r.cols, "cols mismatch n={n} D={}", r.degree);
            assert_eq!(g.rank, r.rank, "rank mismatch n={n} D={}", r.degree);
        }
        assert_eq!(fall, Some(3), "documented FFD=3 not reproduced at n={n}");
        // Also cross-check against the run_sweep row for the same n.
        let refr = reference.iter().find(|r| r.n == n).expect("sweep row");
        assert_eq!(refr.fall_degree, Some(3));
        println!("   n={n}: rows/cols/rank match at D=2..4, FFD=3 ✓");
        checked.push(n);
    }
    println!("   sanity PASS\n");
    json!({
        "seed": format!("0x{:X}", S3_SANITY_SEED),
        "n_checked": checked,
        "match_reference": true,
        "ffd_all": 3,
    })
}

// ── e-elimination ──────────────────────────────────────────────────

/// Packed images of every `e`-variable under the correspondence
/// constraints (`e_{i,d} = σ_{i,d}(X)`, an ANF over the x-bits).
/// Independent of the target `x_R`.
fn build_singles(sys: &S4System) -> Vec<Vec<u64>> {
    let mut singles: Vec<Vec<u64>> = vec![Vec::new(); sys.n_e_vars() as usize];
    for i in 0..3usize {
        for d in 0..e_len(i + 1, L) {
            let idx = sys.e_var(i, d) as usize;
            singles[idx] = sys.correspondence[i][d]
                .monomials()
                .map(|m| pack_mono(m))
                .collect();
        }
    }
    singles
}

/// Substitute every `e`-variable in the descended `S₄` equations by
/// its correspondence image, yielding `n` equations of degree ≤ 6 over
/// the `3l` x-bits.  Pair images are memoized in `pair_cache` (they
/// depend only on `l`, never on `x_R`).
fn eliminate_e(
    sys: &S4System,
    singles: &[Vec<u64>],
    pair_cache: &mut HashMap<u64, Vec<u64>>,
) -> Vec<Vec<Vec<u32>>> {
    let mut eqs = Vec::with_capacity(sys.semaev.len());
    for eq in &sys.semaev {
        let mut acc: HashMap<u64, ()> = HashMap::new();
        for m in eq.monomials() {
            match m.len() {
                0 => toggle(&mut acc, 0),
                1 => {
                    for &p in &singles[m[0] as usize] {
                        toggle(&mut acc, p);
                    }
                }
                2 => {
                    let (a, b) = (m[0] as u64, m[1] as u64);
                    debug_assert!(a < b, "BTreeSet monomials are sorted and squarefree");
                    let key = a | (b << 6);
                    let img = pair_cache.entry(key).or_insert_with(|| {
                        let av: Vec<Vec<u32>> =
                            singles[a as usize].iter().map(|&p| unpack_mono(p)).collect();
                        let bv: Vec<Vec<u32>> =
                            singles[b as usize].iter().map(|&p| unpack_mono(p)).collect();
                        fast_mul_packed(&av, &bv)
                    });
                    for &p in img.iter() {
                        toggle(&mut acc, p);
                    }
                }
                k => panic!("semaev equations are quadratic in e-space; got monomial of len {k}"),
            }
        }
        let mut monos: Vec<Vec<u32>> = acc.into_keys().map(|p| unpack_mono(p)).collect();
        monos.sort();
        eqs.push(monos);
    }
    eqs
}

/// Freeze `X₁ = x1` (bit `j` of `x1` ← variable `j`), collapsing the
/// system to the 16 bits of `X₂` (vars 8..15 → 0..7) and `X₃`
/// (vars 16..23 → 8..15).
fn fold_x1(eqs: &[Vec<Vec<u32>>], x1: u64) -> Vec<Vec<Vec<u32>>> {
    eqs.iter()
        .map(|eq| {
            let mut acc: HashMap<u64, ()> = HashMap::new();
            for m in eq {
                let mut folded: Vec<u32> = Vec::with_capacity(m.len());
                let mut killed = false;
                for &v in m {
                    if v < L {
                        if (x1 >> v) & 1 == 0 {
                            killed = true;
                            break;
                        }
                        // var folds to 1: drops out
                    } else {
                        folded.push(v - L);
                    }
                }
                if killed {
                    continue;
                }
                toggle(&mut acc, pack_mono(&folded));
            }
            let mut monos: Vec<Vec<u32>> = acc.into_keys().map(|p| unpack_mono(p)).collect();
            monos.sort();
            monos
        })
        .collect()
}

// ── Evaluation helpers ─────────────────────────────────────────────

fn eval_eq(eq: &[Vec<u32>], bits: &[bool]) -> bool {
    let mut acc = false;
    for mono in eq {
        if mono.iter().all(|&v| bits[v as usize]) {
            acc = !acc;
        }
    }
    acc
}

fn bits24(x1: u64, x2: u64, x3: u64) -> [bool; 24] {
    let mut b = [false; 24];
    for j in 0..L {
        b[j as usize] = (x1 >> j) & 1 == 1;
        b[(L + j) as usize] = (x2 >> j) & 1 == 1;
        b[(2 * L + j) as usize] = (x3 >> j) & 1 == 1;
    }
    b
}

fn deg_hist(eqs: &[Vec<Vec<u32>>]) -> (HashMap<usize, usize>, usize, usize, usize) {
    let mut hist: HashMap<usize, usize> = HashMap::new();
    let mut zero = 0usize;
    let mut unsat_cert = 0usize;
    let mut dmin = usize::MAX;
    for eq in eqs {
        if eq.is_empty() {
            zero += 1;
            continue;
        }
        let d = eq.iter().map(|m| m.len()).max().unwrap_or(0);
        if eq.len() == 1 && eq[0].is_empty() {
            unsat_cert += 1; // the equation "1 = 0"
        }
        dmin = dmin.min(d);
        *hist.entry(d).or_insert(0) += 1;
    }
    (hist, zero, unsat_cert, dmin)
}

// ── Measurement with a memory cap ──────────────────────────────────

struct VariantResult {
    measured: Vec<Value>,
    skipped_memory: Vec<Value>,
    fall_degree: Option<u32>,
    fall_structural: Option<u32>,
    deg_min: usize,
}

fn measure_variant(
    eqs: &[Vec<Vec<u32>>],
    num_vars: u32,
    d_max: u32,
    cap_bytes: u64,
) -> VariantResult {
    let (_hist, _zero, _unsat, dmin) = deg_hist(eqs);
    // Find the largest contiguous [D_MIN, d_eff] whose dense memory
    // estimate fits under the cap.
    let mut d_eff = D_MIN - 1;
    let mut skipped = Vec::new();
    for d in D_MIN..=d_max {
        let est = macaulay_memory_estimate_monomial(eqs, num_vars, d);
        if est > cap_bytes {
            // The estimate grows with d; once over the cap, stay over.
            for dd in d..=d_max {
                let est_dd = macaulay_memory_estimate_monomial(eqs, num_vars, dd);
                skipped.push(json!({ "degree": dd, "estimate_bytes": est_dd, "reason": "memory_cap" }));
            }
            break;
        }
        d_eff = d;
    }
    let mut measured = Vec::new();
    let mut fall = None;
    if d_eff >= D_MIN {
        let (per_deg, f) = measure_monomial_system(eqs, num_vars, D_MIN, d_eff);
        fall = f;
        for m in &per_deg {
            measured.push(json!({
                "degree": m.degree,
                "rows": m.rows_constructed,
                "cols": m.cols,
                "rank": m.rank,
                "rank_generic": m.rank_generic,
                "fall_signal": m.fall_signal,
            }));
        }
    }
    // A fall at D ≥ 2·δ_min can be a Koszul syzygy f_i·f_j − f_j·f_i
    // (δ_i + δ_j ≥ 2·δ_min); only falls strictly below are structural.
    let fall_structural = fall.filter(|d| (*d as usize) < 2 * dmin);
    VariantResult {
        measured,
        skipped_memory: skipped,
        fall_degree: fall,
        fall_structural,
        deg_min: dmin,
    }
}

fn variant_json(v: &VariantResult) -> Value {
    json!({
        "measured": v.measured,
        "skipped_memory": v.skipped_memory,
        "fall_degree": v.fall_degree,
        "fall_structural": v.fall_structural,
        "deg_min": v.deg_min,
    })
}

// ── Main ───────────────────────────────────────────────────────────

fn main() {
    // ── args ──
    let argv: Vec<String> = std::env::args().collect();
    let mut draws: u32 = 16;
    let mut seed: u64 = 0x54FF_D518;
    let mut dmax_full: u32 = 7;
    let mut dmax_fix: u32 = 9;
    let mut mem_cap_gb: f64 = 1.2;
    let mut out = String::from("ffd_s4_subspace_summary.json");
    let mut i = 1;
    while i < argv.len() {
        match argv[i].as_str() {
            "--draws" => { i += 1; draws = argv[i].parse().expect("--draws"); }
            "--seed" => { i += 1; seed = parse_u64(&argv[i]); }
            "--dmax-full" => { i += 1; dmax_full = argv[i].parse().expect("--dmax-full"); }
            "--dmax-fix" => { i += 1; dmax_fix = argv[i].parse().expect("--dmax-fix"); }
            "--mem-cap-gb" => { i += 1; mem_cap_gb = argv[i].parse().expect("--mem-cap-gb"); }
            "--out" => { i += 1; out = argv[i].clone(); }
            other => panic!("unknown argument {other}"),
        }
        i += 1;
    }
    let cap_bytes = (mem_cap_gb * (1u64 << 30) as f64) as u64;
    assert!(draws >= 1, "--draws must be >= 1");
    let t_start = Instant::now();

    println!("ffd_s4_subspace: l={L} n={N} draws={draws} seed=0x{seed:X} D∈[{D_MIN},{dmax_full}]×[{D_MIN},{dmax_fix}] cap={mem_cap_gb}GB");

    // ── gate 1: sanity ──
    let sanity = sanity_s3();

    // ── field setup (identical to the recorded baseline arm) ──
    let irr = IrreduciblePoly { degree: N, low_terms: LOW_TERMS.to_vec() };
    let gf = Gf2::new(&irr);
    let b = F2mElement::one(N);

    // Correspondence images are target-independent: build once from a
    // probe system, and verify independence on a second probe.
    let t0 = Instant::now();
    let probe1 = weil_descend_s4(N, L, &irr, &b, &gf.to_element(1));
    let singles = build_singles(&probe1);
    let probe2 = weil_descend_s4(N, L, &irr, &b, &gf.to_element(0xABCD));
    let singles2 = build_singles(&probe2);
    assert_eq!(singles, singles2, "correspondence images must not depend on x_R");
    let build_probe_ms = t0.elapsed().as_secs_f64() * 1e3;
    println!("probe build: {:.0} ms each; e-vars: {}; x_R-independence ✓", build_probe_ms, singles.len());

    let mut pair_cache: HashMap<u64, Vec<u64>> = HashMap::new();
    let mut rng_main = SplitMix64::new(seed);
    let mut draw_records = Vec::new();
    let mut decomposable_count = 0usize;
    let mut decompose_times = Vec::new();

    for draw in 0..draws {
        let t_draw = Instant::now();
        // target
        let mut xr = rng_main.next_u64() & ((1u64 << N) - 1);
        while xr == 0 {
            xr = rng_main.next_u64() & ((1u64 << N) - 1);
        }
        let mut rng_draw = SplitMix64::new(rng_main.next_u64());

        // baseline decision on this host (informational; also gives the
        // decomposability attestation and, when present, a witness).
        let t_dec = Instant::now();
        let witness = decompose(xr, L, &gf);
        let decompose_ms = t_dec.elapsed().as_secs_f64() * 1e3;
        decompose_times.push(decompose_ms);
        if witness.is_some() {
            decomposable_count += 1;
        }

        // build + eliminate
        let t_b = Instant::now();
        let x_r_elt = gf.to_element(xr);
        let sys = weil_descend_s4(N, L, &irr, &b, &x_r_elt);
        let build_ms = t_b.elapsed().as_secs_f64() * 1e3;
        let t_e = Instant::now();
        let eqs24 = eliminate_e(&sys, &singles, &mut pair_cache);
        let elim_ms = t_e.elapsed().as_secs_f64() * 1e3;

        // gate 2: elimination exactness on random points of V³
        let t_v = Instant::now();
        let mut mismatches = 0usize;
        for _ in 0..64 {
            let xm1 = rng_draw.next_u64() & 0xFF;
            let xm2 = rng_draw.next_u64() & 0xFF;
            let xm3 = rng_draw.next_u64() & 0xFF;
            let bits = bits24(xm1, xm2, xm3);
            let sys_sat = eqs24.iter().all(|eq| !eval_eq(eq, &bits));
            let (e1, e2, e3) = elementary_symmetric_3(
                &gf.to_element(xm1),
                &gf.to_element(xm2),
                &gf.to_element(xm3),
                &irr,
            );
            let ref_sat = symmetrised_s4_eval(&e1, &e2, &e3, &x_r_elt, &irr).is_zero();
            if sys_sat != ref_sat {
                mismatches += 1;
            }
        }
        let verify_ms = t_v.elapsed().as_secs_f64() * 1e3;
        assert_eq!(
            mismatches, 0,
            "eliminated system disagrees with symmetrised_s4_eval on {mismatches}/64 points (draw {draw}, x_R=0x{xr:X})"
        );

        // X₁-fixed fold
        let x1_fix = rng_draw.next_u64() & 0xFF;
        let eqs16 = fold_x1(&eqs24, x1_fix);

        // gate 3: witness vanishing (both variants)
        let mut witness_checked = false;
        if let Some([w1, w2, w3]) = witness {
            let bits = bits24(w1, w2, w3);
            assert!(
                eqs24.iter().all(|eq| !eval_eq(eq, &bits)),
                "witness decomposition does not satisfy the eliminated system (draw {draw})"
            );
            let eqs16w = fold_x1(&eqs24, w1);
            let mut bits16 = [false; 16];
            for j in 0..L {
                bits16[j as usize] = (w2 >> j) & 1 == 1;
                bits16[(L + j) as usize] = (w3 >> j) & 1 == 1;
            }
            assert!(
                eqs16w.iter().all(|eq| !eval_eq(eq, &bits16)),
                "witness does not satisfy the X₁-fixed fold at its own X₁ (draw {draw})"
            );
            witness_checked = true;
        }

        // measure both variants
        let t_mf = Instant::now();
        let full = measure_variant(&eqs24, 3 * L, dmax_full, cap_bytes);
        let meas_full_ms = t_mf.elapsed().as_secs_f64() * 1e3;
        let t_mx = Instant::now();
        let fix = measure_variant(&eqs16, 2 * L, dmax_fix, cap_bytes);
        let meas_fix_ms = t_mx.elapsed().as_secs_f64() * 1e3;

        let (hist_full, zero_full, unsat_full, _) = deg_hist(&eqs24);
        let (hist_fix, zero_fix, unsat_fix, _) = deg_hist(&eqs16);

        println!(
            "draw {draw:>2}: x_R=0x{xr:06X} decomp={} ({decompose_ms:.0} ms) x1_fix=0x{x1_fix:02X} \
             deg_full={:?} deg_fix={:?} fall_full={:?} fall_fix={:?} ({} ms total)",
            witness.is_some(),
            hist_full,
            hist_fix,
            full.fall_degree,
            fix.fall_degree,
            t_draw.elapsed().as_secs_f64() * 1e3
        );

        draw_records.push(json!({
            "draw": draw,
            "xr": format!("0x{:06X}", xr),
            "decomposable": witness.is_some(),
            "witness": witness.map(|w| [format!("0x{:02X}", w[0]), format!("0x{:02X}", w[1]), format!("0x{:02X}", w[2])]),
            "decompose_ms": decompose_ms,
            "x1_fix": format!("0x{:02X}", x1_fix),
            "build_ms": build_ms,
            "elim_ms": elim_ms,
            "verify": { "points": 64, "mismatches": mismatches, "verify_ms": verify_ms, "witness_checked": witness_checked },
            "deg_hist_full": hist_full.iter().map(|(k, v)| (k.to_string(), json!(v))).collect::<serde_json::Map<_, _>>(),
            "deg_hist_fix": hist_fix.iter().map(|(k, v)| (k.to_string(), json!(v))).collect::<serde_json::Map<_, _>>(),
            "zero_eqs_full": zero_full,
            "zero_eqs_fix": zero_fix,
            "unsat_cert_eqs_full": unsat_full,
            "unsat_cert_eqs_fix": unsat_fix,
            "meas_full_ms": meas_full_ms,
            "meas_fix_ms": meas_fix_ms,
            "full": variant_json(&full),
            "fix": variant_json(&fix),
        }));
    }

    // ── gate 3 coverage: at least one witness-checked decomposable target ──
    let mut extras = Vec::new();
    if decomposable_count == 0 {
        println!("no decomposable draw in {draws} — running extra draws for the witness gate");
        let mut rng_x = SplitMix64::new(seed ^ 0xA5A5_5A5A_A5A5_5A5A);
        for k in 0..16u32 {
            let mut xr = rng_x.next_u64() & ((1u64 << N) - 1);
            while xr == 0 {
                xr = rng_x.next_u64() & ((1u64 << N) - 1);
            }
            if let Some([w1, w2, w3]) = decompose(xr, L, &gf) {
                let sys = weil_descend_s4(N, L, &irr, &b, &gf.to_element(xr));
                let eqs24 = eliminate_e(&sys, &singles, &mut pair_cache);
                let bits = bits24(w1, w2, w3);
                assert!(
                    eqs24.iter().all(|eq| !eval_eq(eq, &bits)),
                    "extra-draw witness does not satisfy the eliminated system"
                );
                let eqs16w = fold_x1(&eqs24, w1);
                let mut bits16 = [false; 16];
                for j in 0..L {
                    bits16[j as usize] = (w2 >> j) & 1 == 1;
                    bits16[(L + j) as usize] = (w3 >> j) & 1 == 1;
                }
                assert!(
                    eqs16w.iter().all(|eq| !eval_eq(eq, &bits16)),
                    "extra-draw witness does not satisfy the X₁-fixed fold"
                );
                println!("   extra draw {k}: x_R=0x{xr:06X} witness-checked ✓");
                extras.push(json!({
                    "extra": k,
                    "xr": format!("0x{:06X}", xr),
                    "witness": [format!("0x{:02X}", w1), format!("0x{:02X}", w2), format!("0x{:02X}", w3)],
                    "witness_checked": true,
                }));
                break;
            }
        }
        assert!(
            !extras.is_empty(),
            "could not find a decomposable target in 16 extra draws — witness gate unmet"
        );
    }

    // ── summary ──
    let falls_full: Vec<Option<u32>> = draw_records
        .iter().map(|r| r["full"]["fall_degree"].as_u64().map(|v| v as u32)).collect();
    let falls_fix: Vec<Option<u32>> = draw_records
        .iter().map(|r| r["fix"]["fall_degree"].as_u64().map(|v| v as u32)).collect();
    let struct_full: Vec<u32> = draw_records
        .iter().filter_map(|r| r["full"]["fall_structural"].as_u64().map(|v| v as u32)).collect();
    let struct_fix: Vec<u32> = draw_records
        .iter().filter_map(|r| r["fix"]["fall_structural"].as_u64().map(|v| v as u32)).collect();

    let ffd_stats = |falls: &[Option<u32>]| -> Value {
        let some: Vec<u32> = falls.iter().filter_map(|f| *f).collect();
        if some.is_empty() {
            json!({ "min": null, "max": null, "mean": null, "draws_with_fall": 0, "draws_censored": falls.len() })
        } else {
            let min = *some.iter().min().unwrap();
            let max = *some.iter().max().unwrap();
            let mean = some.iter().sum::<u32>() as f64 / some.len() as f64;
            json!({ "min": min, "max": max, "mean": mean, "draws_with_fall": some.len(), "draws_censored": falls.len() - some.len() })
        }
    };

    decompose_times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median_decompose_ms = decompose_times[decompose_times.len() / 2];

    let verdict = if !struct_fix.is_empty() {
        format!("L8_S4_SUBSPACE_X1FIXED_STRUCTURAL_FALL_D{}", struct_fix.iter().min().unwrap())
    } else if !struct_full.is_empty() {
        format!("L8_S4_SUBSPACE_FULL_STRUCTURAL_FALL_D{}", struct_full.iter().min().unwrap())
    } else if falls_fix.iter().any(|f| f.is_some()) || falls_full.iter().any(|f| f.is_some()) {
        String::from("L8_S4_SUBSPACE_FALLS_ALL_KOSZUL_AMBIGUOUS")
    } else {
        format!("L8_S4_SUBSPACE_NO_FALL_CENSORED_FULL_D{dmax_full}_FIX_D{dmax_fix}")
    };

    let per_degree_table = |records: &[Value], variant: &str| -> Vec<Value> {
        let mut by_deg: HashMap<u64, Vec<(u64, u64, u64)>> = HashMap::new(); // deg -> (rows, cols, rank)
        for r in records {
            for m in r[variant]["measured"].as_array().unwrap_or(&Vec::new()).clone() {
                by_deg.entry(m["degree"].as_u64().unwrap())
                    .or_default()
                    .push((
                        m["rows"].as_u64().unwrap(),
                        m["cols"].as_u64().unwrap(),
                        m["rank"].as_u64().unwrap(),
                    ));
            }
        }
        let mut degs: Vec<u64> = by_deg.keys().copied().collect();
        degs.sort_unstable();
        degs.iter().map(|d| {
            let v = &by_deg[d];
            let rows_min = v.iter().map(|t| t.0).min().unwrap();
            let rows_max = v.iter().map(|t| t.0).max().unwrap();
            let rank_min = v.iter().map(|t| t.2).min().unwrap();
            let rank_max = v.iter().map(|t| t.2).max().unwrap();
            let rank_mean = v.iter().map(|t| t.2).sum::<u64>() as f64 / v.len() as f64;
            let deficit_rows_mean = rows_max as f64 - rank_mean;
            json!({
                "degree": d, "draws": v.len(),
                "rows_min": rows_min, "rows_max": rows_max, "cols": v[0].1,
                "rank_min": rank_min, "rank_max": rank_max, "rank_mean": rank_mean,
                "mean_deficit_vs_rows": deficit_rows_mean,
            })
        }).collect()
    };

    let created_at = std::process::Command::new("date")
        .args(["-u", "+%Y-%m-%dT%H:%M:%SZ"])
        .output()
        .ok()
        .and_then(|o| if o.status.success() { String::from_utf8(o.stdout).ok() } else { None })
        .map(|s| s.trim().to_string())
        .unwrap_or_else(|| {
            format!("unix_{}", SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs())
        });

    let summary = json!({
        "schema": "ffd_s4_subspace_v1",
        "created_at": created_at,
        "purpose": "FFD min/max/mean over >=16 draws on the subspace-restricted Weil-descended S4 system at l=8 — the 'FFD required' gate of the binary decomposition beat (docs/ic/boundary_targets.json). Measurement only: no oracle, no speedup claim.",
        "config": {
            "n": N, "l": L, "d_min": D_MIN,
            "dmax_full": dmax_full, "dmax_fix": dmax_fix,
            "draws": draws, "seed": format!("0x{seed:X}"),
            "extra_witness_seed": format!("0x{:X}", seed ^ 0xA5A5_5A5A_A5A5_5A5A),
            "mem_cap_gb": mem_cap_gb,
            "irr": { "degree": N, "low_terms": LOW_TERMS, "source": "recorded baseline arm decomp_bench_to_l8.stdout.txt" },
            "curve": "Koblitz b=1 (y^2+xy = x^3+x^2+1)",
            "factor_base": "V = <1,z,...,z^(l-1)>, same as pairs-and-solve baseline",
        },
        "sanity_s3": sanity,
        "correspondence_probe": {
            "build_ms_each": build_probe_ms,
            "x_r_independence_checked": true,
            "e_vars": singles.len(),
        },
        "draws": draw_records,
        "extra_witness_draws": extras,
        "summary": {
            "draws_total": draws,
            "decomposable_draws": decomposable_count,
            "ffd_full": ffd_stats(&falls_full),
            "ffd_fix": ffd_stats(&falls_fix),
            "structural_falls_full": struct_full,
            "structural_falls_fix": struct_fix,
            "median_decompose_ms_pairs_baseline_recheck": median_decompose_ms,
            "per_degree_full": per_degree_table(&draw_records, "full"),
            "per_degree_fix": per_degree_table(&draw_records, "fix"),
            "pair_cache_size": pair_cache.len(),
            "verdict": verdict,
        },
        "timing_s": { "total": t_start.elapsed().as_secs_f64() },
        "honesty": {
            "koszul_note": "falls at D >= 2*deg_min may be Koszul syzygies; fall_structural requires D < 2*deg_min",
            "truncation_note": "Macaulay truncation at D and multipliers <= D-deg per equation, identical convention to the quadratic harness (sanity gate)",
            "no_oracle_claim": true,
            "scope": "l=8 (n=24) only; says nothing about asymptotics",
        },
    });

    let json_str = serde_json::to_string_pretty(&summary).expect("serialize");
    std::fs::write(&out, &json_str).expect("write out");

    println!("\n── summary ──");
    println!("decomposable draws: {decomposable_count}/{draws}");
    println!("ffd_full: {}", summary["summary"]["ffd_full"]);
    println!("ffd_fix:  {}", summary["summary"]["ffd_fix"]);
    println!("structural falls full: {struct_full:?}  fix: {struct_fix:?}");
    println!("median decompose (pairs recheck): {median_decompose_ms:.1} ms");
    println!("pair cache entries: {}", pair_cache.len());
    println!("verdict: {verdict}");
    println!("total: {:.1}s → {out}", t_start.elapsed().as_secs_f64());
}

fn parse_u64(s: &str) -> u64 {
    if let Some(hex) = s.strip_prefix("0x").or_else(|| s.strip_prefix("0X")) {
        u64::from_str_radix(hex, 16).expect("hex seed")
    } else {
        s.parse().expect("decimal seed")
    }
}
