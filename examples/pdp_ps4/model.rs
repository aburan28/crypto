//! The bit-operation model that puts the two fields in one unit, and the
//! `control` subcommand (the subspace-leaf kernel on the PS-0 corpus).
//!
//! **The model is a convention, not a measurement.**  Each field operation is
//! priced at its schoolbook two-input gate count (AND and XOR count one each):
//!
//! | operation | normal basis (type II) | polynomial basis, modulus with `k` low terms |
//! |:--|:--|:--|
//! | product | `n²` AND + `2n² − 2n` XOR | `n²` AND + `(n−1)²` XOR + `(n−1)·k` XOR |
//! | square | `0` (coordinate permutation) | `Σ_j (popcount(row_j) − 1)` of the squaring matrix |
//! | addition | `n` XOR | `n` XOR |
//! | leaf test | popcount `≤ w`: `2n` | `x₃ < 2^ℓ`: `n − ℓ − 1` OR |
//!
//! The normal-basis product has `2n² − n` partial terms (`γ_iγ_j = γ_{i+j} +
//! γ_{i−j}`, minus the `n` terms with `γ_0 = 0`), folded onto `n` outputs.  The
//! polynomial reduction is priced without cascades, which is a lower bound on
//! its cost and so errs against the claim that the normal-basis leaf is no
//! dearer.  Neither side uses Karatsuba or any other sub-quadratic product.
//! Counted field operations come from the instrumented kernels; the model only
//! converts them.

use crate::field::Counts;

#[derive(Clone, Copy, Debug)]
pub struct BitCosts {
    pub mul: f64,
    pub sqr: f64,
    pub add: f64,
    /// The leaf test of one candidate third.
    pub test: f64,
}

pub fn onb_costs(n: usize) -> BitCosts {
    let n = n as f64;
    BitCosts {
        mul: n * n + 2.0 * n * n - 2.0 * n,
        sqr: 0.0,
        add: n,
        test: 2.0 * n,
    }
}

/// `sqr_xors`: the XOR count of the squaring matrix, `k`: low terms of the modulus.
#[cfg(feature = "pdp-opcount")]
pub fn poly_costs(n: usize, k: usize, l: usize, sqr_xors: u64) -> BitCosts {
    let nf = n as f64;
    BitCosts {
        mul: nf * nf + (nf - 1.0) * (nf - 1.0) + (nf - 1.0) * k as f64,
        sqr: sqr_xors as f64,
        add: nf,
        test: (n - l).saturating_sub(1) as f64,
    }
}

/// Model bit operations for a batch of counted operations.
pub fn bitops(c: &Counts, tests: u64, k: &BitCosts) -> f64 {
    (c.mul + c.batch_mul) as f64 * k.mul
        + (c.sqr + c.sqr_k) as f64 * k.sqr
        + c.add as f64 * k.add
        + tests as f64 * k.test
}

#[cfg(feature = "pdp-opcount")]
pub fn cmd_control(args: &[String]) {
    use crate::control::{field_from_modulus, read_corpus, Control};
    use crate::{arg, arg_s};
    use crypto_lib::cryptanalysis::semaev_decomp::opcount;
    use std::io::Write;

    let l: u32 = arg(args, "--l").expect("--l");
    let path = arg_s(args, "--corpus").expect("--corpus");
    let n_uniform: usize = arg(args, "--targets").unwrap_or(6);
    let n_planted: usize = arg(args, "--planted").unwrap_or(2);
    let mut out: Box<dyn Write> = match arg_s(args, "--out") {
        Some(p) => Box::new(
            std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(p)
                .unwrap(),
        ),
        None => Box::new(std::io::stdout()),
    };
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let corpus = read_corpus(&text);
    assert_eq!(corpus.l, l, "corpus is for a different rung");
    let gf = field_from_modulus(corpus.n, corpus.modulus);
    // The squaring matrix of the polynomial basis, for the model.
    let mut sqr_xors = 0u64;
    {
        let cols: Vec<u64> = (0..corpus.n).map(|i| gf.sqr(1u64 << i)).collect();
        for j in 0..corpus.n {
            let pc = cols.iter().filter(|c| (*c >> j) & 1 == 1).count() as u64;
            sqr_xors += pc.saturating_sub(1);
        }
    }
    let kernel = Control::new(&gf, l);
    let chosen: Vec<&crate::control::CorpusTarget> = corpus
        .targets
        .iter()
        .filter(|t| !t.planted)
        .take(n_uniform)
        .chain(corpus.targets.iter().filter(|t| t.planted).take(n_planted))
        .collect();
    let mut tot = opcount::Counts::default();
    let (mut pairs, mut hits, mut agree, mut found_targets) = (0u64, 0u64, 0usize, 0usize);
    let t0 = std::time::Instant::now();
    for t in &chosen {
        opcount::reset();
        let r = kernel.run((t.xr, t.yr));
        let c = opcount::snapshot();
        tot.mul += c.mul;
        tot.sqr += c.sqr;
        tot.sqr_k += c.sqr_k;
        tot.inv += c.inv;
        tot.batch_calls += c.batch_calls;
        tot.batch_elems += c.batch_elems;
        tot.batch_mul += c.batch_mul;
        pairs += r.pairs;
        hits += r.hits;
        found_targets += usize::from(r.found);
        agree += usize::from(r.found == t.has_rational_root);
    }
    let ms = t0.elapsed().as_secs_f64() * 1e3;
    // The kernel's static per-pair tallies, read off control.rs line by line:
    // 16 field additions (2 denominators, 2 recomputed differences, 2 ordinate
    // sums, 2 for the shifted slopes, 4 x 2 for the candidates) plus 4
    // single-bit flips (`^ 1`, the curve constant `a = 1`, not priced), and 4
    // leaf tests.  The weight kernel's `add` counter tallies the same 16.
    let adds_per_pair = 16u64;
    let tests_per_pair = 4u64;
    let per = |v: u64| v as f64 / pairs as f64;
    let bc = poly_costs(
        corpus.n as usize,
        (corpus.modulus.count_ones() - 1) as usize,
        l as usize,
        sqr_xors,
    );
    let mut cnt = crate::field::Counts {
        mul: tot.mul,
        sqr: tot.sqr,
        sqr_k: tot.sqr_k,
        inv: tot.inv,
        batch_calls: tot.batch_calls,
        batch_elems: tot.batch_elems,
        batch_mul: tot.batch_mul,
        add: adds_per_pair * pairs,
        popcnt: 0,
        lookup: 0,
    };
    cnt.popcnt = 0;
    let model = bitops(&cnt, tests_per_pair * pairs, &bc);
    writeln!(
        out,
        "{{\"kind\":\"control\",\"l\":{l},\"n\":{},\"modulus\":\"0x{:x}\",\"targets\":{},\"base_liftable_x\":{},\"pairs\":{pairs},\"hits\":{hits},\"targets_found\":{found_targets},\"verdict_agrees_with_corpus\":{agree},\"counts\":{{\"mul\":{},\"sqr\":{},\"sqr_k\":{},\"inv\":{},\"batch_calls\":{},\"batch_elems\":{},\"batch_mul\":{},\"mulq\":{}}},\"per_pair\":{{\"mul\":{:.4},\"sqr\":{:.4},\"sqr_k\":{:.4},\"inv\":{:.5},\"batch_mul\":{:.4},\"mulq\":{:.4},\"add_static\":{adds_per_pair},\"tests_static\":{tests_per_pair}}},\"model\":{{\"mul_bits\":{:.1},\"sqr_bits\":{:.1},\"add_bits\":{:.1},\"test_bits\":{:.1},\"bitops_total\":{model:.1},\"bitops_per_pair\":{:.3}}},\"ms_unisolated\":{ms:.1}}}",
        corpus.n,
        corpus.modulus,
        chosen.len(),
        kernel.base.len(),
        tot.mul,
        tot.sqr,
        tot.sqr_k,
        tot.inv,
        tot.batch_calls,
        tot.batch_elems,
        tot.batch_mul,
        tot.mul + tot.batch_mul,
        per(tot.mul),
        per(tot.sqr),
        per(tot.sqr_k),
        per(tot.inv),
        per(tot.batch_mul),
        per(tot.mul + tot.batch_mul),
        bc.mul,
        bc.sqr,
        bc.add,
        bc.test,
        model / pairs as f64
    )
    .unwrap();
}

#[cfg(not(feature = "pdp-opcount"))]
pub fn cmd_control(_args: &[String]) {
    eprintln!("control needs the library's counters: build with --features pdp-opcount");
    std::process::exit(2);
}
