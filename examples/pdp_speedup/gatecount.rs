//! PS-1/PS-2: gate counts of the root-test formulations, and execution of the
//! bit-sliced group-law pipeline on 64 lanes against the exact reference.

use crate::bitslice::*;
use crate::bspairs::*;
use crate::common::*;
use crate::group::*;
use std::sync::atomic::{AtomicU64, Ordering};

fn arg<T: std::str::FromStr>(args: &[String], name: &str) -> Option<T> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
}

fn counted<R>(f: impl FnOnce() -> R) -> (R, Gates) {
    gates_reset();
    let r = f();
    (r, gates_get())
}

/// Gate tallies of the bit-sliced primitives and of one task of the
/// group-law pipeline, for the field of rung `(n, l)`.
pub fn cmd_gates(args: &[String]) {
    let n: u32 = arg(args, "--n").unwrap_or(0);
    let l: u32 = arg(args, "--l").expect("--l");
    let n = if n == 0 { 3 * l } else { n };
    let t_batch: usize = arg(args, "--batch").unwrap_or(64);
    let rung = Rung::with(n, l);
    let f = BsField::new(&rung.gf);
    let nn = n as usize;
    let e = zero_elt::<Cnt>();

    let (_, g_mul) = counted(|| f.mul(&e, &e));
    let (_, g_sqr) = counted(|| f.sqr(&e));
    let (_, g_inv) = counted(|| f.inv(&e));
    let (_, g_zero) = counted(|| f.zero_mask(&e));
    let (_, g_in) = counted(|| f.in_subspace_mask(&e, l as usize));
    let row = RowPlanes::<Cnt>::new(nn, 0, 0, 0, 0);
    let (a, g_a) = counted(|| stage_a(&f, &row, &e));
    // One task of the pipeline given the inverses.
    let (_, g_c) = counted(|| stage_c(&f, &row, &e, &e, &a, &e, &e, l as usize));
    // Montgomery across a batch of tasks and two chains (U and V).
    let ds: Vec<E<Cnt>> = vec![e; t_batch];
    let refs: Vec<&E<Cnt>> = ds.iter().collect();
    let (_, g_mont) = counted(|| montgomery(&f, &refs));
    let per_task_mont_fused = 2 * g_mont.fused();
    let g_task_fused = g_a.fused() + g_c.fused() + per_task_mont_fused / t_batch as u64;
    let g_task_unfused = g_a.unfused() + g_c.unfused() + 2 * g_mont.unfused() / t_batch as u64;

    // The baseline (gcd) formulation, derived from the counted field operations.
    let g_prod = (nn * nn) as u64; // n² fused XOR-AND gates, unreduced product
    let g_red = g_mul.fused() - g_prod; // reduction of a 2n−1 product
    let counts = std::fs::read_to_string(
        "research/pdp_oracle_speedup_20261011/baseline_v1/opcount/counts.jsonl",
    )
    .ok()
    .and_then(|t| {
        t.lines()
            .filter_map(|ln| serde_json::from_str::<serde_json::Value>(ln).ok())
            .find(|v| v["l"].as_u64() == Some(l as u64))
    });
    let f1 = counts.map(|v| {
        let pp = |k: &str| v["per_pair"][k].as_f64().unwrap();
        let mul = pp("mul") + pp("batch_mul");
        let products = pp("clmul") + pp("sqr_wide");
        let reduces = pp("reduce");
        let sqrs = pp("sqr") + pp("sqr_k");
        // Standalone products are summed before one reduction: (products − reduces)
        // extra (2n−1)-plane XORs; `sqr_wide` is a wiring (free) squaring.
        let sum_xors = (products - reduces).max(0.0) * (2 * nn - 1) as f64;
        let g = mul * g_mul.fused() as f64
            + sqrs * g_sqr.fused() as f64
            + pp("clmul") * g_prod as f64
            + reduces * g_red as f64
            + sum_xors;
        (g, mul, pp("clmul"), reduces, sqrs, sum_xors)
    });

    println!(
        "{{\"n\":{n},\"l\":{l},\"mul\":{},\"mul_product_part\":{g_prod},\"mul_reduction_part\":{g_red},\"sqr\":{},\"inv\":{},\"zero_mask\":{},\"in_subspace_mask\":{},\"stage_a\":{},\"stage_c\":{},\"montgomery_batch\":{t_batch},\"montgomery_one_chain\":{},\"task_fused\":{g_task_fused},\"task_unfused\":{g_task_unfused},\"f1_gcd_derived\":{}}}",
        g_mul.fused(),
        g_sqr.fused(),
        g_inv.fused(),
        g_zero.fused(),
        g_in.fused(),
        g_a.fused(),
        g_c.fused(),
        g_mont.fused(),
        match f1 {
            Some((g, mul, clm, red, sq, sx)) => format!(
                "{{\"gates\":{:.0},\"mul_eq\":{mul:.2},\"clmul\":{clm:.1},\"reduce\":{red:.1},\"sqr_eq\":{sq:.2},\"sum_xor_gates\":{sx:.0}}}",
                g
            ),
            None => "null".into(),
        }
    );
}

/// Run the 64-lane pipeline over every row against one target and compare
/// each pair's hit with the exact reference.  Returns (pairs, mismatches,
/// exceptional pairs).
fn check_target(rung: &Rung, base: &[(u64, u64)], r: (u64, u64)) -> (u64, u64, u64) {
    let gf = &rung.gf;
    let n = rung.n as usize;
    let l = rung.l as usize;
    let span = rung.span();
    let f = BsField::new(gf);
    let m = base.len();
    let blocks = m.div_ceil(64);
    // Pre-transposed base blocks.
    let mut xb: Vec<E<u64>> = Vec::new();
    let mut yb: Vec<E<u64>> = Vec::new();
    for b in 0..blocks {
        let mut xs = [0u64; 64];
        let mut ys = [0u64; 64];
        for j in 0..64 {
            if let Some(&(x, y)) = base.get(b * 64 + j) {
                xs[j] = x;
                ys[j] = y;
            }
        }
        xb.push(to_planes_u64(&xs, n));
        yb.push(to_planes_u64(&ys, n));
    }
    let rp = Pt::Aff(r.0, r.1);
    let nr = neg(rp);
    let (mut pairs, mut bad, mut exc) = (0u64, 0u64, 0u64);
    for i1 in 0..m {
        let p1 = Pt::Aff(base[i1].0, base[i1].1);
        let (Pt::Aff(xu, yu), Pt::Aff(xv, yv)) = (add(gf, p1, rp), add(gf, p1, nr)) else {
            continue; // x_R in V: excluded by the caller
        };
        let row = RowPlanes::<u64>::new(n, xu, yu, xv, yv);
        let first = i1 / 64;
        let stages: Vec<StageA<u64>> = (first..blocks).map(|b| stage_a(&f, &row, &xb[b])).collect();
        let dus: Vec<&E<u64>> = stages.iter().map(|s| &s.du).collect();
        let dvs: Vec<&E<u64>> = stages.iter().map(|s| &s.dv).collect();
        let ju = montgomery(&f, &dus);
        let jv = montgomery(&f, &dvs);
        for (bi, b) in (first..blocks).enumerate() {
            let hit = stage_c(&f, &row, &xb[b], &yb[b], &stages[bi], &ju[bi], &jv[bi], l);
            for j in 0..64 {
                let i2 = b * 64 + j;
                if i2 < i1 || i2 >= m {
                    continue;
                }
                pairs += 1;
                let reference = thirds(gf, base[i1], base[i2], r)
                    .into_iter()
                    .flatten()
                    .any(|t| t < span);
                let is_exc = (stages[bi].exceptional >> j) & 1 == 1;
                if is_exc {
                    exc += 1; // recomputed exactly by the caller; not a mismatch
                    continue;
                }
                if ((hit >> j) & 1 == 1) != reference {
                    bad += 1;
                }
            }
        }
    }
    (pairs, bad, exc)
}

/// `bs-check --n N --l L [--targets K]`: the bit-sliced pipeline against the
/// exact reference, per pair, for every rational target (or the first `K`).
pub fn cmd_bscheck(args: &[String]) {
    let n: u32 = arg(args, "--n").expect("--n");
    let l: u32 = arg(args, "--l").expect("--l");
    let k: usize = arg(args, "--targets").unwrap_or(usize::MAX);
    let threads: usize = arg(args, "--threads").unwrap_or(3);
    let rung = Rung::with(n, l);
    let gf = &rung.gf;
    let solver = ArtinSchreier::new(gf);
    let base = base_points(gf, &solver, l);
    let targets: Vec<(u64, u64)> = (1..(1u64 << n))
        .filter(|&x| !degenerate_target(gf, l, x))
        .filter_map(|x| lift(gf, &solver, x).map(|y| (x, y)))
        .take(k)
        .collect();
    let next = AtomicU64::new(0);
    let tot = std::sync::Mutex::new((0u64, 0u64, 0u64));
    std::thread::scope(|s| {
        for _ in 0..threads {
            s.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::Relaxed) as usize;
                if i >= targets.len() {
                    break;
                }
                let (p, b, e) = check_target(&rung, &base, targets[i]);
                let mut t = tot.lock().unwrap();
                t.0 += p;
                t.1 += b;
                t.2 += e;
            });
        }
    });
    let t = tot.lock().unwrap();
    println!(
        "{{\"n\":{n},\"l\":{l},\"base_points\":{},\"targets\":{},\"pairs_checked\":{},\"mismatches\":{},\"exceptional_pairs\":{}}}",
        base.len(),
        targets.len(),
        t.0,
        t.1,
        t.2
    );
    if t.1 != 0 {
        std::process::exit(1);
    }
}

/// `export-program --l L --tasks T`: record the straight-line program of a
/// batch of `T` tasks on the recording lane, eliminate dead operations, hash
/// it, and check that interpreting it on 64 lanes equals running the same
/// code directly on `u64` lanes (and that its live gate count equals the
/// counting lane's).
pub fn cmd_export(args: &[String]) {
    let l: u32 = arg(args, "--l").expect("--l");
    let t: usize = arg(args, "--tasks").unwrap_or(8);
    let n = 3 * l;
    let rung = Rung::with(n, l);
    let f = BsField::new(&rung.gf);
    let nn = n as usize;
    // Inputs: per task xu, yu, xv, yv, x2, y2 planes (6n), numbered in order.
    let mut next = 0u32;
    let mut take = |_: usize| -> E<Tape> {
        let mut e = zero_elt::<Tape>();
        for plane in e.iter_mut().take(nn) {
            *plane = Tape::input(next);
            next += 1;
        }
        e
    };
    tape_reset();
    let mut rows: Vec<RowPlanes<Tape>> = Vec::new();
    let mut x2s = Vec::new();
    let mut y2s = Vec::new();
    for k in 0..t {
        rows.push(RowPlanes {
            xu: take(k),
            yu: take(k),
            xv: take(k),
            yv: take(k),
        });
        x2s.push(take(k));
        y2s.push(take(k));
    }
    let n_inputs = next as usize;
    let res = batch::<Tape>(&f, &rows, &x2s, &y2s, l as usize);
    let outputs: Vec<Tape> = res.iter().flat_map(|&(h, e)| [h, e]).collect();
    let ops = tape_take();
    let (live, outs) = eliminate_dead(&ops, &outputs);
    let bytes = program_bytes(&live, &outs);
    let sha = hex(&crypto_lib::hash::sha256::sha256(&bytes));
    if let Some(path) = arg::<String>(args, "--out") {
        // The program as data (`program_bytes`), for the generators of other
        // back ends (the GPU kernel of PS-5), and a sidecar that says what the
        // numbered inputs and the outputs are.
        std::fs::write(&path, &bytes).unwrap();
        std::fs::write(
            format!("{path}.json"),
            format!(
                "{{\"format\":\"pdp-slp/1\",\"n\":{n},\"l\":{l},\"tasks\":{t},\"inputs\":{n_inputs},\"input_order\":\"per task k: xu, yu, xv, yv, x2, y2, each n planes with plane 0 first; input number = k*6n + field*n + plane\",\"outputs\":\"per task k: hit mask, exceptional mask (two wires each, in task order)\",\"encoding\":\"per op a tag byte (0 input,1 zero,2 ones,3 not,4 xor,5 and,6 or,7 xand=a^(b&c)) then u32 LE operands; 255 then the output wires as u32 LE\",\"sha256\":\"{sha}\"}}\n"
            ),
        )
        .unwrap();
    }
    let mut g = Gates::default();
    for op in &live {
        match op {
            TapeOp::Xor(..) => g.xor += 1,
            TapeOp::And(..) => g.and += 1,
            TapeOp::Or(..) => g.or += 1,
            TapeOp::Not(..) => g.not += 1,
            TapeOp::Xand(..) => g.xand += 1,
            _ => {}
        }
    }
    // The counting lane on the same code.
    let rc: Vec<RowPlanes<Cnt>> = (0..t)
        .map(|_| RowPlanes::<Cnt>::new(nn, 0, 0, 0, 0))
        .collect();
    let xc = vec![zero_elt::<Cnt>(); t];
    gates_reset();
    let _ = batch::<Cnt>(&f, &rc, &xc, &xc, l as usize);
    let counted = gates_get();
    // Interpreter against direct execution on random planes.
    let mut rng = SplitMix64(20261011 + l as u64);
    let mut ok = true;
    for _ in 0..4 {
        let inputs: Vec<u64> = (0..n_inputs).map(|_| rng.next()).collect();
        let want = {
            let mut at = 0usize;
            let mut take_u = |_: usize| -> E<u64> {
                let mut e = zero_elt::<u64>();
                for plane in e.iter_mut().take(nn) {
                    *plane = inputs[at];
                    at += 1;
                }
                e
            };
            let mut rows_u = Vec::new();
            let mut x2u = Vec::new();
            let mut y2u = Vec::new();
            for k in 0..t {
                rows_u.push(RowPlanes {
                    xu: take_u(k),
                    yu: take_u(k),
                    xv: take_u(k),
                    yv: take_u(k),
                });
                x2u.push(take_u(k));
                y2u.push(take_u(k));
            }
            batch::<u64>(&f, &rows_u, &x2u, &y2u, l as usize)
                .into_iter()
                .flat_map(|(h, e)| [h, e])
                .collect::<Vec<u64>>()
        };
        let got = eval_program(&live, &outs, &inputs);
        ok &= want == got;
    }
    println!(
        "{{\"n\":{n},\"l\":{l},\"tasks\":{t},\"inputs\":{n_inputs},\"outputs\":{},\"recorded_ops\":{},\"live_ops\":{},\"live_gates_fused\":{},\"live_xor\":{},\"live_and\":{},\"live_or\":{},\"live_not\":{},\"live_xand\":{},\"counted_gates_fused\":{},\"program_bytes\":{},\"sha256\":\"{sha}\",\"interpreter_equals_direct\":{ok}}}",
        outs.len(),
        ops.len(),
        live.len(),
        g.fused(),
        g.xor,
        g.and,
        g.or,
        g.not,
        g.xand,
        counted.fused(),
        bytes.len()
    );
    if !ok {
        std::process::exit(1);
    }
}
