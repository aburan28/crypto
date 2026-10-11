//! **Scalar against batched pair-table construction and three-summand
//! probing**, same process, same inputs, interleaved.
//!
//! ```text
//! ic_batch_bench N A BASE REPS [TARGETS [phases]]
//!   `phases` instead times where a batched build's time goes (additions,
//!   canonicalisation with shift, hash-map inserts) and exits.
//!   BASE = scan:COLUMNS[:NORM]   compact-orbit-scan (raw-x cap 4e8)
//!        | subspace:CAP          Frobenius-invariant subspace of at most CAP points
//! ```
//!
//! For each repeat it builds the Frobenius-folded pair table with
//! [`FrobeniusPairTable::build_scalar`] and with [`FrobeniusPairTable::build`]
//! (order alternating), then runs `TARGETS` three-summand decompositions of
//! the same random subgroup points through [`decompose_mitm_frobenius_scalar`]
//! and [`decompose_mitm_frobenius`].  It asserts that the two tables have the
//! same entries and charged operations and that every decomposition and
//! charged operation count agrees, then prints the **minimum** and median
//! wall time of each over the repeats.  The minimum is the figure to read on
//! a loaded host: contention only ever adds time.
//!
//! This measures wall time, which operation counts do not see: the two paths
//! do the same charged work.  It is a speed measurement of one process on one
//! host, not an isolated-host receipt.
use std::time::Instant;

use crypto_lib::cryptanalysis::ic_boundary::{
    decompose_mitm_frobenius, decompose_mitm_frobenius_scalar, koblitz_factor_base,
    koblitz_instance, BinaryGroup, ColumnFold, CountedGroup, FactorBase, FrobeniusPairTable,
    GroupOps, OracleCounters,
};
use crypto_lib::cryptanalysis::ic_framework::plugins::CompactOrbitScanBase;
use crypto_lib::cryptanalysis::ic_framework::stages::{FactorBaseBuilder, InstanceCtx, Params};
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn main() -> Result<(), String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() < 4 {
        return Err("usage: ic_batch_bench N A BASE REPS [TARGETS]".into());
    }
    let n: u32 = args[0].parse().map_err(|_| "N")?;
    let a: u8 = args[1].parse().map_err(|_| "A")?;
    let reps: usize = args[3].parse().map_err(|_| "REPS")?;
    let targets: usize = args
        .get(4)
        .map_or(Ok(40), |s| s.parse())
        .map_err(|_| "TARGETS")?;
    let inst = koblitz_instance(a, n).ok_or("no such Koblitz instance")?;
    let g = BinaryGroup(&inst.fast);
    let spec: Vec<&str> = args[2].split(':').collect();
    let fb: FactorBase<FastPoint> = match spec[0] {
        "scan" => {
            let columns = spec.get(1).ok_or("scan:COLUMNS")?;
            let ctx = InstanceCtx {
                group: &g,
                generator: inst.generator,
                target: inst.generator,
                r: inst.r,
                cofactor: inst.cofactor,
                group_order: inst.group_order,
                name: inst.name.clone(),
                field_degree: Some(inst.n),
            };
            let mut p = Params::default();
            p.set("columns", *columns);
            p.set("raw_x_cap", "400000000");
            if let Some(norm) = spec.get(2) {
                p.set("norm", *norm);
            }
            CompactOrbitScanBase { instance: &inst }.build(&ctx, &p, &mut GroupOps::default())?
        }
        "subspace" => {
            let cap: usize = spec
                .get(1)
                .ok_or("subspace:CAP")?
                .parse()
                .map_err(|_| "CAP")?;
            let (frob, desc) =
                crypto_lib::cryptanalysis::ic_boundary::choose_koblitz_base(&inst, 1, cap)
                    .ok_or("no invariant subspace within the cap")?;
            koblitz_factor_base(&inst, &frob, ColumnFold::SignedFrobeniusOrbit, desc)
                .ok_or("base")?
        }
        other => return Err(format!("unknown base `{other}`")),
    };
    println!(
        "n = {n} K_{a}: {} points, {} columns; {reps} repeats, {targets} targets (m = 3)",
        fb.points.len(),
        fb.columns
    );

    if args.get(5).map(String::as_str) == Some("phases") {
        // Where the batched build's time goes: the additions, the
        // canonicalisation with its shift, and the hash-map inserts.
        use crypto_lib::cryptanalysis::koblitz_fast::{BatchScratch, FrobeniusCanon};
        use std::collections::HashMap;
        let canon = FrobeniusCanon::new(&inst.gf, inst.n).ok_or("normal basis")?;
        let reps: Vec<usize> = (0..fb.points.len()).step_by(2 * n as usize).collect();
        let mut scratch = BatchScratch::default();
        let (mut t_add, mut t_canon, mut t_map) = (0f64, 0f64, 0f64);
        let mut sums: Vec<FastPoint> = Vec::new();
        let mut keyed: Vec<(u64, u32)> = Vec::new();
        let mut map: HashMap<u64, (u32, u32)> =
            HashMap::with_capacity(reps.len() * fb.points.len() / 2);
        let mut pairs = 0u64;
        for (o, &rep) in reps.iter().enumerate() {
            let qs = &fb.points[..];
            sums.clear();
            let t = Instant::now();
            inst.fast
                .add_many(fb.points[rep], qs, &mut sums, &mut scratch);
            t_add += t.elapsed().as_secs_f64();
            keyed.clear();
            let t = Instant::now();
            for u in sums.iter().filter(|u| !u.infinity) {
                keyed.push(canon.canon_with_shift(u.x));
            }
            t_canon += t.elapsed().as_secs_f64();
            let t = Instant::now();
            for (j, &(key, shift)) in keyed.iter().enumerate() {
                map.entry(key)
                    .or_insert((o as u32, j as u32 | (shift << 26)));
            }
            t_map += t.elapsed().as_secs_f64();
            pairs += keyed.len() as u64;
        }
        let total = t_add + t_canon + t_map;
        println!(
            "build phases over {} representatives, {} pairs ({} entries): add_many {:.0} ms ({:.1} ns/pair, {:.0} %), canon_with_shift {:.0} ms ({:.1} ns/pair, {:.0} %), map insert {:.0} ms ({:.1} ns/pair, {:.0} %)",
            reps.len(), pairs, map.len(),
            t_add * 1e3, t_add * 1e9 / pairs as f64, 100.0 * t_add / total,
            t_canon * 1e3, t_canon * 1e9 / pairs as f64, 100.0 * t_canon / total,
            t_map * 1e3, t_map * 1e9 / pairs as f64, 100.0 * t_map / total,
        );
        return Ok(());
    }

    // Targets: random subgroup points, fixed across repeats.
    let mut rng = StdRng::seed_from_u64(0x1CBA7C4);
    let pts: Vec<FastPoint> = (0..targets)
        .map(|_| {
            g.mul(
                &mut GroupOps::default(),
                inst.generator,
                rng.gen_range(1..inst.r),
            )
        })
        .collect();

    let (mut build_scalar_ms, mut build_batched_ms) = (Vec::new(), Vec::new());
    let (mut probe_scalar_ms, mut probe_batched_ms) = (Vec::new(), Vec::new());
    let mut tables: Option<(FrobeniusPairTable, FrobeniusPairTable)> = None;
    for rep in 0..reps {
        let timed_build = |batched: bool| {
            let t = Instant::now();
            let table = if batched {
                FrobeniusPairTable::build(&inst, &fb)
            } else {
                FrobeniusPairTable::build_scalar(&inst, &fb)
            }
            .expect("table");
            (t.elapsed().as_secs_f64() * 1e3, table)
        };
        let ((ms_s, ts), (ms_b, tb)) = if rep % 2 == 0 {
            let s = timed_build(false);
            let b = timed_build(true);
            (s, b)
        } else {
            let b = timed_build(true);
            let s = timed_build(false);
            (s, b)
        };
        build_scalar_ms.push(ms_s);
        build_batched_ms.push(ms_b);
        assert_eq!(ts.entries, tb.entries, "table entries");
        assert_eq!(ts.build_ops, tb.build_ops, "charged build operations");
        tables = Some((ts, tb));
    }
    let (ts, tb) = tables.expect("at least one repeat");

    for rep in 0..reps {
        let run = |batched: bool| {
            let mut ops = GroupOps::default();
            let mut ctr = OracleCounters::default();
            let t = Instant::now();
            let mut out = Vec::with_capacity(pts.len());
            for &q in &pts {
                out.push(if batched {
                    decompose_mitm_frobenius(&g, &fb, &tb, 3, &mut ops, &mut ctr, q)
                } else {
                    decompose_mitm_frobenius_scalar(&g, &fb, &ts, 3, &mut ops, &mut ctr, q)
                });
            }
            (t.elapsed().as_secs_f64() * 1e3, out, ops, ctr.lookups)
        };
        let (s, b) = if rep % 2 == 0 {
            let s = run(false);
            let b = run(true);
            (s, b)
        } else {
            let b = run(true);
            let s = run(false);
            (s, b)
        };
        assert_eq!(s.1, b.1, "decompositions");
        assert_eq!(s.2, b.2, "charged operations");
        assert_eq!(s.3, b.3, "lookups");
        probe_scalar_ms.push(s.0);
        probe_batched_ms.push(b.0);
    }
    let found = pts
        .iter()
        .filter(|&&q| {
            decompose_mitm_frobenius(
                &g,
                &fb,
                &tb,
                3,
                &mut GroupOps::default(),
                &mut OracleCounters::default(),
                q,
            )
            .is_some()
        })
        .count();

    let line = |what: &str, s: &mut Vec<f64>, b: &mut Vec<f64>| {
        let (smin, bmin) = (
            s.iter().cloned().fold(f64::MAX, f64::min),
            b.iter().cloned().fold(f64::MAX, f64::min),
        );
        let (smed, bmed) = (median(s), median(b));
        println!(
            "{what:<26} scalar min {smin:>10.1} ms  median {smed:>10.1}   batched min {bmin:>10.1} ms  median {bmed:>10.1}   speed-up min {:.2}x  median {:.2}x",
            smin / bmin,
            smed / bmed
        );
    };
    line("table build", &mut build_scalar_ms, &mut build_batched_ms);
    line(
        &format!("{targets} probes (m = 3)"),
        &mut probe_scalar_ms,
        &mut probe_batched_ms,
    );
    println!(
        "identical: {} entries, charged build ops {:?}, {found}/{targets} targets decomposed",
        tb.entries, tb.build_ops
    );
    Ok(())
}
