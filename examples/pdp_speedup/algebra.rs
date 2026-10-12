//! PS-1 checks: bit-degree of the quartic's coefficients, and agreement of
//! the group-law formulation with the summation-polynomial one.

use crate::common::*;
use crate::corpus::*;
use crate::group::*;
use crypto_lib::cryptanalysis::semaev_decomp::{quartic_in_x3, Gf2};
use std::sync::atomic::{AtomicU64, Ordering};

fn arg<T: std::str::FromStr>(args: &[String], name: &str) -> Option<T> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
}

/// Algebraic degree of a Boolean function given by its truth table over `k` bits.
fn anf_degree(table: &[u8], k: u32) -> u32 {
    let mut a = table.to_vec();
    for i in 0..k {
        let step = 1usize << i;
        for j in 0..a.len() {
            if j & step != 0 {
                a[j] ^= a[j ^ step];
            }
        }
    }
    a.iter()
        .enumerate()
        .filter(|(_, &v)| v != 0)
        .map(|(m, _)| (m as u32).count_ones())
        .max()
        .unwrap_or(0)
}

/// `degrees --l L [--n N] [--xr K]`: exact algebraic degree, per coefficient
/// `c₀…c₄` of the quartic in the last abscissa, of every output bit as a
/// Boolean function of the bits of `(X₁, X₂)` (joint), and of `X₂` alone with
/// `X₁` fixed.  The maximum over output bits, over `K` targets (joint) and
/// over `K` fixed `X₁` (per-`X₂`).
pub fn cmd_degrees(args: &[String]) {
    let l: u32 = arg(args, "--l").expect("--l");
    let n: u32 = arg(args, "--n").unwrap_or(3 * l);
    let k: usize = arg(args, "--xr").unwrap_or(6);
    let rung = Rung::with(n, l);
    let gf = &rung.gf;
    let mut rng = SplitMix64(20261011 + n as u64 * 131 + l as u64);
    let span = rung.span() as usize;
    let mut joint = [0u32; 5];
    let mut single = [0u32; 5];
    for _ in 0..k {
        let xr = loop {
            let v = rng.next() & rung.mask();
            if v != 0 {
                break v;
            }
        };
        // Joint: x = x1 | x2 << l.
        let size = span * span;
        for ci in 0..5 {
            for b in 0..n {
                let mut table = vec![0u8; size];
                for x2 in 0..span {
                    for x1 in 0..span {
                        let q = quartic_in_x3(x1 as u64, x2 as u64, xr, gf);
                        table[x1 | (x2 << l)] = ((q.c[ci] >> b) & 1) as u8;
                    }
                }
                joint[ci] = joint[ci].max(anf_degree(&table, 2 * l));
            }
        }
        // Per-x2 with x1 fixed.
        let x1 = rng.next() & (rung.span() - 1);
        for ci in 0..5 {
            for b in 0..n {
                let mut table = vec![0u8; span];
                for x2 in 0..span {
                    let q = quartic_in_x3(x1, x2 as u64, xr, gf);
                    table[x2] = ((q.c[ci] >> b) & 1) as u8;
                }
                single[ci] = single[ci].max(anf_degree(&table, l));
            }
        }
    }
    println!(
        "{{\"n\":{n},\"l\":{l},\"targets\":{k},\"joint_degree\":{:?},\"x2_degree_x1_fixed\":{:?}}}",
        joint, single
    );
}

/// All `(n, ℓ)`-pair roots of the generic summation-polynomial enumerator
/// against the group-law reference, for every rational target.
pub fn cmd_agree(args: &[String]) {
    let n: u32 = arg(args, "--n").expect("--n");
    let l: u32 = arg(args, "--l").expect("--l");
    let threads: usize = arg(args, "--threads").unwrap_or(3);
    let rung = Rung::with(n, l);
    let gf = &rung.gf;
    let solver = ArtinSchreier::new(gf);
    let base = base_points(gf, &solver, l);
    // Every rational target, one y each (±R give the same root set).
    let targets: Vec<(u64, u64)> = (1..(1u64 << n))
        .filter_map(|x| lift(gf, &solver, x).map(|y| (x, y)))
        .collect();
    let next = AtomicU64::new(0);
    let stats = std::sync::Mutex::new((0u64, 0u64, 0u64, 0u64, 0u64, 0u64));
    std::thread::scope(|s| {
        for _ in 0..threads {
            s.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::Relaxed) as usize;
                if i >= targets.len() {
                    break;
                }
                let r = targets[i];
                let s4 = all_roots_serial(&rung, r.0);
                let grp = roots_reference(gf, &base, l, r);
                let degenerate = degenerate_target(gf, l, r.0);
                let mut st = stats.lock().unwrap();
                st.0 += 1; // targets
                st.1 += s4.len() as u64;
                if degenerate {
                    st.2 += 1;
                }
                if s4 == grp {
                    st.3 += 1; // equal
                } else {
                    // The extra roots must all be of the two degenerate families.
                    let extra: Vec<&[u64; 3]> = s4.iter().filter(|t| !grp.contains(t)).collect();
                    let missing = grp.iter().filter(|t| !s4.contains(t)).count();
                    let explained = extra
                        .iter()
                        .all(|t| explained_by_degeneracy(gf, &solver, r.0, t));
                    if degenerate && missing == 0 && explained {
                        st.4 += 1; // differs, explained by the degenerate family
                    } else {
                        st.5 += 1; // unexplained disagreement
                        println!(
                            "UNEXPLAINED n={n} l={l} xr={:#x} degenerate={degenerate} extra={:?} missing={missing}",
                            r.0, extra
                        );
                    }
                }
            });
        }
    });
    let st = stats.lock().unwrap();
    println!(
        "{{\"n\":{n},\"l\":{l},\"base_points\":{},\"targets\":{},\"s4_roots\":{},\"degenerate_targets\":{},\"equal\":{},\"differ_explained\":{},\"unexplained\":{}}}",
        base.len(), st.0, st.1, st.2, st.3, st.4, st.5
    );
    if st.5 != 0 {
        std::process::exit(1);
    }
}

/// The two families of extra roots a degenerate target may carry: two
/// non-liftable abscissae `a, b` and a liftable `c` with `a = b, c = x_R`
/// (F1) or `a·b = 1, c·x_R = 1` (F2).
fn explained_by_degeneracy(gf: &Gf2, solver: &ArtinSchreier, xr: u64, t: &[u64; 3]) -> bool {
    let nl = |x: u64| lift(gf, solver, x).is_none();
    [(t[0], t[1], t[2]), (t[0], t[2], t[1]), (t[1], t[2], t[0])]
        .into_iter()
        .any(|(a, b, c)| {
            nl(a)
                && nl(b)
                && !nl(c)
                && ((a == b && c == xr) || (gf.mul(a, b) == 1 && gf.mul(c, xr) == 1))
        })
}

fn all_roots_serial(rung: &Rung, xr: u64) -> Vec<[u64; 3]> {
    all_roots(rung, xr, 1)
}

/// `sample --n N --l L --pairs M`: `M` random pairs of base points against
/// random rational targets, the four group-law thirds (kept if they lie in
/// `V`) against the summation-polynomial roots of the same pair.
pub fn cmd_sample(args: &[String]) {
    let n: u32 = arg(args, "--n").expect("--n");
    let l: u32 = arg(args, "--l").expect("--l");
    let m: u64 = arg(args, "--pairs").unwrap_or(1_000_000);
    let rung = Rung::with(n, l);
    let gf = &rung.gf;
    let solver = ArtinSchreier::new(gf);
    let base = base_points(gf, &solver, l);
    let lv = crypto_lib::cryptanalysis::semaev_decomp::subspace_poly(l, gf);
    let mut rng = SplitMix64(20261011 + 7 * n as u64 + l as u64);
    let (mut agree, mut with_root, mut bad, mut degenerate) = (0u64, 0u64, 0u64, 0u64);
    for _ in 0..m {
        let r = loop {
            let x = rng.next() & rung.mask();
            if x == 0 {
                continue;
            }
            if let Some(y) = lift(gf, &solver, x) {
                break (x, y);
            }
        };
        let a = base[(rng.next() % base.len() as u64) as usize];
        let b = base[(rng.next() % base.len() as u64) as usize];
        let (a, b) = if a.0 <= b.0 { (a, b) } else { (b, a) };
        let mut s4 = Vec::new();
        pair_roots(&rung, &lv, r.0, a.0, b.0, &mut s4);
        s4.sort_unstable();
        s4.dedup();
        let mut grp: Vec<[u64; 3]> = thirds(gf, a, b, r)
            .into_iter()
            .flatten()
            .filter(|&t| t < rung.span())
            .map(|t| sorted3(a.0, b.0, t))
            .collect();
        grp.sort_unstable();
        grp.dedup();
        if s4 == grp {
            agree += 1;
            if !grp.is_empty() {
                with_root += 1;
            }
        } else if degenerate_target(gf, l, r.0) {
            degenerate += 1; // outside the formulation's domain; counted, not compared
        } else {
            bad += 1;
            println!(
                "DISAGREE xr={:#x} pair=({:#x},{:#x}) s4={:?} grp={:?}",
                r.0, a.0, b.0, s4, grp
            );
        }
    }
    println!(
        "{{\"n\":{n},\"l\":{l},\"pairs\":{m},\"agree\":{agree},\"agree_with_root\":{with_root},\"disagree\":{bad},\"degenerate_skipped\":{degenerate}}}"
    );
    if bad != 0 {
        std::process::exit(1);
    }
}

/// `torsion --n N [--samples M]`: translation by the 2-torsion point
/// `T = (0, 1)` is `x ↦ 1/x` on this curve.  Every liftable `x ≠ 0` for
/// `n ≤ 12`, a sample of `M` otherwise; also that `T` has order 2.
pub fn cmd_torsion(args: &[String]) {
    let n: u32 = arg(args, "--n").expect("--n");
    let samples: u64 = arg(args, "--samples").unwrap_or(200_000);
    let rung = Rung::with(n, n.min(3));
    let gf = &rung.gf;
    let solver = ArtinSchreier::new(gf);
    let t = Pt::Aff(0, 1);
    assert!(on_curve(gf, 0, 1), "T is on the curve");
    assert_eq!(add(gf, t, t), Pt::Inf, "T has order 2");
    let mut rng = SplitMix64(20261011 + n as u64);
    let (mut checked, mut bad) = (0u64, 0u64);
    let mut visit = |x: u64| {
        if x == 0 {
            return;
        }
        if let Some(y) = lift(gf, &solver, x) {
            for yy in [y, y ^ x] {
                checked += 1;
                match add(gf, Pt::Aff(x, yy), t) {
                    Pt::Aff(x2, _) if x2 == gf.inv(x) => {}
                    _ => bad += 1,
                }
            }
        }
    };
    if n <= 12 {
        for x in 1..(1u64 << n) {
            visit(x);
        }
    } else {
        for _ in 0..samples {
            visit(rng.next() & rung.mask());
        }
    }
    println!("{{\"n\":{n},\"points_checked\":{checked},\"mismatches\":{bad}}}");
    if bad != 0 {
        std::process::exit(1);
    }
}

/// `corpus-check --l L --corpus FILE`: the exact group-law reference against
/// the PS-0 corpus's root sets, target by target.
pub fn cmd_corpus_check(args: &[String]) {
    let l: u32 = arg(args, "--l").expect("--l");
    let path: String = arg(args, "--corpus").expect("--corpus");
    let threads: usize = arg(args, "--threads").unwrap_or(3);
    let loaded = crate::corpus::read_corpus(&std::fs::read_to_string(&path).unwrap());
    let rung = Rung::new(l);
    let gf = &rung.gf;
    let solver = ArtinSchreier::new(gf);
    let base = base_points(gf, &solver, l);
    let next = AtomicU64::new(0);
    let res = std::sync::Mutex::new((0u64, 0u64, 0u64));
    std::thread::scope(|s| {
        for _ in 0..threads {
            s.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::Relaxed) as usize;
                if i >= loaded.recs.len() {
                    break;
                }
                let rec = &loaded.recs[i];
                let grp = roots_reference(gf, &base, l, (rec.t.xr, rec.t.yr));
                let deg = degenerate_target(gf, l, rec.t.xr);
                let mut r = res.lock().unwrap();
                r.0 += 1;
                if deg {
                    r.2 += 1;
                }
                if grp != rec.roots {
                    r.1 += 1;
                    println!(
                        "DIFFERENT target {} roots {:?} vs corpus {:?}",
                        rec.t.idx, grp, rec.roots
                    );
                }
            });
        }
    });
    let r = res.lock().unwrap();
    println!(
        "{{\"l\":{l},\"base_points\":{},\"targets\":{},\"differ\":{},\"degenerate_targets\":{}}}",
        base.len(),
        r.0,
        r.1,
        r.2
    );
    if r.1 != 0 {
        std::process::exit(1);
    }
}
