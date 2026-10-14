//! Native instance screen for the frozen quadratic-extension comparison.
//! `screen` prints candidate parameters only; the paired DLP run is a
//! separate, frozen stage in the study protocol.

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::fmt::Write as _;
use std::rc::Rc;
use std::time::Instant;

use crypto_lib::cryptanalysis::ecbench::canonical::{canonical, sha256_hex};
use crypto_lib::cryptanalysis::ext_curve::{jacobi, take_field_counters, ExtDiagonalAutomorphism};
use crypto_lib::cryptanalysis::ext_curve::{random_point, ExtCurve, ExtField, ExtPoint, Fp2};
use crypto_lib::cryptanalysis::glv_invariant_base::{
    factor_u64, rho_reference_folded, Endomorphism, Negation,
};
use crypto_lib::cryptanalysis::ic_boundary::{
    calibrate_group, calibrate_row_ops, rho_reference_negation, Calibration, CountedGroup,
    FactorBase, GroupOps, OracleCounters,
};
use crypto_lib::cryptanalysis::ic_framework::plugins::MitmOracle;
use crypto_lib::cryptanalysis::ic_framework::stages::{
    DecompositionOracle, FactorBaseBuilder, InstanceCtx, Params, Targets,
};
use crypto_lib::cryptanalysis::ic_framework::{run_pipeline, PipelineSpec};
use crypto_lib::cryptanalysis::residual_walk::is_prime_u64;
use crypto_lib::cryptanalysis::residual_walk::pow_mod;
use crypto_lib::hash::sha256::sha256;
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde_json::json;

type Point = ExtPoint<[u64; 2]>;
type Curve = ExtCurve<Fp2>;

fn parse_point(v: &serde_json::Value) -> Point {
    let pair = |name: &str| -> [u64; 2] {
        let a = v[name].as_array().expect("coordinate array");
        [a[0].as_u64().unwrap(), a[1].as_u64().unwrap()]
    };
    if v["infinity"].as_bool().unwrap() {
        Point::infinity([0, 0])
    } else {
        Point::affine(pair("x"), pair("y"))
    }
}

struct LineBase {
    frozen: Vec<Point>,
    projected: Rc<RefCell<Vec<Point>>>,
    columns: usize,
    line: [u64; 2],
}

impl FactorBaseBuilder<Curve> for LineBase {
    fn name(&self) -> &str {
        "extension-line-projected-columns"
    }
    fn describe(&self, _params: &Params) -> String {
        "signed raw points on x=s*t; distinct cofactor projections are columns".into()
    }
    fn build(
        &self,
        ctx: &InstanceCtx<Curve>,
        _params: &Params,
        ops: &mut GroupOps,
    ) -> Result<FactorBase<Point>, String> {
        let curve = ctx.group;
        let f = &curve.f;
        let mut raw = Vec::new();
        let mut projected = Vec::new();
        let mut seen = BTreeSet::new();
        let mut sqrt_calls = 0;
        let mut tested = 0;
        for t in 1..f.p() {
            let mut lifts = curve.lift_x(f.scale(self.line, t));
            sqrt_calls += 1;
            lifts.sort_by_key(|p| f.pack(p.y));
            for p in lifts {
                tested += 1;
                let q = curve.mul(ops, p, ctx.cofactor);
                if q.infinity {
                    continue;
                }
                let orbit = curve.key(&q).min(curve.key(&curve.neg(q)));
                if seen.insert(orbit) {
                    raw.push(p);
                    projected.push(q);
                    if raw.len() == self.columns {
                        break;
                    }
                }
            }
            if raw.len() == self.columns {
                break;
            }
        }
        if raw != self.frozen {
            return Err("built factor base differs from frozen fixture".into());
        }
        *self.projected.borrow_mut() = projected;
        let mut points = Vec::with_capacity(2 * self.columns);
        let mut cols = Vec::with_capacity(2 * self.columns);
        let mut coeffs = Vec::with_capacity(2 * self.columns);
        for (j, p) in raw.iter().enumerate() {
            points.push(*p);
            cols.push(j);
            coeffs.push(1);
            points.push(curve.neg(*p));
            cols.push(j);
            coeffs.push(ctx.r - 1);
        }
        let mut fb = FactorBase::from_column_map(
            format!("{} distinct projected columns", self.columns),
            points,
            cols,
            coeffs,
            self.columns,
            |p| curve.key(p),
            |p| curve.key(&curve.neg(*p)),
            |p| f.pack(p.x),
        )?;
        fb.cost.count("sqrt_solves", sqrt_calls);
        fb.cost.count("abscissae_scanned", sqrt_calls);
        fb.cost.count("lifted_points_tested", tested);
        Ok(fb)
    }
}

#[derive(Clone, Copy)]
enum IndexMode {
    Projected,
    Chi,
    Isogeny,
}

fn orbit_key(curve: &Curve, p: Point) -> Option<[u64; 2]> {
    let mut z = chi(curve, p);
    let mut best = z;
    for _ in 0..2 {
        z = mobius(&curve.f, z);
        if z < best {
            best = z;
        }
    }
    best
}

struct ProjectedOracle {
    mode: IndexMode,
    projected: Rc<RefCell<Vec<Point>>>,
    exact: BTreeMap<u64, Vec<(usize, usize)>>,
    chi_index: BTreeMap<Option<[u64; 2]>, Vec<(u64, usize, usize)>>,
    setup_lookups: u64,
    setup_transforms: u64,
    pair_count: u64,
    false_candidates: u64,
}

impl ProjectedOracle {
    fn new(mode: IndexMode, projected: Rc<RefCell<Vec<Point>>>) -> Self {
        Self {
            mode,
            projected,
            exact: BTreeMap::new(),
            chi_index: BTreeMap::new(),
            setup_lookups: 0,
            setup_transforms: 0,
            pair_count: 0,
            false_candidates: 0,
        }
    }
}

impl DecompositionOracle<Curve> for ProjectedOracle {
    fn name(&self) -> &str {
        match self.mode {
            IndexMode::Projected => "projected-pair",
            IndexMode::Chi => "chi-orbit-index",
            IndexMode::Isogeny => "isogeny-projected-pair",
        }
    }
    fn summands(&self) -> u32 {
        2
    }
    fn describe(&self, _: &Params) -> String {
        match self.mode {
            IndexMode::Projected=>"pair sums after exact cofactor projection".into(),
            IndexMode::Chi=>"pair sums indexed by the chi Mobius orbit, with exact projected-point check; this is not a polynomial-system oracle".into(),
            IndexMode::Isogeny=>"source pair sums mapped by the explicit degree-three isogeny after cofactor projection".into(),
        }
    }
    fn prepare(
        &mut self,
        ctx: &InstanceCtx<Curve>,
        fb: &FactorBase<Point>,
        _: &Params,
        ops: &mut GroupOps,
    ) -> Result<(), String> {
        let src = ctx.group;
        let f = &src.f;
        let dst = Curve {
            f: *f,
            a: f.zero(),
            b: f.from_fp(f.p - 27),
        };
        let base = self.projected.borrow();
        if base.len() != fb.columns {
            return Err("projected base not built".into());
        }
        let mut vals = Vec::with_capacity(fb.points.len());
        for (i, p) in base.iter().enumerate() {
            let _ = i;
            vals.push(*p);
            vals.push(src.neg(*p));
        }
        if let IndexMode::Isogeny = self.mode {
            vals = vals
                .into_iter()
                .map(|p| {
                    self.setup_transforms += 1;
                    phi(src, p)
                })
                .collect();
        }
        for i in 0..vals.len() {
            for j in i..vals.len() {
                let sum = match self.mode {
                    IndexMode::Isogeny => dst.add(ops, vals[i], vals[j]),
                    _ => src.add(ops, vals[i], vals[j]),
                };
                if sum.infinity {
                    continue;
                }
                self.pair_count += 1;
                let key = match self.mode {
                    IndexMode::Isogeny => dst.key(&sum),
                    _ => src.key(&sum),
                };
                match self.mode {
                    IndexMode::Chi => {
                        let orbit = orbit_key(src, sum);
                        self.setup_transforms += 1;
                        self.chi_index.entry(orbit).or_default().push((key, i, j));
                    }
                    _ => {
                        self.exact.entry(key).or_default().push((i, j));
                    }
                }
                self.setup_lookups += 1;
            }
        }
        Ok(())
    }
    fn setup_native(&self) -> BTreeMap<String, u64> {
        BTreeMap::from([
            ("lookups".into(), self.setup_lookups),
            ("canonicalisations".into(), self.setup_transforms),
            ("pair_entries".into(), self.pair_count),
        ])
    }
    fn decompose(
        &mut self,
        ctx: &InstanceCtx<Curve>,
        fb: &FactorBase<Point>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: Point,
    ) -> Option<Vec<usize>> {
        let src = ctx.group;
        let hpoint = src.mul(ops, point, ctx.cofactor % ctx.r);
        let key = match self.mode {
            IndexMode::Projected => src.key(&hpoint),
            IndexMode::Chi => src.key(&hpoint),
            IndexMode::Isogeny => {
                counters.canonicalisations += 1;
                let f = &src.f;
                let dst = Curve {
                    f: *f,
                    a: f.zero(),
                    b: f.from_fp(f.p - 27),
                };
                dst.key(&phi(src, hpoint))
            }
        };
        counters.lookups += 1;
        let candidate = match self.mode {
            IndexMode::Chi => {
                counters.canonicalisations += 1;
                self.chi_index.get(&orbit_key(src, hpoint)).and_then(|xs| {
                    xs.iter()
                        .find(|(k, _, _)| *k == key)
                        .map(|(_, i, j)| (*i, *j))
                })
            }
            _ => self.exact.get(&key).and_then(|xs| xs.first().copied()),
        };
        let (i, j) = candidate?;
        // Independently check the relation in the source subgroup.
        let pi = src.mul(ops, fb.points[i], ctx.cofactor);
        let pj = src.mul(ops, fb.points[j], ctx.cofactor);
        let sum = src.add(ops, pi, pj);
        if sum != hpoint {
            self.false_candidates += 1;
            return None;
        }
        Some(vec![i, j])
    }
}

fn digest(data: &str) -> String {
    sha256_hex(data.as_bytes())
}

fn icv1(p: u64, nu: u64, n2: u64, b: u64) -> serde_json::Value {
    let c0 = (p - nu) % p;
    let field = format!(
        "fpk-{p}-2-{}",
        &digest(&format!("fpk-modulus:{p}:{c0},0"))[..8]
    );
    let model = format!("{{\"a\":[\"0\",\"0\"],\"b\":[\"{b}\",\"0\"],\"field\":\"{field}\",\"form\":\"y^2=x^3+a*x+b\",\"k\":\"2\",\"modulus\":[\"{c0}\",\"0\"],\"p\":\"{p}\",\"v\":\"1\"}}");
    let model_hash = digest(&model);
    let trace = p as i64 * p as i64 + 1 - n2 as i64;
    let id = format!(
        "ICV1:{field}:{trace}:{n2}:0,0:unk:unk:r:{}",
        &model_hash[..12]
    );
    json!({"id":id,"model":model,"field":field,"model_sha256":model_hash})
}

fn chi(curve: &ExtCurve<Fp2>, point: Point) -> Option<[u64; 2]> {
    let f = &curve.f;
    if point.infinity {
        return Some(f.one());
    }
    let den = f.sub(point.y, f.one());
    f.inv(den).map(|inv| f.mul(f.add(point.y, f.one()), inv))
}

fn mobius(f: &Fp2, z: Option<[u64; 2]>) -> Option<[u64; 2]> {
    match z {
        None => Some(f.zero()),
        Some(z) => f.inv(f.sub(f.one(), z)),
    }
}

fn phi(curve: &ExtCurve<Fp2>, point: Point) -> Point {
    let f = &curve.f;
    if point.infinity || point.x == f.zero() {
        return curve.identity();
    }
    let x2 = f.sqr(point.x);
    let x3 = f.mul(x2, point.x);
    let xp = f.mul(f.add(x3, f.from_fp(4)), f.inv(x2).unwrap());
    let yp = f.mul(f.mul(point.y, f.sub(x3, f.from_fp(8))), f.inv(x3).unwrap());
    Point::affine(xp, yp)
}

fn freeze() {
    // Chosen by selector v2, SCREEN_V2.jsonl. The solver reads these
    // parameters, never the answer scalar supplied only to verification.
    for (index, (p, expected_r)) in [(61, 19), (271, 61), (787, 211), (3319, 811)]
        .into_iter()
        .enumerate()
    {
        let n1 = base_order(p);
        let trace = p as i64 + 1 - n1 as i64;
        let n2 = n1 * (p as i64 + 1 + trace) as u64;
        let r = factor_u64(n2 / n1).last().unwrap().0;
        assert_eq!(r, expected_r);
        let h = n2 / r;
        assert_ne!(h % r, 0);
        let f = Fp2::new(p).unwrap();
        let curve = ExtCurve {
            f,
            a: f.zero(),
            b: f.one(),
        };
        let dest = ExtCurve {
            f,
            a: f.zero(),
            b: f.from_fp(p - 27),
        };
        let line = f.s();
        let columns = if index == 0 {
            8
        } else {
            ((r as f64).sqrt().ceil() as usize).min(32)
        };
        let mut raw = Vec::new();
        let mut projected = Vec::new();
        let mut seen = BTreeSet::new();
        let mut ops = GroupOps::default();
        let mut scanned = 0u64;
        for t in 1..p {
            let mut lifts = curve.lift_x(f.scale(line, t));
            lifts.sort_by_key(|point| f.pack(point.y));
            for point in lifts {
                scanned += 1;
                let q = curve.mul(&mut ops, point, h);
                if q.infinity {
                    continue;
                }
                let orbit = curve.key(&q).min(curve.key(&curve.neg(q)));
                if seen.insert(orbit) {
                    raw.push(point);
                    projected.push(q);
                    if raw.len() == columns {
                        break;
                    }
                }
            }
            if raw.len() == columns {
                break;
            }
        }
        assert_eq!(raw.len(), columns, "insufficient projected base");
        let mut rng = StdRng::seed_from_u64(0x2026_1014 ^ p);
        let generator = loop {
            let g = curve.mul(&mut ops, random_point(&curve, &mut rng), h);
            if !g.infinity {
                break g;
            }
        };
        assert!(curve.mul(&mut ops, generator, r).infinity);
        assert!(generator.x[1] != 0 || generator.y[1] != 0);
        let target_seed = format!("torsion-extension/2026-10-14/p={p}/r={r}");
        let hash = sha256(target_seed.as_bytes());
        let k = 1 + u64::from_be_bytes(hash[0..8].try_into().unwrap()) % (r - 1);
        let target = curve.mul(&mut ops, generator, k);
        assert!(curve.is_on_curve(target));
        let torsion = Point::affine(f.zero(), f.one());
        assert!(curve.mul(&mut ops, torsion, 3).infinity);
        let mut chi_checks = 0;
        let mut map_checks = 0;
        for _ in 0..64 {
            let a = random_point(&curve, &mut rng);
            let b = random_point(&curve, &mut rng);
            let sum = curve.add(&mut ops, a, b);
            assert_eq!(
                chi(&curve, curve.add(&mut ops, a, torsion)),
                mobius(&f, chi(&curve, a))
            );
            chi_checks += 1;
            let mapped = phi(&curve, a);
            assert!(dest.is_on_curve(mapped));
            assert_eq!(phi(&curve, sum), dest.add(&mut ops, mapped, phi(&curve, b)));
            map_checks += 1;
        }
        assert_eq!(phi(&curve, torsion), dest.identity());
        assert_eq!(phi(&curve, curve.neg(torsion)), dest.identity());
        println!(
            "{}",
            json!({
                "schema":"torsion_extension_fixture/v1", "p":p,"nu":f.nu,
                "source":icv1(p,f.nu,n2,1),"destination":icv1(p,f.nu,n2,p-27),
                "n1":n1,"trace_fp":trace,"n2":n2,"r":r,"cofactor":h,
                "line":line,"generator":generator,"target":target,"answer_scalar":k,
                "target_seed":target_seed,"base_raw":raw,"base_projected":projected,
                "columns":columns,"raw_points_scanned":scanned,
                "chi_checks":chi_checks,"map_checks":map_checks,"setup_group_ops":ops,
                "embedding_degree":embedding_degree(p,r),
            })
        );
    }
}

fn thread_cpu_ns() -> u64 {
    let mut ts = libc::timespec {
        tv_sec: 0,
        tv_nsec: 0,
    };
    // The experiment runs on macOS and Linux; thread CPU avoids charging
    // scheduler descheduling as cryptographic work. Wall is retained too.
    let rc = unsafe { libc::clock_gettime(libc::CLOCK_THREAD_CPUTIME_ID, &mut ts) };
    assert_eq!(rc, 0, "thread CPU clock unavailable");
    ts.tv_sec as u64 * 1_000_000_000 + ts.tv_nsec as u64
}

fn timed<T>(f: impl FnOnce() -> T) -> (T, u64, u64, (u64, u64)) {
    take_field_counters();
    let wall_start = Instant::now();
    let cpu_start = thread_cpu_ns();
    let value = f();
    let cpu_ns = thread_cpu_ns() - cpu_start;
    let wall_ns = wall_start.elapsed().as_nanos() as u64;
    let field = take_field_counters();
    (value, cpu_ns, wall_ns, field)
}

fn cpu_add_calibration(curve: &Curve, points: &[Point]) -> f64 {
    let mut ops = GroupOps::default();
    let mut acc = points[0];
    let start = thread_cpu_ns();
    for i in 0..100_000 {
        acc = curve.add(&mut ops, acc, points[(i * 7 + 1) % points.len()]);
        if acc.infinity {
            acc = points[(i + 3) % points.len()];
        }
    }
    std::hint::black_box(acc);
    (thread_cpu_ns() - start) as f64 / 100_000.0
}

fn calibrate(curve: &Curve, points: &[Point], r: u64) -> (Calibration, f64) {
    let mut c = Calibration::default();
    calibrate_group(curve, points, &mut c);
    calibrate_row_ops(r, &mut c);
    let add_cpu = cpu_add_calibration(curve, points);
    let mut map = BTreeMap::new();
    for i in 0..1024u64 {
        map.insert(i, i);
    }
    let start = thread_cpu_ns();
    let mut accum = 0;
    for i in 0..200_000u64 {
        accum ^= *map.get(&(i % 2048)).unwrap_or(&0);
    }
    std::hint::black_box(accum);
    c.ns_per_lookup = (thread_cpu_ns() - start) as f64 / 200_000.0;
    let start = thread_cpu_ns();
    for i in 0..2000u64 {
        std::hint::black_box(curve.lift_x(curve.f.scale(curve.f.s(), 1 + i % curve.f.p)));
    }
    c.ns_per_sqrt = Some((thread_cpu_ns() - start) as f64 / 2000.0);
    (c, add_cpu)
}

fn transform_calibration(curve: &Curve, point: Point, mode: IndexMode) -> f64 {
    let n = 4000;
    let start = thread_cpu_ns();
    for _ in 0..n {
        match mode {
            IndexMode::Chi => {
                std::hint::black_box(orbit_key(curve, point));
            }
            IndexMode::Isogeny => {
                std::hint::black_box(phi(curve, point));
            }
            IndexMode::Projected => {}
        }
    }
    (thread_cpu_ns() - start) as f64 / n as f64
}

fn six_class_endomorphism(
    curve: &Curve,
    g: Point,
    r: u64,
) -> (ExtDiagonalAutomorphism<[u64; 2]>, GroupOps) {
    let p = curve.f.p();
    let z = (2..p)
        .map(|a| pow_mod(a, (p - 1) / 3, p))
        .find(|&z| z != 1)
        .unwrap();
    let image = Point::affine(curve.f.mul(curve.f.from_fp(z), g.x), g.y);
    let mut ops = GroupOps::default();
    let lambda = (1..r)
        .filter(|&a| (a * a + a + 1) % r == 0)
        .find(|&a| curve.mul(&mut ops, g, a) == image)
        .expect("zeta eigenvalue");
    (
        ExtDiagonalAutomorphism {
            cx: curve.f.from_fp(z),
            cy: curve.f.one(),
            eigenvalue: lambda,
            order: 3,
            label: "j0-zeta",
        },
        ops,
    )
}

fn solve(path: &str) {
    let fixtures = std::fs::read_to_string(path).expect("fixture JSONL");
    for (index, line) in fixtures.lines().enumerate() {
        let v: serde_json::Value = serde_json::from_str(line).unwrap();
        let p = v["p"].as_u64().unwrap();
        let r = v["r"].as_u64().unwrap();
        let h = v["cofactor"].as_u64().unwrap();
        let n2 = v["n2"].as_u64().unwrap();
        let f = Fp2::new(p).unwrap();
        assert_eq!(f.nu, v["nu"].as_u64().unwrap());
        let curve = Curve {
            f,
            a: f.zero(),
            b: f.one(),
        };
        let g = parse_point(&v["generator"]);
        let q = parse_point(&v["target"]);
        assert!(curve.is_on_curve(g) && curve.is_on_curve(q));
        let mut check_ops = GroupOps::default();
        assert!(curve.mul(&mut check_ops, g, r).infinity);
        assert!(curve.mul(&mut check_ops, q, r).infinity);
        let frozen_raw: Vec<Point> = v["base_raw"]
            .as_array()
            .unwrap()
            .iter()
            .map(parse_point)
            .collect();
        let columns = v["columns"].as_u64().unwrap() as usize;
        assert_eq!(columns, frozen_raw.len());
        let answers = v["answer_scalar"].as_u64().unwrap();
        let (calib, ns_add_cpu) = calibrate(&curve, &[g, q, curve.neg(g), curve.neg(q)], r);
        let seed_base = 0x2026_1014u64 ^ p;
        let cap = (20 * r).min(1_000_000);
        let ctx = InstanceCtx {
            group: &curve,
            generator: g,
            target: q,
            r,
            cofactor: h,
            group_order: n2,
            name: v["source"]["id"].as_str().unwrap().into(),
            field_degree: Some(2),
        };
        println!(
            "{}",
            json!({"kind":"calibration","p":p,"r":r,"index":index,"calibration":calib,"ns_per_add_thread_cpu":ns_add_cpu,"cap":cap,"target":q})
        );
        for round in 0..3u64 {
            let seed = seed_base.wrapping_add(round * 0x9e37_79b9);
            let mut arm_names = vec![
                "raw-pair",
                "projected-pair",
                "chi-orbit-index",
                "isogeny-projected-pair",
                "rho-negation",
                "rho-negation-aa",
                "rho-six-class",
            ];
            if round % 2 == 1 {
                arm_names.reverse();
            }
            for arm in arm_names {
                let (result, extra, cpu_ns, wall_ns, field) = match arm {
                    "rho-negation" | "rho-negation-aa" => {
                        let (res, cpu_ns, wall_ns, field) =
                            timed(|| rho_reference_negation(&curve, g, q, r, seed, 1_000_000));
                        (
                            json!(res),
                            json!({"status":if res.verified{"verified"}else{"exhausted"}}),
                            cpu_ns,
                            wall_ns,
                            field,
                        )
                    }
                    "rho-six-class" => {
                        let ((res, eigen_ops, lambda), cpu_ns, wall_ns, field) = timed(|| {
                            let (aut, eigen_ops) = six_class_endomorphism(&curve, g, r);
                            let neg = Negation { r };
                            let gens: [&dyn Endomorphism<Curve>; 2] = [&neg, &aut];
                            let res = rho_reference_folded(&curve, g, q, r, seed, 1_000_000, &gens)
                                .unwrap();
                            (res, eigen_ops, aut.eigenvalue)
                        });
                        (
                            json!(res),
                            json!({"eigen_setup_group_ops":eigen_ops,"eigenvalue":lambda,"status":if res.verified{"verified"}else{"exhausted"}}),
                            cpu_ns,
                            wall_ns,
                            field,
                        )
                    }
                    _ => {
                        let mut c = calib.clone();
                        if arm != "raw-pair" {
                            let mode = match arm {
                                "projected-pair" => IndexMode::Projected,
                                "chi-orbit-index" => IndexMode::Chi,
                                _ => IndexMode::Isogeny,
                            };
                            c.ns_per_canon = Some(transform_calibration(&curve, g, mode));
                        }
                        let (report, cpu_ns, wall_ns, field) = timed(|| {
                            let projected = Rc::new(RefCell::new(Vec::new()));
                            let base = LineBase {
                                frozen: frozen_raw.clone(),
                                projected: projected.clone(),
                                columns,
                                line: f.s(),
                            };
                            let spec = PipelineSpec {
                                factor_base: base.name().into(),
                                oracle: arm.into(),
                                targets: Targets::Walk,
                                linalg: "incremental-gauss".into(),
                                max_trials: cap,
                                seed,
                                ..Default::default()
                            };
                            if arm == "raw-pair" {
                                let mut oracle = MitmOracle::new(2);
                                run_pipeline(&ctx, &spec, &base, &mut oracle, answers, &c, None)
                                    .unwrap()
                            } else {
                                let mode = match arm {
                                    "projected-pair" => IndexMode::Projected,
                                    "chi-orbit-index" => IndexMode::Chi,
                                    _ => IndexMode::Isogeny,
                                };
                                let mut oracle = ProjectedOracle::new(mode, projected);
                                let rep = run_pipeline(
                                    &ctx,
                                    &spec,
                                    &base,
                                    &mut oracle,
                                    answers,
                                    &c,
                                    None,
                                )
                                .unwrap();
                                assert_eq!(
                                    oracle.false_candidates, 0,
                                    "invalid projected relation"
                                );
                                rep
                            }
                        });
                        (
                            json!(report),
                            json!({"status":if report.verified{"verified"}else if report.exhausted{"exhausted"}else{"wrong_answer"}}),
                            cpu_ns,
                            wall_ns,
                            field,
                        )
                    }
                };
                let s_cpu = (cpu_ns as f64 / ns_add_cpu) / (r as f64).sqrt();
                println!(
                    "{}",
                    json!({
                        "kind":"run","p":p,"r":r,"round":round,"seed":seed,
                        "arm":arm,"cpu_ns":cpu_ns,"wall_ns":wall_ns,
                        "s_cpu_gae":s_cpu,"field_fp_muls":field.0,"field_fp_invs":field.1,
                        "result":result,"extra":extra,
                    })
                );
            }
        }
    }
}

fn verify(path: &str, fixture_path: &str) {
    let fixture: Vec<serde_json::Value> = std::fs::read_to_string(fixture_path)
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    let mut checked = 0u64;
    let mut exhausted = 0u64;
    for line in std::fs::read_to_string(path).unwrap().lines() {
        let v: serde_json::Value = serde_json::from_str(line).unwrap();
        if v["kind"] != "run" {
            continue;
        }
        let p = v["p"].as_u64().unwrap();
        let fix = fixture.iter().find(|f| f["p"].as_u64() == Some(p)).unwrap();
        let f = Fp2::new(p).unwrap();
        let curve = Curve {
            f,
            a: f.zero(),
            b: f.one(),
        };
        let g = parse_point(&fix["generator"]);
        let q = parse_point(&fix["target"]);
        let rec = v["result"]["recovered"].as_u64();
        if let Some(k) = rec {
            let mut ops = GroupOps::default();
            assert_eq!(curve.mul(&mut ops, g, k), q, "scalar replay fails: p={p}");
            assert_eq!(k, fix["answer_scalar"].as_u64().unwrap());
            assert_eq!(v["extra"]["status"], "verified");
            checked += 1;
        } else {
            assert_eq!(v["extra"]["status"], "exhausted");
            exhausted += 1;
        }
    }
    println!(
        "{}",
        json!({"schema":"torsion_extension_replay/v1","verified":checked,"exhausted":exhausted,"fixture_sha256":sha256_hex(std::fs::read(fixture_path).unwrap().as_slice()),"records_sha256":sha256_hex(std::fs::read(path).unwrap().as_slice())})
    );
}

fn support(path: &str) {
    for line in std::fs::read_to_string(path).unwrap().lines() {
        let v: serde_json::Value = serde_json::from_str(line).unwrap();
        let p = v["p"].as_u64().unwrap();
        let r = v["r"].as_u64().unwrap();
        let h = v["cofactor"].as_u64().unwrap();
        let f = Fp2::new(p).unwrap();
        let curve = Curve {
            f,
            a: f.zero(),
            b: f.one(),
        };
        let dest = Curve {
            f,
            a: f.zero(),
            b: f.from_fp(p - 27),
        };
        let g = parse_point(&v["generator"]);
        let target = parse_point(&v["target"]);
        let raw: Vec<Point> = v["base_raw"]
            .as_array()
            .unwrap()
            .iter()
            .map(parse_point)
            .collect();
        let mut points = Vec::with_capacity(raw.len() * 2);
        for point in raw {
            points.push(point);
            points.push(curve.neg(point));
        }
        let mut ops = GroupOps::default();
        let torsion = Point::affine(f.zero(), f.one());
        assert!(curve.mul(&mut ops, torsion, 3).infinity);
        assert!(
            curve.mul(&mut ops, torsion, h).infinity,
            "cofactor must kill T"
        );
        let t2 = curve.add(&mut ops, torsion, torsion);
        let shifts = [curve.identity(), torsion, t2];
        let rshifts: Vec<Point> = shifts
            .iter()
            .map(|&pt| curve.mul(&mut ops, pt, r))
            .collect();
        let mut supports: [BTreeSet<u64>; 3] = std::array::from_fn(|_| BTreeSet::new());
        let mut raw_pair_counts = [0u64; 3];
        let mut projected_support = BTreeSet::new();
        let mut mapped_support = BTreeSet::new();
        let mut chi_keys = BTreeMap::<Option<[u64; 2]>, BTreeSet<u64>>::new();
        let mut chi_checks = 0u64;
        let mut pairs = 0u64;
        for i in 0..points.len() {
            for j in i..points.len() {
                let sum = curve.add(&mut ops, points[i], points[j]);
                if sum.infinity {
                    continue;
                }
                pairs += 1;
                let projected = curve.mul(&mut ops, sum, h);
                if !projected.infinity {
                    let key = curve.key(&projected);
                    projected_support.insert(key);
                    mapped_support.insert(dest.key(&phi(&curve, projected)));
                }
                assert_eq!(
                    chi(&curve, curve.add(&mut ops, sum, torsion)),
                    mobius(&f, chi(&curve, sum))
                );
                chi_checks += 1;
                chi_keys
                    .entry(orbit_key(&curve, projected))
                    .or_default()
                    .insert(curve.key(&projected));
                let rsum = curve.mul(&mut ops, sum, r);
                for a in 0..3 {
                    if rsum == rshifts[a] {
                        let corrected = curve.add(&mut ops, sum, curve.neg(shifts[a]));
                        assert!(curve.mul(&mut ops, corrected, r).infinity);
                        if !corrected.infinity {
                            supports[a].insert(curve.key(&corrected));
                            raw_pair_counts[a] += 1;
                        }
                    }
                }
            }
        }
        let mut union = BTreeSet::new();
        for set in &supports {
            union.extend(set.iter().copied());
        }
        let mut subgroup = BTreeSet::new();
        for k in 1..r {
            subgroup.insert(curve.key(&curve.mul(&mut ops, g, k)));
        }
        assert_eq!(subgroup.len(), (r - 1) as usize);
        assert!(union.is_subset(&subgroup));
        assert_eq!(
            projected_support.len(),
            mapped_support.len(),
            "isogeny injectivity on observed subgroup sums"
        );
        let target_key = curve.key(&target);
        println!(
            "{}",
            json!({
                "schema":"torsion_extension_support/v1","p":p,"r":r,"signed_base_points":points.len(),
                "pair_sums_nonidentity":pairs,"raw_pair_counts_by_shift":raw_pair_counts,
                "raw_support_sizes_by_shift":supports.iter().map(BTreeSet::len).collect::<Vec<_>>(),
                "raw_support_union":union.len(),"raw_target_present":supports[0].contains(&target_key),
                "torsion_target_present":union.contains(&target_key),
                "projected_support":projected_support.len(),"isogeny_mapped_support":mapped_support.len(),
                "chi_orbit_keys":chi_keys.len(),"chi_key_collisions":chi_keys.values().filter(|s|s.len()>1).count(),
                "chi_checks":chi_checks,"group_ops":ops,
            })
        );
    }
}

fn median(mut values: Vec<f64>) -> f64 {
    values.sort_by(f64::total_cmp);
    values[values.len() / 2]
}

fn summarize(path: &str) {
    let mut by_p = BTreeMap::<u64, BTreeMap<u64, BTreeMap<String, serde_json::Value>>>::new();
    for line in std::fs::read_to_string(path).unwrap().lines() {
        let v: serde_json::Value = serde_json::from_str(line).unwrap();
        if v["kind"] != "run" {
            continue;
        }
        let p = v["p"].as_u64().unwrap();
        let round = v["round"].as_u64().unwrap();
        let arm = v["arm"].as_str().unwrap().to_string();
        assert!(by_p
            .entry(p)
            .or_default()
            .entry(round)
            .or_default()
            .insert(arm, v)
            .is_none());
    }
    for (p, rounds) in by_p {
        assert_eq!(rounds.len(), 3);
        let r = rounds[&0]["rho-negation"]["r"].as_u64().unwrap();
        let floor = (std::f64::consts::PI / 4.0).sqrt();
        let mut arms = BTreeSet::new();
        for round in rounds.values() {
            arms.extend(round.keys().cloned());
        }
        for arm in arms {
            let rows: Vec<&serde_json::Value> = rounds.values().map(|row| &row[&arm]).collect();
            let cpu: Vec<f64> = rows.iter().map(|v| v["cpu_ns"].as_f64().unwrap()).collect();
            let s: Vec<f64> = rows
                .iter()
                .map(|v| v["s_cpu_gae"].as_f64().unwrap())
                .collect();
            let counted: Vec<f64> = rows
                .iter()
                .map(|v| v["result"]["s"].as_f64().unwrap())
                .collect();
            let success = rows
                .iter()
                .filter(|v| v["extra"]["status"] == "verified")
                .count();
            let exhausted = rows
                .iter()
                .filter(|v| v["extra"]["status"] == "exhausted")
                .count();
            let paired_neg: Vec<f64> = rounds
                .values()
                .map(|row| {
                    row[&arm]["cpu_ns"].as_f64().unwrap()
                        / row["rho-negation"]["cpu_ns"].as_f64().unwrap()
                })
                .collect();
            let paired_six: Vec<f64> = rounds
                .values()
                .map(|row| {
                    row[&arm]["cpu_ns"].as_f64().unwrap()
                        / row["rho-six-class"]["cpu_ns"].as_f64().unwrap()
                })
                .collect();
            let trials: Vec<f64> = rows
                .iter()
                .map(|v| {
                    v["result"]["decomposition"]["targets_tried"]
                        .as_f64()
                        .unwrap_or(0.0)
                })
                .collect();
            let ranks: Vec<f64> = rows
                .iter()
                .map(|v| {
                    v["result"]["linear_algebra"]["rank"]
                        .as_f64()
                        .unwrap_or(0.0)
                })
                .collect();
            let pairs = rows[0]["result"]["decomposition"]["setup"]["native"]["pair_entries"]
                .as_u64()
                .unwrap_or(0);
            println!(
                "{}",
                json!({
                    "schema":"torsion_extension_summary/v1","source_sha256":sha256_hex(std::fs::read(path).unwrap().as_slice()),
                    "p":p,"r":r,"arm":arm,"n":rows.len(),"verified":success,"exhausted":exhausted,
                    "s_cpu_median":median(s.clone()),"s_cpu_min":s.iter().copied().fold(f64::INFINITY,f64::min),"s_cpu_max":s.iter().copied().fold(0.0,f64::max),
                    "cpu_ns_median":median(cpu),"s_counted_median":median(counted),
                    "ratio_floor":if success==rows.len(){Some(median(s)/floor)}else{None},
                    "paired_ratio_neg_rho":if success==rows.len(){Some(median(paired_neg))}else{None},
                    "paired_ratio_six_rho":if success==rows.len(){Some(median(paired_six))}else{None},
                    "trials_median":median(trials),"rank_median":median(ranks),"pair_entries":pairs,
                    "pair_entry_payload_lower_bytes":pairs*std::mem::size_of::<(u64,usize,usize)>() as u64,
                })
            );
        }
    }
}

fn plot(path: &str) {
    let mut rows = BTreeMap::<(u64, String), serde_json::Value>::new();
    for line in std::fs::read_to_string(path).unwrap().lines() {
        let v: serde_json::Value = serde_json::from_str(line).unwrap();
        rows.insert(
            (v["p"].as_u64().unwrap(), v["arm"].as_str().unwrap().into()),
            v,
        );
    }
    let mut svg = String::new();
    writeln!(svg,"<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 1300 1050\" role=\"img\" aria-label=\"Extension field DLP cold CPU cost by arm and subgroup order\">").unwrap();
    writeln!(svg,"<rect width=\"1300\" height=\"1050\" fill=\"#ffffff\"/><style>text{{font-family:Arial,Helvetica,sans-serif;fill:#17202a}} .title{{font-size:25px;font-weight:bold}} .group{{font-size:19px;font-weight:bold}} .label{{font-size:15px}} .small{{font-size:13px;fill:#4b5563}} .tick{{font-size:15px;fill:#4b5563}}</style>").unwrap();
    writeln!(svg,"<text x=\"30\" y=\"39\" class=\"title\">Cold extension-field DLP cost, one target per curve</text>").unwrap();
    writeln!(svg,"<text x=\"30\" y=\"65\" class=\"small\">Median and range over three paired rounds · S = whole-solve thread CPU / calibrated source addition / sqrt(r) · Darwin arm64 L0</text>").unwrap();
    let x = |s: f64| 235.0 + 255.0 * s.log10();
    for tick in [1.0f64, 10.0, 100.0, 1000.0] {
        let xx = x(tick);
        writeln!(svg,"<line x1=\"{xx:.1}\" y1=\"92\" x2=\"{xx:.1}\" y2=\"941\" stroke=\"#e5e7eb\" stroke-width=\"1\"/>").unwrap();
        writeln!(
            svg,
            "<text x=\"{xx:.1}\" y=\"969\" text-anchor=\"middle\" class=\"tick\">{tick:.0}</text>"
        )
        .unwrap();
    }
    let floor = x((std::f64::consts::PI / 4.0).sqrt());
    writeln!(svg,"<line x1=\"{floor:.1}\" y1=\"92\" x2=\"{floor:.1}\" y2=\"941\" stroke=\"#b91c1c\" stroke-dasharray=\"4 4\"/><text x=\"{floor:.1}\" y=\"987\" text-anchor=\"middle\" class=\"small\">two-class floor</text>").unwrap();
    let arms = [
        ("rho-negation", "rho negation", "#2563eb"),
        ("rho-negation-aa", "rho A/A", "#60a5fa"),
        ("rho-six-class", "rho six class", "#1e3a8a"),
        ("projected-pair", "projected pair", "#0f766e"),
        ("chi-orbit-index", "chi orbit index", "#ea580c"),
        ("isogeny-projected-pair", "isogeny image", "#7c3aed"),
        ("raw-pair", "raw pair", "#6b7280"),
    ];
    for (group, (p, r)) in [(61, 19), (271, 61), (787, 211), (3319, 811)]
        .iter()
        .enumerate()
    {
        let top = 105.0 + group as f64 * 210.0;
        writeln!(
            svg,
            "<text x=\"30\" y=\"{top:.1}\" class=\"group\">F_{}² · r = {}</text>",
            p, r
        )
        .unwrap();
        for (i, (id, label, color)) in arms.iter().enumerate() {
            let v = &rows[&(*p, id.to_string())];
            let y = top + 22.0 + i as f64 * 23.0;
            let lo = x(v["s_cpu_min"].as_f64().unwrap());
            let hi = x(v["s_cpu_max"].as_f64().unwrap());
            let med = x(v["s_cpu_median"].as_f64().unwrap());
            let correct = v["verified"].as_u64().unwrap() == 3;
            let fill = if correct { *color } else { "#ffffff" };
            writeln!(
                svg,
                "<text x=\"205\" y=\"{:.1}\" text-anchor=\"end\" class=\"label\">{label}</text>",
                y + 5.0
            )
            .unwrap();
            writeln!(svg,"<line x1=\"{lo:.1}\" y1=\"{y:.1}\" x2=\"{hi:.1}\" y2=\"{y:.1}\" stroke=\"{color}\" stroke-width=\"3\"/><line x1=\"{lo:.1}\" y1=\"{:.1}\" x2=\"{lo:.1}\" y2=\"{:.1}\" stroke=\"{color}\"/><line x1=\"{hi:.1}\" y1=\"{:.1}\" x2=\"{hi:.1}\" y2=\"{:.1}\" stroke=\"{color}\"/>",y-4.0,y+4.0,y-4.0,y+4.0).unwrap();
            writeln!(svg,"<circle cx=\"{med:.1}\" cy=\"{y:.1}\" r=\"5.5\" fill=\"{fill}\" stroke=\"{color}\" stroke-width=\"2\"/>").unwrap();
            let suffix = if correct { "" } else { " (cap)" };
            writeln!(
                svg,
                "<text x=\"{:.1}\" y=\"{:.1}\" class=\"small\">{:.1}{suffix}</text>",
                hi + 10.0,
                y + 4.0,
                v["s_cpu_median"].as_f64().unwrap()
            )
            .unwrap();
        }
    }
    writeln!(svg,"<text x=\"665\" y=\"1012\" text-anchor=\"middle\" class=\"label\">Whole-solve CPU in calibrated group-addition equivalents per sqrt(r), logarithmic axis</text>").unwrap();
    writeln!(svg,"<text x=\"30\" y=\"1037\" class=\"small\">Source: SUMMARY_V3.jsonl · hollow raw-pair points spent their cap without a scalar · no isolated-host wall claim</text></svg>").unwrap();
    print!("{svg}");
}

fn identities(path: &str) {
    for line in std::fs::read_to_string(path).unwrap().lines() {
        let v: serde_json::Value = serde_json::from_str(line).unwrap();
        let p = v["p"].as_u64().unwrap();
        let nu = v["nu"].as_u64().unwrap();
        let r = v["r"].as_u64().unwrap();
        let h = v["cofactor"].as_u64().unwrap();
        let n2 = v["n2"].as_u64().unwrap();
        let f = Fp2::new(p).unwrap();
        assert_eq!(f.nu, nu);
        assert_eq!(n2, h * r);
        let curve = Curve {
            f,
            a: f.zero(),
            b: f.one(),
        };
        let generator = parse_point(&v["generator"]);
        let dest = Curve {
            f,
            a: f.zero(),
            b: f.from_fp(p - 27),
        };
        let dest_generator = phi(&curve, generator);
        let mut ops = GroupOps::default();
        assert!(curve.is_on_curve(generator));
        assert!(curve.mul(&mut ops, generator, r).infinity);
        assert!(dest.is_on_curve(dest_generator));
        assert!(dest.mul(&mut ops, dest_generator, r).infinity);
        let field = json!({
            "characteristic":p,"degree":2,"representation":"polynomial_basis",
            "defining_polynomial":[p-nu,0,1],"element_encoding":"two_u64_coefficients_le",
        });
        for (name, b, g, id) in [
            ("source", 1, generator, v["source"]["id"].as_str().unwrap()),
            (
                "destination",
                p - 27,
                dest_generator,
                v["destination"]["id"].as_str().unwrap(),
            ),
        ] {
            let model = json!({
                "model":"short_weierstrass","a":[0,0],"b":[b,0],
                "group_order":n2,"subgroup_order":r,"cofactor":h,
                "generator":{"x":g.x,"y":g.y},"target_group":"cyclic_prime_order_subgroup",
            });
            let input = json!({"field":field,"curve":model});
            let sha = sha256_hex(canonical(&input).unwrap().as_bytes());
            let field_sha = sha256_hex(canonical(&field).unwrap().as_bytes());
            let alias = format!("EC1Q{p}D2C{name}h{}", &sha[..12]);
            println!(
                "{}",
                json!({
                    "schema":"torsion_extension_identity/v1","p":p,"name":name,
                    "icv1":id,"ec1_alias":alias,"curve_uid":format!("urn:ec-record:1:sha256:{sha}"),
                    "curve_sha256":sha,"field_sha256":field_sha,"field":field,"curve":model,
                })
            );
        }
    }
}

// A deliberately separate arithmetic path for checking frozen answers. None
// of these operations call the run's Fp2, Curve, GroupOps, or scalar routine.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct CheckF2(u64, u64);

type CheckPoint = Option<(CheckF2, CheckF2)>;

#[derive(Clone, Copy)]
struct CheckField {
    p: u64,
    nu: u64,
}

impl CheckField {
    fn residue(self, x: u128) -> u64 {
        (x % self.p as u128) as u64
    }
    fn add_fp(self, a: u64, b: u64) -> u64 {
        self.residue(a as u128 + b as u128)
    }
    fn sub_fp(self, a: u64, b: u64) -> u64 {
        self.residue(a as u128 + self.p as u128 - b as u128)
    }
    fn mul_fp(self, a: u64, b: u64) -> u64 {
        self.residue(a as u128 * b as u128)
    }
    fn pow_fp(self, mut a: u64, mut n: u64) -> u64 {
        let mut result = 1;
        while n > 0 {
            if n & 1 != 0 {
                result = self.mul_fp(result, a);
            }
            a = self.mul_fp(a, a);
            n >>= 1;
        }
        result
    }
    fn add(self, a: CheckF2, b: CheckF2) -> CheckF2 {
        CheckF2(self.add_fp(a.0, b.0), self.add_fp(a.1, b.1))
    }
    fn sub(self, a: CheckF2, b: CheckF2) -> CheckF2 {
        CheckF2(self.sub_fp(a.0, b.0), self.sub_fp(a.1, b.1))
    }
    fn mul(self, a: CheckF2, b: CheckF2) -> CheckF2 {
        CheckF2(
            self.add_fp(
                self.mul_fp(a.0, b.0),
                self.mul_fp(self.nu, self.mul_fp(a.1, b.1)),
            ),
            self.add_fp(self.mul_fp(a.0, b.1), self.mul_fp(a.1, b.0)),
        )
    }
    fn inv(self, a: CheckF2) -> CheckF2 {
        let norm = self.sub_fp(
            self.mul_fp(a.0, a.0),
            self.mul_fp(self.nu, self.mul_fp(a.1, a.1)),
        );
        assert_ne!(norm, 0, "inverse of zero");
        let scale = self.pow_fp(norm, self.p - 2);
        CheckF2(
            self.mul_fp(a.0, scale),
            self.mul_fp(self.sub_fp(0, a.1), scale),
        )
    }
    fn on_curve(self, point: CheckPoint) -> bool {
        match point {
            None => true,
            Some((x, y)) => self.mul(y, y) == self.add(self.mul(self.mul(x, x), x), CheckF2(1, 0)),
        }
    }
    fn add_points(self, a: CheckPoint, b: CheckPoint) -> CheckPoint {
        let Some((x1, y1)) = a else {
            return b;
        };
        let Some((x2, y2)) = b else {
            return a;
        };
        let slope = if x1 == x2 {
            if y1 != y2 || y1 == CheckF2(0, 0) {
                return None;
            }
            let x_sq = self.mul(x1, x1);
            self.mul(
                self.mul(CheckF2(3, 0), x_sq),
                self.inv(self.mul(CheckF2(2, 0), y1)),
            )
        } else {
            self.mul(self.sub(y2, y1), self.inv(self.sub(x2, x1)))
        };
        let x3 = self.sub(self.sub(self.mul(slope, slope), x1), x2);
        let y3 = self.sub(self.mul(slope, self.sub(x1, x3)), y1);
        Some((x3, y3))
    }
    fn scalar(self, mut point: CheckPoint, mut k: u64) -> CheckPoint {
        let mut result = None;
        while k > 0 {
            if k & 1 != 0 {
                result = self.add_points(result, point);
            }
            point = self.add_points(point, point);
            k >>= 1;
        }
        result
    }
}

fn check_point(v: &serde_json::Value) -> CheckPoint {
    if v["infinity"].as_bool().unwrap() {
        return None;
    }
    Some((
        CheckF2(v["x"][0].as_u64().unwrap(), v["x"][1].as_u64().unwrap()),
        CheckF2(v["y"][0].as_u64().unwrap(), v["y"][1].as_u64().unwrap()),
    ))
}

fn verify_independent(run_path: &str, fixture_path: &str) {
    let fixture_bytes = std::fs::read(fixture_path).expect("fixture bytes");
    let run_bytes = std::fs::read(run_path).expect("run bytes");
    assert_eq!(
        sha256_hex(&fixture_bytes),
        "1a63cc22a22e0fe27d44359325d502038d26d3ac616374fc5558c041939218b7",
        "fixture bytes differ from the frozen input"
    );
    assert_eq!(
        sha256_hex(&run_bytes),
        "ca9a0ec7e31de59c58a78ae5c269350b37963b62ff9de815c805c36c2ec26693",
        "run bytes differ from the frozen V3 result"
    );
    let fixtures: Vec<serde_json::Value> = String::from_utf8(fixture_bytes.clone())
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert_eq!(fixtures.len(), 4);
    let arms = [
        "raw-pair",
        "projected-pair",
        "chi-orbit-index",
        "isogeny-projected-pair",
        "rho-negation",
        "rho-negation-aa",
        "rho-six-class",
    ];
    let mut seen = BTreeSet::new();
    let mut verified = 0;
    let mut exhausted = 0;
    for line in String::from_utf8(run_bytes.clone()).unwrap().lines() {
        let v: serde_json::Value = serde_json::from_str(line).unwrap();
        if v["kind"] == "calibration" {
            continue;
        }
        assert_eq!(v["kind"], "run");
        let p = v["p"].as_u64().unwrap();
        let round = v["round"].as_u64().unwrap();
        let arm = v["arm"].as_str().unwrap();
        assert!(round < 3 && arms.contains(&arm));
        assert!(
            seen.insert((p, round, arm.to_owned())),
            "duplicate run slot"
        );
        let fix = fixtures
            .iter()
            .find(|f| f["p"].as_u64() == Some(p))
            .expect("unknown field");
        let r = fix["r"].as_u64().unwrap();
        assert_eq!(v["r"].as_u64(), Some(r));
        assert_eq!(
            v["seed"].as_u64(),
            Some((0x2026_1014u64 ^ p).wrapping_add(round * 0x9e37_79b9))
        );
        let nu = fix["nu"].as_u64().unwrap();
        let field = CheckField { p, nu };
        assert!(p > 3 && (2..).take_while(|d| d * d <= p).all(|d| p % d != 0));
        assert_eq!(field.pow_fp(nu, (p - 1) / 2), p - 1);
        let g = check_point(&fix["generator"]);
        let q = check_point(&fix["target"]);
        assert!(g.is_some() && q.is_some());
        assert!(field.on_curve(g) && field.on_curve(q));
        assert_eq!(field.scalar(g, r), None, "generator order: p={p}");
        assert_eq!(field.scalar(q, r), None, "target order: p={p}");
        if let Some(k) = v["result"]["recovered"].as_u64() {
            assert!(k < r);
            assert_eq!(field.scalar(g, k), q, "independent scalar replay: p={p}");
            assert_eq!(k, fix["answer_scalar"].as_u64().unwrap());
            assert_eq!(v["extra"]["status"], "verified");
            assert_eq!(v["result"]["verified"], true);
            verified += 1;
        } else {
            assert_eq!(arm, "raw-pair");
            assert_eq!(v["extra"]["status"], "exhausted");
            assert_eq!(v["result"]["exhausted"], true);
            assert_eq!(v["result"]["verified"], false);
            exhausted += 1;
        }
    }
    assert_eq!(seen.len(), fixtures.len() * 3 * arms.len());
    for fix in &fixtures {
        let p = fix["p"].as_u64().unwrap();
        for round in 0..3 {
            for arm in &arms {
                assert!(seen.contains(&(p, round, (*arm).to_owned())));
            }
        }
    }
    assert_eq!((verified, exhausted), (75, 9));
    println!(
        "{}",
        json!({
            "schema":"torsion_extension_independent_replay/v2",
            "arithmetic":"independent_u128_fp2_chord_tangent",
            "fixture_sha256":sha256_hex(&fixture_bytes),
            "records_sha256":sha256_hex(&run_bytes),
            "slots":seen.len(),"verified":verified,"exhausted":exhausted,
        })
    );
}

fn iso_scale(curve: &Curve, u: u64, point: Point) -> Point {
    if point.infinity {
        return curve.identity();
    }
    let f = &curve.f;
    let u2 = (u as u128 * u as u128 % f.p as u128) as u64;
    let u3 = (u2 as u128 * u as u128 % f.p as u128) as u64;
    Point::affine(f.scale(point.x, u2), f.scale(point.y, u3))
}

fn iso_check(path: &str) {
    let bytes = std::fs::read(path).expect("frozen fixture bytes");
    let fixture_sha256 = sha256_hex(&bytes);
    for line in String::from_utf8(bytes.clone()).unwrap().lines() {
        let v: serde_json::Value = serde_json::from_str(line).unwrap();
        let p = v["p"].as_u64().unwrap();
        let r = v["r"].as_u64().unwrap();
        let f = Fp2::new(p).unwrap();
        let src = Curve {
            f,
            a: f.zero(),
            b: f.one(),
        };
        let dst = Curve {
            f,
            a: f.zero(),
            b: f.from_fp(p - 27),
        };
        let u = (1..p)
            .find(|u| ((*u as u128 * *u as u128) % p as u128) as u64 == p - 3)
            .expect("no base-field scaling isomorphism");
        let u2 = (u as u128 * u as u128 % p as u128) as u64;
        let u3 = (u2 as u128 * u as u128 % p as u128) as u64;
        let u6 = (u3 as u128 * u3 as u128 % p as u128) as u64;
        assert_eq!(u6, (p - 27) % p);
        let g = parse_point(&v["generator"]);
        let q = parse_point(&v["target"]);
        let mut ops = GroupOps::default();
        for point in std::iter::once(g)
            .chain(std::iter::once(q))
            .chain(v["base_raw"].as_array().unwrap().iter().map(parse_point))
        {
            assert!(src.is_on_curve(point));
            assert!(dst.is_on_curve(iso_scale(&src, u, point)));
        }
        assert!(dst.mul(&mut ops, iso_scale(&src, u, g), r).infinity);
        assert!(dst.mul(&mut ops, iso_scale(&src, u, q), r).infinity);
        let seed = v["target_seed"].as_str().unwrap();
        let seed_sha = sha256_hex(seed.as_bytes());
        let mut rng = StdRng::seed_from_u64(u64::from_str_radix(&seed_sha[..16], 16).unwrap());
        for _ in 0..64 {
            let a = random_point(&src, &mut rng);
            let b = random_point(&src, &mut rng);
            let lhs = iso_scale(&src, u, src.add(&mut ops, a, b));
            let rhs = dst.add(&mut ops, iso_scale(&src, u, a), iso_scale(&src, u, b));
            assert_eq!(lhs, rhs, "scaling homomorphism p={p}");
        }
        println!(
            "{}",
            json!({
                "schema":"torsion_extension_iso_check/v1", "p":p,"r":r,"u":u,
                "u2":u2,"u3":u3,"u6":u6,
                "source_icv1":v["source"]["id"],"destination_icv1":v["destination"]["id"],
                "fixture_sha256":fixture_sha256,"pair_checks":64,
                "fixed_point_checks":v["base_raw"].as_array().unwrap().len()+2,
                "group_ops":ops,
            })
        );
    }
}

fn base_order(p: u64) -> u64 {
    let mut n = 1;
    for x in 0..p {
        let rhs = ((x as u128 * x as u128 % p as u128 * x as u128 + 1) % p as u128) as u64;
        n += (1 + jacobi(rhs, p)) as u64;
    }
    n
}

fn embedding_degree(p: u64, r: u64) -> u64 {
    let mut power = p % r;
    for k in 1..=r {
        if power == 1 {
            return k;
        }
        power = ((power as u128 * p as u128) % r as u128) as u64;
    }
    unreachable!("r must be prime and coprime to p")
}

fn screen() {
    let intervals = [
        (31, 127, 17),
        (128, 511, 53),
        (512, 2047, 211),
        (2048, 4093, 811),
    ];
    for (lo, hi, min_r) in intervals {
        let mut selected = None;
        for p in lo..=hi {
            if p % 3 != 1 || !is_prime_u64(p) {
                continue;
            }
            let n1 = base_order(p);
            let t = p as i64 + 1 - n1 as i64;
            let anti = (p as i64 + 1 + t) as u64;
            let n2 = n1 * anti;
            let r = factor_u64(anti).last().map(|&(v, _)| v).unwrap_or(0);
            if r < min_r || r == 3 || n1 % r == 0 || (n2 / r) % r == 0 {
                continue;
            }
            let f = Fp2::new(p).unwrap();
            let curve = ExtCurve {
                f,
                a: f.zero(),
                b: f.one(),
            };
            let h = n2 / r;
            let mut rng = StdRng::seed_from_u64(0x2026_1014 ^ p);
            let mut ops = GroupOps::default();
            let mut generator = None;
            for _ in 0..64 {
                let raw = random_point(&curve, &mut rng);
                let g = curve.mul(&mut ops, raw, h);
                if !g.infinity && (g.x[1] != 0 || g.y[1] != 0) && curve.mul(&mut ops, g, r).infinity
                {
                    generator = Some(g);
                    break;
                }
            }
            let Some(g) = generator else {
                continue;
            };
            let mut line = None;
            for a in 0..p {
                let c = [a, 1];
                let mut projected = BTreeSet::new();
                let mut raw_points = 0u64;
                for s in 1..p {
                    let x = f.scale(c, s);
                    for point in curve.lift_x(x) {
                        raw_points += 1;
                        let q = curve.mul(&mut ops, point, h);
                        if !q.infinity {
                            projected.insert(curve.key(&q));
                        }
                    }
                    if projected.len() >= 16 {
                        break;
                    }
                }
                if projected.len() >= 16 {
                    line = Some((c, raw_points, projected.len()));
                    break;
                }
            }
            if let Some((c, raw_points, projected_points)) = line {
                selected = Some(json!({
                    "p":p,"nu":f.nu,"n1":n1,"trace":t,"n2":n2,
                    "r":r,"h":h,"embedding_degree":embedding_degree(p,r),
                    "generator":{"x":g.x,"y":g.y},
                    "line":c,"raw_points_inspected":raw_points,
                    "projected_signed_points_seen":projected_points,
                    "group_ops_screen":ops,
                }));
                break;
            }
        }
        println!(
            "{}",
            json!({"interval":[lo,hi],"min_r":min_r,"selected":selected})
        );
    }
}

fn main() {
    match std::env::args().nth(1).as_deref() {
        Some("screen") => screen(),
        Some("freeze") => freeze(),
        Some("solve") => solve(&std::env::args().nth(2).expect("fixture path")),
        Some("verify") => verify(
            &std::env::args().nth(2).expect("records path"),
            &std::env::args().nth(3).expect("fixture path"),
        ),
        Some("support") => support(&std::env::args().nth(2).expect("fixture path")),
        Some("summarize") => summarize(&std::env::args().nth(2).expect("run path")),
        Some("plot") => plot(&std::env::args().nth(2).expect("summary path")),
        Some("identities") => identities(&std::env::args().nth(2).expect("fixture path")),
        Some("verify-independent") => verify_independent(
            &std::env::args().nth(2).expect("records path"),
            &std::env::args().nth(3).expect("fixture path"),
        ),
        Some("iso-check") => iso_check(&std::env::args().nth(2).expect("fixture path")),
        _ => panic!("usage: isogeny_torsion_extension screen"),
    }
}
