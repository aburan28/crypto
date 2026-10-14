//! Native instance screen for the frozen quadratic-extension comparison.
//! `screen` prints candidate parameters only; the paired DLP run is a
//! separate, frozen stage in the study protocol.

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::rc::Rc;
use std::time::Instant;

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
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
                take_field_counters();
                let started = Instant::now();
                let cpu_start = thread_cpu_ns();
                let (result, extra) = match arm {
                    "rho-negation" | "rho-negation-aa" => {
                        let res = rho_reference_negation(&curve, g, q, r, seed, 1_000_000);
                        (
                            json!(res),
                            json!({"status":if res.verified{"verified"}else{"exhausted"}}),
                        )
                    }
                    "rho-six-class" => {
                        let (aut, eigen_ops) = six_class_endomorphism(&curve, g, r);
                        let neg = Negation { r };
                        let gens: [&dyn Endomorphism<Curve>; 2] = [&neg, &aut];
                        let res =
                            rho_reference_folded(&curve, g, q, r, seed, 1_000_000, &gens).unwrap();
                        (
                            json!(res),
                            json!({"eigen_setup_group_ops":eigen_ops,"eigenvalue":aut.eigenvalue,"status":if res.verified{"verified"}else{"exhausted"}}),
                        )
                    }
                    _ => {
                        let projected = Rc::new(RefCell::new(Vec::new()));
                        let base = LineBase {
                            frozen: frozen_raw.clone(),
                            projected: projected.clone(),
                            columns,
                            line: f.s(),
                        };
                        let mut spec = PipelineSpec {
                            factor_base: base.name().into(),
                            oracle: arm.into(),
                            targets: Targets::Walk,
                            linalg: "incremental-gauss".into(),
                            max_trials: cap,
                            seed,
                            ..Default::default()
                        };
                        // All arms receive exactly the same seed, factor base and target.
                        spec.factor_base_params = Params::default();
                        let mut c = calib.clone();
                        let report = if arm == "raw-pair" {
                            let mut oracle = MitmOracle::new(2);
                            run_pipeline(&ctx, &spec, &base, &mut oracle, answers, &c, None)
                                .unwrap()
                        } else {
                            let mode = match arm {
                                "projected-pair" => IndexMode::Projected,
                                "chi-orbit-index" => IndexMode::Chi,
                                _ => IndexMode::Isogeny,
                            };
                            c.ns_per_canon = Some(transform_calibration(&curve, g, mode));
                            let mut oracle = ProjectedOracle::new(mode, projected);
                            let rep =
                                run_pipeline(&ctx, &spec, &base, &mut oracle, answers, &c, None)
                                    .unwrap();
                            assert_eq!(oracle.false_candidates, 0, "invalid projected relation");
                            rep
                        };
                        (
                            json!(report),
                            json!({"status":if report.verified{"verified"}else if report.exhausted{"exhausted"}else{"wrong_answer"}}),
                        )
                    }
                };
                let cpu_ns = thread_cpu_ns() - cpu_start;
                let wall_ns = started.elapsed().as_nanos() as u64;
                let field = take_field_counters();
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
        let rec = if v["arm"].as_str().unwrap().starts_with("rho") {
            v["result"]["recovered"].as_u64()
        } else {
            v["result"]["recovered"].as_u64()
        };
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
        _ => panic!("usage: isogeny_torsion_extension screen"),
    }
}
