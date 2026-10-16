//! One-hop prime61 isogeny transport and matched capped factor-base runs.
//! Usage: prime61_isogeny_experiment ROUTES.json [trials=100000] [runs=12]
//! The route degree bound is coverage, not an exhaustive class certificate.
use std::{collections::BTreeMap, fs};

use crypto_lib::cryptanalysis::{
    isogeny_degree_search::class_number::{factor_u64, is_prime_u64},
    isogeny_eval::{
        counted_fp::{NoCount, PrimeField},
        small_degree::{kernel_x_coordinates, ShortWeierstrass},
    },
    isogeny_walk::{
        curve::Model,
        field::{Fe, Field},
        kernel,
        poly::{self, Poly},
        record::{ec1_identity, ec1_records, sha256_hex},
    },
    prime_fast::{collect_relations, solve_relations, FastCurve, FastFactorBase, Pt},
};
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde_json::{json, Value};

const P: u64 = (1u64 << 61) - 1;
const R: u64 = 2305843010818082053;
const SECRET: u64 = 539365399;
const G: (u64, u64) = (1, 1796364306280342617);
const Q: (u64, u64) = (992726311114892094, 558260110083040833);

fn number(v: &Value) -> u64 {
    v.as_u64()
        .unwrap_or_else(|| v.as_str().unwrap().parse().unwrap())
}

fn word(f: &Field, x: &Fe) -> u64 {
    f.to_big(x).to_u64().unwrap()
}

/// Known Vélu formula, with sums evaluated as traces in F_p[T]/h(T).
/// This handles Frobenius-stable kernels without rational kernel points.
struct KernelMap<'a> {
    f: &'a Field,
    h: Poly,
    v: Poly,
    u: Poly,
    sums: Vec<Fe>,
    u2: Fe,
    u3: Fe,
}

impl<'a> KernelMap<'a> {
    fn new(f: &'a Field, source: &Model, target: &Model, ell: u64, h: Poly) -> Self {
        let codomain = kernel::verify_kernel(f, source, ell, &h).expect("kernel certificate");
        let u2 = kernel::link_to_target(f, &codomain, target).expect("codomain isomorphism");
        let u_root = f.sqrt(&u2).expect("u2 is square");
        let u3 = f.mul(&u2, &u_root);
        let v = poly::rem(
            f,
            &vec![f.mul(&f.from_u64(2), &source.a), f.zero(), f.from_u64(6)],
            &h,
        );
        let u = poly::rem(
            f,
            &vec![
                f.mul(&f.from_u64(4), &source.b),
                f.mul(&f.from_u64(4), &source.a),
                f.zero(),
                f.from_u64(4),
            ],
            &h,
        );
        let d = h.len() - 1;
        let mut sums = vec![f.from_u64(d as u64)];
        sums.extend(poly::power_sums(f, &h, d.saturating_sub(1)));
        Self {
            f,
            h,
            v,
            u,
            sums,
            u2,
            u3,
        }
    }

    fn trace(&self, a: &Poly) -> Fe {
        a.iter().zip(&self.sums).fold(self.f.zero(), |acc, (a, s)| {
            self.f.add(&acc, &self.f.mul(a, s))
        })
    }

    fn eval(&self, point: Option<(u64, u64)>) -> Option<(u64, u64)> {
        let (x, y) = point?;
        let f = self.f;
        let x = f.from_u64(x);
        let y = f.from_u64(y);
        let z = poly::invmod(f, &vec![x, f.neg(&f.one())], &self.h)?;
        let z2 = poly::mulmod(f, &z, &z, &self.h);
        let z3 = poly::mulmod(f, &z2, &z, &self.h);
        let dx = poly::add(
            f,
            &poly::mulmod(f, &self.v, &z, &self.h),
            &poly::mulmod(f, &self.u, &z2, &self.h),
        );
        let dy = poly::add(
            f,
            &poly::mulmod(f, &self.v, &z2, &self.h),
            &poly::scale(f, &poly::mulmod(f, &self.u, &z3, &self.h), &f.from_u64(2)),
        );
        let xx = f.mul(&self.u2, &f.add(&x, &self.trace(&dx)));
        let yy = f.mul(&self.u3, &f.mul(&y, &f.sub(&f.one(), &self.trace(&dy))));
        Some((word(f, &xx), word(f, &yy)))
    }
}

fn fast(a: u64, b: u64) -> FastCurve {
    FastCurve::new("", P, a, b, R)
}

fn pair_support(c: &FastCurve, pts: &[Pt]) -> Vec<Option<(u64, u64)>> {
    let mut signed: Vec<_> = pts.iter().flat_map(|&p| [p, c.neg(p)]).collect();
    signed.sort_by_key(|&p| c.canonical(p));
    signed
        .iter()
        .flat_map(|&a| {
            signed
                .iter()
                .map(move |&b| c.add(Some(a), Some(b)).map(|p| c.canonical(p)))
        })
        .collect()
}

fn run_curve(
    id: &str,
    ell: u64,
    c: &FastCurve,
    g: (u64, u64),
    q: (u64, u64),
    trials: u64,
    runs: u64,
) {
    let (g, q) = (c.point(g.0, g.1), c.point(q.0, q.1));
    assert!(c.is_on_curve(g) && c.is_on_curve(q));
    assert!(c.mul(Some(g), R).is_none() && c.mul(Some(q), R).is_none());
    assert_eq!(c.mul(Some(g), SECRET), Some(q));
    let base = FastFactorBase::build(c, 3);
    for &p in &base.pts {
        assert!(c.is_on_curve(p) && c.mul(Some(p), R).is_none());
    }
    let (gx, gy) = c.canonical(g);
    let (field_record, curve_record) = ec1_records(
        &BigUint::from(P),
        &BigUint::from(c.f.from_mont(c.a)),
        &BigUint::from(c.f.from_mont(c.b)),
        &BigUint::from(R),
        &BigUint::from(1u8),
        (&BigUint::from(gx), &BigUint::from(gy)),
    );
    let (ec1, uid) = ec1_identity(&BigUint::from(P), &field_record, &curve_record, "fp");
    println!(
        "{}",
        json!({"kind":"curve", "edge_id":id, "degree":ell,
        "measured_ec1":ec1,"measured_curve_uid":uid,"generator_policy":"image of frozen source G under the recorded map",
        "a":c.f.from_mont(c.a).to_string(),"b":c.f.from_mont(c.b).to_string(),
        "transported_G":([g.x,g.y].map(|x|c.f.from_mont(x).to_string())),
        "transported_Q":([q.x,q.y].map(|x|c.f.from_mont(x).to_string())),
        "base":base.pts.iter().map(|&p| {let (x,y)=c.canonical(p); [x.to_string(),y.to_string()]}).collect::<Vec<_>>() })
    );
    for run in 0..runs {
        let seed = 2026101600 + run;
        let (rels, stats) = collect_relations(c, g, q, &base, None, 20, trials, seed);
        for rel in &rels {
            let mut rhs = None;
            for &(i, sign) in &rel.entries {
                let pt = base.pts[i as usize];
                let pt = if sign < 0 { c.neg(pt) } else { pt };
                rhs = c.add(rhs, c.mul(Some(pt), sign.unsigned_abs() as u64));
            }
            assert_eq!(c.combine(g, q, rel.a, rel.b), rhs);
        }
        assert_eq!(stats.rejected, 0);
        let recovered = solve_relations(c, 3, &rels);
        let verified = recovered.is_some_and(|k| k == SECRET && c.mul(Some(g), k) == Some(q));
        if recovered.is_some() {
            assert!(verified);
        }
        println!(
            "{}",
            json!({"kind":"run","edge_id":id,"degree":ell,"run":run,"seed":seed,
            "targets":stats.targets,"relations":stats.relations,"rejected":stats.rejected,
            "replayed_relations":rels.len(),"scalar":recovered.map(|k|k.to_string()),"verified":verified,
            "status":if verified {"VERIFIED_SOLVE"} else if stats.targets==trials {"EXHAUSTED_TRIAL_CAP"} else {"UNRESOLVED_RANK"}})
        );
    }
}

fn main() {
    let argv: Vec<_> = std::env::args().collect();
    let bytes = fs::read(argv.get(1).expect("ROUTES.json [trials] [runs]")).unwrap();
    let routes: Value = serde_json::from_slice(&bytes).unwrap();
    let trials: u64 = argv.get(2).map(|x| x.parse().unwrap()).unwrap_or(100000);
    let runs: u64 = argv.get(3).map(|x| x.parse().unwrap()).unwrap_or(12);
    assert!(trials > 0 && runs > 0);
    let d = (-1604388101i128).pow(2) - 4 * P as i128;
    let factors = factor_u64((-d).try_into().unwrap());
    let mut product = 1u128;
    for &(p, e) in &factors {
        assert!(is_prime_u64(p));
        product *= (p as u128).pow(e);
    }
    assert_eq!(product, (-d) as u128);
    assert!(factors.iter().all(|&(_, e)| e == 1) && d.rem_euclid(4) == 1);
    assert!(is_prime_u64(P) && is_prime_u64(R));
    let cf = PrimeField::<NoCount>::new(P).unwrap();
    let sw = ShortWeierstrass::new(&cf, P - 3, 111).unwrap();
    assert!(kernel_x_coordinates(&cf, &sw, 2).unwrap().is_empty());
    println!(
        "{}",
        json!({"kind":"manifest","schema":"prime61.horizontal-factor-base/v1",
        "routes_sha256":sha256_hex(std::str::from_utf8(&bytes).unwrap()),"p":P.to_string(),"r":R.to_string(),"trace":"-1604388101",
        "discriminant":d.to_string(),"factorization":factors,"all_factors_deterministically_prime":true,
        "frobenius_conductor":1,"endomorphism_conductor":1,"vertical_Fp_edges_exist":false,
        "degree2_Fp_edges":0,"coverage":"one-hop edges in supplied route file; not the entire horizontal class",
        "runs_per_curve":runs,"max_targets_per_run":trials,"inverse_pairs":3,
        "full_operation_cost":null,"speedup":null})
    );
    let f = Field::new(&BigUint::from(P)).unwrap();
    let source = fast(P - 3, 111);
    let model = Model {
        a: f.from_u64(P - 3),
        b: f.from_u64(111),
    };
    let root = &routes["curve_nodes"][0];
    assert_eq!(
        root["coefficients_a1_a2_a3_a4_a6"][3],
        json!((P - 3).to_string())
    );
    assert_eq!(root["coefficients_a1_a2_a3_a4_a6"][4], json!("111"));
    let root_ref = root["ref"].as_str().unwrap();
    let nodes: BTreeMap<_, _> = routes["curve_nodes"]
        .as_array()
        .unwrap()
        .iter()
        .map(|n| (n["ref"].as_str().unwrap(), n))
        .collect();
    run_curve("identity", 1, &source, G, Q, trials, runs);
    let base = FastFactorBase::build(&source, 3);
    for edge in routes["edges"].as_array().unwrap() {
        assert_eq!(
            edge["source_curve_ref"],
            json!(root_ref),
            "one-hop input required"
        );
        assert_eq!(edge["direction"], json!("h"));
        let ell = number(&edge["degree"]);
        let node = nodes[edge["target_curve_ref"].as_str().unwrap()];
        assert_eq!(number(&node["field_characteristic"]), P);
        assert_eq!(number(&node["group_order"]), R);
        let (a, b) = (
            number(&node["coefficients_a1_a2_a3_a4_a6"][3]),
            number(&node["coefficients_a1_a2_a3_a4_a6"][4]),
        );
        let target = fast(a, b);
        let tm = Model {
            a: f.from_u64(a),
            b: f.from_u64(b),
        };
        let h: Poly = edge["kernel_polynomial"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| f.from_u64(number(x)))
            .collect();
        let map = KernelMap::new(&f, &model, &tm, ell, h);
        let (mg, mq) = (map.eval(Some(G)).unwrap(), map.eval(Some(Q)).unwrap());
        for k in 0..64 {
            let p = source
                .mul(Some(source.point(G.0, G.1)), k)
                .map(|p| source.canonical(p));
            let want = target
                .mul(Some(target.point(mg.0, mg.1)), k)
                .map(|p| target.canonical(p));
            assert_eq!(map.eval(p), want, "scalar transport");
        }
        let transported: Vec<Pt> = base
            .pts
            .iter()
            .map(|&p| {
                let q = map.eval(Some(source.canonical(p))).unwrap();
                target.point(q.0, q.1)
            })
            .collect();
        let src_sums = pair_support(&source, &base.pts);
        let mut images: Vec<_> = src_sums.into_iter().map(|p| map.eval(p)).collect();
        let mut dst_sums = pair_support(&target, &transported);
        images.sort();
        dst_sums.sort();
        assert_eq!(images, dst_sums, "transported pair-support control");
        println!(
            "{}",
            json!({"kind":"transport_check","edge_id":edge["id"],"degree":ell,
            "kernel_verified":true,"scalar_checks":64,"transported_pair_sums_checked":36,
            "map_u2":word(&f,&map.u2).to_string(),"map_u3":word(&f,&map.u3).to_string(),"status":"PASS"})
        );
        run_curve(
            edge["id"].as_str().unwrap(),
            ell,
            &target,
            mg,
            mq,
            trials,
            runs,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_lib::cryptanalysis::isogeny_eval::small_degree::SmallIsogeny;
    #[test]
    fn quotient_trace_evaluation_matches_existing_degree3_map_including_pole() {
        let f = Field::new(&BigUint::from(P)).unwrap();
        let pf = PrimeField::<NoCount>::new(P).unwrap();
        let e = ShortWeierstrass::new(&pf, P - 3, 111).unwrap();
        let x0 = kernel_x_coordinates(&pf, &e, 3).unwrap()[0];
        let small = SmallIsogeny::new(&pf, &e, 3, x0).unwrap();
        let m = Model {
            a: f.from_u64(P - 3),
            b: f.from_u64(111),
        };
        let cm = small.codomain();
        let target = Model {
            a: f.from_u64(cm.a),
            b: f.from_u64(cm.b),
        };
        let map = KernelMap::new(&f, &m, &target, 3, vec![f.neg(&f.from_u64(x0)), f.one()]);
        for x in 1..128 {
            if let Some(y) = pf.sqrt(e.rhs(&pf, x)) {
                let raw = small.eval(&pf, Some((x, y)));
                let got = map.eval(Some((x, y)));
                assert_eq!(got.map(|q| q.0), raw.map(|q| q.0));
                assert!(got == raw || got == raw.map(|(x, y)| (x, pf.neg(y))));
            }
        }
        assert_eq!(map.eval(None), None);
        assert_eq!(map.eval(Some((x0, 1))), None);
    }
}
