//! Certified property sidecar for one immutable synthetic-model inventory.
mod field;
mod identities;
mod presentation;
use crypto_lib::hash::sha256::sha256;
use field::{divisors, Element, Field};
use num_bigint::{BigInt, BigUint};
use num_traits::{One, ToPrimitive, Zero};
use serde_json::{json, Value};
use std::{
    collections::BTreeMap,
    fs,
    path::Path,
    process::{Command, Stdio},
    time::{Duration, Instant},
};
const INPUT: &str = "docs/curves/weak-families/models.json";
const INPUT_SHA: &str = "9738f6bd6c016732f88e4e7979ac09e9aa6fd16eba96220ec8781c0017b6055e";
const STUDY: &str = "research/weak_curve_properties_20261010";
const OUTPUT: &str = "docs/curves/weak-families/properties.json";
fn hash(b: &[u8]) -> String {
    hex::encode(sha256(b))
}
fn load(repo: &Path) -> Value {
    let b = fs::read(repo.join(INPUT)).unwrap();
    assert_eq!(hash(&b), INPUT_SHA, "fixed model input changed");
    let v: Value = serde_json::from_slice(&b).unwrap();
    assert_eq!(v["models"].as_array().unwrap().len(), 814);
    v
}
fn read_json(p: &Path) -> Value {
    serde_json::from_slice(&fs::read(p).unwrap()).unwrap()
}
fn numbers(v: &Value) -> Element {
    v.as_array()
        .unwrap()
        .iter()
        .map(|c| {
            c.as_str()
                .map(|s| s.parse().unwrap())
                .unwrap_or_else(|| c.as_u64().unwrap())
        })
        .collect()
}
fn encoded(m: &Value, f: &Field, k: &str) -> Element {
    let a = numbers(&m[k]);
    assert_eq!(a.len(), f.n);
    assert!(a.iter().all(|a| *a < f.p));
    a
}
fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}
fn valuation(n: &BigUint, p: u32) -> u32 {
    assert!(!n.is_zero());
    let mut n = n.clone();
    let mut k = 0;
    while (&n % p).is_zero() {
        n /= p;
        k += 1;
    }
    k
}
fn local_two(delta: &BigInt) -> (u32, u32, &'static str) {
    assert!(*delta < BigInt::zero());
    let abs = (-delta).to_biguint().unwrap();
    let s = valuation(&abs, 2);
    let odd = &abs >> s as usize;
    let residue = (8 - (&odd % 8u32).to_u64().unwrap()) % 8;
    let dk2 = if s % 2 == 1 {
        3
    } else if residue % 4 == 1 {
        0
    } else {
        2
    };
    assert!(
        s >= dk2 && (s - dk2) % 2 == 0,
        "invalid quadratic-order discriminant"
    );
    let split = if dk2 != 0 {
        "ramified"
    } else if residue == 1 {
        "split"
    } else {
        assert_eq!(residue, 5);
        "inert"
    };
    ((s - dk2) / 2, dk2, split)
}
fn endomorphism_two(depth: u64, r2: usize, full4: bool) -> Value {
    assert!(!full4 || r2 == 3);
    let (lo, hi, status) = if r2 == 1 {
        (depth, depth, "exact_floor_from_one_2_point")
    } else {
        let gain = if full4 {
            2
        } else if r2 == 3 {
            1
        } else {
            assert_eq!(depth, 0);
            0
        };
        let hi = depth.checked_sub(gain).unwrap();
        (
            0,
            hi,
            if hi == 0 {
                "exact_from_torsion_bound"
            } else if full4 {
                "bounded_by_full_4_torsion"
            } else {
                "bounded_by_full_2_torsion"
            },
        )
    };
    json!({"value":if lo==hi{json!(lo)}else{Value::Null},"bounds":[lo,hi],"status":status})
}
fn descent(f: &Field, a: &Element, b: &Element, j: &Element) -> (usize, &'static str, Value) {
    let (m, kind, parameter) = if *j == f.zero() {
        (6, "sextic", b.clone())
    } else if *j == f.constant(1728) {
        (4, "quartic", a.clone())
    } else {
        (
            2,
            "quadratic",
            f.mul(&f.mul(&f.constant(3), b), &f.inv(&f.mul(&f.constant(2), a))),
        )
    };
    let mut certificates = vec![];
    let mut minimum = None;
    for d in divisors(f.n) {
        let mut fj = j.clone();
        for _ in 0..d {
            fj = f.frob(&fj);
        }
        if fj != *j {
            continue;
        }
        let index = (&f.q - 1u32) / (BigUint::from(f.p).pow(d as u32) - 1u32);
        let h = gcd(m, (&index % m).to_u64().unwrap());
        let result = f.pow(&parameter, &((&f.q - 1u32) / h));
        let accepted = result == f.one();
        if accepted && minimum.is_none() {
            minimum = Some(d);
        }
        certificates.push(json!({"candidate_degree":d,"unit_group_index":index.to_string(),"image_divisor":h,"power_value":result,"accepted":accepted}));
    }
    (
        minimum.unwrap(),
        kind,
        json!({"twist_parameter":parameter,"candidates":certificates}),
    )
}
fn certify_roots(
    f: &Field,
    a: &Element,
    b: &Element,
    n: &BigUint,
    roots: &[Element],
) -> Result<(usize, bool), String> {
    for (i, r) in roots.iter().enumerate() {
        if r.len() != f.n || r.iter().any(|c| *c >= f.p) || roots[..i].contains(r) {
            return Err("malformed root witness".into());
        }
        let value = f.add(&f.add(&f.mul(&f.mul(r, r), r), &f.mul(a, r)), b);
        if value != f.zero() {
            return Err("false rational 2-root".into());
        }
    }
    let denominator = f.add(
        &f.mul(&f.constant(4), &f.mul(&f.mul(a, a), a)),
        &f.mul(&f.constant(27), &f.mul(b, b)),
    );
    if denominator == f.zero() {
        return Err("singular curve".into());
    }
    let square = f.square(&f.neg(&denominator));
    let even = (n % 2u32).is_zero();
    let expected = if !even {
        if !square {
            return Err("order/discriminant parity mismatch".into());
        }
        0
    } else if square {
        3
    } else {
        1
    };
    if roots.len() != expected {
        return Err("incomplete rational 2-root witness".into());
    }
    let mut full4 = expected == 3 && f.square(&f.neg(&f.one()));
    if full4 {
        for i in 0..3 {
            for j in i + 1..3 {
                full4 &= f.square(&f.sub(&roots[i], &roots[j]));
            }
        }
    }
    if expected == 3 && !(n % 4u32).is_zero() {
        return Err("full 2-torsion contradicts order".into());
    }
    if full4 && !(n % 16u32).is_zero() {
        return Err("full 4-torsion contradicts order".into());
    }
    Ok((expected, full4))
}
fn audit_map(repo: &Path) -> BTreeMap<usize, BTreeMap<String, String>> {
    let at = repo.join(STUDY).join("audit_run2");
    let receipt = read_json(&at.join("receipt.json"));
    assert_eq!(receipt["status"], "PROCESS_COMPLETE");
    assert_eq!(
        receipt["input_sha256"],
        hash(&fs::read(repo.join(STUDY).join("inputs.gp")).unwrap())
    );
    assert_eq!(
        receipt["source_sha256"],
        hash(&fs::read(repo.join(STUDY).join("replay.gp")).unwrap())
    );
    assert_eq!(
        receipt["source_sha256"],
        hash(&fs::read(at.join("source/replay.gp")).unwrap())
    );
    assert_eq!(
        receipt["wrapper_source_sha256"],
        hash(&fs::read(at.join("source/main.rs")).unwrap())
    );
    assert!(fs::read(at.join("stderr.txt")).unwrap().is_empty());
    let s = fs::read_to_string(at.join("stdout.txt")).unwrap();
    assert_eq!(
        s.lines().last().unwrap(),
        "PROPERTY_REPLAY_COMPLETE|models=814"
    );
    let mut rows = BTreeMap::new();
    for line in s.lines().filter(|l| l.starts_with("PROPERTY|")) {
        let kv: BTreeMap<String, String> = line
            .split('|')
            .skip(1)
            .map(|s| {
                let (k, v) = s.split_once('=').unwrap();
                (k.into(), v.into())
            })
            .collect();
        let i = kv["index"].parse().unwrap();
        assert!(rows.insert(i, kv).is_none());
    }
    assert_eq!(rows.len(), 814);
    rows
}
fn build(repo: &Path) -> Value {
    let input = load(repo);
    let audits = audit_map(repo);
    let mut fields = BTreeMap::<String, Field>::new();
    let mut rows = vec![];
    let mut counts = BTreeMap::<String, BTreeMap<String, u64>>::new();
    for (i, c) in input["models"].as_array().unwrap().iter().enumerate() {
        let model: Value = serde_json::from_str(c["model_json"].as_str().unwrap()).unwrap();
        assert_eq!(
            hash(c["model_json"].as_str().unwrap().as_bytes()),
            c["model_sha256"].as_str().unwrap()
        );
        let key = model["field"].as_str().unwrap().to_string();
        let p: u64 = model["p"].as_str().unwrap().parse().unwrap();
        let f = fields
            .entry(key)
            .or_insert_with(|| Field::new(p, numbers(&model["modulus"])));
        assert_eq!(f.p, p);
        assert_eq!(f.modulus, numbers(&model["modulus"]));
        let a = encoded(&model, f, "a");
        let b = encoded(&model, f, "b");
        let j = encoded(c, f, "j_coefficients");
        let numerator = f.mul(&f.constant(6912), &f.mul(&f.mul(&a, &a), &a));
        let denominator = f.add(
            &f.mul(&f.constant(4), &f.mul(&f.mul(&a, &a), &a)),
            &f.mul(&f.constant(27), &f.mul(&b, &b)),
        );
        assert_ne!(denominator, f.zero());
        assert_eq!(
            f.mul(&numerator, &f.inv(&denominator)),
            j,
            "recorded j failed"
        );
        let n: BigUint = c["order"].as_str().unwrap().parse().unwrap();
        let t: BigInt = c["trace"].as_str().unwrap().parse().unwrap();
        assert_eq!(
            BigInt::from(n.clone()),
            BigInt::from(f.q.clone()) + 1u32 - &t
        );
        assert!(&t * &t <= BigInt::from(f.q.clone() * 4u32));
        let ordinary = !(&t % BigInt::from(p)).is_zero();
        let delta = &t * &t - BigInt::from(f.q.clone() * 4u32);
        let mut cm = if ordinary {
            let (depth, dk2, split) = local_two(&delta);
            json!({"frobenius_discriminant":delta.to_string(),"frobenius_conductor_v2":depth,"fundamental_discriminant_v2":dk2,"two_splitting":split,"endomorphism_conductor":null,"endomorphism_conductor_v2":null,"endomorphism_conductor_v2_bounds":[0,depth],"full_frobenius_conductor":null,"full_fundamental_discriminant":null})
        } else {
            json!({"frobenius_discriminant":delta.to_string(),"frobenius_conductor_v2":null,"fundamental_discriminant_v2":null,"two_splitting":"not_applicable","endomorphism_conductor":null,"endomorphism_conductor_v2":null,"endomorphism_conductor_v2_bounds":null,"full_frobenius_conductor":null,"full_fundamental_discriminant":null})
        };
        let jd = f.degree(&[j.clone()]);
        let cd = f.degree(&[a.clone(), b.clone()]);
        let (iso, kind, descent) = descent(f, &a, &b, &j);
        let audit = &audits[&i];
        let roots: Vec<Element> = serde_json::from_str(&audit["roots"]).unwrap();
        let (r2, full4) = certify_roots(f, &a, &b, &n, &roots).unwrap();
        if ordinary {
            let level = endomorphism_two(cm["frobenius_conductor_v2"].as_u64().unwrap(), r2, full4);
            cm["endomorphism_conductor_v2"] = level["value"].clone();
            cm["endomorphism_conductor_v2_bounds"] = level["bounds"].clone();
            cm["endomorphism_conductor_v2_status"] = level["status"].clone();
            cm["ordinary_two_volcano_level"] = level["value"].clone();
        } else {
            cm["endomorphism_conductor_v2_status"] = json!("not_applicable");
            cm["ordinary_two_volcano_level"] = Value::Null;
        }

        assert_eq!(jd.to_string(), audit["j_degree"]);
        assert_eq!(cd.to_string(), audit["coefficient_degree"]);
        assert_eq!(iso.to_string(), audit["isomorphic_degree"]);
        assert_eq!(kind, audit["twist_kind"]);
        assert_eq!(full4, audit["full4"] == "1");
        if ordinary {
            assert_eq!(
                cm["frobenius_conductor_v2"].as_u64().unwrap().to_string(),
                audit["f2"]
            );
            assert_eq!(
                cm["fundamental_discriminant_v2"]
                    .as_u64()
                    .unwrap()
                    .to_string(),
                audit["dk2"]
            );
            let sp = match cm["two_splitting"].as_str().unwrap() {
                "split" => "1",
                "inert" => "-1",
                _ => "0",
            };
            assert_eq!(sp, audit["split2"]);
        } else {
            assert_eq!(audit["f2"], "-1");
        }
        if c["labels"]["model_g3_cubic"] == "yes" {
            assert_eq!(r2, 3);
            assert!(cm["frobenius_conductor_v2"].as_u64().unwrap() >= 2);
        }
        if c["labels"]["model_g3_quadratic"] == "yes" {
            assert_eq!(r2, 1);
            assert_eq!(jd, 6);
        }
        if c["inventory"] == "independent_model_law_large_source" {
            assert!(ordinary);
            assert_eq!(
                cm["frobenius_conductor_v2"].as_u64().unwrap() >= 2,
                c["historical_cubic_admission"].as_bool().unwrap()
            );
        }
        let mut small = BTreeMap::new();
        let mut remainder = n.clone();
        for prime in [2u32, 3, 5, 7, 11, 13] {
            let v = valuation(&n, prime);
            small.insert(prime.to_string(), v);
            for _ in 0..v {
                remainder /= prime;
            }
        }
        let smooth = remainder.is_one();
        let cohort = if c["inventory"] == "independent_model_law_large_source" {
            "independent_large"
        } else if p == 7 && f.n == 6 {
            "p7_replay_control"
        } else {
            "higher_degree_control"
        };
        let property = json!({"status":"native_and_pari_verified","ordinary":ordinary,"supersingular":!ordinary,"cm":cm,"rational_2_nonidentity_points":r2,"full_rational_4_torsion":full4,"ambient_order_small_prime_valuations":small,"ambient_order_13_smooth":smooth,"ambient_order_unfactored_remainder":remainder.to_string(),"j_minimum_degree":jd,"coefficient_minimum_degree":cd,"isomorphic_model_minimum_degree":iso,"twist_kind":kind,"j_in_proper_subfield":jd<f.n,"displayed_model_in_proper_subfield":cd<f.n,"isomorphic_model_in_proper_subfield":iso<f.n,"selected_subgroup_order":null,"selected_point_subfield_containment":"unspecified","subgroup_advantage":"unmeasured"});
        let checks = json!({
            "WCF1:supersingular":{"scope":"curve","status":if ordinary{"no"}else{"yes"}},
            "WCF1:anomalous-prime":{"scope":"prime_field_model","status":if f.n!=1{"outside_domain"}else if n==BigUint::from(p){"yes"}else{"no"}},
            "WCF1:smooth-selected-order":{"scope":"ambient_order_bound_13","ambient_order_13_smooth":smooth,"selected_subgroup_order":"unspecified"},
            "WCF1:proper-subfield-points":{"scope":"model_availability","isomorphic_model_available":iso<f.n,"selected_point_containment":"unspecified"},
            "WCF1:g3-cubic":{"scope":"ordinary_class_necessary_condition","status":if f.n!=6||!ordinary{"outside_recorded_domain"}else if property["cm"]["frobenius_conductor_v2"].as_u64().unwrap()>=2{"necessary_condition_passed"}else{"excluded_by_depth_condition"}},
            "WCF1:g3-quadratic":{"scope":"ordinary_model_necessary_conditions","status":if f.n!=6||!ordinary{"outside_recorded_domain"}else if r2==1&&jd==6{"necessary_conditions_passed"}else{"excluded_by_model_conditions"}}
        });
        let count = counts.entry(cohort.into()).or_default();
        *count.entry("models".into()).or_default() += 1;
        for (k, v) in [
            ("ordinary", ordinary),
            ("supersingular", !ordinary),
            ("full_2_torsion", r2 == 3),
            ("full_4_torsion", full4),
            ("j_proper_subfield", jd < f.n),
            ("literal_proper_subfield", cd < f.n),
            ("isomorphic_proper_subfield", iso < f.n),
            ("ambient_13_smooth", smooth),
        ] {
            if v {
                *count.entry(k.into()).or_default() += 1;
            }
        }
        if ordinary {
            *count
                .entry(format!(
                    "f2_depth:{}",
                    property["cm"]["frobenius_conductor_v2"].as_u64().unwrap()
                ))
                .or_default() += 1;
            *count
                .entry(format!(
                    "two_splitting:{}",
                    property["cm"]["two_splitting"].as_str().unwrap()
                ))
                .or_default() += 1;
        }
        if ordinary {
            if let Some(level) = property["cm"]["endomorphism_conductor_v2"].as_u64() {
                *count.entry(format!("endomorphism_v2:{level}")).or_default() += 1;
                if level == 0 {
                    *count.entry("certified_2_surface".into()).or_default() += 1;
                }
                if level == property["cm"]["frobenius_conductor_v2"].as_u64().unwrap() {
                    *count.entry("certified_2_floor".into()).or_default() += 1;
                }
            } else {
                *count.entry("endomorphism_v2_unknown".into()).or_default() += 1;
            }
        }
        rows.push(json!({"icv1":c["icv1"],"slug":c["slug"],"model_sha256":c["model_sha256"],"inventory":c["inventory"],"cohort":cohort,"field":c["field"],"trace":c["trace"],"order":c["order"],"family_labels":c["labels"],"family_property_checks":checks,"properties":property,"witnesses":{"rational_2_roots":roots,"isomorphic_descent":descent,"pari_index":i}}));
    }
    json!({"schema":"weak-curve-structural-properties/v1","model_input":{"path":INPUT,"sha256":INPUT_SHA},"pari_audit":{"stdout_sha256":hash(&fs::read(repo.join(STUDY).join("audit_run2/stdout.txt")).unwrap()),"replay_source_sha256":hash(&fs::read(repo.join(STUDY).join("replay.gp")).unwrap())},"scope":"fixed synthetic-model inventory; class labels preserved; no subgroup-cost inference","summary":{"models":rows.len(),"fields":fields.len(),"cohorts":counts},"records":rows})
}
fn emit(repo: &Path) {
    let input = load(repo);
    let mut text = String::new();
    for (i, c) in input["models"].as_array().unwrap().iter().enumerate() {
        let m: Value = serde_json::from_str(c["model_json"].as_str().unwrap()).unwrap();
        text.push_str(&format!(
            "audit({i},{},{},{},{},{},{},{},{});\n",
            m["p"].as_str().unwrap(),
            m["k"].as_str().unwrap(),
            serde_json::to_string(&numbers(&m["modulus"])).unwrap(),
            serde_json::to_string(&numbers(&m["a"])).unwrap(),
            serde_json::to_string(&numbers(&m["b"])).unwrap(),
            serde_json::to_string(&numbers(&c["j_coefficients"])).unwrap(),
            c["trace"].as_str().unwrap(),
            c["order"].as_str().unwrap()
        ));
    }
    fs::write(repo.join(STUDY).join("inputs.gp"), text).unwrap();
}
fn audit(repo: &Path) {
    let base = repo.join(STUDY);
    let run = base.join("audit_run2");
    assert!(!run.exists(), "audit run exists; preserve it");
    fs::create_dir(&run).unwrap();
    fs::create_dir(run.join("source")).unwrap();
    for name in [
        "main.rs",
        "field.rs",
        "identities.rs",
        "presentation.rs",
        "replay.gp",
        "inputs.gp",
    ] {
        fs::copy(base.join(name), run.join("source").join(name)).unwrap();
    }
    let gp = ["/opt/homebrew/bin/gp", "/usr/bin/gp"]
        .into_iter()
        .find(|p| Path::new(p).is_file())
        .unwrap();
    let mut child = Command::new(gp)
        .args(["-q", "-f", "-s", "128M"])
        .arg(base.join("replay.gp"))
        .env("WEAK_PROPERTIES_INPUTS", base.join("inputs.gp"))
        .stdout(Stdio::from(
            fs::File::create(run.join("stdout.txt")).unwrap(),
        ))
        .stderr(Stdio::from(
            fs::File::create(run.join("stderr.txt")).unwrap(),
        ))
        .spawn()
        .unwrap();
    let start = Instant::now();
    let mut capped = false;
    let status = loop {
        if let Some(s) = child.try_wait().unwrap() {
            break s;
        }
        if start.elapsed() > Duration::from_secs(300) {
            capped = true;
            child.kill().unwrap();
            break child.wait().unwrap();
        }
        std::thread::sleep(Duration::from_millis(100));
    };
    let receipt = json!({"schema":"weak-properties-pari-run/v1","gp_path":gp,"status":if capped{"TIME_CAP"}else if status.success()&&fs::read(run.join("stderr.txt")).unwrap().is_empty(){"PROCESS_COMPLETE"}else{"PROCESS_ERROR"},"exit_code":status.code(),"cap_seconds":300,"wall_ms":start.elapsed().as_millis().to_string(),"accounting":"correctness diagnostic including process launch; shared host","input_sha256":hash(&fs::read(base.join("inputs.gp")).unwrap()),"source_sha256":hash(&fs::read(base.join("replay.gp")).unwrap()),"wrapper_source_sha256":hash(&fs::read(base.join("main.rs")).unwrap()),"wrapper_binary_sha256":hash(&fs::read(std::env::current_exe().unwrap()).unwrap())});
    fs::write(
        run.join("receipt.json"),
        serde_json::to_vec_pretty(&receipt).unwrap(),
    )
    .unwrap();
    assert_eq!(receipt["status"], "PROCESS_COMPLETE");
}
fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(
        args.len(),
        2,
        "emit|audit|identities|build|check|render REPO"
    );
    let repo = Path::new(&args[1]);
    match args[0].as_str() {
        "emit" => emit(repo),
        "audit" => audit(repo),
        "identities" => identities::issue_all(&repo.join(STUDY)),
        "build" => {
            let data = build(repo);
            identities::check_all(&repo.join(STUDY));
            fs::write(repo.join(OUTPUT), serde_json::to_vec_pretty(&data).unwrap()).unwrap();
            fs::write(
                repo.join(STUDY).join("summary.json"),
                serde_json::to_vec_pretty(&data["summary"]).unwrap(),
            )
            .unwrap();
        }
        "check" => {
            let data = build(repo);
            identities::check_all(&repo.join(STUDY));
            assert_eq!(data, read_json(&repo.join(OUTPUT)));
            assert_eq!(
                data["summary"],
                read_json(&repo.join(STUDY).join("summary.json"))
            );
            presentation::check(repo, &data);
        }
        "render" => {
            let data = build(repo);
            presentation::render(repo, &data);
        }
        _ => panic!("unknown operation"),
    }
    println!("WEAK_PROPERTIES_COMPLETE|operation={}|models=814", args[0]);
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn torsion_certifies_local_endomorphism_level() {
        assert_eq!(endomorphism_two(2, 1, false)["value"], 2);
        assert_eq!(endomorphism_two(1, 3, false)["value"], 0);
        assert_eq!(endomorphism_two(2, 3, true)["value"], 0);
        assert_eq!(endomorphism_two(4, 3, true)["bounds"], json!([0, 2]));
        assert_eq!(endomorphism_two(4, 3, false)["value"], Value::Null);
    }
    #[test]
    fn discriminant_two_cases() {
        for (d, f, s) in [
            (-3, 0, "inert"),
            (-4, 0, "ramified"),
            (-7, 0, "split"),
            (-8, 0, "ramified"),
            (-12, 1, "inert"),
            (-28, 1, "split"),
            (-32, 1, "ramified"),
            (-48, 2, "inert"),
        ] {
            let (a, _, b) = local_two(&BigInt::from(d));
            assert_eq!((a, b), (f, s));
        }
    }
    #[test]
    fn descent_distinguishes_j_and_model() {
        let f = Field::new(5, vec![2, 0]);
        let z = f.generator();
        let zero = f.zero();
        assert_eq!(descent(&f, &zero, &z, &zero).0, 2);
        assert_eq!(f.degree(&[zero.clone()]), 1);
        assert_eq!(descent(&f, &zero, &f.one(), &zero).0, 1);
        let f = Field::new(7, vec![1, 0]);
        let z = f.generator();
        assert_eq!(descent(&f, &z, &f.zero(), &f.constant(1728)).0, 1);
        assert_eq!(
            descent(&f, &f.add(&z, &f.one()), &f.zero(), &f.constant(1728)).0,
            2
        );
    }
    #[test]
    fn root_mutation_is_rejected() {
        let f = Field::new(7, vec![1, 0]);
        let a = f.constant(6);
        let b = f.zero();
        let roots = vec![f.zero(), f.one(), f.constant(6)];
        let order = BigUint::from(64u32);
        assert!(certify_roots(&f, &a, &b, &order, &roots).is_ok());
        let mut bad = roots;
        bad[1] = f.constant(2);
        assert!(certify_roots(&f, &a, &b, &order, &bad).is_err());
    }
}
