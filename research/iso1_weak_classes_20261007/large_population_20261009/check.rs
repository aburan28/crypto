//! Verify raw receipts, construct canonical curve identities, and report censored labels.
use crypto_lib::hash::sha256::sha256;
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet},
    env, fs,
    path::Path,
};
fn sha(bytes: &[u8]) -> String {
    hex::encode(sha256(bytes))
}
fn digits(s: &str) -> Vec<u8> {
    s.bytes()
        .map(|b| {
            assert!(b.is_ascii_digit());
            b - b'0'
        })
        .collect()
}
fn add(a: &str, b: &str) -> String {
    let mut a = digits(a);
    let mut b = digits(b);
    a.reverse();
    b.reverse();
    let mut carry = 0u8;
    let mut r = Vec::new();
    for i in 0..a.len().max(b.len()) {
        let x = a.get(i).copied().unwrap_or(0) + b.get(i).copied().unwrap_or(0) + carry;
        r.push(x % 10 + b'0');
        carry = x / 10;
    }
    if carry > 0 {
        r.push(carry + b'0');
    }
    r.reverse();
    String::from_utf8(r).unwrap()
}
fn subtract(a: &str, b: &str) -> String {
    assert!(a.len() > b.len() || (a.len() == b.len() && a >= b));
    let mut a = digits(a);
    let mut b = digits(b);
    a.reverse();
    b.reverse();
    let mut borrow = 0i16;
    let mut r = Vec::new();
    for (i, &v) in a.iter().enumerate() {
        let mut x = i16::from(v) - i16::from(b.get(i).copied().unwrap_or(0)) - borrow;
        borrow = 0;
        if x < 0 {
            x += 10;
            borrow = 1;
        }
        r.push(x as u8 + b'0');
    }
    assert_eq!(borrow, 0);
    while r.len() > 1 && r.last() == Some(&b'0') {
        r.pop();
    }
    r.reverse();
    String::from_utf8(r).unwrap()
}
fn neg(s: &str) -> String {
    if s == "0" {
        s.into()
    } else if let Some(t) = s.strip_prefix('-') {
        t.into()
    } else {
        format!("-{s}")
    }
}
fn rem(s: &str, m: u64) -> u64 {
    let mut r = 0;
    for b in s.trim_start_matches('-').bytes() {
        r = (r * 10 + u64::from(b - b'0')) % m;
    }
    if s.starts_with('-') && r > 0 {
        m - r
    } else {
        r
    }
}
fn bitlen(s: &str) -> usize {
    let mut a = digits(s);
    let mut n = 0;
    while a.iter().any(|&v| v != 0) {
        let mut carry = 0u8;
        for d in &mut a {
            let x = carry * 10 + *d;
            *d = x / 2;
            carry = x % 2;
        }
        n += 1;
    }
    n
}
fn wilson(k: usize, n: usize) -> (f64, f64) {
    if n == 0 {
        return (0.0, 1.0);
    }
    let z = 1.959963984540054;
    let n = n as f64;
    let p = k as f64 / n;
    let den = 1.0 + z * z / n;
    let mid = (p + z * z / (2.0 * n)) / den;
    let h = z * (p * (1.0 - p) / n + z * z / (4.0 * n * n)).sqrt() / den;
    (mid - h, mid + h)
}
#[derive(Default)]
struct Stats {
    p: u64,
    n: u64,
    bits: usize,
    attempts: usize,
    complete: usize,
    ordinary: usize,
    admitted: usize,
    direct: usize,
    witness: usize,
    full4: usize,
    already_full4: usize,
    tested: usize,
    expanded: usize,
    edges2: usize,
    edges3: usize,
    statuses: BTreeMap<String, usize>,
    depths: BTreeMap<u64, usize>,
    small_weak: usize,
    small_zero: usize,
    small_missed: usize,
}
fn canon(
    p: u64,
    k: u64,
    modulus: &Value,
    a: &Value,
    b: &Value,
    j: &Value,
    t: &str,
    order: &str,
) -> Value {
    let coefficients = modulus
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_str().unwrap())
        .collect::<Vec<_>>()
        .join(",");
    let field = format!(
        "fpk-{p}-{k}-{}",
        &sha(format!("fpk-modulus:{p}:{coefficients}").as_bytes())[..8]
    );
    let model = json!({"a":a,"b":b,"field":field,"form":"y^2=x^3+a*x+b","k":k.to_string(),"modulus":modulus,"p":p.to_string(),"v":"1"});
    let model_json = serde_json::to_string(&model).unwrap();
    let h = sha(model_json.as_bytes());
    let jj = j
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_str().unwrap())
        .collect::<Vec<_>>()
        .join(",");
    let signed = if let Some(u) = t.strip_prefix('-') {
        format!("tm{u}")
    } else {
        format!("t{t}")
    };
    json!({"slug":format!("icv1-fp{}k{k}-{signed}-{}",64-p.leading_zeros(),&h[..8]),"icv1":format!("ICV1:{field}:{t}:{order}:{jj}:unk:unk:r:{}",&h[..12]),"field_id":field,"model_json":model_json,"model_sha256":h,"trace":t,"group_order":order,"j_coefficients":j,"subgroup_order":null,"endomorphism_order_conductor":null,"volcano_level":null})
}
fn run(study: &Path, dir: &Path, write: bool) {
    let receipt = fs::read_to_string(dir.join("receipts.tsv")).unwrap();
    let freeze = fs::read_to_string(dir.join("source_freeze.txt")).unwrap();
    for (file, key) in [
        ("sample.gp", "gp_script_sha256"),
        ("runner.rs", "native_driver_sha256"),
        ("PROTOCOL.md", "protocol_sha256"),
    ] {
        let expected = freeze
            .lines()
            .find_map(|l| l.strip_prefix(&format!("{key}=")))
            .unwrap();
        assert_eq!(
            sha(&fs::read(study.join(file)).unwrap()),
            expected,
            "source freeze {file}"
        );
    }
    let exact_csv =
        fs::read_to_string(study.parent().unwrap().join("p7_exact_absolute.csv")).unwrap();
    let exact: BTreeMap<String, usize> = exact_csv
        .lines()
        .skip(1)
        .map(|l| {
            let v = l.split(',').collect::<Vec<_>>();
            (v[2].into(), v[3].parse().unwrap())
        })
        .collect();
    let mut stats: BTreeMap<(usize, String), Stats> = BTreeMap::new();
    let mut records = Vec::new();
    let mut names = BTreeSet::new();
    let mut ids = BTreeSet::new();
    for row in receipt.lines().skip(1) {
        let c = row.split('\t').collect::<Vec<_>>();
        assert_eq!(c.len(), 12);
        assert!(names.insert(c[11]));
        let f = c[0].parse::<usize>().unwrap();
        let p = c[1].parse::<u64>().unwrap();
        let n = c[2].parse::<u64>().unwrap();
        let mode = c[5];
        let s = stats.entry((f, mode.into())).or_default();
        s.p = p;
        s.n = n;
        s.attempts += 1;
        let stdout = dir.join(c[11]);
        let stderr = dir.join(c[11].replace(".stdout", ".stderr"));
        let raw = fs::read_to_string(&stdout).unwrap();
        assert_eq!(sha(raw.as_bytes()), c[9]);
        assert_eq!(sha(&fs::read(&stderr).unwrap()), c[10]);
        if c[6] != "COMPLETE" {
            *s.statuses.entry(c[6].into()).or_default() += 1;
            continue;
        }
        assert_eq!(c[7], "0");
        assert_eq!(fs::metadata(&stderr).unwrap().len(), 0);
        assert!(raw
            .lines()
            .any(|l| l == format!("COMPLETE|{p}|{n}|{}|{mode}", c[4])));
        let lines = raw
            .lines()
            .map(|l| l.split('|').collect::<Vec<_>>())
            .collect::<Vec<_>>();
        let field = lines.iter().find(|v| v[0] == "FIELD").unwrap();
        assert_eq!(field[1], c[1]);
        assert_eq!(field[2], c[2]);
        let q = field[3];
        s.bits = bitlen(q);
        let modulus: Value = serde_json::from_str(field[4]).unwrap();
        assert_eq!(modulus.as_array().unwrap().len(), (2 * n) as usize);
        let sample = lines.iter().find(|v| v[0] == "SAMPLE").unwrap();
        assert_eq!(sample.len(), 13);
        assert_eq!(sample[1], c[4]);
        assert_eq!(sample[2], mode);
        let t = sample[3];
        let order = sample[4];
        let ordinary = sample[5] == "1";
        let admitted = sample[6] == "1";
        let direct = sample[7] == "1";
        let epsilon = sample[10];
        assert_eq!(ordinary, rem(t, p) != 0);
        assert_eq!(rem(order, 4), 0);
        assert_eq!(
            order,
            if let Some(u) = t.strip_prefix('-') {
                add(&add(q, "1"), u)
            } else {
                subtract(&add(q, "1"), t)
            }
        );
        assert_eq!(
            admitted,
            rem(t, 16) == rem(&add(q, "1"), 16) || rem(t, 16) == rem(&neg(&add(q, "1")), 16)
        );
        if ordinary {
            let depth = sample[8].parse::<u64>().unwrap();
            assert_eq!(admitted, depth >= 2);
            *s.depths.entry(depth).or_default() += 1;
        }
        assert!(!ordinary || !direct || admitted);
        assert_eq!(epsilon, if rem(order, 16) == 0 { "1" } else { "-1" });
        let search = lines.iter().find(|v| v[0] == "SEARCH").unwrap();
        assert_eq!(search.len(), 8);
        let status = search[1];
        assert_eq!(status == "NOT_ADMITTED", !ordinary || !admitted);
        assert!(search[2].parse::<usize>().unwrap() <= 256);
        let witness = status == "WITNESS";
        assert_eq!(
            witness,
            lines.iter().any(|v| v[0] == "CURVE" && v[1] == "witness")
        );
        s.complete += 1;
        s.ordinary += usize::from(ordinary);
        s.admitted += usize::from(ordinary && admitted);
        s.direct += usize::from(direct);
        s.witness += usize::from(witness);
        s.tested += search[2].parse::<usize>().unwrap();
        s.expanded += search[3].parse::<usize>().unwrap();
        s.edges2 += search[4].parse::<usize>().unwrap();
        s.edges3 += search[5].parse::<usize>().unwrap();
        *s.statuses.entry(status.into()).or_default() += 1;
        if ordinary && admitted {
            let torsion = lines.iter().find(|v| v[0] == "TORSION").unwrap();
            assert_eq!(torsion[1], "1");
            s.full4 += 1;
            s.already_full4 += usize::from(torsion[2] == "0");
        }
        if mode == "control" {
            assert!(
                !direct && ordinary && admitted && witness,
                "known-positive control must be found"
            );
        }
        if p == 7 && ordinary {
            let count = exact[t];
            if count > 0 {
                s.small_weak += 1;
                if admitted && !witness {
                    s.small_missed += 1;
                }
            } else {
                s.small_zero += 1;
                assert!(!witness, "false witness in exact zero class");
            }
        }
        let selected_t = if epsilon == "-1" { neg(t) } else { t.into() };
        let selected_order = if epsilon == "-1" {
            subtract(&add(&add(q, q), "2"), order)
        } else {
            order.into()
        };
        let mut curves = Vec::new();
        for v in lines.iter().filter(|v| v[0] == "CURVE") {
            assert_eq!(v.len(), 6);
            let a: Value = serde_json::from_str(v[2]).unwrap();
            let b: Value = serde_json::from_str(v[3]).unwrap();
            let j: Value = serde_json::from_str(v[4]).unwrap();
            let roots: Value = serde_json::from_str(v[5]).unwrap();
            for x in [&a, &b, &j] {
                assert_eq!(x.as_array().unwrap().len(), 2 * n as usize);
                for d in x.as_array().unwrap() {
                    assert!(d.as_str().unwrap().parse::<u64>().unwrap() < p);
                }
            }
            let original = v[1] == "source" || v[1] == "control_known_weak";
            let mut curve = canon(
                p,
                2 * n,
                &modulus,
                &a,
                &b,
                &j,
                if original { t } else { &selected_t },
                if original { order } else { &selected_order },
            );
            curve["role"] = json!(v[1]);
            curve["original_root_coordinates"] = roots;
            curve["original_model"] = json!("y^2=(x-r1)*(x-r2)*(x-r3)");
            curve["short_model_coordinate_change"] = json!("X=x+a2/3; Y=y; a2=-(r1+r2+r3)");
            ids.insert(curve["icv1"].as_str().unwrap().to_string());
            curves.push(curve);
        }
        let routes=lines.iter().filter(|v|v[0]=="ROUTE").map(|v|json!({"parent":v[1],"child":v[2],"degree":v[3],"kernel_x":serde_json::from_str::<Value>(v[4]).unwrap()})).collect::<Vec<_>>();
        records.push(json!({"field_number":f,"characteristic":p.to_string(),"odd_degree":n,"field_degree":2*n,"field_cardinality":q,"bit_length":s.bits,"modulus_low_coefficients":modulus,"sample_number":c[3],"seed":c[4],"mode":mode,"raw_output":c[11],"ordinary":ordinary,"admitted":admitted,"direct_weak":direct,"frobenius_conductor_2_depth":sample[8],"search_status":status,"tested_vertices":search[2],"selected_trace_sign":epsilon,"curves":curves,"route":routes,"class_label":if witness{"VERIFIED_CUBIC_NORM_ONE_WEAK"}else if ordinary&&!admitted{"PROVED_CUBIC_NORM_ONE_EXCLUDED"}else{"CUBIC_NORM_ONE_UNRESOLVED"}}));
    }
    let mut summary = Vec::new();
    let mut csv=String::from("field,mode,p,odd_degree,field_degree,bit_length,attempts,completed,ordinary,admitted,excluded,direct_weak,verified_weak,unresolved,full4_representative,already_full4_selected,tested_vertices,expanded_vertices,degree2_edges,degree3_edges,admission_rate,wilson95_low,wilson95_high,precision_lower,precision_upper\n");
    for ((f, mode), s) in &stats {
        let (lo, hi) = wilson(s.admitted, s.ordinary);
        let rate = if s.ordinary > 0 {
            s.admitted as f64 / s.ordinary as f64
        } else {
            0.0
        };
        let lower = if s.admitted > 0 {
            s.witness as f64 / s.admitted as f64
        } else {
            0.0
        };
        csv.push_str(&format!("{f},{mode},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{rate:.9},{lo:.9},{hi:.9},{lower:.9},1.000000000\n",s.p,s.n,2*s.n,s.bits,s.attempts,s.complete,s.ordinary,s.admitted,s.ordinary-s.admitted,s.direct,s.witness,s.admitted-s.witness,s.full4,s.already_full4,s.tested,s.expanded,s.edges2,s.edges3));
        summary.push(json!({"field":f,"mode":mode,"p":s.p.to_string(),"odd_degree":s.n,"bit_length":s.bits,"attempts":s.attempts,"completed":s.complete,"ordinary":s.ordinary,"admitted":s.admitted,"excluded":s.ordinary-s.admitted,"direct_weak":s.direct,"verified_weak":s.witness,"unresolved":s.admitted-s.witness,"full4_representative":s.full4,"already_full4_selected":s.already_full4,"tested_vertices":s.tested,"expanded_vertices":s.expanded,"degree2_edges":s.edges2,"degree3_edges":s.edges3,"statuses":s.statuses,"depths":s.depths,"admission_wilson95":[lo,hi],"precision_identified_interval":[lower,1.0],"p7_exact_weak_classes":s.small_weak,"p7_exact_zero_classes":s.small_zero,"p7_weak_classes_without_witness":s.small_missed}));
        println!("field={f} {mode} bits={} complete={}/{} admitted={}/{} witnesses={} unresolved={} vertices={} {:?}",s.bits,s.complete,s.attempts,s.admitted,s.ordinary,s.witness,s.admitted-s.witness,s.tested,s.statuses);
    }
    let artifact = json!({"schema":"iso1.population/v1","source_receipts_sha256":sha(receipt.as_bytes()),"source_freeze_sha256":sha(freeze.as_bytes()),"unique_curve_models":ids.len(),"summary":summary,"records":records});
    let serialized = serde_json::to_string_pretty(&artifact).unwrap() + "\n";
    if write {
        fs::write(dir.join("curve_records.json"), serialized).unwrap();
        fs::write(dir.join("summary.csv"), csv).unwrap();
    } else {
        assert_eq!(
            fs::read_to_string(dir.join("curve_records.json")).unwrap(),
            serialized
        );
        assert_eq!(fs::read_to_string(dir.join("summary.csv")).unwrap(), csv);
    }
    println!(
        "PASS: raw hashes, source freeze, class-label gates, {} exact model identities",
        ids.len()
    );
}
fn main() {
    let a = env::args().skip(1).collect::<Vec<_>>();
    assert!(
        a.len() == 2 || a.len() == 3,
        "iso1_population_check STUDY_DIR RUN_DIR [--write]"
    );
    run(
        Path::new(&a[0]),
        Path::new(&a[1]),
        a.get(2).map(String::as_str) == Some("--write"),
    );
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn retained_receipts_and_labels() {
        let s = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("large_population_20261009");
        for r in [
            "validation_p7_run1",
            "pilot_run1",
            "positive_controls_run1",
            "population_run1",
        ] {
            run(&s, &s.join(r), false);
        }
    }
}
