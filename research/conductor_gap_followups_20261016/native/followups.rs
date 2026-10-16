//! Frozen native follow-up preparation, execution and independent replay.
mod analysis;
#[cfg(test)]
mod tests;
mod tracking;
use crypto_lib::cryptanalysis::{
    ecbench::canonical::sha256_hex,
    ic_boundary::{
        rho_walk_with, CountedGroup, GroupOps, IncrementalGauss, NegationClasses, PrimeCurve,
        PrimeInstance, PrimePoint, RhoClasses, RhoWalk, RowStatus,
    },
    ic_framework::cm_orbit::{action, derive_map},
};
use rand::rngs::StdRng;
use rand::{seq::SliceRandom, SeedableRng};
use serde_json::{json, Value};
use std::{
    cell::Cell,
    collections::BTreeSet,
    fs,
    path::{Path, PathBuf},
    process::Command,
    time::Instant,
};
const ROOT: &str = "research/conductor_gap_followups_20261016";
const OLD: &str = "research/conductor_gap_cm_orbits_20261015";
fn read(p: impl AsRef<Path>) -> Value {
    serde_json::from_slice(&fs::read(p).unwrap()).unwrap()
}
fn write(p: impl AsRef<Path>, v: &Value) {
    let p = p.as_ref();
    let temporary = p.with_extension(format!("tmp-{}", std::process::id()));
    fs::write(&temporary, serde_json::to_string_pretty(v).unwrap() + "\n").unwrap();
    fs::File::open(&temporary).unwrap().sync_all().unwrap();
    fs::rename(temporary, p).unwrap();
}
fn u(v: &Value, k: &str) -> u64 {
    v[k].as_u64().unwrap()
}
fn product(a: u64, b: u64, m: u64) -> u64 {
    (a as u128 * b as u128 % m as u128) as u64
}
fn power(mut a: u64, mut e: u64, m: u64) -> u64 {
    let mut b = 1;
    while e > 0 {
        if e & 1 != 0 {
            b = product(b, a, m)
        }
        a = product(a, a, m);
        e >>= 1;
    }
    b
}
fn inst(v: &Value) -> PrimeInstance {
    PrimeInstance {
        name: "public conductor control".into(),
        curve: PrimeCurve {
            p: u(v, "p"),
            a: u(v, "a"),
            b: u(v, "b"),
        },
        group_order: u(v, "N"),
        r: u(v, "N"),
        cofactor: 1,
        generator: (u(v, "Gx"), u(v, "Gy")),
    }
}
fn target(c: &PrimeCurve, r: u64, label: &str) -> PrimePoint {
    for attempt in 0.. {
        let h = sha256_hex(
            format!("conductor-followups/public-target/v1/{label}/{attempt}").as_bytes(),
        );
        let x = u64::from_str_radix(&h[..16], 16).unwrap() % c.p;
        if let Some(y) = c.sqrt(c.rhs(x)) {
            if y != 0 {
                let p = PrimePoint::affine(x, y.min(c.p - y));
                assert_eq!(c.mul(&mut GroupOps::default(), p, r), PrimePoint::INFINITY);
                return p;
            }
        }
    }
    unreachable!()
}
fn seed(label: &str) -> u64 {
    u64::from_str_radix(&sha256_hex(label.as_bytes())[..16], 16).unwrap()
}
fn endpoint(pair: &Value, role: u64) -> Value {
    json!({"p":pair["p"],"N":pair["N"],"a":pair[format!("a{role}")],"b":pair[format!("b{role}")],"Gx":pair[format!("G{role}")][0],"Gy":pair[format!("G{role}")][1],"conductor":if role==0{json!(1)}else{pair["ell"].clone()},"D":-19})
}
fn job(
    cohort: &str,
    pair: u64,
    rep: u64,
    arm: &str,
    v: &Value,
    q: PrimePoint,
    qs: u64,
    partition: u64,
) -> Value {
    json!({"cohort":cohort,"pair":pair,"rep":rep,"arm":arm,"curve":v,"Qx":q.x,"Qy":q.y,"query_seed":qs,"partition_seed":partition,"B":8192,"query_limit":300000000,"wall_cap_seconds":180,"rho_max_steps":1000000000})
}
fn prepare() {
    let root = Path::new(ROOT);
    let pairs = read(Path::new(OLD).join("inputs/pairs.json"));
    let pairs = pairs.as_array().unwrap();
    assert_eq!(pairs.len(), 100);
    let mut rng = StdRng::seed_from_u64(2026101601);
    let mut selected: Vec<usize> = (0..100).collect();
    selected.shuffle(&mut rng);
    selected.truncate(20);
    write(root.join("prepared/selected-gaps.json"), &json!(selected));
    let mut cross = Vec::new();
    let mut yields = Vec::new();
    let mut rho = Vec::new();
    let mut params = String::new();
    let mut orbit = Vec::new();
    for id in 0..100 {
        let pair = &pairs[id];
        for role in 0..2 {
            let v = endpoint(pair, role);
            let i = inst(&v);
            let m = if role == 0 {
                Some(derive_map(&i, 19, &mut GroupOps::default()).unwrap())
            } else {
                None
            };
            let (_, lam, _) = action(&i, 19).unwrap();
            let mut ps = vec![
                i.curve.p,
                i.r,
                2 * id as u64 + role,
                lam,
                1,
                1,
                0,
                0,
                0,
                0,
                0,
                0,
            ];
            if let Some(m) = m {
                ps[4] = m.inv_c;
                ps[5] = m.inv_cu;
                ps[6] = m.kernel_x[0];
                ps[7] = m.kernel_x[1];
                ps[8] = m.kernel_t[0];
                ps[9] = m.kernel_t[1];
                ps[10] = m.kernel_rhs[0];
                ps[11] = m.kernel_rhs[1];
            }
            params += &(ps.iter().map(u64::to_string).collect::<Vec<_>>().join(" ") + "\n");
            for length in [256, 512, 1024, 2048] {
                let mut j = job(
                    "orbit",
                    id as u64,
                    length,
                    if role == 0 {
                        "crater-recurrence"
                    } else {
                        "floor-recurrence"
                    },
                    &v,
                    i.generator_point(),
                    0,
                    0,
                );
                j["B"] = json!(4 * length);
                j["query_limit"] = json!(0);
                orbit.push(j);
            }
            if !selected.contains(&id) {
                continue;
            }
            for rep in 0..10 {
                let q = target(&i.curve, i.r, &format!("yield/{id}/{role}/{rep}"));
                yields.push(job(
                    "yield",
                    id as u64,
                    rep,
                    if role == 0 {
                        "crater-recurrence"
                    } else {
                        "floor-recurrence"
                    },
                    &v,
                    q,
                    seed(&format!("yield-query/{id}/{rep}")),
                    0,
                ));
            }
            for rep in 0..5 {
                let q = target(&i.curve, i.r, &format!("rho/{id}/{role}/{rep}"));
                for part in 0..2 {
                    rho.push(job(
                        "rho",
                        id as u64,
                        rep,
                        if role == 0 { "crater-rho" } else { "floor-rho" },
                        &v,
                        q,
                        seed(&format!("rho-jumps/{id}/{rep}")),
                        seed(&format!("rho-partition/{id}/{part}")),
                    ));
                }
            }
        }
        if !selected.contains(&id) {
            continue;
        }
        for rep in 0..4 {
            let mut block = Vec::new();
            for (role, arm) in [
                (0, "crater-cm"),
                (0, "crater-recurrence"),
                (1, "floor-recurrence"),
            ] {
                let v = endpoint(pair, role);
                let i = inst(&v);
                let q = target(&i.curve, i.r, &format!("crossover/{id}/{role}/{rep}"));
                let mut j = job(
                    "crossover",
                    id as u64,
                    rep,
                    arm,
                    &v,
                    q,
                    seed(&format!("cross-query/{id}/{rep}")),
                    0,
                );
                j["query_limit"] = json!(1000000000);
                j["wall_cap_seconds"] = json!(240);
                block.push(j)
            }
            let perms = [
                [0, 1, 2],
                [0, 2, 1],
                [1, 0, 2],
                [1, 2, 0],
                [2, 0, 1],
                [2, 1, 0],
            ];
            let order = perms[cross.len() % 6];
            cross.push(order.map(|k| block[k].clone()).to_vec());
        }
    }
    cross.shuffle(&mut rng);
    let cross: Vec<_> = cross.into_iter().flatten().collect();
    yields.shuffle(&mut rng);
    rho.shuffle(&mut rng);
    let mut coordinates = Vec::new();
    for old in &rho {
        if u(old, "rep") != 0 || !selected[..4].contains(&(u(old, "pair") as usize)) {
            continue;
        }
        let mut j = old.clone();
        j["cohort"] = json!("rho-coordinate");
        let p = u(&j["curve"], "p");
        for (key, factor) in [("a", 16), ("b", 64), ("Gx", 4), ("Gy", 8)] {
            j["curve"][key] = json!(product(u(&old["curve"], key), factor, p));
        }
        j["Qx"] = json!(product(u(old, "Qx"), 4, p));
        j["Qy"] = json!(product(u(old, "Qy"), 8, p));
        let i = inst(&j["curve"]);
        assert!(i.curve.is_on_curve(i.generator_point()));
        assert!(i
            .curve
            .is_on_curve(PrimePoint::affine(u(&j, "Qx"), u(&j, "Qy"))));
        coordinates.push(j);
    }
    assert_eq!(coordinates.len(), 16);
    for (name, jobs) in [
        ("crossover", cross),
        ("yield", yields),
        ("rho", rho),
        ("orbit", orbit),
        ("rho-coordinate", coordinates),
    ] {
        write(root.join(format!("prepared/{name}.json")), &json!(jobs));
    }
    fs::write(root.join("prepared/cm-params.txt"), params).unwrap();
    if root.join("results/dose-construction-v3.jsonl").exists() {
        let text = fs::read_to_string(root.join("results/dose-construction-v3.jsonl")).unwrap();
        let mut jobs = Vec::new();
        let mut par = String::new();
        let models: Vec<Value> = text
            .lines()
            .map(|s| serde_json::from_str(s).unwrap())
            .collect();
        assert_eq!(models.len(), 32);
        for v in &models {
            let i = inst(v);
            let id = 1000 + u(v, "class") * 10 + u(v, "stratum");
            let (_, lam, qdegree) = action(&i, (-v["D"].as_i64().unwrap()) as u64).unwrap();
            par += &format!("{} {} {id} {lam} 1 1 0 0 0 0 0 0\n", i.curve.p, i.r);
            for rep in 0..5 {
                let q = target(&i.curve, i.r, &format!("dose/{id}/{rep}"));
                for arm in ["dose-recurrence", "dose-rho"] {
                    let mut j = job(
                        "dose",
                        id,
                        rep,
                        arm,
                        v,
                        q,
                        seed(&format!("dose-query/{}/{rep}", v["class"])),
                        0,
                    );
                    j["degree"] = json!(qdegree);
                    j["B"] = json!(8192);
                    j["query_limit"] = json!(1000000000);
                    j["wall_cap_seconds"] = json!(240);
                    jobs.push(j)
                }
            }
        }
        jobs.shuffle(&mut rng);
        write(root.join("prepared/dose.json"), &json!(jobs));
        fs::write(root.join("prepared/dose-params.txt"), par).unwrap();
    }
}
fn canonical_base(
    c: &PrimeCurve,
    g: PrimePoint,
    r: u64,
    lambda: u64,
    b: usize,
) -> Vec<(PrimePoint, usize, u64)> {
    let mut base = Vec::new();
    let mut x = 1;
    let mut ops = GroupOps::default();
    for col in 0..4 {
        let s = loop {
            let a = x;
            x += 1;
            if a != g.x {
                if let Some(y) = c.sqrt(c.rhs(a)) {
                    if y != 0 {
                        break PrimePoint::affine(a, y.min(c.p - y));
                    }
                }
            }
        };
        let mut p = s;
        let mut k = 1;
        for _ in 0..b / 4 {
            base.push((p, col, k));
            p = c.mul(&mut ops, p, lambda);
            k = product(k, lambda, r);
            if p.y > c.p - p.y {
                p = c.neg(p);
                k = r - k
            }
        }
    }
    base
}
fn base_text(dir: &Path) -> String {
    let local = dir.join("factor-bases/job-0.csv");
    if local.exists() {
        return fs::read_to_string(local).unwrap();
    }
    let reference = read(dir.join("base-reference.json"));
    let hash = reference["sha256"].as_str().unwrap();
    assert_eq!(hash.len(), 64);
    assert!(hash.bytes().all(|b| b.is_ascii_hexdigit()));
    let path = Path::new(ROOT).join(format!("bases/{hash}.csv.gz"));
    let compressed = fs::read(path).unwrap();
    assert_eq!(sha256_hex(&compressed), reference["compressed_sha256"]);
    let mut text = String::new();
    std::io::Read::read_to_string(
        &mut flate2::read::GzDecoder::new(compressed.as_slice()),
        &mut text,
    )
    .unwrap();
    assert_eq!(sha256_hex(text.as_bytes()), hash);
    text
}
fn archive_base(dir: &Path) {
    use std::io::Write;
    let csv = fs::read(dir.join("factor-bases/job-0.csv")).unwrap();
    let hash = sha256_hex(&csv);
    let root = Path::new(ROOT).join("bases");
    fs::create_dir_all(&root).unwrap();
    let p = root.join(format!("{hash}.csv.gz"));
    if !p.exists() {
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        enc.write_all(&csv).unwrap();
        fs::write(&p, enc.finish().unwrap()).unwrap();
    }
    write(
        dir.join("base-reference.json"),
        &json!({"sha256":hash,"compressed_sha256":sha256_hex(&fs::read(p).unwrap()),"scope":"Lossless content-addressed base artifact; identical bytes shared across runs."}),
    );
    fs::remove_file(dir.join("factor-bases/job-0.csv")).unwrap();
}
fn replay(j: &Value, dir: &Path, raw: &Value) -> Value {
    let i = inst(&j["curve"]);
    let c = &i.curve;
    let (_, lambda, _) = action(&i, (-j["curve"]["D"].as_i64().unwrap()) as u64).unwrap();
    let base = canonical_base(c, i.generator_point(), i.r, lambda, u(j, "B") as usize);
    if j["cohort"] == "orbit" && j["curve"]["conductor"] == 1 {
        let map = derive_map(&i, 19, &mut GroupOps::default()).unwrap();
        let mut fields = Default::default();
        for k in 0..base.len() - 1 {
            if base[k].1 == base[k + 1].1 {
                let p = map.apply(c, base[k].0, &mut fields).unwrap();
                let canonical = if p.y > c.p - p.y { c.neg(p) } else { p };
                assert_eq!(canonical, base[k + 1].0);
            }
        }
    }
    let mut csv = String::from("point,column,coefficient,x,y\n");
    let mut xs = BTreeSet::new();
    for (idx, (p, col, k)) in base.iter().enumerate() {
        assert!(c.is_on_curve(*p));
        assert!(xs.insert(p.x));
        csv += &format!("{idx},{col},{k},{},{}\n", p.x, p.y)
    }
    assert_eq!(base_text(dir), csv);
    let q = PrimePoint::affine(u(j, "Qx"), u(j, "Qy"));
    let queries = u(raw, "queries");
    assert!(queries <= u(j, "query_limit"));
    let mut gauss = IncrementalGauss::new(5, i.r);
    let mut profiles = BTreeSet::new();
    let mut ops = GroupOps::default();
    let events = fs::read_to_string(dir.join("events.jsonl")).unwrap();
    let retained: Vec<Value> = fs::read_to_string(dir.join("relations/job-0.jsonl"))
        .unwrap()
        .lines()
        .map(|s| serde_json::from_str(s).unwrap())
        .collect();
    let mut checked_retained = 0;
    let mut hits = 0;
    let mut duplicates = 0;
    let mut last = None;
    for line in events.lines() {
        let e: Value = serde_json::from_str(line).unwrap();
        let offset = u(&e, "query");
        assert!(offset < queries && last.is_none_or(|l| offset >= l));
        last = Some(offset);
        assert_eq!(
            u(&e, "alpha"),
            ((u(j, "query_seed") % i.r) as u128 + offset as u128).rem_euclid(i.r as u128) as u64
        );
        let mut row = vec![0u64; 5];
        let mut sum = PrimePoint::INFINITY;
        for (k, sk) in [("i", "si"), ("j", "sj")] {
            let (p, col, coef) = base[u(&e, k) as usize];
            let sign = e[sk].as_i64().unwrap();
            assert!(sign == 1 || sign == -1);
            sum = c.add(&mut ops, sum, if sign == 1 { p } else { c.neg(p) });
            row[col] = ((row[col] as u128 + if sign == 1 { coef } else { i.r - coef } as u128)
                % i.r as u128) as u64;
        }
        let ag = c.mul(&mut ops, i.generator_point(), u(&e, "alpha"));
        assert_eq!(sum, c.add(&mut ops, ag, q));
        assert_eq!(sum, PrimePoint::affine(u(&e, "Rx"), u(&e, "Ry")));
        let nonzero = *row[..4].iter().find(|&&x| x != 0).unwrap();
        let inv = power(nonzero, i.r - 2, i.r);
        let profile: Vec<_> = row[..4].iter().map(|&x| product(x, inv, i.r)).collect();
        let fresh = profiles.insert(profile);
        hits += 1;
        if !fresh {
            duplicates += 1;
            assert_eq!(e["kind"], "proportional")
        } else {
            assert_eq!(e["kind"], "retained");
            let mut claimed = e.clone();
            claimed.as_object_mut().unwrap().remove("kind");
            claimed["beta"] = json!(1);
            assert_eq!(claimed, retained[checked_retained]);
            checked_retained += 1;
            row[4] = i.r - 1;
            let before = gauss.rank();
            assert!(!matches!(
                gauss.add_row(row, u(&e, "alpha")),
                RowStatus::Inconsistent
            ));
            assert_eq!(e["independent"], gauss.rank() > before);
        }
        assert_eq!(gauss.rank() as u64, u(&e, "rank"));
    }
    assert_eq!(hits, u(raw, "hits"));
    assert_eq!(checked_retained, retained.len());
    assert_eq!(duplicates, u(raw, "rejected_duplicate_columns"));
    assert_eq!(gauss.rank() as u64, u(raw, "rank"));
    if raw["verified"] == true {
        let d = u(raw, "recovered");
        assert_eq!(gauss.pinned(4), Some(d));
        assert_eq!(c.mul(&mut ops, i.generator_point(), d), q);
    } else {
        assert!(gauss.pinned(4).is_none())
    }
    json!({"ok":true,"base_sha256":sha256_hex(csv.as_bytes()),"events_sha256":sha256_hex(events.as_bytes()),"base_points":base.len(),"hits":hits,"proportional":duplicates,"rank":gauss.rank(),"recovery_verified":raw["verified"]})
}
struct Salted<'a> {
    c: &'a PrimeCurve,
    salt: u64,
    keys: Cell<u64>,
    negs: Cell<u64>,
}
impl CountedGroup for Salted<'_> {
    type Elt = PrimePoint;
    fn identity(&self) -> PrimePoint {
        self.c.identity()
    }
    fn is_identity(&self, p: &PrimePoint) -> bool {
        self.c.is_identity(p)
    }
    fn add(&self, o: &mut GroupOps, p: PrimePoint, q: PrimePoint) -> PrimePoint {
        self.c.add(o, p, q)
    }
    fn double(&self, o: &mut GroupOps, p: PrimePoint) -> PrimePoint {
        self.c.double(o, p)
    }
    fn neg(&self, p: PrimePoint) -> PrimePoint {
        self.negs.set(self.negs.get() + 1);
        self.c.neg(p)
    }
    fn key(&self, p: &PrimePoint) -> u64 {
        self.keys.set(self.keys.get() + 1);
        let mut x = self.c.key(p) ^ self.salt;
        x ^= x >> 30;
        x = x.wrapping_mul(0xbf58476d1ce4e5b9);
        x ^= x >> 27;
        x = x.wrapping_mul(0x94d049bb133111eb);
        x ^ (x >> 31)
    }
    fn counts_field_ops(&self) -> bool {
        true
    }
}
struct FixedNeg {
    r: u64,
}
impl RhoClasses<Salted<'_>> for FixedNeg {
    fn automorphisms(&self) -> u32 {
        2
    }
    fn canon(&self, g: &Salted<'_>, o: &mut GroupOps, p: PrimePoint) -> (PrimePoint, u64) {
        NegationClasses { r: self.r }.canon(g.c, o, p)
    }
    fn method(&self) -> &'static str {
        "fixed-canonical-negation/salted-partition"
    }
}
fn rho_run(j: &Value) -> Value {
    let i = inst(&j["curve"]);
    let c = &i.curve;
    let g = Salted {
        c,
        salt: u(j, "partition_seed"),
        keys: Cell::new(0),
        negs: Cell::new(0),
    };
    let q = PrimePoint::affine(u(j, "Qx"), u(j, "Qy"));
    let out = rho_walk_with(
        &g,
        &FixedNeg { r: i.r },
        i.generator_point(),
        q,
        i.r,
        u(j, "query_seed"),
        u(j, "rho_max_steps"),
        RhoWalk::negation(),
    );
    let mut check = GroupOps::default();
    if let Some(d) = out.recovered {
        assert!(out.verified);
        assert_eq!(c.mul(&mut check, i.generator_point(), d), q);
        assert_ne!(c.mul(&mut check, i.generator_point(), (d + 1) % i.r), q);
    }
    let field = out.group_ops.field_ops;
    json!({"result":out,"field_ops":field,"key_calls":g.keys.get(),"explicit_neg_calls":g.negs.get(),"sqrt_r":(i.r as f64).sqrt(),"hash_canonicalization_cost":"included in process CPU; key calls and canonicalizations counted, not assigned field-operation equivalents"})
}
fn run_one(cohort: &str, index: usize) {
    let root = Path::new(ROOT);
    let jobs = read(root.join(format!("prepared/{cohort}.json")));
    let j = &jobs[index];
    let dir = root.join(format!("runs/{cohort}-{index:04}"));
    assert!(!dir.exists(), "never overwrite a run; inspect its receipt");
    fs::create_dir_all(&dir).unwrap();
    write(dir.join("input.json"), j);
    let start = Instant::now();
    let result = if j["arm"].as_str().unwrap().contains("rho") {
        rho_run(j)
    } else {
        for name in ["factor-bases", "relations"] {
            fs::create_dir(dir.join(name)).unwrap()
        }
        let i = inst(&j["curve"]);
        let role = if j["arm"] == "crater-cm" { 0 } else { 4 };
        let pid = if cohort == "dose" {
            u(j, "pair")
        } else {
            2 * u(j, "pair") + u(&j["curve"], "conductor").ne(&1) as u64
        };
        let csv=format!("id,pair,role,p,r,a,b,gx,gy,qx,qy,seed,repeat\n0,{pid},{role},{},{},{},{},{},{},{},{},{},0\n",i.curve.p,i.r,i.curve.a,i.curve.b,i.generator.0,i.generator.1,j["Qx"],j["Qy"],j["query_seed"]);
        fs::write(dir.join("work.csv"), csv).unwrap();
        let q = if j["degree"] == 11 { 11 } else { 5 };
        let binary = format!("/workspace/native-deps/cm-folded-events-q{q}");
        let params = root.join(if cohort == "dose" {
            "prepared/dose-params.txt"
        } else if cohort == "control" {
            "prepared/control-params.txt"
        } else {
            "prepared/cm-params.txt"
        });
        let out = Command::new(binary)
            .arg(dir.join("work.csv"))
            .arg(j["B"].to_string())
            .arg(j["query_limit"].to_string())
            .arg(j["wall_cap_seconds"].to_string())
            .arg("1")
            .arg(&dir)
            .arg(params)
            .output()
            .unwrap();
        fs::write(dir.join("stdout.json"), &out.stdout).unwrap();
        fs::write(dir.join("stderr.log"), &out.stderr).unwrap();
        assert!(out.status.success(), "backend failed {:?}", out.status);
        let raw: Value = serde_json::from_slice(&out.stdout).unwrap();
        let replay = replay(j, &dir, &raw);
        if cohort != "control" {
            archive_base(&dir);
        }
        json!({"native":raw,"replay":replay})
    };
    write(
        dir.join("result.json"),
        &json!({"schema":"conductor-followup.run/v1","input_sha256":sha256_hex(&fs::read(dir.join("input.json")).unwrap()),"result":result,"driver_wall_seconds":start.elapsed().as_secs_f64()}),
    );
}
fn controls() {
    let root = Path::new(ROOT);
    let mut jobs = Vec::new();
    let mut params = String::new();
    for (idx, name, d) in [
        (0, "control-fixture-fix", 19),
        (1, "control-q11-crater", 43),
        (2, "control-q11-floor", 43),
    ] {
        let spec = read(Path::new(OLD).join(format!("sessions/{name}/spec.json")));
        let cs = &spec["workloads"]["curves"][0];
        let v = json!({"p":cs["p"],"N":cs["group_order"],"a":cs["a"],"b":cs["b"],"Gx":cs["gx"],"Gy":cs["gy"],"conductor":1,"D":-(d as i64)});
        let i = inst(&v);
        let (_, lambda, degree) = action(&i, d).unwrap();
        params += &format!(
            "{} {} {} {lambda} 1 1 0 0 0 0 0 0\n",
            i.curve.p,
            i.r,
            2 * idx
        );
        for arm in ["crater-recurrence", "crater-rho"] {
            let q = target(&i.curve, i.r, &format!("small-control/{idx}"));
            let mut j = job("control", idx, 0, arm, &v, q, 12345, 6789);
            j["B"] = json!(64);
            j["query_limit"] = json!(300000);
            j["degree"] = json!(degree);
            jobs.push(j);
        }
    }
    write(root.join("prepared/control.json"), &json!(jobs));
    fs::write(root.join("prepared/control-params.txt"), params).unwrap();
}
fn certificate_inputs() {
    let root = Path::new(ROOT);
    let text = fs::read_to_string(root.join("results/dose-construction-v3.jsonl")).unwrap();
    let mut out = String::new();
    for line in text.lines() {
        let v: Value = serde_json::from_str(line).unwrap();
        let vals = [
            v["class"].clone(),
            v["stratum"].clone(),
            v["p"].clone(),
            v["N"].clone(),
            v["t"].clone(),
            v["fpi"].clone(),
            v["conductor"].clone(),
            v["D"].clone(),
            v["a"].clone(),
            v["b"].clone(),
            v["Gx"].clone(),
            v["Gy"].clone(),
        ];
        out += &(serde_json::to_string(&vals).unwrap() + "\n");
    }
    fs::write(root.join("prepared/dose-cert.gpvec"), out).unwrap();
    let text = fs::read_to_string(Path::new(OLD).join("prepared/cm-params.txt")).unwrap();
    let mut out = String::new();
    for line in text.lines() {
        let v: Vec<u64> = line
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        out += &format!("[{},{},{}]\n", v[2], v[1], v[3]);
    }
    fs::write(root.join("prepared/orbit-order.gpvec"), out).unwrap();
}
fn main() {
    let a: Vec<_> = std::env::args().collect();
    match a[1].as_str() {
        "prepare" => prepare(),
        "freeze" => tracking::freeze(),
        "campaign" => tracking::campaign(),
        "progress" => analysis::progress(),
        "analyze" => analysis::analyze(),
        "audit-all" => analysis::audit(),
        "instrumentation-controls" => tracking::instrumentation(),
        "certificate-inputs" => certificate_inputs(),
        "controls" => controls(),
        "run-one" => run_one(&a[2], a[3].parse().unwrap()),
        "replay" => {
            let dir = PathBuf::from(&a[2]);
            let j = read(dir.join("input.json"));
            let v = read(dir.join("result.json"));
            if let Some(raw) = v["result"].get("native") {
                write(dir.join("audit.json"), &replay(&j, &dir, raw));
            } else {
                let i = inst(&j["curve"]);
                let r = &v["result"]["result"];
                if let Some(d) = r["recovered"].as_u64() {
                    assert_eq!(
                        i.curve
                            .mul(&mut GroupOps::default(), i.generator_point(), d),
                        PrimePoint::affine(u(&j, "Qx"), u(&j, "Qy"))
                    );
                }
                write(
                    dir.join("audit.json"),
                    &json!({"ok":true,"recovery_checked":r["verified"]}),
                );
            }
        }
        _ => panic!("prepare | controls | run-one COHORT INDEX | replay DIR"),
    }
}
