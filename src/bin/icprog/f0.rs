//! Track B's F0 instance runs and their replays, natively: B3's gate rho
//! (its measurement 5), B3b's two-word instances (measurement 5) and B4's
//! three-word instances (measurement 5, `rounds/B4-multi-word/run.py` and
//! `analyse.py`, ported).
//!
//! Each document runs as one `ic price` at F0 through the programme's
//! runner: isolated, one core, `RAYON_NUM_THREADS=1`, under `timeout`.  A
//! contended or failed attempt is kept and the run made again, at most
//! twice; the first clean, complete attempt is the run's figure.  Nothing
//! is overwritten, so each command resumes where the last stopped.
//!
//! The replay is `[d]G = Q` on the document's curve, for every arm's
//! certificate, in this module's own arithmetic over GF(2^n): the
//! generator's (`conformance/v2/make_cases.py`'s `Field` and `Curve`),
//! ported, which shares nothing with the tool.  A run passes when its
//! report is `complete`, every arm's scalar replays, the arms agree, and a
//! known logarithm equals them.

use std::path::{Path, PathBuf};

use num_bigint::BigUint;
use num_traits::{One, Zero};

use super::b5a;
use super::bench::Bench;
use super::json::{self, obj, J};
use super::runs;
use super::suite;

// ── GF(2^n) and the curve y² + xy = x³ + ax² + b ───────────────────

/// GF(2)[z] / (f), elements as integers below 2^n.
pub struct Field {
    n: u64,
    f: BigUint,
}

impl Field {
    pub fn new(n: u64, f: BigUint) -> Result<Field, String> {
        if f.bits() != n + 1 || !f.bit(0) {
            return Err(format!(
                "the modulus {f:#x} is not of degree {n} with a constant term"
            ));
        }
        Ok(Field { n, f })
    }

    pub fn mul(&self, a: &BigUint, b: &BigUint) -> BigUint {
        let (mut r, mut a) = (BigUint::zero(), a.clone());
        for i in 0..b.bits() {
            if b.bit(i) {
                r ^= &a;
            }
            a <<= 1u32;
            if a.bit(self.n) {
                a ^= &self.f;
            }
        }
        r
    }

    pub fn sqr(&self, a: &BigUint) -> BigUint {
        self.mul(a, a)
    }

    /// Algorithm 2.48 of Hankerson, Menezes and Vanstone.
    pub fn inv(&self, a: &BigUint) -> Option<BigUint> {
        if a.is_zero() {
            return None;
        }
        let (mut u, mut v) = (a.clone(), self.f.clone());
        let (mut g1, mut g2) = (BigUint::one(), BigUint::zero());
        while !u.is_one() {
            if u.bits() < v.bits() {
                std::mem::swap(&mut u, &mut v);
                std::mem::swap(&mut g1, &mut g2);
            }
            let j = u.bits() - v.bits();
            u ^= &v << j;
            g1 ^= &g2 << j;
        }
        Some(g1)
    }
}

pub type Point = Option<(BigUint, BigUint)>;

pub struct Curve {
    k: Field,
    a: BigUint,
    b: BigUint,
}

impl Curve {
    pub fn new(k: Field, a: BigUint, b: BigUint) -> Result<Curve, String> {
        if b.is_zero() {
            return Err("b = 0: the curve is singular".into());
        }
        Ok(Curve { k, a, b })
    }

    pub fn on_curve(&self, p: &Point) -> bool {
        let Some((x, y)) = p else { return true };
        let k = &self.k;
        let x2 = k.sqr(x);
        k.sqr(y) ^ k.mul(x, y) == k.mul(&x2, x) ^ k.mul(&self.a, &x2) ^ &self.b
    }

    pub fn add(&self, p: &Point, q: &Point) -> Point {
        let (Some((x1, y1)), Some((x2, y2))) = (p, q) else {
            return if p.is_none() { q.clone() } else { p.clone() };
        };
        let k = &self.k;
        if x1 == x2 {
            if (y1 ^ y2) == *x1 {
                return None;
            }
            return self.dbl(p);
        }
        let lam = k.mul(&(y1 ^ y2), &k.inv(&(x1 ^ x2))?);
        let x3 = k.sqr(&lam) ^ &lam ^ x1 ^ x2 ^ &self.a;
        let y3 = k.mul(&lam, &(x1 ^ &x3)) ^ &x3 ^ y1;
        Some((x3, y3))
    }

    pub fn dbl(&self, p: &Point) -> Point {
        let (x, y) = p.as_ref()?;
        let k = &self.k;
        let lam = x ^ k.mul(y, &k.inv(x)?);
        let x3 = k.sqr(&lam) ^ &lam ^ &self.a;
        let y3 = k.sqr(x) ^ k.mul(&(lam ^ BigUint::one()), &x3);
        Some((x3, y3))
    }

    pub fn mul(&self, e: &BigUint, p: &Point) -> Point {
        let mut r: Point = None;
        for i in (0..e.bits()).rev() {
            r = self.dbl(&r);
            if e.bit(i) {
                r = self.add(&r, p);
            }
        }
        r
    }
}

/// A document's integer: `0x`-hex or decimal, as a string or a number.
fn as_int(v: &J) -> Result<BigUint, String> {
    match v {
        J::Int(i) if *i >= 0 => Ok(BigUint::from(*i as u128)),
        J::Str(s) => match s.strip_prefix("0x") {
            Some(h) => BigUint::parse_bytes(h.as_bytes(), 16),
            None => BigUint::parse_bytes(s.as_bytes(), 10),
        }
        .ok_or_else(|| format!("`{s}` is not an integer")),
        other => Err(format!(
            "{} is not an integer",
            json::dumps_line(other, false)
        )),
    }
}

fn point(v: &J) -> Result<Point, String> {
    Ok(Some((as_int(v.at("x")?)?, as_int(v.at("y")?)?)))
}

/// `[scalar]G = Q` on the document's curve, and the known logarithm: a
/// binary document here, a `prime_extension` one in B5a's arithmetic.
pub fn replay(doc: &J, scalar: Option<&BigUint>) -> Result<J, String> {
    let field = doc.at("field")?;
    if field.get("kind").and_then(J::as_str) == Some("prime_extension") {
        return b5a::replay(doc, scalar);
    }
    let n = field
        .at("degree")?
        .as_i128()
        .and_then(|n| u64::try_from(n).ok())
        .ok_or("`field.degree` is not a degree")?;
    let curve_doc = doc.at("curve")?;
    let b = match curve_doc.get("b") {
        Some(b) => as_int(b)?,
        None => BigUint::one(),
    };
    let curve = Curve::new(
        Field::new(n, as_int(field.at("modulus")?)?)?,
        as_int(curve_doc.at("a")?)?,
        b,
    )?;
    let sub = doc.at("subgroup")?;
    let g = point(sub.at("generator")?)?;
    let order = as_int(sub.at("order")?)?;
    let target = doc.at("target")?;
    let known = target.get("known_log").map(as_int).transpose()?;
    let q = match &known {
        Some(d) => curve.mul(d, &g),
        None => point(target.at("point")?)?,
    };
    // A point off the curve replays nothing: the document is not what it says.
    let on_curve = curve.on_curve(&g) && curve.on_curve(&q);
    let replays = on_curve && scalar.is_some_and(|d| curve.mul(d, &g) == q);
    let known_matches = match (&known, scalar) {
        (None, _) => J::Null,
        (Some(d), Some(s)) => J::Bool(*s == d % &order),
        (Some(_), None) => J::Bool(false),
    };
    Ok(obj([
        ("replays", J::Bool(replays)),
        ("points_on_curve", J::Bool(on_curve)),
        (
            "known_log",
            known.map_or(J::Null, |d| J::Str(d.to_str_radix(10))),
        ),
        ("known_matches", known_matches),
    ]))
}

// ── the instance sets ──────────────────────────────────────────────

/// One run: its id in the run tree, its document, the method put in its
/// place if any, and its time limit.
pub struct Instance {
    pub run_id: String,
    pub document: PathBuf,
    pub method: Option<J>,
    pub timeout_s: u64,
    /// Whether the run passes `--repeats 1 --repeats-fast 1`: every set's
    /// runs but B5a's `solve: rho` documents, as its `run.py` declared.
    pub repeats_one: bool,
}

/// B4's `HOUR_S`; the tool's own default budget for the others.
const HOUR_S: u64 = 3600;
const DAY_S: u64 = 86_400;

fn slug_of(doc: &Path) -> Result<String, String> {
    let name = json::read(doc)?
        .at("name")?
        .as_str()
        .ok_or("`name` is not a string")?
        .to_string();
    let start = name
        .find("icv1-")
        .ok_or_else(|| format!("{}: its name has no ICV1 slug", doc.display()))?;
    Ok(name[start..]
        .chars()
        .take_while(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || *c == '-')
        .collect())
}

/// The runs of a set: `b4` (C088–C093), `b3b` (C078–C082), each instance's
/// known-logarithm case and its public point `T001`; or `b3-gate`, the
/// gate's public point under rho alone.
pub fn instances(programme: &Path, set: &str) -> Result<Vec<Instance>, String> {
    let conf = programme.join("conformance");
    let pairs = |dir: &str, word: &str, list: &[(&str, u32)], timeout_s: u64| {
        let mut out = Vec::new();
        for (cid, n) in list {
            let params = conf.join(dir).join("params");
            let known = params.join(format!("{cid}-kic-{word}-n{n}.json"));
            let slug = slug_of(&known)?;
            out.push(Instance {
                run_id: format!("{slug}/{cid}-known"),
                document: known,
                method: None,
                timeout_s,
                repeats_one: true,
            });
            out.push(Instance {
                run_id: format!("{slug}/{cid}-T001"),
                document: params.join(format!("{cid}-T001-n{n}.json")),
                method: None,
                timeout_s,
                repeats_one: true,
            });
        }
        Ok::<_, String>(out)
    };
    match set {
        "b4" => pairs(
            "v2-b4",
            "three-word",
            &[
                ("C088", 127),
                ("C089", 137),
                ("C090", 151),
                ("C091", 157),
                ("C092", 173),
                ("C093", 179),
            ],
            HOUR_S,
        ),
        "b3b" => pairs(
            "v2-b3b",
            "two-word",
            &[
                ("C078", 67),
                ("C079", 67),
                ("C080", 79),
                ("C081", 71),
                ("C082", 83),
            ],
            DAY_S,
        ),
        "b3-gate" => {
            // B3's measurement 5: the gate's document with rho alone at F0
            // and the document's own seed.
            let document = conf.join("v2").join("params").join("gate-m83-T001.json");
            let seed = json::read(&document)?
                .at("method")?
                .at("rho")?
                .at("seed")?
                .clone();
            Ok(vec![Instance {
                run_id: format!("{}/gate-T001-rho", suite::curve_slug(0, 83)?),
                document,
                method: Some(obj([
                    ("solve", J::Str("rho".into())),
                    ("fidelity", J::Str("F0".into())),
                    (
                        "rho",
                        obj([("pipeline", J::Str("auto".into())), ("seed", seed)]),
                    ),
                ])),
                timeout_s: DAY_S,
                repeats_one: true,
            }])
        }
        "b5a" => {
            // B5a's measurement 5: each instance's two documents, keyed by
            // its slug, then C050's on its own known-answer target.  Only
            // the paired documents pass the repetition flags.
            let params = conf.join("v2-b5a").join("params");
            let mut out = Vec::new();
            for id in ["G1", "G2", "G3", "E2", "E5", "E11", "B2"] {
                for which in ["known", "T001"] {
                    let document = params.join(format!("{id}-{which}.json"));
                    let paired = json::read(&document)?.at("method")?.at("solve")?.as_str()
                        == Some("paired");
                    out.push(Instance {
                        run_id: format!("{}/{id}-{which}", slug_of(&document)?),
                        document,
                        method: None,
                        timeout_s: HOUR_S,
                        repeats_one: paired,
                    });
                }
            }
            out.push(Instance {
                run_id: "icv1-fp10k3-t41822-dfa3991d/C050-known".into(),
                document: conf
                    .join("v2-b2")
                    .join("params")
                    .join("C050-cubic-extension.json"),
                method: None,
                timeout_s: HOUR_S,
                repeats_one: false,
            });
            Ok(out)
        }
        other => Err(format!(
            "unknown instance set {other:?}; the sets are b3-gate, b3b, b4 and b5a"
        )),
    }
}

/// The document a run prices: the instance's own, or a copy with its
/// method replaced, written beside the run once.
fn priced_document(tree: &Path, inst: &Instance) -> Result<PathBuf, String> {
    let Some(method) = &inst.method else {
        return Ok(inst.document.clone());
    };
    let path = tree.join(format!("{}.params.json", inst.run_id));
    if !path.exists() {
        let J::Obj(mut kv) = json::read(&inst.document)? else {
            return Err(format!("{} is not an object", inst.document.display()));
        };
        for (k, v) in kv.iter_mut() {
            if k == "method" {
                *v = method.clone();
            }
        }
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
        }
        std::fs::write(&path, json::dumps(&J::Obj(kv), 1) + "\n")
            .map_err(|e| format!("{}: {e}", path.display()))?;
    }
    Ok(path)
}

/// Every run of the set, each at most three attempts.
pub fn run(b: &Bench, ic: &Path, list: &[Instance], runs_dir: &Path) -> Result<J, String> {
    let tree = runs_dir.join("f0");
    for inst in list {
        let doc = priced_document(&tree, inst)?;
        let out = tree.join(format!("{}.price.json", inst.run_id));
        let mut rep = J::Null;
        for k in 0..=runs::RETRIES {
            let attempt = runs::attempt(&out, k);
            if !attempt.exists() {
                if let Some(parent) = attempt.parent() {
                    std::fs::create_dir_all(parent)
                        .map_err(|e| format!("{}: {e}", parent.display()))?;
                }
                let mut cmd: Vec<String> = vec![
                    "timeout".into(),
                    inst.timeout_s.to_string(),
                    ic.to_string_lossy().into_owned(),
                    "price".into(),
                    "--params".into(),
                    doc.to_string_lossy().into_owned(),
                    "--json".into(),
                    "--out".into(),
                    attempt.to_string_lossy().into_owned(),
                ];
                if inst.repeats_one {
                    cmd.extend(["--repeats", "1", "--repeats-fast", "1"].map(String::from));
                }
                b.launch(&cmd, &attempt, &tree)?;
            }
            rep = runs::load(&attempt);
            if runs::clean(&attempt)? && rep.get("status").and_then(J::as_str) == Some("complete") {
                break;
            }
        }
        let result = rep.get("result");
        let shown = |k: &str| {
            result
                .and_then(|r| r.get(k))
                .map_or("None".to_string(), |v| json::dumps_line(v, false))
        };
        println!(
            "{}: {} scalar {} verified {}",
            inst.run_id,
            rep.get("status").and_then(J::as_str).unwrap_or("None"),
            shown("scalar"),
            shown("verified")
        );
    }
    Ok(J::Null)
}

fn scalar_of(rep: &J, arm: &str) -> Result<Option<BigUint>, String> {
    rep.get("certificates")
        .and_then(|c| c.get(arm))
        .and_then(|c| c.get("scalar"))
        .filter(|s| !matches!(s, J::Null))
        .map(as_int)
        .transpose()
}

/// One run's row: its figure, both arms' replays, its costs and counts.
fn row(programme: &Path, inst: &Instance, tree: &Path) -> Result<J, String> {
    let out = tree.join(format!("{}.price.json", inst.run_id));
    let attempts: Vec<PathBuf> = (0..=runs::RETRIES)
        .map(|k| runs::attempt(&out, k))
        .collect();
    let mut figure = out.clone();
    for a in &attempts {
        if a.exists()
            && runs::clean(a)?
            && runs::load(a).get("status").and_then(J::as_str) == Some("complete")
        {
            figure = a.clone();
            break;
        }
    }
    let rep = runs::load(&figure);
    let doc = json::read(&priced_document(tree, inst)?)?;
    let rho_only = doc
        .get("method")
        .and_then(|m| m.get("solve"))
        .and_then(J::as_str)
        == Some("rho");
    let ic_scalar = scalar_of(&rep, "ic")?;
    // A `solve: rho` report may name its scalar only in its result.
    let rho_scalar = match scalar_of(&rep, "rho")? {
        None if rho_only => rep
            .get("result")
            .and_then(|r| r.get("scalar"))
            .filter(|s| !matches!(s, J::Null))
            .map(as_int)
            .transpose()?,
        s => s,
    };
    let ic_replay = if rho_only {
        J::Null
    } else {
        replay(&doc, ic_scalar.as_ref())?
    };
    let rho_replay = replay(&doc, rho_scalar.as_ref())?;
    let replays = |r: &J| r.get("replays").is_some_and(J::truthy);
    let known_ok = |r: &J| !matches!(r.get("known_matches"), Some(J::Bool(false)));
    let complete = rep.get("status").and_then(J::as_str) == Some("complete");
    let passed = complete
        && replays(&rho_replay)
        && known_ok(&rho_replay)
        && (rho_only || (replays(&ic_replay) && ic_scalar == rho_scalar && known_ok(&ic_replay)));
    let record = runs::run_record(&figure)?.unwrap_or(J::Null);
    let median = rep.get("median").cloned().unwrap_or(J::Null);
    let first = rep
        .get("repetitions")
        .and_then(J::as_arr)
        .and_then(|r| r.first())
        .cloned()
        .unwrap_or(J::Null);
    let get = |v: &J, k: &str| v.get(k).cloned().unwrap_or(J::Null);
    let document = inst
        .document
        .strip_prefix(programme)
        .unwrap_or(&inst.document)
        .to_string_lossy()
        .into_owned();
    Ok(obj([
        ("run", J::Str(inst.run_id.clone())),
        ("document", J::Str(document)),
        (
            "figure",
            J::Str(
                figure
                    .file_name()
                    .map(|f| f.to_string_lossy().into_owned())
                    .unwrap_or_default(),
            ),
        ),
        (
            "attempts",
            J::Int(attempts.iter().filter(|a| a.exists()).count() as i128),
        ),
        ("exit_status", get(&record, "exit_status")),
        ("contended", get(&record, "contended")),
        ("status", get(&rep, "status")),
        (
            "words",
            rep.get("ic").map_or(J::Null, |ic| get(ic, "words")),
        ),
        (
            "scalar",
            rho_scalar
                .as_ref()
                .or(ic_scalar.as_ref())
                .map_or(J::Null, |s| J::Str(s.to_str_radix(10))),
        ),
        (
            "verified_in_run",
            rep.get("result").map_or(J::Null, |r| get(r, "verified")),
        ),
        ("ic_and_rho_agree", get(&rep, "ic_and_rho_agree")),
        ("replay_ic", ic_replay),
        ("replay_rho", rho_replay),
        ("s_ic_cold", get(&median, "s_ic_cold")),
        ("s_rho_cold", get(&median, "s_rho_cold")),
        (
            "cold_ratio_ic_over_rho",
            get(&median, "cold_ratio_ic_over_rho"),
        ),
        ("online_speedup", get(&median, "online_speedup")),
        ("setup_phases_ns", get(&first, "setup_phases_ns")),
        ("ic_online_phases_ns", get(&median, "ic_online_phases_ns")),
        ("counts", get(&rep, "counts")),
        ("rho_counts", get(&rep, "rho_counts")),
        ("pass", J::Bool(passed)),
    ]))
}

pub fn analyse(
    programme: &Path,
    what: &str,
    list: &[Instance],
    runs_dir: &Path,
) -> Result<J, String> {
    let tree = runs_dir.join("f0");
    let rows = list
        .iter()
        .map(|inst| row(programme, inst, &tree))
        .collect::<Result<Vec<J>, String>>()?;
    let passed = rows
        .iter()
        .filter(|r| r.get("pass").is_some_and(J::truthy))
        .count();
    Ok(obj([
        ("measurement", J::Str(what.into())),
        ("runs", J::Int(rows.len() as i128)),
        ("passed", J::Int(passed as i128)),
        ("pass", J::Bool(passed == rows.len())),
        ("rows", J::Arr(rows)),
    ]))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn programme() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("research/ic_tool_program")
    }

    #[test]
    fn the_gates_frozen_points_are_on_the_curve_and_of_the_subgroups_order() {
        let doc =
            json::read(&programme().join("conformance/v2/params/gate-m83-T001.json")).unwrap();
        let field = doc.at("field").unwrap();
        let curve = Curve::new(
            Field::new(83, as_int(field.at("modulus").unwrap()).unwrap()).unwrap(),
            BigUint::zero(),
            BigUint::one(),
        )
        .unwrap();
        let sub = doc.at("subgroup").unwrap();
        let order = as_int(sub.at("order").unwrap()).unwrap();
        for p in [
            point(sub.at("generator").unwrap()).unwrap(),
            point(doc.at("target").unwrap().at("point").unwrap()).unwrap(),
        ] {
            assert!(curve.on_curve(&p));
            assert!(curve.mul(&order, &p).is_none());
            assert!(curve.mul(&(&order - 1u32), &p).is_some());
        }
    }

    #[test]
    fn a_known_logarithm_replays_and_a_wrong_one_does_not() {
        for (dir, file) in [
            ("v2-b3b", "C078-kic-two-word-n67.json"),
            ("v2-b4", "C088-kic-three-word-n127.json"),
        ] {
            let doc = json::read(
                &programme()
                    .join("conformance")
                    .join(dir)
                    .join("params")
                    .join(file),
            )
            .unwrap();
            let known = as_int(doc.at("target").unwrap().at("known_log").unwrap()).unwrap();
            let order = as_int(doc.at("subgroup").unwrap().at("order").unwrap()).unwrap();
            let good = replay(&doc, Some(&known)).unwrap();
            assert_eq!(good.get("replays"), Some(&J::Bool(true)), "{file}");
            assert_eq!(good.get("known_matches"), Some(&J::Bool(true)), "{file}");
            // d + r replays too (it names the same point), and is reduced
            // before it is compared with the known logarithm.
            let wrapped = replay(&doc, Some(&(&known + &order))).unwrap();
            assert_eq!(wrapped.get("replays"), Some(&J::Bool(true)), "{file}");
            let bad = replay(&doc, Some(&(&known + 1u32))).unwrap();
            assert_eq!(bad.get("replays"), Some(&J::Bool(false)), "{file}");
            assert_eq!(bad.get("known_matches"), Some(&J::Bool(false)), "{file}");
            let none = replay(&doc, None).unwrap();
            assert_eq!(none.get("replays"), Some(&J::Bool(false)), "{file}");
        }
    }

    #[test]
    fn the_field_inverts_and_the_group_law_is_consistent() {
        // GF(2^7) under z^7 + z + 1, and y² + xy = x³ + x² + 1 over it.
        let k = Field::new(7, BigUint::from(0b1000_0011u32)).unwrap();
        for a in 1u32..128 {
            let a = BigUint::from(a);
            assert!(k.mul(&a, &k.inv(&a).unwrap()).is_one());
        }
        let curve = Curve::new(k, BigUint::one(), BigUint::one()).unwrap();
        let mut points = vec![None];
        for x in 0u32..128 {
            for y in 0u32..128 {
                let p = Some((BigUint::from(x), BigUint::from(y)));
                if curve.on_curve(&p) {
                    points.push(p);
                }
            }
        }
        let order = BigUint::from(points.len());
        for p in points.iter().take(12) {
            assert!(curve.mul(&order, p).is_none());
            for q in points.iter().skip(5).take(6) {
                let sum = curve.add(p, q);
                assert!(curve.on_curve(&sum));
                assert_eq!(sum, curve.add(q, p));
            }
        }
    }

    #[test]
    fn the_sets_name_their_runs_by_slug() {
        let b4 = instances(&programme(), "b4").unwrap();
        assert_eq!(b4.len(), 12);
        assert!(b4[0].run_id.starts_with("icv1-f2m127-") && b4[0].run_id.ends_with("/C088-known"));
        assert!(b4
            .iter()
            .all(|i| i.document.exists() && i.timeout_s == HOUR_S));
        let b3b = instances(&programme(), "b3b").unwrap();
        assert_eq!(b3b.len(), 10);
        assert!(b3b.iter().all(|i| i.document.exists()));
        let gate = instances(&programme(), "b3-gate").unwrap();
        assert_eq!(gate.len(), 1);
        assert_eq!(
            gate[0].run_id,
            "icv1-f2m83-tm6151469093347-debefd74/gate-T001-rho"
        );
        assert!(instances(&programme(), "b5").is_err());
        // B5a: two runs an instance and C050's, only the paired ones with
        // the repetition flags, every document on file.
        let b5a = instances(&programme(), "b5a").unwrap();
        assert_eq!(b5a.len(), 15);
        assert!(b5a
            .iter()
            .all(|i| i.document.exists() && i.timeout_s == HOUR_S));
        assert_eq!(b5a.iter().filter(|i| i.repeats_one).count(), 6);
        assert_eq!(b5a[0].run_id, "icv1-fp9k3-tm4575-aba104ba/G1-known");
        assert_eq!(b5a[14].run_id, "icv1-fp10k3-t41822-dfa3991d/C050-known");
    }
}
