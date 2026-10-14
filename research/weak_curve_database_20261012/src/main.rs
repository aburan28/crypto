//! Fixed public research-corpus importer. It accepts no curve, target or input path.
use num_bigint::{BigInt, BigUint};
use rusqlite::{params, types::ValueRef, Connection};
use serde_json::Value;
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet, HashMap, HashSet},
    env,
    error::Error,
    fs, io,
    path::{Path, PathBuf},
};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
type ModelIdentity = HashMap<String, (String, String, String)>;
const EXACT_PRIMES: [i64; 7] = [7, 11, 13, 17, 19, 29, 37];
const DIAGNOSTIC_PRIMES: [i64; 3] = [41, 43, 47];

fn bad(message: impl Into<String>) -> Box<dyn Error> {
    Box::new(io::Error::new(io::ErrorKind::InvalidData, message.into()))
}
fn ensure(ok: bool, message: impl Into<String>) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(bad(message))
    }
}
fn s<'a>(v: &'a Value, key: &str) -> Result<&'a str> {
    v.get(key)
        .and_then(Value::as_str)
        .ok_or_else(|| bad(format!("missing string {key}")))
}
fn a<'a>(v: &'a Value, key: &str) -> Result<&'a [Value]> {
    v.get(key)
        .and_then(Value::as_array)
        .map(Vec::as_slice)
        .ok_or_else(|| bad(format!("missing array {key}")))
}
fn str_array(v: &Value, key: &str) -> Result<Vec<String>> {
    a(v, key)?
        .iter()
        .map(|x| {
            x.as_str()
                .map(str::to_owned)
                .ok_or_else(|| bad(format!("non-string {key}")))
        })
        .collect()
}
fn sha(bytes: &[u8]) -> String {
    hex::encode(Sha256::digest(bytes))
}
fn biguint(x: &str) -> Result<BigUint> {
    BigUint::parse_bytes(x.as_bytes(), 10).ok_or_else(|| bad(format!("bad unsigned integer {x}")))
}
fn bigint(x: &str) -> Result<BigInt> {
    BigInt::parse_bytes(x.as_bytes(), 10).ok_or_else(|| bad(format!("bad signed integer {x}")))
}
fn json(v: &Value) -> String {
    serde_json::to_string(v).expect("JSON value serializes")
}
fn nullable(v: &Value, key: &str) -> Result<Option<String>> {
    match v.get(key) {
        Some(Value::Null) | None => Ok(None),
        Some(Value::String(x)) => Ok(Some(x.clone())),
        _ => Err(bad(format!("bad optional string {key}"))),
    }
}
fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf()
}
fn study() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}
fn database_path() -> PathBuf {
    study().join("weak_curves.sqlite")
}

struct Importer<'a> {
    db: &'a Connection,
    root: PathBuf,
    lock: BTreeMap<String, String>,
    seen: BTreeSet<String>,
    sources: HashSet<String>,
}
impl<'a> Importer<'a> {
    fn new(db: &'a Connection) -> Result<Self> {
        let bytes = fs::read(study().join("INPUT_LOCK.json"))?;
        let lock: BTreeMap<String, String> = serde_json::from_slice(&bytes)?;
        Ok(Self {
            db,
            root: root(),
            lock,
            seen: BTreeSet::new(),
            sources: HashSet::new(),
        })
    }
    fn input(&mut self, id: &str, path: &str, grade: &str, scope: &str) -> Result<Vec<u8>> {
        ensure(!self.seen.contains(path), format!("duplicate input {path}"))?;
        let bytes = fs::read(self.root.join(path))?;
        let digest = sha(&bytes);
        let expected = self
            .lock
            .get(path)
            .ok_or_else(|| bad(format!("unlocked input {path}")))?;
        ensure(expected == &digest, format!("input hash mismatch {path}"))?;
        self.db.execute(
            "INSERT INTO input_file VALUES (?1,?2,?3,?4,?5,?6)",
            params![id, path, digest, bytes.len() as i64, grade, scope],
        )?;
        self.seen.insert(path.to_owned());
        Ok(bytes)
    }
    fn local_json(&mut self, id: &str, name: &str) -> Result<Value> {
        let path = format!("research/weak_curve_database_20261012/{name}");
        Ok(serde_json::from_slice(&self.input(
            id,
            &path,
            "authored_registry",
            "metadata",
        )?)?)
    }
    fn catalog(&mut self) -> Result<Value> {
        Ok(serde_json::from_slice(&self.input(
            "families",
            "docs/curves/weak-families/catalog.json",
            "published_and_parent_registry",
            "family_and_condition",
        )?)?)
    }
    fn source_refs(&mut self, catalog: &Value) -> Result<()> {
        let extra = self.local_json("extra_sources", "sources.json")?;
        ensure(
            s(&extra, "schema")? == "weak_curve.extra_sources/v1",
            "extra source schema",
        )?;
        for source in a(catalog, "sources")?.iter().chain(a(&extra, "sources")?) {
            let id = s(source, "id")?;
            ensure(
                self.sources.insert(id.to_owned()),
                format!("duplicate source id {id}"),
            )?;
            let title = s(source, "title")?;
            let url = nullable(source, "url")?;
            let path = nullable(source, "path")?;
            ensure(
                url.is_some() ^ path.is_some(),
                format!("source {id} has wrong link"),
            )?;
            self.db.execute(
                "INSERT INTO source_ref VALUES (?1,?2,?3,?4)",
                params![id, title, url, path],
            )?;
        }
        Ok(())
    }
    fn refs(&self, record: &Value, key: &str) -> Result<String> {
        let ids = str_array(record, key)?;
        ensure(!ids.is_empty(), format!("empty {key}"))?;
        for id in &ids {
            ensure(self.sources.contains(id), format!("unknown source id {id}"))?;
        }
        Ok(json(&serde_json::to_value(ids)?))
    }
    fn families(&self, catalog: &Value) -> Result<usize> {
        ensure(
            s(catalog, "schema")? == "weak-curve-family-catalog/v1",
            "family schema",
        )?;
        let rows = a(catalog, "families")?;
        for f in rows {
            self.db.execute(
                "INSERT INTO family VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11,?12,?13,?14,?15)",
                params![
                    s(f, "id")?,
                    s(f, "name")?,
                    s(f, "category")?,
                    s(f, "scope")?,
                    s(f, "evidence")?,
                    s(f, "local_status")?,
                    s(f, "domain")?,
                    s(f, "recognition")?,
                    s(f, "find_members")?,
                    s(f, "exhaustive_scope")?,
                    s(f, "subgroup_gate")?,
                    s(f, "cost_gate")?,
                    self.refs(f, "sources")?,
                    json(f),
                    "families"
                ],
            )?;
        }
        ensure(rows.len() == 21, "family count changed")?;
        Ok(rows.len())
    }
    fn models(&mut self) -> Result<(Value, ModelIdentity)> {
        let path = "docs/curves/weak-families/models.json";
        let bytes = self.input("models", path, "parent_model_inventory", "exact_model")?;
        let data: Value = serde_json::from_slice(&bytes)?;
        ensure(
            s(&data, "schema")? == "weak-curve-model-index/v1",
            "model schema",
        )?;
        let catalog_hash: String = self.db.query_row(
            "SELECT sha256 FROM input_file WHERE id='families'",
            [],
            |r| r.get(0),
        )?;
        ensure(
            s(&data, "family_catalog_sha256")? == catalog_hash,
            "family hash binding",
        )?;
        let mut identity = HashMap::new();
        for m in a(&data, "models")? {
            let icv1 = s(m, "icv1")?;
            let model_sha = s(m, "model_sha256")?;
            let trace = s(m, "trace")?;
            let order = s(m, "order")?;
            let field = m.get("field").ok_or_else(|| bad("missing model field"))?;
            let p = s(field, "characteristic")?;
            let degree = field
                .get("degree")
                .and_then(Value::as_u64)
                .ok_or_else(|| bad("missing degree"))? as u32;
            let bits = field
                .get("bit_length")
                .and_then(Value::as_u64)
                .ok_or_else(|| bad("missing bit length"))? as i64;
            let cardinality = s(field, "cardinality")?;
            ensure(
                biguint(p)?.pow(degree) == biguint(cardinality)?,
                format!("bad field size {icv1}"),
            )?;
            let expected_order = BigInt::from(biguint(cardinality)?) + 1 - bigint(trace)?;
            ensure(
                expected_order == bigint(order)?,
                format!("trace/order mismatch {icv1}"),
            )?;
            ensure(
                model_sha.len() == 64 && icv1.ends_with(&model_sha[..12]),
                format!("model identity mismatch {icv1}"),
            )?;
            let labels = m.get("labels").ok_or_else(|| bad("missing labels"))?;
            let model_json = s(m, "model_json")?;
            let _: Value = serde_json::from_str(model_json)?;
            let subgroup_order = nullable(m, "subgroup_order")?;
            self.db.execute(
                "INSERT INTO curve_model VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11,?12,?13,?14,?15,?16,?17,?18)",
                params![icv1,model_sha,s(m,"slug")?,p,degree as i64,bits,cardinality,
                    trace,order,subgroup_order,s(m,"inventory")?,s(labels,"combined_class")?,
                    s(labels,"model_g3_cubic")?,s(labels,"model_g3_quadratic")?,
                    s(labels,"subgroup_advantage")?,model_json,json(m),"models"],
            )?;
            ensure(
                identity
                    .insert(
                        icv1.to_owned(),
                        (model_sha.to_owned(), trace.to_owned(), order.to_owned()),
                    )
                    .is_none(),
                format!("duplicate ICV1 {icv1}"),
            )?;
        }
        ensure(identity.len() == 814, "model count changed")?;
        Ok((data, identity))
    }
    fn properties(&mut self, identity: &ModelIdentity) -> Result<()> {
        let bytes = self.input(
            "properties",
            "docs/curves/weak-families/properties.json",
            "native_and_pari_property_inventory",
            "exact_model_property",
        )?;
        let data: Value = serde_json::from_slice(&bytes)?;
        ensure(
            s(&data, "schema")? == "weak-curve-structural-properties/v1",
            "property schema",
        )?;
        let models_hash: String =
            self.db
                .query_row("SELECT sha256 FROM input_file WHERE id='models'", [], |r| {
                    r.get(0)
                })?;
        ensure(
            s(
                data.get("model_input").ok_or_else(|| bad("model_input"))?,
                "sha256",
            )? == models_hash,
            "model/property hash binding",
        )?;
        let mut seen = HashSet::new();
        for p in a(&data, "records")? {
            let id = s(p, "icv1")?;
            let parent = identity
                .get(id)
                .ok_or_else(|| bad(format!("orphan property {id}")))?;
            ensure(
                s(p, "model_sha256")? == parent.0
                    && s(p, "trace")? == parent.1
                    && s(p, "order")? == parent.2,
                format!("property/model mismatch {id}"),
            )?;
            ensure(seen.insert(id), format!("duplicate property {id}"))?;
            self.db.execute(
                "INSERT INTO curve_property VALUES (?1,?2,?3,?4,?5,?6)",
                params![
                    id,
                    json(
                        p.get("properties")
                            .ok_or_else(|| bad("missing properties"))?
                    ),
                    json(
                        p.get("family_property_checks")
                            .ok_or_else(|| bad("missing checks"))?
                    ),
                    json(p.get("witnesses").ok_or_else(|| bad("missing witnesses"))?),
                    json(p),
                    "properties"
                ],
            )?;
        }
        ensure(
            seen.len() == 814 && seen.len() == identity.len(),
            "property coverage changed",
        )?;
        Ok(())
    }
    fn exact_classes(&mut self, models: &Value) -> Result<()> {
        let mut p7 = HashMap::new();
        for row in a(models, "class_oracle")? {
            let trace: i64 = s(row, "trace")?.parse()?;
            p7.insert(trace, row.clone());
        }
        ensure(p7.len() == 294, "p7 parent oracle count")?;
        let mut total = 0;
        for p in EXACT_PRIMES {
            let path =
                format!("research/iso1_weak_classes_20261007/two_branch_20261009/p{p}_classes.csv");
            let id = format!("exact_p{p}");
            let bytes = self.input(&id, &path, "complete_exact", "two_branch_union_ordinary")?;
            let text = std::str::from_utf8(&bytes)?;
            let mut lines = text.lines();
            ensure(
                lines.next() == Some("p,trace,cubic,quadratic,combined"),
                format!("bad header {path}"),
            )?;
            let mut traces = BTreeSet::new();
            for (i, line) in lines.enumerate() {
                let fields: Vec<_> = line.split(',').collect();
                ensure(
                    fields.len() == 5,
                    format!("bad exact CSV line {path}:{}", i + 2),
                )?;
                ensure(fields[0].parse::<i64>()? == p, "wrong prime in exact CSV")?;
                let t: i64 = fields[1].parse()?;
                let c = bool01(fields[2])?;
                let q = bool01(fields[3])?;
                let u = bool01(fields[4])?;
                ensure(u == (c | q), format!("bad union {p}:{t}"))?;
                ensure(
                    is_ordinary(p, t),
                    format!("out-of-domain exact trace {p}:{t}"),
                )?;
                ensure(traces.insert(t), format!("duplicate exact trace {p}:{t}"))?;
                if p == 7 {
                    let parent = p7
                        .get(&t)
                        .ok_or_else(|| bad(format!("p7 oracle missing {t}")))?;
                    ensure(
                        parent["cubic_support"] == Value::Bool(c == 1)
                            && parent["quadratic_support"] == Value::Bool(q == 1)
                            && parent["combined_support"] == Value::Bool(u == 1),
                        format!("p7 oracle mismatch {t}"),
                    )?;
                }
                self.db.execute(
                    "INSERT INTO class_observation VALUES (?1,?2,'two_branch_union','ordinary','complete_exact',?3,?4,?5,NULL,NULL,NULL,NULL,NULL,NULL,?6,?7)",
                    params![p,t,c,q,u,(i+2) as i64,id],
                )?;
            }
            verify_exact_domain(p, &traces)?;
            total += traces.len();
        }
        ensure(total == 87_486, "exact class total changed")?;
        Ok(())
    }
    fn diagnostic_classes(&mut self) -> Result<()> {
        let mut total = 0;
        for p in DIAGNOSTIC_PRIMES {
            let suffix = if p == 47 { "orbit" } else { "twist_derived" };
            let path = format!("research/iso1_weak_classes_20261007/p{p}_{suffix}.csv");
            let id = format!("diagnostic_p{p}");
            let bytes = self.input(
                &id,
                &path,
                "probabilistic_diagnostic",
                "cubic_norm_one_hasse_candidates",
            )?;
            let text = std::str::from_utf8(&bytes)?;
            let mut lines = text.lines();
            let header = lines.next().ok_or_else(|| bad("empty diagnostic CSV"))?;
            ensure(header=="p,q,trace,weak_representatives,fundamental_discriminant,frobenius_2_depth,two_split,v2_delta_quarter,class_number_odd,trace_status",
                format!("bad diagnostic header {p}"))?;
            let names: Vec<_> = header.split(',').collect();
            let mut traces = BTreeSet::new();
            for (i, line) in lines.enumerate() {
                let fields: Vec<_> = line.split(',').collect();
                ensure(
                    fields.len() == names.len(),
                    format!("bad diagnostic line {p}:{}", i + 2),
                )?;
                ensure(
                    fields[0].parse::<i64>()? == p && fields[1].parse::<i64>()? == p * p,
                    format!("diagnostic p/q mismatch {p}:{}", i + 2),
                )?;
                let t: i64 = fields[2].parse()?;
                let reps: i64 = fields[3].parse()?;
                ensure(reps >= 0, "negative representative count")?;
                let status = fields[9];
                ensure(
                    status
                        == if is_ordinary(p, t) {
                            "ordinary"
                        } else {
                            "nonordinary_or_unrealized"
                        },
                    format!("diagnostic trace status mismatch {p}:{t}"),
                )?;
                ensure(
                    traces.insert(t),
                    format!("duplicate diagnostic trace {p}:{t}"),
                )?;
                let depth = if fields[5].is_empty() {
                    None
                } else {
                    Some(fields[5].parse::<i64>()?)
                };
                let fundamental = if fields[4].is_empty() {
                    None
                } else {
                    Some(fields[4].to_owned())
                };
                let split = if fields[6].is_empty() {
                    None
                } else {
                    Some(fields[6].parse::<i64>()?)
                };
                let v2_delta = if fields[7].is_empty() {
                    None
                } else {
                    Some(fields[7].parse::<i64>()?)
                };
                let class_odd = if fields[8].is_empty() {
                    None
                } else {
                    Some(fields[8].parse::<i64>()?)
                };
                self.db.execute(
                    "INSERT INTO class_observation VALUES (?1,?2,'cubic_norm_one',?3,'probabilistic_diagnostic',NULL,NULL,NULL,?4,?5,?6,?7,?8,?9,?10,?11)",
                    params![p,t,status,reps,depth,fundamental,split,v2_delta,class_odd,(i+2) as i64,id],
                )?;
            }
            verify_diagnostic_domain(p, &traces)?;
            total += traces.len();
        }
        ensure(total == 252_254, "diagnostic row total changed")?;
        Ok(())
    }
    fn conjectures(&mut self) -> Result<()> {
        let bytes = self.input(
            "conjectures",
            "research/weak_curve_conjectures_20261012/conjectures.json",
            "authored_hypothesis_ledger",
            "unproved_hypotheses",
        )?;
        let data: Value = serde_json::from_slice(&bytes)?;
        ensure(
            s(&data, "schema")? == "weak_curve.conjectures/v1",
            "conjecture schema",
        )?;
        let mut open = 0;
        let mut refuted = 0;
        for c in a(&data, "entries")? {
            match s(c, "status")? {
                "open" => open += 1,
                "refuted" => refuted += 1,
                _ => return Err(bad("unknown conjecture status")),
            }
            self.db.execute(
                "INSERT INTO conjecture VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9)",
                params![
                    s(c, "id")?,
                    s(c, "title")?,
                    s(c, "status")?,
                    s(c, "kind")?,
                    s(c, "domain")?,
                    s(c, "falsifier")?,
                    s(c, "lean_status")?,
                    json(c),
                    "conjectures"
                ],
            )?;
        }
        ensure(open == 7 && refuted == 1, "conjecture counts changed")?;
        Ok(())
    }
    fn theorems(&mut self) -> Result<()> {
        let data = self.local_json("theorems", "theorems.json")?;
        ensure(
            s(&data, "schema")? == "weak_curve.theorem_registry/v1",
            "theorem schema",
        )?;
        for t in a(&data, "theorems")? {
            self.db.execute(
                "INSERT INTO theorem VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11)",
                params![
                    s(t, "id")?,
                    s(t, "name")?,
                    s(t, "exact_statement")?,
                    s(t, "hypotheses")?,
                    s(t, "result_status")?,
                    nullable(t, "lean_declaration")?,
                    s(t, "lean_status")?,
                    nullable(t, "open_bridge")?,
                    self.refs(t, "source_ids")?,
                    json(t),
                    "theorems"
                ],
            )?;
        }
        ensure(a(&data, "theorems")?.len() == 10, "theorem count changed")?;
        Ok(())
    }
    fn heuristics(&mut self) -> Result<()> {
        let data = self.local_json("heuristics", "heuristics.json")?;
        ensure(
            s(&data, "schema")? == "weak_curve.heuristic_registry/v1",
            "heuristic schema",
        )?;
        for h in a(&data, "heuristics")? {
            self.db.execute(
                "INSERT INTO heuristic VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9,?10)",
                params![
                    s(h, "id")?,
                    s(h, "statement")?,
                    s(h, "scope")?,
                    s(h, "status")?,
                    s(h, "calibration")?,
                    s(h, "falsifier")?,
                    nullable(h, "linked_conjecture")?,
                    self.refs(h, "source_ids")?,
                    json(h),
                    "heuristics"
                ],
            )?;
        }
        ensure(
            a(&data, "heuristics")?.len() == 4,
            "heuristic count changed",
        )?;
        Ok(())
    }
    fn mechanisms(&mut self) -> Result<()> {
        let data = self.local_json("mechanisms", "mechanisms.json")?;
        ensure(
            s(&data, "schema")? == "weak_curve.mechanism_recipes/v1",
            "mechanism schema",
        )?;
        for m in a(&data, "mechanisms")? {
            self.db.execute(
                "INSERT INTO mechanism VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11,?12)",
                params![
                    s(m, "id")?,
                    s(m, "family_id")?,
                    s(m, "name")?,
                    s(m, "construction")?,
                    s(m, "recognition")?,
                    s(m, "witness_requirements")?,
                    s(m, "selected_subgroup_gate")?,
                    s(m, "measured_cost_status")?,
                    s(m, "validation_status")?,
                    self.refs(m, "source_ids")?,
                    json(m),
                    "mechanisms"
                ],
            )?;
        }
        ensure(
            a(&data, "mechanisms")?.len() == 12,
            "mechanism count changed",
        )?;
        Ok(())
    }
    fn proof_obligations(&mut self) -> Result<()> {
        let data = self.local_json("proof_obligations", "proof_obligations.json")?;
        ensure(
            s(&data, "schema")? == "weak_curve.family_proof_obligations/v1",
            "proof obligation schema",
        )?;
        let rows = a(&data, "rows")?;
        let mut covered = BTreeSet::new();
        for row in rows {
            let family_id = s(row, "family_id")?;
            ensure(
                covered.insert(family_id.to_owned()),
                format!("duplicate proof obligation {family_id}"),
            )?;
            self.db.execute(
                "INSERT INTO family_proof_obligation VALUES (?1,?2,?3,?4,?5,?6,?7,?8)",
                params![
                    family_id,
                    nullable(row, "theorem_id")?,
                    nullable(row, "heuristic_id")?,
                    s(row, "principal_claim")?,
                    s(row, "proof_state")?,
                    s(row, "next_bridge")?,
                    json(row),
                    "proof_obligations"
                ],
            )?;
        }
        let mut stmt = self.db.prepare("SELECT id FROM family ORDER BY id")?;
        let families: BTreeSet<String> = stmt
            .query_map([], |r| r.get(0))?
            .collect::<rusqlite::Result<_>>()?;
        ensure(covered == families, "proof obligation family coverage")?;
        Ok(())
    }
    fn finish(self) -> Result<()> {
        ensure(
            self.seen.len() == self.lock.len(),
            "input lock has unconsumed entries",
        )?;
        Ok(())
    }
}
fn bool01(x: &str) -> Result<i64> {
    match x {
        "0" => Ok(0),
        "1" => Ok(1),
        _ => Err(bad(format!("bad Boolean {x}"))),
    }
}
fn is_ordinary(p: i64, t: i64) -> bool {
    let bound = 2 * p * p * p;
    t.abs() < bound && t.rem_euclid(4) == 2 && t % p != 0
}
fn verify_exact_domain(p: i64, traces: &BTreeSet<i64>) -> Result<()> {
    let bound = 2 * p * p * p;
    let expected: BTreeSet<_> = (-bound..=bound).filter(|&t| is_ordinary(p, t)).collect();
    ensure(
        *traces == expected,
        format!("incomplete exact trace domain p={p}"),
    )?;
    ensure(
        traces.len() == (p * p * (p - 1)) as usize,
        format!("exact count p={p}"),
    )
}
fn verify_diagnostic_domain(p: i64, traces: &BTreeSet<i64>) -> Result<()> {
    let bound = 2 * p * p * p;
    let expected: BTreeSet<_> = (-bound..=bound).filter(|&t| t.rem_euclid(4) == 2).collect();
    ensure(
        *traces == expected,
        format!("incomplete Hasse candidate domain p={p}"),
    )?;
    ensure(
        traces.len() == (p * p * p + 1) as usize,
        format!("diagnostic count p={p}"),
    )
}

fn build_at(path: &Path) -> Result<()> {
    if path.exists() {
        fs::remove_file(path)?;
    }
    let mut db = Connection::open(path)?;
    db.execute_batch("PRAGMA journal_mode=DELETE; PRAGMA page_size=4096;")?;
    db.execute_batch(include_str!("../schema.sql"))?;
    let transaction = db.transaction()?;
    {
        let mut importer = Importer::new(&transaction)?;
        let catalog = importer.catalog()?;
        importer.source_refs(&catalog)?;
        importer.families(&catalog)?;
        let (models, identity) = importer.models()?;
        importer.properties(&identity)?;
        importer.exact_classes(&models)?;
        importer.diagnostic_classes()?;
        importer.conjectures()?;
        importer.heuristics()?;
        importer.theorems()?;
        importer.mechanisms()?;
        importer.proof_obligations()?;
        importer.finish()?;
    }
    transaction.commit()?;
    db.execute_batch("VACUUM;")?;
    integrity(&db)?;
    Ok(())
}
fn integrity(db: &Connection) -> Result<()> {
    let result: String = db.query_row("PRAGMA integrity_check", [], |r| r.get(0))?;
    ensure(result == "ok", format!("SQLite integrity: {result}"))?;
    let foreign: i64 = db.query_row("SELECT COUNT(*) FROM pragma_foreign_key_check", [], |r| {
        r.get(0)
    })?;
    ensure(foreign == 0, "foreign key errors")?;
    for (table, count) in [
        ("family", 21),
        ("curve_model", 814),
        ("curve_property", 814),
        ("class_observation", 339_740),
        ("conjecture", 8),
        ("heuristic", 4),
        ("theorem", 10),
        ("mechanism", 12),
        ("family_proof_obligation", 21),
    ] {
        let n: i64 = db.query_row(&format!("SELECT COUNT(*) FROM {table}"), [], |r| r.get(0))?;
        ensure(n == count, format!("{table}: expected {count}, got {n}"))?;
    }
    let exact: i64 = db.query_row(
        "SELECT COUNT(*) FROM class_observation WHERE label_grade='complete_exact'",
        [],
        |r| r.get(0),
    )?;
    let diag: i64 = db.query_row(
        "SELECT COUNT(*) FROM class_observation WHERE label_grade='probabilistic_diagnostic'",
        [],
        |r| r.get(0),
    )?;
    ensure(exact == 87_486 && diag == 252_254, "class grade totals")?;
    Ok(())
}
fn canonical_digest(db: &Connection) -> Result<String> {
    let mut hasher = Sha256::new();
    let tables = [
        ("input_file", "id"),
        ("source_ref", "id"),
        ("family", "id"),
        ("curve_model", "icv1"),
        ("curve_property", "icv1"),
        ("class_observation", "characteristic,trace,mechanism_scope"),
        ("conjecture", "id"),
        ("heuristic", "id"),
        ("theorem", "id"),
        ("mechanism", "id"),
        ("family_proof_obligation", "family_id"),
    ];
    for (table, key) in tables {
        hasher.update(table.as_bytes());
        let mut stmt = db.prepare(&format!("SELECT * FROM {table} ORDER BY {key}"))?;
        let cols = stmt.column_count();
        let mut rows = stmt.query([])?;
        while let Some(row) = rows.next()? {
            for i in 0..cols {
                match row.get_ref(i)? {
                    ValueRef::Null => hasher.update([0]),
                    ValueRef::Integer(x) => {
                        hasher.update([1]);
                        hasher.update(x.to_le_bytes());
                    }
                    ValueRef::Real(x) => {
                        hasher.update([2]);
                        hasher.update(x.to_le_bytes());
                    }
                    ValueRef::Text(x) => {
                        hasher.update([3]);
                        hasher.update((x.len() as u64).to_le_bytes());
                        hasher.update(x);
                    }
                    ValueRef::Blob(x) => {
                        hasher.update([4]);
                        hasher.update((x.len() as u64).to_le_bytes());
                        hasher.update(x);
                    }
                }
            }
        }
    }
    Ok(hex::encode(hasher.finalize()))
}
type SummaryRow = (i64, String, String, i64, i64, i64, i64);
fn summary_rows(db: &Connection) -> Result<Vec<SummaryRow>> {
    let mut stmt=db.prepare("SELECT characteristic,mechanism_scope,label_grade,rows,ordinary_rows,positive_ordinary_rows,zero_ordinary_rows FROM class_summary ORDER BY characteristic,mechanism_scope")?;
    let rows = stmt.query_map([], |r| {
        Ok((
            r.get(0)?,
            r.get(1)?,
            r.get(2)?,
            r.get(3)?,
            r.get(4)?,
            r.get(5)?,
            r.get(6)?,
        ))
    })?;
    Ok(rows.collect::<std::result::Result<Vec<_>, _>>()?)
}
fn products(db: &Connection) -> Result<(String, String)> {
    let rows = summary_rows(db)?;
    let mut tsv =
        String::from("p\tscope\tgrade\trows\tordinary\tpositive_ordinary\tzero_ordinary\n");
    let mut svg=String::from(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"900\" height=\"640\" viewBox=\"0 0 900 640\">\n\
         <rect width=\"900\" height=\"640\" fill=\"#f8fafc\"/>\n\
         <text x=\"24\" y=\"32\" font-family=\"sans-serif\" font-size=\"21\" font-weight=\"bold\">Class support in the fixed weak-curve censuses</text>\n\
         <rect x=\"24\" y=\"47\" width=\"16\" height=\"12\" fill=\"#075985\"/><text x=\"46\" y=\"58\" font-family=\"sans-serif\" font-size=\"13\">support recorded</text>\n\
         <rect x=\"190\" y=\"47\" width=\"16\" height=\"12\" fill=\"#bfcbd6\"/><text x=\"212\" y=\"58\" font-family=\"sans-serif\" font-size=\"13\">zero row in that mechanism</text>\n\
         <text x=\"24\" y=\"85\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"bold\">Exact two-branch ordinary trace domain</text>\n",
    );
    for (i, (p, scope, grade, n, ordinary, pos, zero)) in rows.iter().enumerate() {
        ensure(pos + zero == *ordinary, format!("summary partition p={p}"))?;
        tsv.push_str(&format!(
            "{p}\t{scope}\t{grade}\t{n}\t{ordinary}\t{pos}\t{zero}\n"
        ));
        if i == 7 {
            svg.push_str("<text x=\"24\" y=\"417\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"bold\">Cubic norm-one diagnostic; zero labels depend on randomized counting</text>\n");
        }
        let y = if i < 7 {
            103 + i as i64 * 41
        } else {
            435 + (i as i64 - 7) * 41
        };
        let width = if *ordinary == 0 {
            0
        } else {
            (*pos * 510 + ordinary / 2) / ordinary
        };
        svg.push_str(&format!(
            "<text x=\"24\" y=\"{y}\" font-family=\"sans-serif\" font-size=\"14\">p={p}</text>\
             <rect x=\"92\" y=\"{}\" width=\"510\" height=\"22\" fill=\"#bfcbd6\"/>\
             <rect x=\"92\" y=\"{}\" width=\"{width}\" height=\"22\" fill=\"#075985\"/>\
             <text x=\"620\" y=\"{y}\" font-family=\"sans-serif\" font-size=\"13\">{pos}/{ordinary} support</text>\n",
            y-16,y-16,
        ));
    }
    svg.push_str("<text x=\"24\" y=\"605\" font-family=\"sans-serif\" font-size=\"12\">Population: signed ordinary traces t=2 mod 4; counts are mechanism-specific. Source: CLASS_SUMMARY.tsv and hash-locked input CSVs.</text>\n</svg>\n");
    Ok((tsv, svg))
}
fn main() -> Result<()> {
    let command = env::args().nth(1).unwrap_or_else(|| "summary".into());
    ensure(
        env::args().count() <= 2,
        "usage: weak-curve-database [build|check|summary]",
    )?;
    let final_path = database_path();
    match command.as_str() {
        "build" => {
            let temp = study().join("weak_curves.new.sqlite");
            build_at(&temp)?;
            fs::rename(&temp, &final_path)?;
            let db = Connection::open(&final_path)?;
            let (tsv, svg) = products(&db)?;
            fs::write(study().join("CLASS_SUMMARY.tsv"), tsv)?;
            fs::write(study().join("class_counts.svg"), svg)?;
            println!(
                "BUILT\t{}\t{}",
                final_path.display(),
                canonical_digest(&db)?
            );
        }
        "check" => {
            let db = Connection::open_with_flags(
                &final_path,
                rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY,
            )?;
            integrity(&db)?;
            let temp = study().join("weak_curves.check.sqlite");
            build_at(&temp)?;
            let fresh =
                Connection::open_with_flags(&temp, rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY)?;
            let actual = canonical_digest(&db)?;
            let expected = canonical_digest(&fresh)?;
            let (tsv, svg) = products(&fresh)?;
            fs::remove_file(&temp)?;
            ensure(
                actual == expected,
                format!("database content differs: {actual} != {expected}"),
            )?;
            ensure(
                fs::read_to_string(study().join("CLASS_SUMMARY.tsv"))? == tsv,
                "summary TSV differs",
            )?;
            ensure(
                fs::read_to_string(study().join("class_counts.svg"))? == svg,
                "graph differs",
            )?;
            println!("CHECKED\t{actual}");
        }
        "summary" => {
            let db = Connection::open_with_flags(
                &final_path,
                rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY,
            )?;
            integrity(&db)?;
            let (tsv, _) = products(&db)?;
            print!("{tsv}");
        }
        _ => return Err(bad("usage: weak-curve-database [build|check|summary]")),
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_domain_rejects_missing_trace() {
        let p = 7;
        let bound = 2 * p * p * p;
        let mut all: BTreeSet<_> = (-bound..=bound).filter(|&t| is_ordinary(p, t)).collect();
        assert!(verify_exact_domain(p, &all).is_ok());
        all.remove(&22);
        assert!(verify_exact_domain(p, &all).is_err());
    }
    #[test]
    fn diagnostic_domain_keeps_boundary_and_rejects_missing_row() {
        let p: i64 = 47;
        let bound = 2 * p * p * p;
        let mut all: BTreeSet<_> = (-bound..=bound).filter(|&t| t.rem_euclid(4) == 2).collect();
        assert!(verify_diagnostic_domain(p, &all).is_ok());
        assert!(!is_ordinary(p, bound));
        all.remove(&bound);
        assert!(verify_diagnostic_domain(p, &all).is_err());
    }
    #[test]
    fn boolean_parser_is_strict() {
        assert_eq!(bool01("0").unwrap(), 0);
        assert_eq!(bool01("1").unwrap(), 1);
        assert!(bool01("2").is_err());
        assert!(bool01("true").is_err());
    }
}
