//! Native, public-data-only factor-base design, export and cross-backend replay.
//! This program does not consume external targets or choose a runtime winner.
use crypto_lib::binary_ecc::curve::{point_add, point_neg, scalar_mul};
use crypto_lib::binary_ecc::{BinaryCurve, BinaryPoint, F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::koblitz_fast_arith::{FastBinaryCurve128, FastPoint128};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    frobenius_eigenvalue, koblitz_point_count,
};
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::{BTreeMap, HashSet};
use std::fs::{self, OpenOptions};
use std::io::{Read, Write};
use std::path::Path;
use std::process::Command;
use std::time::{Duration, Instant};

#[path = "../research/koblitz_n83_factor_base_sweep_20261008/compact_cold.rs"]
mod compact_cold;
#[path = "../research/koblitz_n83_factor_base_sweep_20261008/primary_adapter.rs"]
mod primary_adapter;

const STUDY: &str = "koblitz_n83_factor_base_sweep_20261008";
const PREFIX: &str = "s3://crypto-autoresearcher/factor-bases/icv1/etc";
const SEEDS: [u64; 3] = [2026100801, 2026100802, 2026100803];
const SIZES: [usize; 6] = [32, 64, 128, 256, 600, 900];
const DIMS: [usize; 8] = [4, 8, 12, 16, 20, 28, 41, 82];
const BACKENDS: [&str; 8] = [
    "compact_s3_four_sum",
    "mitm",
    "f4",
    "sat_cdcl",
    "wdsat",
    "fes_gray",
    "fes_moebius",
    "fes_monica",
];
const SYMMETRIES: [&str; 5] = [
    "none",
    "ordered_summands",
    "signed_frobenius",
    "two_torsion",
    "four_torsion",
];
const LA: [&str; 3] = ["dense_u64", "sparse_u64", "modular_biguint"];
const RADICES: [usize; 8] = [5, 4, 5, 8, 2, 3, 3, 3];

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
struct BaseSpec {
    a: u8,
    family: String,
    parameter: usize,
    parameter_unit: String,
    seed: u64,
    closure: String,
}

fn base_specs() -> Vec<BaseSpec> {
    let mut out = Vec::new();
    for a in [0, 1] {
        for family in [
            "public_x_sequential",
            "public_x_hash",
            "public_x_gray_prefix",
            "polynomial_subspace",
            "random_linear_subspace",
            "affine_subspace",
            "two_torsion_u_subspace",
            "four_torsion_v_subspace",
            "frobenius_divisor",
        ] {
            let orbit = family.starts_with("public_x");
            let parameters: Vec<usize> = if orbit {
                SIZES.to_vec()
            } else if family == "frobenius_divisor" {
                vec![1, 82, 83]
            } else {
                DIMS.to_vec()
            };
            let seeds: Vec<u64> = if family == "frobenius_divisor" {
                vec![0]
            } else {
                SEEDS.to_vec()
            };
            let closures: Vec<&str> = if family == "frobenius_divisor" {
                vec!["signed_frobenius"]
            } else {
                vec!["none", "sign", "signed_frobenius"]
            };
            for parameter in parameters {
                for &seed in &seeds {
                    for &closure in &closures {
                        out.push(BaseSpec {
                            a,
                            family: family.into(),
                            parameter,
                            parameter_unit: if orbit {
                                "orbit_columns"
                            } else {
                                "binary_dimension"
                            }
                            .into(),
                            seed,
                            closure: closure.into(),
                        });
                    }
                }
            }
        }
    }
    out
}

fn solver_count() -> usize {
    RADICES.iter().product()
}

fn decode(mut ordinal: usize) -> [usize; 8] {
    let mut digits = [0; 8];
    for (i, radix) in RADICES.iter().enumerate().rev() {
        digits[i] = ordinal % radix;
        ordinal /= radix;
    }
    digits
}

// Every tuple gets exactly one primary disposition. Additional prerequisites
// remain in protocol.json; a disposition never means the backend was executed.
fn disposition(base: &BaseSpec, d: [usize; 8]) -> &'static str {
    let m = d[0] + 2;
    let backend = BACKENDS[d[3]];
    if base.parameter_unit == "binary_dimension" && base.parameter > 20 {
        return "resource_gate_dimension_above_20";
    }
    if base.family == "frobenius_divisor" {
        return "dimension_one_torsion_only";
    }
    if SYMMETRIES[d[2]] == "signed_frobenius" && base.closure != "signed_frobenius" {
        return "incompatible_symmetry_and_domain";
    }
    if backend == "compact_s3_four_sum"
        && (m != 4 || !base.family.starts_with("public_x") || base.closure != "signed_frobenius")
    {
        return "incompatible_compact_four_sum_domain";
    }
    if d[1] != 0 && !["f4", "wdsat"].contains(&backend) {
        return "split_adapter_not_implemented_for_backend";
    }
    if d[5] != 0 {
        return "large_prime_adapter_is_single_word";
    }
    if backend.starts_with("fes_") {
        return "wide_anf_filter_and_full_verifier_required";
    }
    if base.a == 0 && LA[d[7]] != "modular_biguint" {
        return "subgroup_order_exceeds_u64";
    }
    "unexecuted_recipe_requires_n83_pipeline_adapter"
}

fn case(ordinal: usize) -> Result<Value> {
    let bases = base_specs();
    let total = bases.len() * solver_count();
    if ordinal >= total {
        return Err(format!("ordinal {ordinal} outside 0..{total}").into());
    }
    let base = &bases[ordinal / solver_count()];
    let d = decode(ordinal % solver_count());
    let split_bits = [0, 4, 8, 12][d[1]];
    let enumeration = ["binary", "gray"][d[4]];
    let threads = [1, 4, 12][d[6]];
    Ok(
        json!({"ordinal":ordinal,"base":base,"summands":d[0]+2,"split_bits":split_bits,"symmetry":SYMMETRIES[d[2]],"backend":BACKENDS[d[3]],"enumeration":enumeration,"max_large_primes":d[5],"threads":threads,"linear_algebra":LA[d[7]],"disposition":disposition(base,d),"executed":false}),
    )
}

fn write_new_json(path: &Path, value: &Value) -> Result<()> {
    let mut file = OpenOptions::new().write(true).create_new(true).open(path)?;
    serde_json::to_writer_pretty(&mut file, value)?;
    file.write_all(b"\n")?;
    Ok(())
}

fn design() -> Value {
    let bases = base_specs();
    let mut counts = BTreeMap::<&str, usize>::new();
    for base in &bases {
        for ordinal in 0..solver_count() {
            *counts
                .entry(disposition(base, decode(ordinal)))
                .or_default() += 1;
        }
    }
    json!({"schema":"n83.factor-base-design/v1","study":STUDY,"curve_models":["icv1-f2m83-tm6151469093347-debefd74","icv1-f2m83-t6151469093347-cdcc5432"],"primary_curve_a":0,"modulus_low_terms":[0,1,2,45],"base_specs":bases,"solver_axes":{"summands":[2,3,4,5,6],"split_bits":[0,4,8,12],"symmetry":SYMMETRIES,"backend":BACKENDS,"enumeration":["binary","gray"],"max_large_primes":[0,1,2],"threads":[1,4,12],"linear_algebra":LA},"solver_radices":RADICES,"base_spec_count":bases.len(),"solver_tuple_count":solver_count(),"total_tuple_count":bases.len()*solver_count(),"disposition_counts":counts,"coverage":"complete finite Cartesian product, addressable by mixed-radix ordinal; this is design coverage, not experimental coverage","outside_grid":"all other bases, dimensions, split heuristics, polynomial bases, solver versions, restart schedules and parameter values remain unsearched","objective":{"primary":"fully charged single-target cold index-calculus runtime","exclusive_phases":["curve_setup","factor_base_construction","orbit_and_index_construction","relation_generation_and_encoding","relation_solving_and_lifting","partial_relation_graph_and_filtering","linear_algebra","individual_log","verification","artifact_io"],"timeout_handling":"retain failure class and cap; unknown completion time, never rank as a win","selected_winner":null},"pilot_budget_seconds":3600,"s3_prefix":PREFIX,"runtime_claim_gate":"matched complete baseline/candidate runs, A/A controls, five interleaved rounds, disjoint holdouts, L2 isolation, paired 95% CI, independent replay; no total-time winner from stage diagnostics"})
}

struct PinnedCurve {
    curve: BinaryCurve,
    subgroup_order: BigUint,
    cofactor: BigUint,
    lambda: BigUint,
}

impl PinnedCurve {
    fn generator(&self) -> &BinaryPoint {
        &self.curve.generator
    }
    fn mul(&self, p: &BinaryPoint, k: &BigUint) -> BinaryPoint {
        scalar_mul(&self.curve, p, k)
    }
    fn frobenius(&self, p: &BinaryPoint) -> BinaryPoint {
        match p {
            BinaryPoint::Infinity => BinaryPoint::Infinity,
            BinaryPoint::Affine { x, y } => BinaryPoint::Affine {
                x: x.square(&self.curve.irreducible),
                y: y.square(&self.curve.irreducible),
            },
        }
    }
}

fn curve(a: u8) -> Result<PinnedCurve> {
    let (order, cofactor, gx, gy) = match a {
        0 => (
            "2417851639230796216685689",
            "4",
            "477f77103dfad59850800",
            "2fa5e737d542c4e4fd5c3",
        ),
        1 => (
            "8569786107849059",
            "1128547018",
            "68a212cfe19a809fe0598",
            "244a245ea0b17d8cc8297",
        ),
        _ => return Err("curve a must be 0 or 1".into()),
    };
    let subgroup_order = BigUint::parse_bytes(order.as_bytes(), 10).ok_or("order")?;
    let cofactor = BigUint::parse_bytes(cofactor.as_bytes(), 10).ok_or("cofactor")?;
    let generator = BinaryPoint::Affine {
        x: F2mElement::from_hex(gx, 83),
        y: F2mElement::from_hex(gy, 83),
    };
    let curve = BinaryCurve {
        m: 83,
        irreducible: IrreduciblePoly {
            degree: 83,
            low_terms: vec![0, 1, 2, 45],
        },
        a: if a == 0 {
            F2mElement::zero(83)
        } else {
            F2mElement::one(83)
        },
        b: F2mElement::one(83),
        generator,
        order: subgroup_order.clone(),
        cofactor: cofactor.clone(),
    };
    // Rabin criterion for prime degree 83: X^(2^83)=X modulo P and
    // gcd(P, X^2+X)=1. P has constant 1 and P(1)=1, so neither X nor
    // X+1 divides it. This check uses the generic polynomial arithmetic.
    let indeterminate = F2mElement::from_hex("2", 83);
    if indeterminate.square_k_times(83, &curve.irreducible) != indeterminate {
        return Err("pinned modulus failed the prime-degree irreducibility check".into());
    }
    if koblitz_point_count(a, 83) != &subgroup_order * &cofactor
        || !curve.is_on_curve(&curve.generator)
        || scalar_mul(&curve, &curve.generator, &subgroup_order) != BinaryPoint::Infinity
    {
        return Err("pinned curve/group validation failed".into());
    }
    let lambda = frobenius_eigenvalue(&curve, if a == 0 { -1 } else { 1 }, &subgroup_order)
        .ok_or("Frobenius eigenvalue")?;
    Ok(PinnedCurve {
        curve,
        subgroup_order,
        cofactor,
        lambda,
    })
}

fn general(p: FastPoint128) -> BinaryPoint {
    match p {
        None => BinaryPoint::Infinity,
        Some((x, y)) => BinaryPoint::Affine {
            x: F2mElement::from_hex(&format!("{x:x}"), 83),
            y: F2mElement::from_hex(&format!("{y:x}"), 83),
        },
    }
}

fn words(p: &BinaryPoint) -> FastPoint128 {
    match p {
        BinaryPoint::Infinity => None,
        BinaryPoint::Affine { x, y } => Some((
            x.to_biguint().to_u128().unwrap(),
            y.to_biguint().to_u128().unwrap(),
        )),
    }
}

fn public_x(policy: &str, seed: u64, index: u64) -> u128 {
    match policy {
        "public_x_sequential" => u128::from(index + seed % 1024),
        "public_x_gray_prefix" => u128::from((index ^ (index >> 1)) ^ (seed % 1024)),
        "public_x_hash" => {
            let hash = blake3::hash(format!("{STUDY}:public-x:{seed}:{index}").as_bytes());
            u128::from_le_bytes(hash.as_bytes()[..16].try_into().unwrap()) & ((1u128 << 83) - 1)
        }
        _ => panic!("unknown frozen policy"),
    }
}

fn canonical_x(fast: &FastBinaryCurve128, mut x: u128) -> u128 {
    let mut key = x;
    for _ in 1..83 {
        x = fast.gf.sqr(x);
        key = key.min(x);
    }
    key
}

fn point_json(p: FastPoint128) -> Value {
    match p {
        None => Value::Null,
        Some((x, y)) => json!([x.to_string(), y.to_string()]),
    }
}

fn parse_point(value: &Value) -> Result<FastPoint128> {
    if value.is_null() {
        return Ok(None);
    }
    let row = value.as_array().ok_or("point must be an array")?;
    if row.len() != 2 {
        return Err("point arity".into());
    }
    let x: u128 = row[0].as_str().ok_or("x encoding")?.parse()?;
    let y: u128 = row[1].as_str().ok_or("y encoding")?.parse()?;
    if x >= 1u128 << 83 || y >= 1u128 << 83 {
        return Err("out-of-field point".into());
    }
    Ok(Some((x, y)))
}

// Independent of artifact metadata, row order and orbit representative choice.
fn point_set_hash(points: impl IntoIterator<Item = FastPoint128>) -> Result<String> {
    let mut points = points
        .into_iter()
        .map(|p| p.ok_or("identity in point set"))
        .collect::<std::result::Result<Vec<_>, _>>()?;
    points.sort_unstable();
    let mut hash = blake3::Hasher::new();
    hash.update(b"n83.sorted-polynomial-points/v1");
    for (x, y) in points {
        hash.update(&x.to_le_bytes());
        hash.update(&y.to_le_bytes());
    }
    Ok(hash.finalize().to_hex().to_string())
}

fn construct(
    a: u8,
    policy: &str,
    columns: usize,
    seed: u64,
    deadline: Instant,
) -> Result<(Vec<u8>, Value)> {
    let kc = curve(a)?;
    let fast =
        FastBinaryCurve128::new(&kc.curve.irreducible, u128::from(a)).ok_or("wide backend")?;
    let start = Instant::now();
    let mut seen = HashSet::new();
    let mut reps = Vec::new();
    let mut sources = Vec::new();
    let mut source_indices = Vec::new();
    let mut scanned = 0u64;
    while reps.len() < columns {
        if Instant::now() >= deadline {
            return Err("BUDGET_EXHAUSTED".into());
        }
        scanned += 1;
        let raw_x = public_x(policy, seed, scanned);
        let Some(source) = fast.points_with_x(1, raw_x).first().copied() else {
            continue;
        };
        let Some((x, y)) = fast.scalar_mul(source, &kc.cofactor) else {
            continue;
        };
        if x <= 1 || !seen.insert(canonical_x(&fast, x)) {
            continue;
        }
        reps.push(Some((x, y)));
        sources.push(source);
        source_indices.push(scanned);
    }
    let header = json!({"schema":"n83.factor-base/v1","study":STUDY,"curve_slug":if a==0{"icv1-f2m83-tm6151469093347-debefd74"}else{"icv1-f2m83-t6151469093347-cdcc5432"},"n":83,"curve_a":a,"curve_b":"1","modulus_low_terms":[0,1,2,45],"field_encoding":"decimal u128 polynomial word; bit i is coefficient of z^i","subgroup_order":kc.subgroup_order.to_string(),"cofactor":kc.cofactor.to_string(),"frobenius_lambda":kc.lambda.to_string(),"generator":point_json(words(kc.generator())),"policy":policy,"seed":seed,"orbit_columns":columns,"point_count":columns*166,"closure":"signed_frobenius","subspace_membership_after_projection":false,"representatives":reps.iter().copied().map(point_json).collect::<Vec<_>>(),"source_points_before_cofactor_clearing":sources.iter().copied().map(point_json).collect::<Vec<_>>(),"accepted_candidate_indices":source_indices});
    let mut bytes = serde_json::to_vec(&header)?;
    bytes.push(b'\n');
    let r = kc.subgroup_order.to_u128().ok_or("order packing")?;
    let mut point_count = 0usize;
    let mut point_set = Vec::with_capacity(columns * 166);
    for (column, rep) in reps.iter().copied().enumerate() {
        let (mut x, mut y) = rep.ok_or("identity representative")?;
        let mut coeff = BigUint::from(1u8);
        for phase in 0..83 {
            for negative in [false, true] {
                let p = Some((x, if negative { x ^ y } else { y }));
                let c = coeff.to_u128().ok_or("coefficient packing")?;
                let label = if negative { r - c } else { c };
                serde_json::to_writer(
                    &mut bytes,
                    &json!({"point":point_json(p),"column":column,"phase":phase,"negative":negative,"coefficient":label.to_string()}),
                )?;
                bytes.push(b'\n');
                point_set.push(p);
                point_count += 1;
            }
            x = fast.gf.sqr(x);
            y = fast.gf.sqr(y);
            coeff = (&coeff * &kc.lambda) % &kc.subgroup_order;
        }
        if Some((x, y)) != rep {
            return Err("Frobenius orbit did not close".into());
        }
    }
    let row = json!({"a":a,"policy":policy,"columns":columns,"seed":seed,"raw_x_candidates":scanned,"points":point_count,"point_set_blake3":point_set_hash(point_set)?,"construction_elapsed_ms_informational_L0":start.elapsed().as_secs_f64()*1000.0,"runtime_rank_eligible":false,"total_index_calculus_runtime_ms":null,"relation_yield":null,"matrix_rank":null,"verified_dlp":null});
    Ok((bytes, row))
}

fn pilot(root: &Path, budget: u64) -> Result<()> {
    if budget == 0 || budget > 3600 {
        return Err("pilot budget must be 1..3600 seconds".into());
    }
    fs::create_dir(root)?;
    fs::create_dir(root.join("objects"))?;
    let start = Instant::now();
    let deadline = start + Duration::from_secs(budget);
    let mut rows = Vec::new();
    let mut stop = "completed_factor_base_panel";
    'panel: for a in [0, 1] {
        for policy in [
            "public_x_sequential",
            "public_x_hash",
            "public_x_gray_prefix",
        ] {
            for seed in SEEDS {
                for columns in [64, 256, 600] {
                    match construct(a, policy, columns, seed, deadline) {
                        Ok((bytes, mut row)) => {
                            let plain_hash = blake3::hash(&bytes).to_hex().to_string();
                            let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
                            encoder.write_all(&bytes)?;
                            let compressed = encoder.finish()?;
                            let hash = blake3::hash(&compressed).to_hex().to_string();
                            let object = format!("objects/{hash}.jsonl.gz");
                            let mut file = OpenOptions::new()
                                .write(true)
                                .create_new(true)
                                .open(root.join(&object))?;
                            file.write_all(&compressed)?;
                            row["object"] = json!(object);
                            row["compressed_blake3"] = json!(hash);
                            row["plain_blake3"] = json!(plain_hash);
                            row["bytes"] = json!(compressed.len());
                            row["s3_uri"] = json!(format!("{PREFIX}/{STUDY}/a{a}/{object}"));
                            rows.push(row);
                            println!("base a={a} policy={policy} seed={seed} K={columns} retained");
                        }
                        Err(error) if error.to_string() == "BUDGET_EXHAUSTED" => {
                            stop = "budget_exhausted";
                            break 'panel;
                        }
                        Err(error) => return Err(error),
                    }
                }
            }
        }
    }
    let git = Command::new("git").args(["rev-parse", "HEAD"]).output()?;
    write_new_json(
        &root.join("manifest.json"),
        &json!({"schema":"n83.factor-base-panel/v1","study":STUDY,"status":stop,"source_commit":String::from_utf8_lossy(&git.stdout).trim(),"source_blake3":blake3::hash(include_bytes!("koblitz_n83_factor_base_export.rs")).to_hex().to_string(),"architecture":std::env::consts::ARCH,"os":std::env::consts::OS,"budget_seconds":budget,"elapsed_seconds_informational":start.elapsed().as_secs_f64(),"expected_base_count":54,"completed_base_count":rows.len(),"bases":rows,"selected_best_total_runtime":null,"measurement_scope":"factor-base construction and exact artifacts only; no relation collection, linear algebra, individual log, or runtime winner"}),
    )?;
    Ok(())
}

fn load_object(root: &Path, row: &Value) -> Result<Vec<u8>> {
    let object = row["object"].as_str().ok_or("missing object")?;
    let hash = row["compressed_blake3"]
        .as_str()
        .ok_or("missing compressed hash")?;
    if object != format!("objects/{hash}.jsonl.gz")
        || hash.len() != 64
        || !hash.bytes().all(|b| b.is_ascii_hexdigit())
    {
        return Err("object path/hash binding".into());
    }
    let bytes = fs::read(root.join(object))?;
    if blake3::hash(&bytes).to_hex().as_str() != hash {
        return Err("compressed hash mismatch".into());
    }
    let mut decoded = Vec::new();
    GzDecoder::new(&bytes[..]).read_to_end(&mut decoded)?;
    if blake3::hash(&decoded).to_hex().as_str()
        != row["plain_blake3"].as_str().ok_or("missing plain hash")?
    {
        return Err("plain hash mismatch".into());
    }
    Ok(decoded)
}

fn replay(root: &Path) -> Result<()> {
    let manifest: Value = serde_json::from_slice(&fs::read(root.join("manifest.json"))?)?;
    let rows = manifest["bases"].as_array().ok_or("manifest bases")?;
    if manifest["schema"] != "n83.factor-base-panel/v1"
        || manifest["study"] != STUDY
        || rows.is_empty()
        || manifest["completed_base_count"].as_u64() != Some(rows.len() as u64)
        || manifest["expected_base_count"] != 54
        || match manifest["status"].as_str() {
            Some("completed_factor_base_panel") => rows.len() != 54,
            Some("budget_exhausted") => rows.len() >= 54,
            _ => true,
        }
    {
        return Err("incomplete or invalid panel manifest".into());
    }
    let mut checks = Vec::new();
    let mut objects = HashSet::new();
    for row in rows {
        if !objects.insert(row["object"].as_str().ok_or("object")?) {
            return Err("duplicate manifest object".into());
        }
        let bytes = load_object(root, row)?;
        let text = std::str::from_utf8(&bytes)?;
        let mut lines = text.lines();
        let header: Value = serde_json::from_str(lines.next().ok_or("header")?)?;
        let a: u8 = header["curve_a"].as_u64().ok_or("curve a")?.try_into()?;
        let kc = curve(a)?;
        let slug = if a == 0 {
            "icv1-f2m83-tm6151469093347-debefd74"
        } else {
            "icv1-f2m83-t6151469093347-cdcc5432"
        };
        if header["schema"] != "n83.factor-base/v1"
            || header["study"] != STUDY
            || header["n"] != 83
            || header["curve_b"] != "1"
            || header["curve_slug"] != slug
            || header["curve_a"] != row["a"]
            || header["policy"] != row["policy"]
            || header["seed"] != row["seed"]
            || header["orbit_columns"] != row["columns"]
            || header["closure"] != "signed_frobenius"
            || header["modulus_low_terms"] != json!([0, 1, 2, 45])
            || header["subgroup_order"] != kc.subgroup_order.to_string()
            || header["cofactor"] != kc.cofactor.to_string()
            || header["frobenius_lambda"] != kc.lambda.to_string()
            || header["generator"] != point_json(words(kc.generator()))
        {
            return Err("curve/subgroup binding mismatch".into());
        }
        let reps: Vec<_> = header["representatives"]
            .as_array()
            .ok_or("reps")?
            .iter()
            .map(parse_point)
            .collect::<Result<_>>()?;
        let sources: Vec<_> = header["source_points_before_cofactor_clearing"]
            .as_array()
            .ok_or("sources")?
            .iter()
            .map(parse_point)
            .collect::<Result<_>>()?;
        let indices: Vec<_> = header["accepted_candidate_indices"]
            .as_array()
            .ok_or("candidate indices")?
            .iter()
            .map(|v| v.as_u64().ok_or("candidate index"))
            .collect::<std::result::Result<_, _>>()?;
        let policy = header["policy"].as_str().ok_or("policy")?;
        if ![
            "public_x_sequential",
            "public_x_hash",
            "public_x_gray_prefix",
        ]
        .contains(&policy)
            || indices.first().copied().unwrap_or(0) == 0
            || indices.windows(2).any(|w| w[0] >= w[1])
            || indices.last().copied() != row["raw_x_candidates"].as_u64()
        {
            return Err("public candidate policy/index binding".into());
        }
        if reps.len() != sources.len()
            || reps.len() != indices.len()
            || reps.len() != row["columns"].as_u64().ok_or("column count")? as usize
        {
            return Err("representative count".into());
        }
        // Separate implementation: generic multi-limb arithmetic, not the
        // producer's Gf2_128 kernel. Same host and repository are disclosed.
        for ((rep, source), index) in reps.iter().zip(&sources).zip(&indices) {
            let gp = general(*rep);
            let sp = general(*source);
            if *rep == None
                || source.map(|p| p.0)
                    != Some(public_x(
                        policy,
                        header["seed"].as_u64().ok_or("seed")?,
                        *index,
                    ))
                || !kc.curve.is_on_curve(&gp)
                || !kc.curve.is_on_curve(&sp)
                || kc.mul(&sp, &kc.cofactor) != gp
                || kc.mul(&gp, &kc.subgroup_order) != BinaryPoint::Infinity
                || kc.mul(&gp, &kc.lambda) != kc.frobenius(&gp)
            {
                return Err("reference subgroup/projection/Frobenius check".into());
            }
        }
        let mut seen = HashSet::new();
        let r = kc.subgroup_order.to_u128().ok_or("r")?;
        let mut count = 0usize;
        for rep in &reps {
            let mut expected = general(*rep);
            let mut coefficient = BigUint::from(1u8);
            for phase in 0..83 {
                for negative in [false, true] {
                    let line = lines.next().ok_or("missing point row")?;
                    let entry: Value = serde_json::from_str(line)?;
                    let point = parse_point(&entry["point"])?;
                    let want = if negative {
                        point_neg(&expected)
                    } else {
                        expected.clone()
                    };
                    let c = coefficient.to_u128().ok_or("coefficient")?;
                    if general(point) != want
                        || !kc.curve.is_on_curve(&want)
                        || !seen.insert(point)
                        || entry["column"].as_u64() != Some((count / 166) as u64)
                        || entry["phase"].as_u64() != Some(phase)
                        || entry["negative"].as_bool() != Some(negative)
                        || entry["coefficient"]
                            != if negative {
                                (r - c).to_string()
                            } else {
                                c.to_string()
                            }
                    {
                        return Err("point, orbit, duplicate or label check".into());
                    }
                    count += 1;
                }
                expected = kc.frobenius(&expected);
                coefficient = (&coefficient * &kc.lambda) % &kc.subgroup_order;
            }
            if expected != general(*rep) {
                return Err("orbit closure".into());
            }
        }
        if lines.next().is_some()
            || count != row["points"].as_u64().ok_or("point count")? as usize
            || header["point_count"].as_u64() != Some(count as u64)
        {
            return Err("point count/trailing data".into());
        }
        if row["point_set_blake3"] != point_set_hash(seen)? {
            return Err("point-set hash mismatch".into());
        }
        checks.push(json!({"object":row["object"],"points_checked":count,"representatives_checked":reps.len(),"status":"PASS"}));
        println!("replayed {}", row["object"].as_str().unwrap());
    }
    write_new_json(
        &root.join("replay.json"),
        &json!({"schema":"n83.factor-base-replay/v1","status":"PASS","panel_manifest_blake3":blake3::hash(&fs::read(root.join("manifest.json"))?).to_hex().to_string(),"backend":"generic multi-limb BinaryCurve versus producer Gf2_128","independence":"same host and repository; external independent replay remains pending","checks":checks,"best_total_runtime":null}),
    )?;
    Ok(())
}

struct TwoSumOutcome {
    relation: Option<(FastPoint128, FastPoint128)>,
    lookups: usize,
    censored: bool,
}

// The oracle receives points only. Fixture scalars stay in a validation
// sidecar and do not enter the lookup procedure.
fn two_sum_probe(
    fast: &FastBinaryCurve128,
    points: &[FastPoint128],
    set: &HashSet<FastPoint128>,
    target: FastPoint128,
    deadline: Instant,
) -> TwoSumOutcome {
    for (i, point) in points.iter().enumerate() {
        if i % 1024 == 0 && Instant::now() >= deadline {
            return TwoSumOutcome {
                relation: None,
                lookups: i,
                censored: true,
            };
        }
        let complement = fast.add(target, FastBinaryCurve128::neg(*point));
        if set.contains(&complement) {
            return TwoSumOutcome {
                relation: Some((*point, complement)),
                lookups: i + 1,
                censored: false,
            };
        }
    }
    TwoSumOutcome {
        relation: None,
        lookups: points.len(),
        censored: false,
    }
}

fn probes(root: &Path, budget: u64) -> Result<()> {
    if budget == 0 || budget > 3600 {
        return Err("probe budget must be 1..3600 seconds".into());
    }
    let start = Instant::now();
    let deadline = start + Duration::from_secs(budget);
    let manifest_bytes = fs::read(root.join("manifest.json"))?;
    let manifest_hash = blake3::hash(&manifest_bytes).to_hex().to_string();
    let manifest: Value = serde_json::from_slice(&manifest_bytes)?;
    let receipt: Value = serde_json::from_slice(&fs::read(root.join("replay.json"))?)?;
    if receipt["status"] != "PASS" || receipt["panel_manifest_blake3"] != manifest_hash {
        return Err("matching replay required for relation probes".into());
    }
    let mut targets = Vec::new();
    let mut public_rows = Vec::new();
    let mut validation_rows = Vec::new();
    for a in [0, 1] {
        let kc = curve(a)?;
        let mut arm = Vec::new();
        for i in 0..32 {
            let digest =
                blake3::hash(format!("{STUDY}:independent-public-probe:a{a}:{i}").as_bytes());
            let scalar = BigUint::from_bytes_le(digest.as_bytes())
                % (&kc.subgroup_order - BigUint::from(1u8))
                + BigUint::from(1u8);
            let target = kc.mul(kc.generator(), &scalar);
            public_rows.push(json!({"a":a,"fixture":i,"point":point_json(words(&target))}));
            validation_rows
                .push(json!({"a":a,"fixture":i,"known_answer_scalar":scalar.to_string()}));
            arm.push(target);
        }
        targets.push(arm);
    }
    let fixture_setup_ms = start.elapsed().as_secs_f64() * 1000.0;
    let corpus = json!({"schema":"n83.public-probe-corpus/v1","study":STUDY,"inputs":"public synthetic fixtures only","targets":public_rows});
    let corpus_hash = blake3::hash(&serde_json::to_vec(&corpus)?)
        .to_hex()
        .to_string();
    write_new_json(&root.join("probe-corpus.json"), &corpus)?;
    write_new_json(
        &root.join("probe-validation.json"),
        &json!({"schema":"n83.public-probe-validation/v1","corpus_canonical_json_blake3":corpus_hash,"oracle_receives_known_answer_scalars":false,"validation":validation_rows}),
    )?;
    let mut results = Vec::new();
    let mut status = "completed_relation_stage_panel";
    for base in manifest["bases"].as_array().ok_or("bases")? {
        if Instant::now() >= deadline {
            status = "budget_exhausted";
            break;
        }
        let preparation_start = Instant::now();
        let decoded = load_object(root, base)?;
        let mut lines = std::str::from_utf8(&decoded)?.lines();
        lines.next().ok_or("header")?;
        let points: Vec<_> = lines
            .map(|line| -> Result<_> {
                let entry: Value = serde_json::from_str(line)?;
                parse_point(&entry["point"])
            })
            .collect::<Result<_>>()?;
        let set: HashSet<_> = points.iter().copied().collect();
        if set.len() != points.len()
            || points.len() != base["points"].as_u64().ok_or("points")? as usize
        {
            return Err("probe point count/distinctness".into());
        }
        let a: u8 = base["a"].as_u64().ok_or("a")?.try_into()?;
        let kc = curve(a)?;
        let fast =
            FastBinaryCurve128::new(&kc.curve.irreducible, u128::from(a)).ok_or("wide backend")?;
        let preparation_ms = preparation_start.elapsed().as_secs_f64() * 1000.0;
        let mut outcomes = Vec::new();
        let mut total_lookups = 0usize;
        let mut relations = 0usize;
        let mut solve_ms = 0.0;
        let mut verify_ms = 0.0;
        for (i, target) in targets
            .get(a as usize)
            .ok_or("target arm")?
            .iter()
            .enumerate()
        {
            let solve_start = Instant::now();
            let outcome = two_sum_probe(&fast, &points, &set, words(target), deadline);
            solve_ms += solve_start.elapsed().as_secs_f64() * 1000.0;
            total_lookups += outcome.lookups;
            let verify_start = Instant::now();
            let relation = if let Some((p, q)) = outcome.relation {
                if point_add(&kc.curve, &general(p), &general(q)) != *target {
                    return Err("reference relation group verification".into());
                }
                relations += 1;
                json!([point_json(p), point_json(q)])
            } else {
                Value::Null
            };
            verify_ms += verify_start.elapsed().as_secs_f64() * 1000.0;
            outcomes.push(json!({"fixture":i,"complement_lookups":outcome.lookups,"relation":relation,"status":if outcome.censored{"UNKNOWN_budget"}else if outcome.relation.is_some(){"verified_relation"}else{"no_two_sum_in_complete_base"}}));
            if outcome.censored {
                status = "budget_exhausted";
                break;
            }
        }
        results.push(json!({"object":base["object"],"point_set_blake3":base["point_set_blake3"],"a":a,"policy":base["policy"],"seed":base["seed"],"columns":base["columns"],"summands":2,"backend":"exact_point_complement_lookup","max_large_primes":0,"threads":1,"targets_attempted":outcomes.len(),"verified_relations":relations,"complement_lookups":total_lookups,"preparation_ms_informational_L0":preparation_ms,"solve_ms_informational_L0":solve_ms,"verification_ms_informational_L0":verify_ms,"outcomes":outcomes,"matrix_rank":null,"total_index_calculus_runtime_ms":null,"runtime_rank_eligible":false}));
        println!(
            "probed a={a} K={} policy={} relations={relations} lookups={total_lookups}",
            base["columns"], base["policy"]
        );
        if status == "budget_exhausted" {
            break;
        }
    }
    let git = Command::new("git").args(["rev-parse", "HEAD"]).output()?;
    write_new_json(
        &root.join("probes.json"),
        &json!({"schema":"n83.relation-stage-probes/v1","study":STUDY,"status":status,"source_commit":String::from_utf8_lossy(&git.stdout).trim(),"source_blake3":blake3::hash(include_bytes!("koblitz_n83_factor_base_export.rs")).to_hex().to_string(),"panel_manifest_blake3":manifest_hash,"corpus_canonical_json_blake3":corpus_hash,"shared_fixture_setup_ms_informational_L0":fixture_setup_ms,"budget_seconds":budget,"elapsed_seconds_informational_L0":start.elapsed().as_secs_f64(),"completed_base_rows":results.len(),"results":results,"inference_scope":"fixed public two-summand relation-stage diagnostics; shared corpora and duplicate bases are dependent observations; no higher-arity or total-runtime inference","selected_best_total_runtime":null}),
    )?;
    Ok(())
}

fn upload(root: &Path) -> Result<()> {
    let manifest: Value = serde_json::from_slice(&fs::read(root.join("manifest.json"))?)?;
    let replay: Value = serde_json::from_slice(&fs::read(root.join("replay.json"))?)?;
    if replay["status"] != "PASS"
        || replay["panel_manifest_blake3"]
            != blake3::hash(&fs::read(root.join("manifest.json"))?)
                .to_hex()
                .to_string()
    {
        return Err("matching replay receipt required".into());
    }
    let mut receipts = Vec::new();
    for row in manifest["bases"].as_array().ok_or("bases")? {
        load_object(root, row)?;
        let object = row["object"].as_str().ok_or("object")?;
        let uri = row["s3_uri"].as_str().ok_or("uri")?;
        let expected = format!(
            "{PREFIX}/{STUDY}/a{}/{object}",
            row["a"].as_u64().ok_or("a")?
        );
        if uri != expected {
            return Err("storage destination binding".into());
        }
        let status = Command::new("aws")
            .args(["s3", "cp"])
            .arg(root.join(object))
            .arg(uri)
            .args(["--only-show-errors", "--metadata"])
            .arg(format!(
                "blake3={}",
                row["compressed_blake3"].as_str().unwrap()
            ))
            .status()?;
        if !status.success() {
            return Err(format!("upload failed: {uri}").into());
        }
        let download = root.join("download-check.jsonl.gz");
        let status = Command::new("aws")
            .args(["s3", "cp", uri])
            .arg(&download)
            .arg("--only-show-errors")
            .status()?;
        if !status.success() {
            return Err(format!("download replay failed: {uri}").into());
        }
        let got = blake3::hash(&fs::read(&download)?).to_hex().to_string();
        fs::remove_file(&download)?;
        if got != row["compressed_blake3"].as_str().ok_or("hash")? {
            return Err("S3 round-trip hash mismatch".into());
        }
        receipts.push(json!({"s3_uri":uri,"compressed_blake3":got,"downloaded_hash_matches":true}));
        println!("uploaded and downloaded {uri}");
    }
    write_new_json(
        &root.join("upload-receipt.json"),
        &json!({"schema":"n83.factor-base-storage/v1","status":"PASS","bucket":"crypto-autoresearcher","prefix":PREFIX,"objects":receipts,"manifest_blake3":blake3::hash(&fs::read(root.join("manifest.json"))?).to_hex().to_string()}),
    )?;
    for name in [
        "manifest.json",
        "replay.json",
        "upload-receipt.json",
        "duplicate-point-sets.json",
        "probes.json",
        "probe-corpus.json",
        "probe-validation.json",
    ] {
        if !root.join(name).exists() {
            continue;
        }
        let uri = format!(
            "{PREFIX}/{STUDY}/panels/{}/{name}",
            blake3::hash(&fs::read(root.join("manifest.json"))?).to_hex()
        );
        if !Command::new("aws")
            .args(["s3", "cp"])
            .arg(root.join(name))
            .arg(&uri)
            .arg("--only-show-errors")
            .status()?
            .success()
        {
            return Err(format!("receipt upload failed: {uri}").into());
        }
    }
    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("plan")=>{ let value=design(); write_new_json(Path::new(args.get(2).ok_or("plan output.json")?),&value)?; println!("{}",json!({"base_specs":value["base_spec_count"],"tuples":value["total_tuple_count"],"dispositions":value["disposition_counts"]})); },
        Some("case")=>println!("{}",serde_json::to_string_pretty(&case(args.get(2).ok_or("case ordinal")?.parse()?)?)?),
        Some("pilot")=>pilot(Path::new(args.get(2).ok_or("pilot new-directory budget-seconds")?),args.get(3).ok_or("budget")?.parse()?)?,
        Some("replay")=>replay(Path::new(args.get(2).ok_or("replay directory")?))?,
        Some("probes")=>probes(Path::new(args.get(2).ok_or("probes directory budget_seconds")?),args.get(3).ok_or("budget")?.parse()?)?,
        Some("upload")=>upload(Path::new(args.get(2).ok_or("upload directory")?))?,
        Some("cold")=>compact_cold::run_cli(std::iter::once(args[0].clone()).chain(args[2..].iter().cloned()).collect()),
        Some("primary-adapter-check")=>{
            let panel=Path::new(args.get(2).ok_or("primary-adapter-check PANEL_DIR K NEW_OUTPUT_JSON")?);
            let columns=args.get(3).ok_or("primary orbit columns")?.parse()?;
            let output=Path::new(args.get(4).ok_or("new output JSON path")?);
            let receipt=primary_adapter::check_panel(panel,columns)?;
            write_new_json(output,&receipt)?;
            println!("{receipt}");
        },
        Some("primary-cold")=>{
            if args.len()!=9 { return Err("primary-cold PANEL_DIR K M STRATEGY MAX_TRIALS BUDGET_SECONDS NEW_RUN_DIR".into()); }
            let summary=primary_adapter::run_cli(
                Path::new(&args[2]),args[3].parse()?,args[4].parse()?,&args[5],
                args[6].parse()?,args[7].parse()?,Path::new(&args[8]),
            )?;
            println!("{summary}");
        },
        _=>return Err("usage: plan output.json | case ordinal | pilot NEW_DIRECTORY budget_seconds | replay directory | probes directory budget_seconds | upload directory | cold PANEL_DIR K MAX_SECONDS [unordered] | primary-adapter-check PANEL_DIR K NEW_OUTPUT_JSON | primary-cold PANEL_DIR K M STRATEGY MAX_TRIALS BUDGET_SECONDS NEW_RUN_DIR".into()),
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn finite_grid_is_addressable_and_accounts_for_every_tuple() {
        let bases = base_specs();
        assert_eq!(bases.len(), 1050);
        assert_eq!(solver_count(), 43200);
        let total = bases.len() * solver_count();
        assert_eq!(total, 45360000);
        assert!(case(total).is_err());
        assert_eq!(case(total - 1).unwrap()["ordinal"], total - 1);
        let mut counts = BTreeMap::new();
        for b in &bases {
            for o in 0..solver_count() {
                *counts.entry(disposition(b, decode(o))).or_insert(0usize) += 1;
            }
        }
        assert_eq!(counts.values().sum::<usize>(), total);
    }
    #[test]
    fn gray_order_visits_the_same_complete_subspace() {
        let binary: HashSet<_> = (0..4096)
            .map(|i| public_x("public_x_sequential", 0, i))
            .collect();
        let gray: HashSet<_> = (0..4096)
            .map(|i| public_x("public_x_gray_prefix", 0, i))
            .collect();
        assert_eq!(binary, gray);
    }
    #[test]
    fn invalid_point_and_out_of_range_case_fail_closed() {
        assert!(parse_point(&json!(["9671406556917033397649408", "0"])).is_err());
        assert!(parse_point(&json!(["1"])).is_err());
        assert!(case(usize::MAX).is_err());
    }
    #[test]
    fn degree_83_cyclotomic_block_order_is_82() {
        let mut v = 1usize;
        let mut first = 0;
        for i in 1..=82 {
            v = v * 2 % 83;
            if v == 1 {
                first = i;
                break;
            }
        }
        assert_eq!(first, 82);
    }
    #[test]
    fn subgroup_moduli_are_kept_separate() {
        assert!(curve(0).unwrap().subgroup_order.to_u64().is_none());
        assert!(curve(1).unwrap().subgroup_order.to_u64().is_some());
    }
    #[test]
    fn point_only_oracle_finds_a_verified_relation_and_retains_censoring() {
        let kc = curve(0).unwrap();
        let fast = FastBinaryCurve128::new(&kc.curve.irreducible, 0).unwrap();
        let p = words(kc.generator());
        let q = fast.scalar_mul(p, &BigUint::from(7u8));
        let target = fast.add(p, q);
        let points = vec![p, q];
        let set = points.iter().copied().collect();
        let result = two_sum_probe(
            &fast,
            &points,
            &set,
            target,
            Instant::now() + Duration::from_secs(10),
        );
        let (left, right) = result.relation.unwrap();
        assert_eq!(
            point_add(&kc.curve, &general(left), &general(right)),
            general(target)
        );
        assert!(!result.censored);
        let absent_target = fast.scalar_mul(p, &BigUint::from(9u8));
        let absent = two_sum_probe(
            &fast,
            &points,
            &set,
            absent_target,
            Instant::now() + Duration::from_secs(10),
        );
        assert!(absent.relation.is_none());
        assert_eq!(absent.lookups, points.len());
        assert!(!absent.censored);
        for left in &points {
            for right in &points {
                assert_ne!(
                    point_add(&kc.curve, &general(*left), &general(*right)),
                    general(absent_target)
                );
            }
        }
        let capped = two_sum_probe(&fast, &points, &set, target, Instant::now());
        assert!(capped.censored);
        assert_eq!(capped.lookups, 0);
        assert!(capped.relation.is_none());
    }
    #[test]
    fn small_export_replays_with_generic_arithmetic() {
        let dir = std::env::temp_dir().join(format!("n83-fb-test-{}", std::process::id()));
        fs::create_dir(&dir).unwrap();
        fs::create_dir(dir.join("objects")).unwrap();
        let (bytes, mut row) = construct(
            0,
            "public_x_hash",
            2,
            17,
            Instant::now() + Duration::from_secs(60),
        )
        .unwrap();
        let mut gz = GzEncoder::new(Vec::new(), Compression::default());
        gz.write_all(&bytes).unwrap();
        let compressed = gz.finish().unwrap();
        let hash = blake3::hash(&compressed).to_hex().to_string();
        let object = format!("objects/{hash}.jsonl.gz");
        row["object"] = json!(object);
        row["compressed_blake3"] = json!(hash);
        row["plain_blake3"] = json!(blake3::hash(&bytes).to_hex().to_string());
        fs::write(dir.join(&object), compressed).unwrap();
        write_new_json(&dir.join("manifest.json"), &json!({"schema":"n83.factor-base-panel/v1","study":STUDY,"status":"budget_exhausted","completed_base_count":1,"expected_base_count":54,"bases":[row.clone()]})).unwrap();
        replay(&dir).unwrap();
        // Rehash a semantically invalid coefficient: byte hashes alone must
        // not allow it through the separate arithmetic verifier.
        let text = std::str::from_utf8(&bytes).unwrap();
        let mut entries: Vec<Value> = text
            .lines()
            .map(|line| serde_json::from_str(line).unwrap())
            .collect();
        entries[1]["coefficient"] = json!("0");
        let mut bad_bytes = Vec::new();
        for entry in entries {
            serde_json::to_writer(&mut bad_bytes, &entry).unwrap();
            bad_bytes.push(b'\n');
        }
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(&bad_bytes).unwrap();
        let bad_compressed = encoder.finish().unwrap();
        let bad_hash = blake3::hash(&bad_compressed).to_hex().to_string();
        let mut bad_row = row.clone();
        bad_row["object"] = json!(format!("objects/{bad_hash}.jsonl.gz"));
        bad_row["compressed_blake3"] = json!(bad_hash);
        bad_row["plain_blake3"] = json!(blake3::hash(&bad_bytes).to_hex().to_string());
        fs::write(
            dir.join(bad_row["object"].as_str().unwrap()),
            bad_compressed,
        )
        .unwrap();
        fs::write(dir.join("manifest.json"), serde_json::to_vec(&json!({"schema":"n83.factor-base-panel/v1","study":STUDY,"status":"budget_exhausted","completed_base_count":1,"expected_base_count":54,"bases":[bad_row]})).unwrap()).unwrap();
        assert!(replay(&dir)
            .unwrap_err()
            .to_string()
            .contains("label check"));
        fs::write(dir.join(object), b"corrupt").unwrap();
        assert!(load_object(&dir, &row).is_err());
        fs::remove_dir_all(dir).unwrap();
    }
}
