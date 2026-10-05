//! Track B's measurements (IC_TOOL_PROGRAM.md §9): `harness/bround.py`,
//! natively, for the steps' protocols and their amendments.
//!
//! The steps, each resumable, every output kept on disk:
//! - `manifest`: the host and both arms' binaries (`host.json`), and where
//!   the A/A bands come from (`aa-source.json`): the named round's if this
//!   host is that round's, the run's own otherwise.
//! - `aa`: the run's own A/A, the base against a byte-identical copy on
//!   `M1`'s 22 rows, five rounds, for a host the named round did not use.
//! - `conformance`: the conformance suite at the steps named, on the
//!   candidate with its build commit and on the base, whose failures record
//!   what the base could not do (`conformance/{candidate,base}.json`).
//! - `pin`: the candidate on all 90 suite rows, untimed, against v0's
//!   outputs from R01's profile pass (`pin/pin.json`).
//! - `translate` (B1 on): `ic check --translate` on each suite row must
//!   equal design §8's translation, built here as C009's generator
//!   (`conformance/v2/make_cases.py`) builds it, and `ic price` on it must
//!   give the row's v0 outputs (`translate/translate.json`).
//! - `timing`: the base against the candidate on `M1`'s 22 rows, five
//!   rounds ABAB, isolated.
//! - `v2timing` (B1's measurement 6): the candidate on `M1`'s rows, the v1
//!   file against its v2 translation, five rounds ABAB.
//! - `chain` (the steps' amendment of 2026-10-01): one interleave over the
//!   base and each step's arm in the queue's order, five rounds, the order
//!   reversed every other round; each step's figure is its paired ratio
//!   against the arm before it.
//! - `analyse`: the verdicts, as JSON.
//!
//! Where the script read R01's A/A bands, this reads the bands of the round
//! named by `--aa-from` (the newest baseline's, R07's on v3), by slug, or
//! the run's own when it has them.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use super::bench::{self, Arm, Bench};
use super::conformance;
use super::json::{self, obj, J};
use super::pin;
use super::rounds::{self, speed, Ctx};
use super::runs;
use super::suite::{self, Row};

pub const ROUNDS: u32 = 5;

/// The chain's arms in the queue's order: the baseline, then each step.
pub const CHAIN_ORDER: [&str; 9] = ["base", "B0", "B1", "B3", "B2", "B2b", "B7a", "B3b", "B4"];

/// A run's settings.
pub struct Run {
    pub ctx: Ctx,
    /// The conformance steps accepted so far and the one judged.
    pub steps: Vec<String>,
    /// The round whose `analysis.json` gives the A/A bands (its `aa`
    /// section, by slug or by size), and whose `runs/host.json` is its host.
    pub aa_from: PathBuf,
}

fn write(path: &Path, doc: &J) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
    }
    std::fs::write(path, json::dumps(doc, 1) + "\n").map_err(|e| format!("{}: {e}", path.display()))
}

fn arm<'a>(arms: &'a [Arm], name: &str) -> Result<&'a Arm, String> {
    arms.iter()
        .find(|a| a.name == name)
        .ok_or_else(|| format!("no {name} arm"))
}

// ── manifest and A/A ───────────────────────────────────────────────

const HOST_KEYS: [&str; 6] = [
    "cpu_model",
    "cpu_flags_relevant",
    "logical_cores",
    "memory",
    "transparent_hugepage",
    "os",
];

/// `host.json`, and `aa-source.json`: whether this host is the A/A
/// round's, and so whose bands the timing is read against.
pub fn manifest(
    r: &Run,
    b: &Bench,
    arms: &[Arm],
    root: &Path,
    commits: &[Option<String>],
) -> Result<J, String> {
    let doc = bench::host_manifest(
        &r.ctx.runs.join("host.json"),
        root,
        speed::binaries(arms, commits)?,
        &b.isolate,
        "",
    )?;
    let out = r.ctx.runs.join("aa-source.json");
    if out.exists() {
        return json::read(&out);
    }
    // A round's analysis carries its host; R01's predates that, and its
    // host is in its run tree.
    let theirs = match json::read(&r.aa_from.join("analysis.json"))?.get("host") {
        Some(h) if h.truthy() => h.clone(),
        _ => json::read(&pin::r01_runs(&r.ctx.programme)?.join("host.json"))?,
    };
    let same = HOST_KEYS.iter().all(|k| match (doc.get(k), theirs.get(k)) {
        (Some(x), Some(y)) => json::py_eq(x, y),
        _ => false,
    });
    let round = r
        .aa_from
        .file_name()
        .map(|f| f.to_string_lossy().into_owned())
        .unwrap_or_default();
    let record = obj([
        ("aa_from", J::Str(round.clone())),
        ("host_matches", J::Bool(same)),
        (
            "aa",
            J::Str(if same {
                format!("{round}'s")
            } else {
                "the run's own".into()
            }),
        ),
    ]);
    write(&out, &record)?;
    Ok(record)
}

/// The A/A bands by slug: the run's own if it ran one, else the named
/// round's.
fn aa_bands(r: &Run) -> Result<(HashMap<String, J>, String), String> {
    let own = r.ctx.runs.join("aa");
    if own.exists() {
        let mut bands = HashMap::new();
        for ((a, n), rs) in suite::by_size(&speed::m1_rows(&r.ctx)?) {
            let cold = rounds::paired_of(&own, &rs, ("A", "A2"), runs::ic_cold_ns, false)?;
            if cold.get("lo").is_some() {
                bands.insert(suite::curve_slug(a, n)?, cold);
            }
        }
        return Ok((bands, "the run's own".into()));
    }
    let doc = json::read(&r.aa_from.join("analysis.json"))?;
    let mut bands = HashMap::new();
    for row in doc.at("aa")?.as_arr().ok_or("`aa` is not a list")? {
        let slug = match (
            row.get("slug").and_then(J::as_str),
            row.get("size").and_then(J::as_str),
        ) {
            (Some(s), _) => s.to_string(),
            (None, Some(size)) => {
                // R01's `k<a>n<n>`.
                let a: u8 = size[1..2].parse().map_err(|_| format!("size `{size}`"))?;
                let n: u32 = size[3..].parse().map_err(|_| format!("size `{size}`"))?;
                suite::curve_slug(a, n)?
            }
            _ => return Err("an A/A row has neither `slug` nor `size`".into()),
        };
        bands.insert(slug, row.at("cold")?.clone());
    }
    let round = r
        .aa_from
        .file_name()
        .map(|f| f.to_string_lossy().into_owned())
        .unwrap_or_default();
    Ok((bands, format!("{round}'s")))
}

// ── conformance ────────────────────────────────────────────────────

/// The suite on the candidate, with its build commit, and on the base; each
/// report written once.
pub fn conformance(r: &Run, arms: &[Arm], cand_commit: Option<&str>) -> Result<J, String> {
    let mut kv = Vec::new();
    for (name, arm_name, commit) in [("candidate", "cand", cand_commit), ("base", "base", None)] {
        let path = r.ctx.runs.join("conformance").join(format!("{name}.json"));
        if !path.exists() {
            let ic = &arm(arms, arm_name)?.binary;
            let (doc, _) = conformance::run_steps(&r.ctx.programme, ic, &r.steps, commit)?;
            write(&path, &doc)?;
        }
        kv.push((name.to_string(), conformance_summary(&json::read(&path)?)?));
    }
    Ok(J::Obj(kv))
}

fn conformance_summary(doc: &J) -> Result<J, String> {
    let results = doc
        .at("results")?
        .as_arr()
        .ok_or("`results` is not a list")?;
    let failed: Vec<J> = results
        .iter()
        .filter(|c| !c.get("pass").is_some_and(J::truthy))
        .filter_map(|c| c.get("id").cloned())
        .collect();
    Ok(obj([
        ("passed", doc.at("passed")?.clone()),
        ("cases", J::Int(results.len() as i128)),
        ("failed", J::Arr(failed)),
    ]))
}

// ── the translation reference (design §8, as C009's generator) ────

/// `#E(GF(2^n))` for `y² + xy = x³ + ax² + 1`, by the Lucas recurrence of
/// the trace of Frobenius.
pub fn koblitz_order(a: u8, n: u32) -> Result<u128, String> {
    if !(1..=125).contains(&n) || a > 1 {
        return Err(format!("no Koblitz order for a = {a}, n = {n} here"));
    }
    let t: i128 = if a == 0 { -1 } else { 1 };
    let (mut prev, mut cur) = (2i128, t);
    for _ in 1..n {
        (prev, cur) = (cur, t * cur - 2 * prev);
    }
    Ok(((1i128 << n) + 1 - cur) as u128)
}

fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

fn powmod(mut b: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1 % m;
    b %= m;
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, b, m);
        }
        b = mulmod(b, b, m);
        e >>= 1;
    }
    r
}

/// Miller–Rabin with the first twelve primes as bases: exact below 2^64.
fn is_prime(n: u64) -> bool {
    const BASES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    if n < 2 {
        return false;
    }
    for p in BASES {
        if n.is_multiple_of(p) {
            return n == p;
        }
    }
    let (mut d, mut s) = (n - 1, 0);
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    'bases: for b in BASES {
        let mut x = powmod(b, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulmod(x, x, n);
            if x == n - 1 {
                continue 'bases;
            }
        }
        return false;
    }
    true
}

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

/// A nontrivial factor of an odd composite, by Brent's variant of
/// Pollard's rho.
fn split(n: u64) -> u64 {
    for c in 1u64.. {
        let f = |x: u64| (mulmod(x, x, n) + c) % n;
        let (mut x, mut y, mut g) = (2u64, 2u64, 1u64);
        while g == 1 {
            x = f(x);
            y = f(f(y));
            g = gcd(x.abs_diff(y), n);
        }
        if g != n {
            return g;
        }
    }
    unreachable!("some constant splits a composite")
}

fn largest_prime_factor_u64(mut m: u64) -> u64 {
    if m < 2 {
        return 1;
    }
    let mut last = 1u64;
    let mut d = 2u64;
    while d <= 1 << 20 && d.saturating_mul(d) <= m {
        while m.is_multiple_of(d) {
            m /= d;
            last = d;
        }
        d += if d == 2 { 1 } else { 2 };
    }
    if m == 1 {
        return last;
    }
    if is_prime(m) {
        return m.max(last);
    }
    let f = split(m);
    largest_prime_factor_u64(f)
        .max(largest_prime_factor_u64(m / f))
        .max(last)
}

/// The largest prime factor, as `make_cases.py`'s `largest_prime_factor`
/// defines it (its trial division is the definition; this reaches the same
/// number by trial division to 2^20, then an exact primality test and
/// Pollard's rho).  Below 2^64 only: the suite's orders are.
pub fn largest_prime_factor(m: u128) -> Result<u128, String> {
    let m = u64::try_from(m).map_err(|_| format!("{m} is past this reference's range"))?;
    Ok(largest_prime_factor_u64(m).into())
}

/// `a·b mod f` over GF(2), for `deg f ≤ 63`.
fn clmul_mod(a: u64, b: u64, f: u64, n: u32) -> u64 {
    let mut p = 0u128;
    for i in 0..64 {
        if b >> i & 1 == 1 {
            p ^= (a as u128) << i;
        }
    }
    for i in (n..128).rev() {
        if p >> i & 1 == 1 {
            p ^= (f as u128) << (i - n);
        }
    }
    p as u64
}

fn poly_gcd(mut a: u64, mut b: u64) -> u64 {
    let deg = |x: u64| 63 - x.leading_zeros();
    while b != 0 {
        while a != 0 && deg(a) >= deg(b) {
            a ^= b << (deg(a) - deg(b));
        }
        (a, b) = (b, a);
    }
    a
}

/// Rabin's test over GF(2), as `curve_id.is_irreducible_f2`.
fn is_irreducible_f2(f: u64) -> bool {
    let n = 63 - f.leading_zeros();
    if n < 1 || f & 1 == 0 {
        return n == 1;
    }
    let x_pow_2k = |k: u32| {
        let mut y = 2u64;
        for _ in 0..k {
            y = clmul_mod(y, y, f, n);
        }
        y
    };
    let x_mod_f = if n == 1 { 2 ^ f } else { 2 };
    if x_pow_2k(n) != x_mod_f {
        return false;
    }
    let primes = (2..=n).filter(|&p| n.is_multiple_of(p) && (2..p).all(|q| p % q != 0));
    for p in primes {
        if poly_gcd(f, x_pow_2k(n / p) ^ 2) != 1 {
            return false;
        }
    }
    true
}

/// The modulus `curve_id.find_irreducible_sparse` gives: the numerically
/// least `x^n + low` with `low` odd and of weight at most four that is
/// irreducible, the lows enumerated as `_sparse_lows` enumerates them.
pub fn find_irreducible_sparse(n: u32) -> Result<u64, String> {
    if n == 0 || n > 62 {
        return Err(format!("no sparse modulus of degree {n} here"));
    }
    if n == 1 {
        return Ok(0b11);
    }
    let mut lows = vec![1u64];
    for top in 1..n {
        let mut below = vec![0u64];
        for i in 1..top {
            below.push(1 << i);
            for j in i + 1..top {
                below.push(1 << i | 1 << j);
            }
        }
        below.sort_unstable();
        lows.extend(below.into_iter().map(|rest| 1 | rest | 1 << top));
    }
    lows.into_iter()
        .map(|low| 1u64 << n | low)
        .find(|&f| is_irreducible_f2(f))
        .ok_or_else(|| format!("no sparse irreducible of degree {n}"))
}

const RECIPE_KEYS: [&str; 12] = [
    "summands",
    "descent_summands",
    "collection_window",
    "collection_aim",
    "pair_table_bytes",
    "pair_table_tier",
    "solver",
    "seed",
    "max_trials",
    "linear_algebra",
    "collection",
    "factor_base",
];

fn subset(doc: &J, keys: &[&str]) -> J {
    J::Obj(
        keys.iter()
            .filter_map(|k| doc.get(k).map(|v| (k.to_string(), v.clone())))
            .collect(),
    )
}

/// Design §8's translation of a suite row's v1 file, as C009's generator
/// writes it (`make_cases.document`, with the Koblitz form, the search
/// rule's generator and the row's rho seed).
pub fn reference_translation(v1: &J, a: u8, n: u32, rho_seed: i128) -> Result<J, String> {
    let f = find_irreducible_sparse(n)?;
    let order = koblitz_order(a, n)?;
    let r = largest_prime_factor(order)?;
    let target = v1
        .at("targets")?
        .as_arr()
        .and_then(|t| t.first())
        .ok_or("the v1 file has no target")?;
    Ok(obj([
        ("schema_version", J::Int(2)),
        ("name", v1.at("name")?.clone()),
        (
            "field",
            obj([
                ("kind", J::Str("binary".into())),
                ("degree", J::Int(n.into())),
                ("modulus", J::Str(format!("0x{f:x}"))),
            ]),
        ),
        (
            "curve",
            obj([("form", J::Str("koblitz".into())), ("a", J::Int(a.into()))]),
        ),
        (
            "subgroup",
            obj([
                ("order", J::Str(r.to_string())),
                ("cofactor", J::Str((order / r).to_string())),
                (
                    "generator",
                    obj([("rule", J::Str("koblitz_search_v1".into()))]),
                ),
            ]),
        ),
        (
            "target",
            subset(target, &["public_hash_seed", "known_log", "random_seed"]),
        ),
        (
            "method",
            obj([
                ("solve", J::Str("paired".into())),
                (
                    "index_calculus",
                    obj([
                        ("pipeline", J::Str("kic".into())),
                        ("recipe", subset(v1, &RECIPE_KEYS)),
                    ]),
                ),
                (
                    "rho",
                    obj([
                        ("pipeline", J::Str("rho-koblitz".into())),
                        ("seed", J::Int(rho_seed)),
                    ]),
                ),
            ]),
        ),
    ]))
}

// ── pin and translate ──────────────────────────────────────────────

pub fn pin(r: &Run, arms: &[Arm]) -> Result<J, String> {
    pin::pin(&r.ctx.programme, &r.ctx.runs, &arm(arms, "cand")?.binary)
}

fn translation_path(r: &Run, row: &Row) -> PathBuf {
    r.ctx.runs.join("translate").join(&row.id).join("v2.json")
}

/// `ic check --translate` on the row's v1 file, written beside the run
/// once.
pub(crate) fn translate_row(r: &Run, cand: &Path, row: &Row) -> Result<J, String> {
    let path = translation_path(r, row);
    if !path.exists() {
        let params = row.params.as_ref().ok_or("a suite row has no parameters")?;
        let rho_seed = row.rho_seed.ok_or("a suite row has no rho seed")?;
        let out = Command::new(cand)
            .arg("check")
            .arg("--params")
            .arg(params)
            .args(["--translate", "--rho-seed", &rho_seed.to_string(), "--json"])
            .env("RAYON_NUM_THREADS", "1")
            .stderr(Stdio::null())
            .output()
            .map_err(|e| format!("{}: {e}", cand.display()))?;
        let doc = json::parse(&String::from_utf8_lossy(&out.stdout))?;
        let docs = doc
            .get("documents")
            .and_then(J::as_arr)
            .map(<[J]>::to_vec)
            .unwrap_or_default();
        let one = if docs.len() == 1 {
            docs[0].clone()
        } else {
            J::Arr(docs)
        };
        write(&path, &one)?;
    }
    json::read(&path)
}

pub fn translate(r: &Run, arms: &[Arm]) -> Result<J, String> {
    let out = r.ctx.runs.join("translate").join("translate.json");
    if out.exists() {
        return json::read(&out);
    }
    let cand = &arm(arms, "cand")?.binary;
    let r01 = pin::r01_runs(&r.ctx.programme)?;
    let mut rows = Vec::new();
    for (tier, row) in pin::suite_rows(&r.ctx.programme)? {
        let doc = translate_row(r, cand, &row)?;
        let params = row.params.as_ref().ok_or("a suite row has no parameters")?;
        let reference = reference_translation(
            &json::read(params)?,
            row.a,
            row.n,
            row.rho_seed.ok_or("a suite row has no rho seed")?,
        )?;
        let reference_path = translation_path(r, &row).with_file_name("reference.json");
        if !reference_path.exists() {
            write(&reference_path, &reference)?;
        }
        let report = translation_path(r, &row).with_file_name("v2.price.json");
        if !report.exists() {
            Command::new("taskset")
                .args(["-c", bench::CPUS])
                .arg(cand)
                .arg("price")
                .arg("--params")
                .arg(translation_path(r, &row))
                .args(["--json", "--out"])
                .arg(&report)
                .env("RAYON_NUM_THREADS", "1")
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .status()
                .map_err(|e| format!("taskset: {e}"))?;
        }
        let new = runs::load(&report);
        let old = pin::v0_report(&r01, tier, &row)?;
        let (a, b) = (pin::outputs(&new), pin::outputs(&old));
        let mut differs: Vec<&str> = a
            .iter()
            .zip(&b)
            .filter(|((_, x), (_, y))| !json::py_eq(x, y))
            .map(|((k, _), _)| *k)
            .collect();
        differs.sort_unstable();
        rows.push(obj([
            ("row", J::Str(row.id.clone())),
            (
                "same_json_as_reference",
                J::Bool(json::py_eq(&doc, &reference)),
            ),
            ("same_outputs", J::Bool(differs.is_empty())),
            (
                "differs_in",
                J::Arr(differs.into_iter().map(|k| J::Str(k.into())).collect()),
            ),
        ]));
    }
    let held = rows.iter().all(|e| {
        e.get("same_json_as_reference").is_some_and(J::truthy)
            && e.get("same_outputs").is_some_and(J::truthy)
    });
    let doc = obj([("rows", J::Arr(rows)), ("held", J::Bool(held))]);
    write(&out, &doc)?;
    Ok(doc)
}

// ── the timed steps ────────────────────────────────────────────────

pub fn timing(r: &Run, b: &Bench, arms: &[Arm]) -> Result<J, String> {
    let pair = [arm(arms, "base")?.clone(), arm(arms, "cand")?.clone()];
    b.interleave(
        &pair,
        &speed::m1_rows(&r.ctx)?,
        ROUNDS,
        &r.ctx.runs.join("timing"),
        &[],
    )?;
    Ok(J::Null)
}

/// The candidate on each `M1` row, its v1 file against its translation,
/// the order alternating by round.
pub fn v2timing(r: &Run, b: &Bench, arms: &[Arm]) -> Result<J, String> {
    let cand = arm(arms, "cand")?.binary.clone();
    let d = r.ctx.runs.join("v2timing");
    for k in 1..=ROUNDS {
        for row in speed::m1_rows(&r.ctx)? {
            translate_row(r, &cand, &row)?;
            let names = if k % 2 == 1 {
                ["v1", "v2"]
            } else {
                ["v2", "v1"]
            };
            for name in names {
                let out = d.join(name).join(&row.id).join(format!("r{k}.price.json"));
                let priced = if name == "v1" {
                    row.clone()
                } else {
                    Row {
                        params: Some(translation_path(r, &row)),
                        ..row.clone()
                    }
                };
                // A v2 document carries its own seed, and its price is
                // single-target by definition: the script's command.
                let extra: &[String] = &[];
                let rep = if name == "v1" {
                    b.price(&cand, &priced, &out, &d, extra)?
                } else {
                    price_v2(b, &cand, &priced, &out, &d)?
                };
                let status = rep.get("status").and_then(J::as_str).unwrap_or("None");
                println!("r{k} {} {name}: {status}", row.id);
            }
        }
    }
    Ok(J::Null)
}

/// `ic price --params <v2> --json --out <out>`, isolated, with the runner's
/// retries.
fn price_v2(b: &Bench, cand: &Path, row: &Row, out: &Path, log_dir: &Path) -> Result<J, String> {
    let params = row.params.as_ref().ok_or("no translation")?;
    let mut rep = J::Null;
    for k in 0..=runs::RETRIES {
        let attempt = runs::attempt(out, k);
        if !attempt.exists() {
            if let Some(parent) = attempt.parent() {
                std::fs::create_dir_all(parent)
                    .map_err(|e| format!("{}: {e}", parent.display()))?;
            }
            let cmd = vec![
                cand.to_string_lossy().into_owned(),
                "price".into(),
                "--params".into(),
                params.to_string_lossy().into_owned(),
                "--json".into(),
                "--out".into(),
                attempt.to_string_lossy().into_owned(),
            ];
            b.launch(&cmd, &attempt, log_dir)?;
        }
        rep = runs::load(&attempt);
        if runs::clean(&attempt)? && rep.get("status").and_then(J::as_str) == Some("complete") {
            return Ok(rep);
        }
    }
    Ok(rep)
}

/// The chain's arms, checked against the queue's order: the base first,
/// then the steps present, in order.
pub fn chain_arms(given: Vec<Arm>) -> Result<Vec<Arm>, String> {
    if !given.iter().any(|a| a.name == "base") {
        return Err("the chain needs a base arm".into());
    }
    if let Some(bad) = given
        .iter()
        .find(|a| !CHAIN_ORDER.contains(&a.name.as_str()))
    {
        return Err(format!(
            "{} is not a chain arm; the arms are {}",
            bad.name,
            CHAIN_ORDER.join(", ")
        ));
    }
    let mut out = Vec::new();
    for name in CHAIN_ORDER {
        let mut these = given.iter().filter(|a| a.name == name);
        if let Some(a) = these.next() {
            if these.next().is_some() {
                return Err(format!("the arm {name} is given twice"));
            }
            if !a.binary.exists() {
                return Err(format!(
                    "the {name} arm's binary does not exist: {}",
                    a.binary.display()
                ));
            }
            out.push(a.clone());
        }
    }
    Ok(out)
}

pub fn chain(
    r: &Run,
    b: &Bench,
    arms: &[Arm],
    root: &Path,
    commits: &[Option<String>],
) -> Result<J, String> {
    bench::host_manifest(
        &r.ctx.runs.join("host.json"),
        root,
        speed::binaries(arms, commits)?,
        &b.isolate,
        "",
    )?;
    b.interleave(
        arms,
        &speed::m1_rows(&r.ctx)?,
        ROUNDS,
        &r.ctx.runs.join("chain"),
        &[],
    )?;
    Ok(obj([(
        "arms",
        J::Arr(arms.iter().map(|a| J::Str(a.name.clone())).collect()),
    )]))
}

// ── the analysis ───────────────────────────────────────────────────

/// Per size, the paired cold ratio of one arm over another, with the A/A
/// band and whether the size regresses beyond it.
fn sizes(r: &Run, d: &Path, arms: (&str, &str), bands: &HashMap<String, J>) -> Result<J, String> {
    let mut out = Vec::new();
    for ((a, n), rs) in suite::by_size(&speed::m1_rows(&r.ctx)?) {
        let slug = suite::curve_slug(a, n)?;
        let mut kv = rounds::size_head(&slug, a, n, &rs);
        kv.push((
            "cold".into(),
            rounds::paired_of(d, &rs, arms, runs::ic_cold_ns, false)?,
        ));
        rounds::aa_fields(&mut kv, bands.get(&slug))?;
        out.push(J::Obj(kv));
    }
    Ok(J::Arr(out))
}

fn any_regression(sizes: &J) -> bool {
    sizes
        .as_arr()
        .unwrap_or(&[])
        .iter()
        .any(|s| s.get("regresses_beyond_aa").is_some_and(J::truthy))
}

fn timed(
    r: &Run,
    step: &str,
    arms: (&str, &str),
    what: &str,
    bands: &HashMap<String, J>,
) -> Result<J, String> {
    let d = r.ctx.runs.join(step);
    let s = sizes(r, &d, arms, bands)?;
    Ok(obj([
        ("what", J::Str(what.into())),
        ("accounting", runs::accounting(&d)?),
        ("any_regression", J::Bool(any_regression(&s))),
        ("sizes", s),
    ]))
}

pub fn analyse(r: &Run) -> Result<J, String> {
    let runs_dir = &r.ctx.runs;
    let mut doc = vec![
        (
            "steps".to_string(),
            J::Arr(r.steps.iter().map(|s| J::Str(s.clone())).collect()),
        ),
        (
            "host".into(),
            rounds::opt(json::read_opt(&runs_dir.join("host.json"))?),
        ),
        (
            "aa_source".into(),
            rounds::opt(json::read_opt(&runs_dir.join("aa-source.json"))?),
        ),
    ];
    let conf = runs_dir.join("conformance");
    if conf.exists() {
        let mut kv = Vec::new();
        for name in ["candidate", "base"] {
            if let Some(d) = json::read_opt(&conf.join(format!("{name}.json")))? {
                kv.push((name.to_string(), conformance_summary(&d)?));
            }
        }
        doc.push(("conformance".into(), J::Obj(kv)));
    }
    if let Some(p) = json::read_opt(&runs_dir.join("pin").join("pin.json"))? {
        let failing: Vec<J> = p
            .at("rows")?
            .as_arr()
            .unwrap_or(&[])
            .iter()
            .filter(|e| !e.get("equal").is_some_and(J::truthy))
            .filter_map(|e| e.get("row").cloned())
            .collect();
        let mut kv = match rounds::pin_summary(&Some(p))? {
            J::Obj(kv) => kv,
            _ => Vec::new(),
        };
        kv.push(("failing_rows".into(), J::Arr(failing)));
        doc.push(("pin".into(), J::Obj(kv)));
    }
    if let Some(t) = json::read_opt(&runs_dir.join("translate").join("translate.json"))? {
        let failing: Vec<J> = t
            .at("rows")?
            .as_arr()
            .unwrap_or(&[])
            .iter()
            .filter(|e| {
                !(e.get("same_outputs").is_some_and(J::truthy)
                    && e.get("same_json_as_reference").is_some_and(J::truthy))
            })
            .filter_map(|e| e.get("row").cloned())
            .collect();
        doc.push((
            "translate".into(),
            obj([
                ("held", t.at("held")?.clone()),
                ("failing_rows", J::Arr(failing)),
            ]),
        ));
    }
    let timed_steps = ["timing", "chain", "v2timing"];
    if timed_steps.iter().any(|s| runs_dir.join(s).exists()) {
        let (bands, source) = aa_bands(r)?;
        doc.push(("aa_bands".into(), J::Str(source)));
        if runs_dir.join("timing").exists() {
            doc.push((
                "timing".into(),
                timed(
                    r,
                    "timing",
                    ("base", "cand"),
                    "paired cold-time ratio, base over candidate, per size; above 1 is faster",
                    &bands,
                )?,
            ));
        }
        let chain_dir = runs_dir.join("chain");
        if chain_dir.exists() {
            let names: Vec<&str> = CHAIN_ORDER
                .iter()
                .copied()
                .filter(|n| chain_dir.join(n).exists())
                .collect();
            let mut steps = Vec::new();
            for w in names.windows(2) {
                let s = sizes(r, &chain_dir, (w[0], w[1]), &bands)?;
                steps.push((
                    w[1].to_string(),
                    obj([
                        ("against", J::Str(w[0].into())),
                        ("any_regression", J::Bool(any_regression(&s))),
                        ("sizes", s),
                    ]),
                ));
            }
            doc.push((
                "chain".into(),
                obj([
                    (
                        "what",
                        J::Str(
                            "per step, the paired cold-time ratio of the arm before it over the \
                             step's arm, per size; above 1 is faster"
                                .into(),
                        ),
                    ),
                    (
                        "arms",
                        J::Arr(names.iter().map(|n| J::Str((*n).into())).collect()),
                    ),
                    ("accounting", runs::accounting(&chain_dir)?),
                    ("steps", J::Obj(steps)),
                ]),
            ));
        }
        if runs_dir.join("v2timing").exists() {
            doc.push((
                "v2timing".into(),
                timed(
                    r,
                    "v2timing",
                    ("v1", "v2"),
                    "paired cold-time ratio, v1 file over v2 translation, per size",
                    &bands,
                )?,
            ));
        }
    }
    Ok(J::Obj(doc))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_koblitz_orders_are_the_point_counts() {
        // Every point of y² + xy = x³ + ax² + 1 over GF(2^n), counted.
        for n in 2..=9 {
            let f = find_irreducible_sparse(n).unwrap();
            let mul = |x: u64, y: u64| clmul_mod(x, y, f, n);
            for a in 0..=1u8 {
                let mut count = 1u128; // the point at infinity
                for x in 0..1u64 << n {
                    let rhs = mul(mul(x, x), x) ^ mul(a.into(), mul(x, x)) ^ 1;
                    count += (0..1u64 << n)
                        .filter(|&y| mul(y, y) ^ mul(x, y) == rhs)
                        .count() as u128;
                }
                assert_eq!(koblitz_order(a, n).unwrap(), count, "a = {a}, n = {n}");
            }
        }
    }

    #[test]
    fn largest_prime_factors_agree_with_trial_division() {
        let naive = |mut m: u64| {
            let (mut d, mut last) = (2u64, 1u64);
            while d * d <= m {
                while m.is_multiple_of(d) {
                    m /= d;
                    last = d;
                }
                d += if d == 2 { 1 } else { 2 };
            }
            if m > 1 {
                m.max(last)
            } else {
                last
            }
        };
        for m in (2u64..5000).chain([1 << 40, 999_999_000_001, 1_000_003 * 1_000_033]) {
            assert_eq!(
                largest_prime_factor(m.into()).unwrap(),
                naive(m).into(),
                "m = {m}"
            );
        }
        // Two factors past the trial bound: Pollard's rho splits them.
        let (p, q) = (1_048_583u64, 2_097_169u64);
        assert!(is_prime(p) && is_prime(q));
        assert_eq!(largest_prime_factor((p * q).into()).unwrap(), q.into());
    }

    #[test]
    fn the_sparse_moduli_are_the_librarys_at_every_one_word_degree() {
        assert_eq!(find_irreducible_sparse(2).unwrap(), 0b111);
        assert_eq!(find_irreducible_sparse(3).unwrap(), 0b1011);
        for n in 1..=62 {
            let f = find_irreducible_sparse(n).unwrap();
            let lib = crypto_lib::cryptanalysis::koblitz_index_calculus::find_irreducible_sparse(n)
                .unwrap();
            let full = lib
                .low_terms
                .iter()
                .fold(1u64 << lib.degree, |acc, &t| acc | 1 << t);
            assert_eq!(f, full, "n = {n}");
        }
    }

    /// C009's frozen document is the generator's translation of the smoke
    /// row `M1-T01` at `icv1-f2m31-tm90707-c95f16f5`, under a name of its
    /// own: the reference must write it again, byte for byte.
    #[test]
    fn the_reference_reproduces_c009s_frozen_translation() {
        let programme = Path::new(env!("CARGO_MANIFEST_DIR")).join("research/ic_tool_program");
        let frozen = programme.join("conformance/v2/params/C009-translation-of-smoke-row.json");
        let v1 = json::read(&programme.join("suite/v1/params/smoke/k0n31/M1-T01.json")).unwrap();
        let name = json::read(&frozen).unwrap().at("name").unwrap().clone();
        let J::Obj(mut kv) = reference_translation(&v1, 0, 31, 0x230001).unwrap() else {
            panic!("a translation is an object")
        };
        for (k, v) in kv.iter_mut() {
            if k == "name" {
                *v = name.clone();
            }
        }
        assert_eq!(
            json::dumps(&J::Obj(kv), 1) + "\n",
            std::fs::read_to_string(&frozen).unwrap()
        );
    }

    #[test]
    fn irreducibility_agrees_with_trial_division_at_small_degrees() {
        // f is irreducible iff no polynomial of degree 1..=deg/2 divides it.
        let divides = |g: u64, mut f: u64| {
            let dg = 63 - g.leading_zeros();
            while f != 0 && 63 - f.leading_zeros() >= dg {
                f ^= g << (63 - f.leading_zeros() - dg);
            }
            f == 0
        };
        for f in 2u64..1 << 11 {
            let n = 63 - f.leading_zeros();
            let naive = n >= 1 && (2u64..1 << (n / 2 + 1)).all(|g| g == f || !divides(g, f));
            assert_eq!(is_irreducible_f2(f), naive, "f = {f:#b}");
        }
    }
}
