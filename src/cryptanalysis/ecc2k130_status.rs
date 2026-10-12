//! The ECC2K-130 campaign dashboard, ported from `ecc2k130/aws/status.py`.
//!
//! It reads the slot table (a rehearsal store's `slots.json`, the DynamoDB
//! table, or the bucket's `slots/` objects) and sums what the walkers
//! report. The text and the `--json` summary are `status.py`'s byte for
//! byte, followed by what it did not show: each slot's steps per walk.
//!
//! Progress is measured against 2^60.9 expected iterations (Bailey et al.),
//! an expectation and not a deadline: the probability that a collision has
//! not happened after c times the expected work is about exp(-pi c^2 / 4),
//! 17% at c = 1.5 and 4% at c = 2. Iterations are each slot's checkpointed
//! per-walk iteration base times that slot's own walk count, so they survive
//! worker restarts, and a slot's live but unuploaded stretch is not counted.
//! The walk count is not a campaign constant: an Ada slot lets the client
//! size its grid and runs tens of times fewer walks than the 385,024 x 16
//! Blackwell preset, so its base climbs that much faster. A slot that does
//! not report its walks (an old worker) is counted at `--walks` and listed
//! as assumed.
//!
//! Every distinguished point ends a walk and the walk restarts from the next
//! seed, so a slot's iterations over its points is its mean walk length.
//! At the campaign's weight 32 that is 2^28.41 steps
//! (`dp_ingest.CAMPAIGN_ITER_PER_DP_LOG2`). Below 50,000 points the ratio is
//! mostly the lag between the checkpoint and the uploads, so it is not
//! judged; a slot more than 2x away is off weight, the rule the ingest host
//! applies.

use std::fmt::Write as _;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::process::Command;

use super::ecc2k130_pyjson::{
    dumps_floats, format_g, int_literal, parse_json, py_float, py_int, py_repr, py_str, py_strip,
    py_sum_floats, Json, Obj, Style,
};

/// `2.0 ** 60.9` as CPython computed it.
pub const EXPECTED_ITERS: f64 = f64::from_bits(0x43bd_db68_0117_ab0a);
/// The Blackwell preset, 385,024 workers x 16 batch.
pub const PRESET_WALKS: i128 = 385_024 * 16;
/// `dp_ingest.CAMPAIGN_ITER_PER_DP_LOG2`: steps per point at weight 32.
pub const STEPS_PER_WALK_LOG2: f64 = 28.41;
/// `dp_ingest.DP_RATIO_MIN_RECORDS`.
pub const MIN_POINTS: i128 = 50_000;
/// `dp_ingest.DP_RATIO_TOLERANCE_LOG2`: more than 2x away is off weight.
pub const TOLERANCE_LOG2: f64 = 1.0;

// ── Reading the slot table ─────────────────────────────────────────

/// Where the slots live, in the order `status.py` preferred them.
pub enum Source {
    Local(PathBuf),
    Table(String),
    Bucket(String),
}

/// `status.scan`. `aws` is the CLI to run, `aws` itself outside tests.
pub fn scan(source: &Source, aws: &str) -> Result<Vec<Obj>, String> {
    match source {
        Source::Local(dir) => scan_local(dir),
        Source::Table(table) => scan_table(table, aws),
        Source::Bucket(bucket) => scan_bucket(bucket, aws),
    }
}

/// `[dict(v, slot=int(k)) for k, v in json.load(slots.json).items()]`.
fn scan_local(dir: &Path) -> Result<Vec<Obj>, String> {
    let path = dir.join("slots.json");
    let text = fs::read_to_string(&path).map_err(|e| os_error(&e, &path))?;
    let Json::Obj(slots) = parse_json(&text)? else {
        return Err("'list' object has no attribute 'items'".into());
    };
    let mut items = Vec::new();
    for (k, v) in slots.iter() {
        let Json::Obj(item) = v else {
            return Err("dictionary update sequence element #0 has length 1; 2 is required".into());
        };
        let mut item = item.clone();
        item.set("slot", Json::Int(py_int(&Json::str(k))?));
        items.push(item);
    }
    Ok(items)
}

/// The DynamoDB table, each attribute through `fromDv`.
fn scan_table(table: &str, aws: &str) -> Result<Vec<Obj>, String> {
    let args = [
        "dynamodb",
        "scan",
        "--table-name",
        table,
        "--output",
        "json",
    ];
    let out = run(aws, &args)?;
    if !out.status.success() {
        let mut argv = vec![aws.to_string()];
        argv.extend(args.iter().map(|a| a.to_string()));
        let shown: Vec<String> = argv.iter().map(|a| py_repr(&Json::str(a))).collect();
        return Err(format!(
            "Command '[{}]' returned non-zero exit status {}.",
            shown.join(", "),
            out.status.code().unwrap_or(-1)
        ));
    }
    let reply = parse_json(&String::from_utf8_lossy(&out.stdout))?;
    let items = match reply.as_obj().and_then(|o| o.get("Items")) {
        None => return Ok(Vec::new()),
        Some(Json::List(items)) => items,
        Some(_) => return Err("DynamoDB Items is not a list".into()),
    };
    let mut slots = Vec::new();
    for it in items {
        let Json::Obj(it) = it else {
            return Err("a DynamoDB item is not an object".into());
        };
        let mut item = Obj::new();
        for (k, v) in it.iter() {
            item.set(k, from_dv(v)?);
        }
        slots.push(item);
    }
    Ok(slots)
}

/// `status.fromDv`: a number, a string or a bool; anything else is None.
fn from_dv(v: &Json) -> Result<Json, String> {
    let Json::Obj(v) = v else {
        return Err("a DynamoDB attribute is not an object".into());
    };
    if let Some(n) = v.get("N") {
        let Json::Str(n) = n else {
            return Err("a DynamoDB number is not a string".into());
        };
        return if n.contains('.') || n.contains('e') {
            py_float(&Json::str(n)).map(Json::Float)
        } else {
            py_int(&Json::str(n)).map(Json::Int)
        };
    }
    if let Some(s) = v.get("S") {
        return Ok(s.clone());
    }
    if let Some(b) = v.get("BOOL") {
        return Ok(b.clone());
    }
    Ok(Json::Null)
}

/// `worker.S3Slots.scan`: every `slots/slot-NNNNN.json`, with its slot
/// number and ETag. The object read is the one the slot number names.
fn scan_bucket(bucket: &str, aws: &str) -> Result<Vec<Obj>, String> {
    let listing = run(
        aws,
        &[
            "s3api",
            "list-objects-v2",
            "--bucket",
            bucket,
            "--prefix",
            "slots/",
            "--query",
            "Contents[].Key",
            "--output",
            "json",
        ],
    )?;
    if !listing.status.success() {
        return Err(format!(
            "s3 list failed: {}",
            py_strip(&String::from_utf8_lossy(&listing.stderr))
        ));
    }
    let keys = match parse_json(&String::from_utf8_lossy(&listing.stdout))? {
        Json::List(keys) => keys,
        k if !k.truthy() => Vec::new(),
        _ => return Err("s3 listing is not a list of keys".into()),
    };
    let tmp = std::env::var("TMPDIR")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("/tmp"))
        .join(format!("ecc-slot-{}.json", std::process::id()));
    let mut items = Vec::new();
    for key in &keys {
        let Json::Str(key) = key else {
            return Err("s3 listing holds a key that is not a string".into());
        };
        let name = key.rsplit('/').next().unwrap_or(key);
        if !(name.starts_with("slot-") && name.ends_with(".json")) {
            continue;
        }
        let slot = py_int(&Json::str(&name[5..name.len() - 5]))?;
        let object = format!("slots/slot-{slot:05}.json");
        let tmp_arg = tmp.to_string_lossy();
        let got = run(
            aws,
            &[
                "s3api",
                "get-object",
                "--bucket",
                bucket,
                "--key",
                &object,
                &tmp_arg,
                "--output",
                "json",
            ],
        )?;
        if !got.status.success() {
            let stderr = String::from_utf8_lossy(&got.stderr);
            if ["NoSuchKey", "404", "Not Found"]
                .iter()
                .any(|m| stderr.contains(m))
            {
                continue;
            }
            return Err(format!("s3 get-object failed: {}", py_strip(&stderr)));
        }
        let meta = parse_json(&String::from_utf8_lossy(&got.stdout))?;
        let etag = meta
            .as_obj()
            .and_then(|m| m.get("ETag"))
            .cloned()
            .ok_or("'ETag'")?;
        let text = fs::read_to_string(&tmp).map_err(|e| os_error(&e, &tmp))?;
        let Json::Obj(mut item) = parse_json(&text)? else {
            return Err(format!("{object} is not a JSON object"));
        };
        item.set("slot", Json::Int(slot));
        item.set("_etag", etag);
        items.push(item);
    }
    let _ = fs::remove_file(&tmp);
    Ok(items)
}

fn run(program: &str, args: &[&str]) -> Result<std::process::Output, String> {
    Command::new(program)
        .args(args)
        .output()
        .map_err(|e| os_error(&e, Path::new(program)))
}

/// An `OSError` as Python prints it: `[Errno 2] No such file or directory: 'x'`.
fn os_error(e: &io::Error, path: &Path) -> String {
    let shown = py_repr(&Json::str(&path.to_string_lossy()));
    match e.raw_os_error() {
        Some(code) => {
            let text = io::Error::from_raw_os_error(code).to_string();
            let text = text.split(" (os error").next().unwrap_or(&text);
            format!("[Errno {code}] {text}: {shown}")
        }
        None => format!("{e}: {shown}"),
    }
}

// ── The snapshot ───────────────────────────────────────────────────

/// `int(it.get(key) or 0)`.
fn int_field(it: &Obj, key: &str) -> Result<i128, String> {
    match it.get(key) {
        Some(v) if v.truthy() => py_int(v),
        _ => Ok(0),
    }
}

/// `float(it.get(key) or 0)`.
fn float_field(it: &Obj, key: &str) -> Result<f64, String> {
    match it.get(key) {
        Some(v) if v.truthy() => py_float(v),
        _ => Ok(0.0),
    }
}

/// `status.slotWalks`: what the slot reported, else the preset.
fn slot_walks(it: &Obj, fallback: i128) -> Result<i128, String> {
    let walks = int_field(it, "walks")?;
    Ok(if walks > 0 { walks } else { fallback })
}

/// `status.slotIters`.
fn slot_iters(it: &Obj, fallback: i128) -> Result<i128, String> {
    let base = int_field(it, "ckptIter")?.max(0);
    let walks = slot_walks(it, fallback)?;
    base.checked_mul(walks)
        .ok_or_else(|| "a slot's iteration count does not fit in 128 bits".into())
}

fn total(terms: impl IntoIterator<Item = Result<i128, String>>) -> Result<i128, String> {
    let mut sum = 0i128;
    for term in terms {
        sum = sum
            .checked_add(term?)
            .ok_or("an iteration sum does not fit in 128 bits")?;
    }
    Ok(sum)
}

/// `int_value >= now`, exactly, as Python compares an int with a float.
fn at_least(value: i128, now: f64) -> bool {
    !now.is_nan()
        && (now == f64::NEG_INFINITY
            || (now.ceil() < 2f64.powi(127) && value >= now.ceil() as i128))
}

/// A slot's steps per walk, as `(log2, walks assumed)`, once it has the
/// points to judge it by.
fn steps_per_walk(it: &Obj, fallback: i128) -> Result<Option<(f64, bool)>, String> {
    let iters = slot_iters(it, fallback)?;
    let points = int_field(it, "dpUploaded")?;
    if points < MIN_POINTS || iters <= 0 {
        return Ok(None);
    }
    let assumed = int_field(it, "walks")? <= 0;
    Ok(Some(((iters as f64 / points as f64).log2(), assumed)))
}

fn off_weight(log2: f64) -> bool {
    (log2 - STEPS_PER_WALK_LOG2).abs() > TOLERANCE_LOG2
}

/// `round(x, 3)`.
fn round3(x: f64) -> f64 {
    format!("{x:.3}").parse().unwrap_or(x)
}

/// `status.report`, written into `out` as it goes: on an error `out` holds
/// what `status.py` had printed before raising.
pub fn report(
    out: &mut String,
    items: &[Obj],
    walks_per_slot: i128,
    now: f64,
    as_json: bool,
) -> Result<(), String> {
    let mut alive = Vec::new();
    for it in items {
        if at_least(int_field(it, "leaseUntil")?, now) {
            alive.push(it);
        }
    }
    let mut rates = Vec::new();
    for it in &alive {
        rates.push(float_field(it, "rate")?);
    }
    let rate_json = py_sum_floats(&rates);
    let rate = match rate_json {
        Json::Float(r) => r,
        _ => 0.0,
    };
    let ckpt_iters = total(items.iter().map(|it| slot_iters(it, walks_per_slot)))?;
    let mut assumed = Vec::new();
    for it in items {
        if int_field(it, "walks")? <= 0 && int_field(it, "ckptIter")? > 0 {
            assumed.push(it);
        }
    }
    let dps = total(items.iter().map(|it| int_field(it, "dpUploaded")))?;
    let state_is = |it: &&Obj, s: &str| it.get("state") == Some(&Json::str(s));
    let solved: Vec<&Obj> = items.iter().filter(|it| state_is(it, "solved")).collect();
    let retired = items.iter().filter(|it| state_is(it, "retired")).count();
    let errors = items.iter().filter(|it| state_is(it, "error")).count();
    let frac = ckpt_iters as f64 / EXPECTED_ITERS;
    let remaining = if rate > 0.0 {
        (EXPECTED_ITERS - ckpt_iters as f64) / rate
    } else {
        f64::INFINITY
    };
    let assumed_iters = total(assumed.iter().map(|it| slot_iters(it, walks_per_slot)))?;

    let mut judged = Vec::new();
    for it in items {
        if let Some((log2, false)) = steps_per_walk(it, walks_per_slot)? {
            judged.push((it, log2));
        }
    }
    let mean_steps = |slots: &[(&Obj, f64)]| -> Result<Option<f64>, String> {
        let iters = total(slots.iter().map(|(it, _)| slot_iters(it, walks_per_slot)))?;
        let points = total(slots.iter().map(|(it, _)| int_field(it, "dpUploaded")))?;
        Ok((!slots.is_empty()).then(|| iters as f64 / points as f64))
    };
    let (off, on): (Vec<(&Obj, f64)>, Vec<(&Obj, f64)>) =
        judged.iter().partition(|(_, log2)| off_weight(*log2));
    let pooled = mean_steps(&judged)?;
    let on_weight = mean_steps(&on)?;

    let count = |n: usize| Json::Int(n as i128);
    let summary = Obj::new()
        .with("slots", count(items.len()))
        .with("alive", count(alive.len()))
        .with("retired", count(retired))
        .with("errors", count(errors))
        .with("rateBps", rate_json)
        .with("checkpointedIterations", Json::Int(ckpt_iters))
        .with("fractionOfExpected", Json::Float(frac))
        .with("dpUploaded", Json::Int(dps))
        .with("etaSecondsAtCurrentRate", Json::Float(remaining))
        .with(
            "solved",
            Json::List(solved.iter().map(|it| Json::Obj((*it).clone())).collect()),
        )
        .with("walksAssumedSlots", count(assumed.len()))
        .with("assumedWalksPerSlot", Json::Int(walks_per_slot))
        .with("checkpointedIterationsAssumed", Json::Int(assumed_iters))
        .with("stepsPerWalk", pooled.map_or(Json::Null, Json::Float))
        .with(
            "stepsPerWalkLog2",
            pooled.map_or(Json::Null, |p| Json::Float(round3(p.log2()))),
        )
        .with("stepsPerWalkSlots", count(judged.len()))
        .with("stepsPerWalkExpectedLog2", Json::Float(STEPS_PER_WALK_LOG2))
        .with(
            "offWeightSlots",
            Json::List(
                off.iter()
                    .map(|(it, _)| it.get("slot").cloned().unwrap_or(Json::Null))
                    .collect(),
            ),
        )
        .with(
            "stepsPerWalkOnWeightLog2",
            on_weight.map_or(Json::Null, |p| Json::Float(round3(p.log2()))),
        );
    if as_json {
        out.push_str(&dumps_floats(&Json::Obj(summary), Style::Indent(1), false));
        out.push('\n');
        return Ok(());
    }

    let _ = writeln!(
        out,
        "slots {}, alive {}, retired {}, errors {}",
        items.len(),
        alive.len(),
        retired,
        errors
    );
    let _ = writeln!(
        out,
        "aggregate {} B it/s   checkpointed {} iterations = {}% of 2^60.9   {} points uploaded",
        fixed(rate / 1e9, 3),
        format_g(ckpt_iters as f64, 4),
        fixed(100.0 * frac, 3),
        dps
    );
    if !assumed.is_empty() {
        let _ = writeln!(
            out,
            "  {} slot(s) do not report their walk count; counted at the {}-walk preset \
             ({} of the iterations above are that guess)",
            assumed.len(),
            walks_per_slot,
            format_g(assumed_iters as f64, 4)
        );
    }
    if rate > 0.0 {
        let _ = writeln!(
            out,
            "at this rate the expected remaining work takes {} (17% chance it needs 1.5x, 4% chance 2x)",
            human_time(remaining)?
        );
    }
    match pooled {
        Some(p) => {
            let _ = writeln!(
                out,
                "steps per walk 2^{:.2} over {} slot(s) with {}+ points (2^{STEPS_PER_WALK_LOG2} expected at weight 32)",
                p.log2(),
                judged.len(),
                MIN_POINTS
            );
        }
        None => {
            let _ = writeln!(
                out,
                "steps per walk: no slot reporting its walks has {MIN_POINTS} points yet \
                 (2^{STEPS_PER_WALK_LOG2} expected at weight 32)"
            );
        }
    }
    if !off.is_empty() {
        let shown: Vec<String> = off
            .iter()
            .map(|(it, log2)| {
                let slot = it.get("slot").map_or("?".into(), py_str);
                format!("slot {slot} (2^{log2:.2})")
            })
            .collect();
        let _ = writeln!(
            out,
            "  off weight, more than 2x from 2^{STEPS_PER_WALK_LOG2}: {}",
            shown.join(", ")
        );
        if let Some(p) = on_weight {
            let _ = writeln!(
                out,
                "  without them, 2^{:.2} over {} slot(s)",
                p.log2(),
                on.len()
            );
        }
    }
    if let Some(first) = solved.first() {
        let solution = first.get("solution").cloned().unwrap_or(Json::Null);
        let _ = writeln!(out, "SOLVED: {}", py_str(&solution));
    }
    let _ = writeln!(
        out,
        "{:>5} {:<28} {:<24} {:>7} {:>10} {:>12} {:>10} {:>8} {:>10}",
        "slot", "owner", "gpu", "B it/s", "dp", "ckpt iter", "walks", "lease", "steps/walk"
    );
    let mut rows: Vec<&Obj> = Vec::with_capacity(items.len());
    for it in items {
        match it.get("slot") {
            Some(Json::Int(_) | Json::Bool(_) | Json::Float(_)) => rows.push(it),
            Some(other) => {
                return Err(format!(
                    "%d format: a real number is required, not {}",
                    python_type(other)
                ))
            }
            None => return Err("'slot'".into()),
        }
    }
    rows.sort_by(|a, b| {
        let key = |it: &Obj| match it.get("slot") {
            Some(Json::Int(i)) => *i as f64,
            Some(Json::Bool(b)) => f64::from(u8::from(*b)),
            Some(Json::Float(f)) => *f,
            _ => f64::NAN,
        };
        key(a).total_cmp(&key(b))
    });
    for it in rows {
        let lease = int_field(it, "leaseUntil")? as f64 - now;
        let state = it.get("state").cloned().unwrap_or(Json::str("?"));
        let status = if lease > 0.0 {
            format!("{}s", percent_d(&Json::Float(lease))?)
        } else if state == Json::str("active") {
            "expired".into()
        } else {
            py_str(&state)
        };
        let walks = int_field(it, "walks")?;
        let slot = percent_d(it.get("slot").expect("checked above"))?;
        let owner: String = py_str(it.get("owner").unwrap_or(&Json::str("-")))
            .chars()
            .take(28)
            .collect();
        let gpu: String = py_str(it.get("gpuName").unwrap_or(&Json::str("-")))
            .chars()
            .take(24)
            .collect();
        let rate = fixed(float_field(it, "rate")? / 1e9, 3);
        let dp = int_field(it, "dpUploaded")?;
        let ckpt = int_field(it, "ckptIter")?;
        let walks = if walks > 0 {
            walks.to_string()
        } else {
            format!("{walks_per_slot}?")
        };
        let steps = match steps_per_walk(it, walks_per_slot)? {
            None => "-".to_string(),
            Some((log2, true)) => format!("2^{log2:.2}?"),
            Some((log2, false)) if off_weight(log2) => format!("2^{log2:.2}!"),
            Some((log2, false)) => format!("2^{log2:.2}"),
        };
        let _ = writeln!(
            out,
            "{slot:>5} {owner:<28} {gpu:<24} {rate:>7} {dp:>10} {ckpt:>12} {walks:>10} {status:>8} {steps:>10}"
        );
    }
    Ok(())
}

fn python_type(v: &Json) -> &'static str {
    match v {
        Json::Null => "NoneType",
        Json::Bool(_) => "bool",
        Json::Int(_) => "int",
        Json::Float(_) => "float",
        Json::Str(_) => "str",
        Json::List(_) => "list",
        Json::Obj(_) => "dict",
    }
}

/// `'%.{decimals}f' % x`.
fn fixed(x: f64, decimals: usize) -> String {
    if x.is_nan() {
        "nan".into()
    } else if x.is_infinite() {
        if x > 0.0 { "inf" } else { "-inf" }.into()
    } else {
        format!("{x:.decimals$}")
    }
}

/// `'%d' % v`: a float is truncated toward zero first.
fn percent_d(v: &Json) -> Result<String, String> {
    match v {
        Json::Int(i) => Ok(i.to_string()),
        Json::Bool(b) => Ok(u8::from(*b).to_string()),
        Json::Float(f) if f.is_nan() => Err("cannot convert float NaN to integer".into()),
        Json::Float(f) if f.is_infinite() => Err("cannot convert float infinity to integer".into()),
        Json::Float(f) => {
            let t = f.trunc();
            Ok(if t == 0.0 {
                "0".into()
            } else {
                format!("{t:.0}")
            })
        }
        other => Err(format!(
            "%d format: a real number is required, not {}",
            python_type(other)
        )),
    }
}

/// `vx // wx` on floats, as CPython rounds it.
fn floor_div(vx: f64, wx: f64) -> f64 {
    let mut m = vx % wx;
    let mut div = (vx - m) / wx;
    if m != 0.0 {
        if (wx < 0.0) != (m < 0.0) {
            m += wx;
            div -= 1.0;
        }
    } else {
        m = 0f64.copysign(wx);
    }
    let _ = m;
    if div != 0.0 {
        let mut floor = div.floor();
        if div - floor > 0.5 {
            floor += 1.0;
        }
        floor
    } else {
        0f64.copysign(vx / wx)
    }
}

/// `status.humanTime`.
fn human_time(sec: f64) -> Result<String, String> {
    Ok(if sec < 3600.0 {
        format!("{}m", percent_d(&Json::Float(floor_div(sec, 60.0)))?)
    } else if sec < 86400.0 {
        format!("{}h", fixed(sec / 3600.0, 1))
    } else {
        format!("{}d", fixed(sec / 86400.0, 1))
    })
}

/// `--walks` as `type=int` reads it.
pub fn parse_walks(s: &str) -> Result<i128, String> {
    int_literal(s).ok_or_else(|| format!("invalid int value: {}", py_repr(&Json::str(s))))
}

#[cfg(test)]
mod tests {
    use super::*;

    const ADA_WALKS: i128 = 14848 * 16;
    const NOW: f64 = 1_790_000_000.25;

    fn slot(n: i128, ckpt: i128, walks: Option<i128>, rate: f64, dp: i128) -> Obj {
        let mut it = Obj::new()
            .with("slot", Json::Int(n))
            .with("ckptIter", Json::Int(ckpt))
            .with("rate", Json::Float(rate))
            .with("dpUploaded", Json::Int(dp))
            .with("state", Json::str("active"))
            .with("leaseUntil", Json::Int(NOW as i128 + 180));
        if let Some(w) = walks {
            it.set("walks", Json::Int(w));
        }
        it
    }

    fn summary(items: &[Obj], walks: i128) -> Obj {
        let mut out = String::new();
        report(&mut out, items, walks, NOW, true).unwrap();
        match parse_json(&out).unwrap() {
            Json::Obj(o) => o,
            other => panic!("{other:?}"),
        }
    }

    fn int(o: &Obj, key: &str) -> i128 {
        o.get(key).and_then(Json::as_int).unwrap()
    }

    fn float(o: &Obj, key: &str) -> f64 {
        match o.get(key) {
            Some(Json::Float(f)) => *f,
            other => panic!("{key}: {other:?}"),
        }
    }

    // `test_status_walks.Walks`, which these replace.

    #[test]
    fn reported_walks_beat_the_preset() {
        let items = [
            slot(0, 1000, Some(PRESET_WALKS), 0.0, 0),
            slot(1, 1000, Some(ADA_WALKS), 0.0, 0),
        ];
        assert_eq!(
            int(&summary(&items, PRESET_WALKS), "checkpointedIterations"),
            1000 * PRESET_WALKS + 1000 * ADA_WALKS
        );
    }

    #[test]
    fn an_ada_slot_is_not_counted_at_the_blackwell_preset() {
        let honest = int(
            &summary(&[slot(1, 1000, Some(ADA_WALKS), 0.0, 0)], PRESET_WALKS),
            "checkpointedIterations",
        );
        assert_eq!(honest, 1000 * ADA_WALKS);
        assert!((1000 * PRESET_WALKS) as f64 / honest as f64 > 25.0);
    }

    #[test]
    fn the_same_group_operations_count_the_same_whatever_the_grid() {
        let big = summary(
            &[slot(0, ADA_WALKS, Some(PRESET_WALKS), 0.0, 0)],
            PRESET_WALKS,
        );
        let small = summary(
            &[slot(1, PRESET_WALKS, Some(ADA_WALKS), 0.0, 0)],
            PRESET_WALKS,
        );
        assert_eq!(
            int(&big, "checkpointedIterations"),
            int(&small, "checkpointedIterations")
        );
    }

    #[test]
    fn a_slot_without_walks_is_assumed_and_said_so() {
        let got = summary(
            &[
                slot(0, 1000, None, 0.0, 0),
                slot(1, 1000, Some(ADA_WALKS), 0.0, 0),
            ],
            PRESET_WALKS,
        );
        assert_eq!(int(&got, "walksAssumedSlots"), 1);
        assert_eq!(int(&got, "assumedWalksPerSlot"), PRESET_WALKS);
        assert_eq!(
            int(&got, "checkpointedIterationsAssumed"),
            1000 * PRESET_WALKS
        );
        assert_eq!(
            int(&got, "checkpointedIterations"),
            1000 * PRESET_WALKS + 1000 * ADA_WALKS
        );
    }

    #[test]
    fn an_idle_slot_is_neither_counted_nor_flagged() {
        let got = summary(&[slot(0, -1, None, 0.0, 0)], PRESET_WALKS);
        assert_eq!(int(&got, "checkpointedIterations"), 0);
        assert_eq!(int(&got, "walksAssumedSlots"), 0);
    }

    #[test]
    fn eta_and_fraction_follow_the_honest_count() {
        let got = summary(&[slot(0, 1000, Some(ADA_WALKS), 1e9, 0)], PRESET_WALKS);
        let iters = (1000 * ADA_WALKS) as f64;
        assert_eq!(float(&got, "fractionOfExpected"), iters / EXPECTED_ITERS);
        assert_eq!(
            float(&got, "etaSecondsAtCurrentRate"),
            (EXPECTED_ITERS - iters) / 1e9
        );
    }

    #[test]
    fn expected_iterations_are_two_to_the_60_9() {
        assert!((EXPECTED_ITERS.log2() - 60.9).abs() < 1e-12);
    }

    // What status.py did not show.

    #[test]
    fn steps_per_walk_pool_the_slots_that_report_their_walks() {
        let walks = 6_160_384i128;
        let at_weight = |n, points: i128| {
            let steps = 2f64.powf(28.41) * points as f64;
            slot(n, (steps / walks as f64) as i128, Some(walks), 14e9, points)
        };
        let items = [
            at_weight(0, 80_000),
            at_weight(1, 120_000),
            at_weight(2, 49_999),
            slot(3, 1 << 20, Some(walks), 14e9, 3_000_000),
            slot(4, 1 << 20, None, 14e9, 3_000_000),
        ];
        let got = summary(&items, PRESET_WALKS);
        assert_eq!(int(&got, "stepsPerWalkSlots"), 3);
        assert_eq!(
            got.get("offWeightSlots"),
            Some(&Json::List(vec![Json::Int(3)]))
        );
        let pooled = float(&got, "stepsPerWalk");
        let want = (2f64.powf(28.41) * 200_000.0 + ((1i128 << 20) * walks) as f64) / 3_200_000.0;
        assert!((pooled / want - 1.0).abs() < 1e-6, "{pooled} vs {want}");
        assert_eq!(float(&got, "stepsPerWalkLog2"), round3(want.log2()));
        assert_eq!(float(&got, "stepsPerWalkExpectedLog2"), 28.41);
        assert_eq!(float(&got, "stepsPerWalkOnWeightLog2"), 28.41);

        let mut text = String::new();
        report(&mut text, &items, PRESET_WALKS, NOW, false).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        let column = |n: usize| {
            lines[lines.len() - 5 + n]
                .split_whitespace()
                .last()
                .unwrap()
        };
        assert_eq!(column(0), "2^28.41");
        assert_eq!(column(1), "2^28.41");
        assert_eq!(column(2), "-");
        assert!(column(3).ends_with('!'), "{}", column(3));
        assert!(column(4).ends_with('?'), "{}", column(4));
        assert!(text.contains("  off weight, more than 2x from 2^28.41: slot 3 (2^"));
        assert!(
            text.contains("\n  without them, 2^28.41 over 2 slot(s)\n"),
            "{text}"
        );
    }

    #[test]
    fn nothing_to_judge_says_so() {
        let mut text = String::new();
        report(
            &mut text,
            &[slot(0, 10, Some(16), 0.0, 3)],
            PRESET_WALKS,
            NOW,
            false,
        )
        .unwrap();
        assert!(text.contains("steps per walk: no slot reporting its walks has 50000 points yet"));
        let got = summary(&[slot(0, 10, Some(16), 0.0, 3)], PRESET_WALKS);
        assert_eq!(got.get("stepsPerWalk"), Some(&Json::Null));
        assert_eq!(int(&got, "stepsPerWalkSlots"), 0);
    }

    // The bytes status.py printed, for slot tables CPython 3.12 rendered on
    // 2026-10-06 (`now` pinned to NOW).

    fn table() -> Vec<Obj> {
        let mut retired = slot(7, 5, Some(16), 2.5e9, 12);
        retired.set("state", Json::str("retired"));
        retired.set("leaseUntil", Json::Int(0));
        retired.set("owner", Json::str("i-0123456789abcdef0-a-very-long-owner"));
        let mut expired = slot(2, 900, None, 13.9e9, 4);
        expired.set("leaseUntil", Json::Int(NOW as i128 - 5));
        expired.set(
            "gpuName",
            Json::str("NVIDIA RTX PRO 6000 Blackwell Server Edition"),
        );
        let mut live = slot(0, 1000, Some(ADA_WALKS), 14.1e9, 77);
        live.set("owner", Json::Null);
        vec![retired, expired, live]
    }

    #[test]
    fn text_is_status_py_text() {
        let mut text = String::new();
        report(&mut text, &table(), PRESET_WALKS, NOW, false).unwrap();
        let python = "slots 3, alive 1, retired 1, errors 0\naggregate 14.100 B it/s   checkpointed 5.782e+09 iterations = 0.000% of 2^60.9   93 points uploaded\n  1 slot(s) do not report their walk count; counted at the 6160384-walk preset (5.544e+09 of the iterations above are that guess)\nat this rate the expected remaining work takes 1766.0d (17% chance it needs 1.5x, 4% chance 2x)\n slot owner                        gpu                       B it/s         dp    ckpt iter      walks    lease\n    0 None                         -                         14.100         77         1000     237568     179s\n    2 -                            NVIDIA RTX PRO 6000 Blac  13.900          4          900   6160384?  expired\n    7 i-0123456789abcdef0-a-very-l -                          2.500         12            5         16  retired\n";
        let ours: Vec<&str> = text.lines().collect();
        let theirs: Vec<&str> = python.lines().collect();
        assert_eq!(ours.len(), theirs.len() + 1, "{text}");
        for (n, want) in theirs.iter().enumerate() {
            let got = if n < 4 { ours[n] } else { ours[n + 1] };
            let got = if n >= 4 { &got[..got.len() - 11] } else { got };
            assert_eq!(got, *want, "line {n}");
        }
        assert!(ours[4].starts_with("steps per walk: no slot"));
    }

    #[test]
    fn json_is_status_py_json() {
        let mut out = String::new();
        report(&mut out, &table(), PRESET_WALKS, NOW, true).unwrap();
        let python = "{\n \"slots\": 3,\n \"alive\": 1,\n \"retired\": 1,\n \"errors\": 0,\n \"rateBps\": 14100000000.0,\n \"checkpointedIterations\": 5781913680,\n \"fractionOfExpected\": 2.6874776904316424e-09,\n \"dpUploaded\": 93,\n \"etaSecondsAtCurrentRate\": 152583517.38432422,\n \"solved\": [],\n \"walksAssumedSlots\": 1,\n \"assumedWalksPerSlot\": 6160384,\n \"checkpointedIterationsAssumed\": 5544345600\n}\n";
        let head = &python[..python.len() - 3];
        assert!(out.starts_with(head), "{out}");
        assert!(out[head.len()..].starts_with(",\n \"stepsPerWalk\": null,"));
    }

    #[test]
    fn no_live_slot_sums_to_the_int_zero() {
        let mut out = String::new();
        report(&mut out, &[], PRESET_WALKS, NOW, true).unwrap();
        let python = "{\n \"slots\": 0,\n \"alive\": 0,\n \"retired\": 0,\n \"errors\": 0,\n \"rateBps\": 0,\n \"checkpointedIterations\": 0,\n \"fractionOfExpected\": 0.0,\n \"dpUploaded\": 0,\n \"etaSecondsAtCurrentRate\": Infinity,\n \"solved\": [],\n \"walksAssumedSlots\": 0,\n \"assumedWalksPerSlot\": 6160384,\n \"checkpointedIterationsAssumed\": 0\n}\n";
        assert!(out.starts_with(&python[..python.len() - 3]), "{out}");
    }

    #[test]
    fn python_formatting_helpers() {
        assert_eq!(floor_div(-30.0, 60.0), -1.0);
        assert_eq!(floor_div(119.9, 60.0), 1.0);
        assert_eq!(human_time(59.0).unwrap(), "0m");
        assert_eq!(human_time(-30.0).unwrap(), "-1m");
        assert_eq!(human_time(5400.0).unwrap(), "1.5h");
        assert_eq!(human_time(1.5 * 86400.0).unwrap(), "1.5d");
        assert_eq!(percent_d(&Json::Float(179.9)).unwrap(), "179");
        assert_eq!(percent_d(&Json::Float(-0.5)).unwrap(), "0");
        assert_eq!(
            percent_d(&Json::Float(1e20)).unwrap(),
            "100000000000000000000"
        );
        assert!(percent_d(&Json::str("x")).is_err());
        assert_eq!(fixed(f64::NAN, 3), "nan");
        assert_eq!(parse_walks(" 1_024 ").unwrap(), 1024);
        assert!(parse_walks("1.5").is_err());
        assert!(at_least(10, 9.5) && at_least(10, 10.0) && !at_least(10, 10.25));
    }

    #[test]
    fn errors_name_what_python_named() {
        let mut bad = slot(0, 1, Some(1), 0.0, 0);
        bad.set("leaseUntil", Json::str("soon"));
        let mut out = String::new();
        assert_eq!(
            report(&mut out, &[bad], PRESET_WALKS, NOW, false).unwrap_err(),
            "invalid literal for int() with base 10: 'soon'"
        );
        let nameless = Obj::new().with("ckptIter", Json::Int(1));
        let mut out = String::new();
        assert_eq!(
            report(&mut out, &[nameless], PRESET_WALKS, NOW, false).unwrap_err(),
            "'slot'"
        );
        assert!(
            out.ends_with("steps/walk\n"),
            "the header comes first: {out}"
        );
    }
}
