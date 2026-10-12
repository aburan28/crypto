//! Independent replay of a compact-orbit IC relation certificate.
mod reference;

use reference::{
    array, ensure, integer, mod_add, mod_mul, mod_pow, mod_sub, one, pair, point, sha256_file,
    BinaryCurve, Check, Point,
};
use serde_json::{json, Value};
use std::collections::{BTreeMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

struct Base {
    curve: BinaryCurve,
    order: u64,
    generator: Point,
    points: Vec<Point>,
    labels: Vec<(usize, u64)>,
    reps: Vec<Point>,
    hash: String,
}

fn numbers(v: &Value, name: &str) -> Check<Vec<u64>> {
    array(v, name)?
        .iter()
        .map(|x| x.as_u64().ok_or_else(|| format!("noninteger in {name}")))
        .collect()
}

fn same(a: &Value, b: &Value, label: &str) -> Check<()> {
    ensure(a == b, format!("{label} mismatch"))
}

fn check_base(base: &Value, header: &Value) -> Check<Base> {
    let n = integer(base, "n")? as u32;
    let a = integer(base, "a")?;
    let curve = BinaryCurve::new(n, a, array(base, "field_modulus_low_terms")?)?;
    let order = integer(base, "subgroup_order")?;
    ensure(order > 2, "invalid subgroup order")?;
    let points: Vec<_> = array(base, "factor_base_point_coordinates")?
        .iter()
        .map(point)
        .collect::<Check<_>>()?;
    let labels: Vec<_> = array(base, "factor_base_point_labels")?
        .iter()
        .map(|v| {
            let p = pair(v)?;
            Ok((p.0 as usize, p.1))
        })
        .collect::<Check<_>>()?;
    let reps: Vec<_> = array(base, "factor_base_representatives")?
        .iter()
        .map(point)
        .collect::<Check<_>>()?;
    let columns = integer(base, "orbit_columns")? as usize;
    ensure(
        reps.len() == columns
            && points.len() == labels.len()
            && points.len() == 2 * n as usize * columns,
        "factor-base count/orbit shape mismatch",
    )?;
    ensure(
        points.iter().all(Option::is_some),
        "identity in factor base",
    )?;
    ensure(
        points.iter().copied().collect::<HashSet<_>>().len() == points.len(),
        "duplicate factor-base point",
    )?;
    for name in [
        "base_hash",
        "factor_base_points",
        "orbit_columns",
        "subgroup_order",
        "n",
        "a",
    ] {
        let rhs = if name == "factor_base_points" {
            json!(points.len())
        } else {
            base[name].clone()
        };
        same(&header[name], &rhs, name)?;
    }
    let generator = point(&header["generator"])?;
    ensure(
        curve.on_curve(generator) && curve.scale(order, generator)?.is_none(),
        "generator subgroup check failed",
    )?;
    ensure(labels.len() > 2, "missing Frobenius label")?;
    let lambda = labels[2].1;
    ensure(lambda > 0 && lambda < order, "invalid Frobenius eigenvalue")?;
    for (col, &rep) in reps.iter().enumerate() {
        ensure(
            curve.on_curve(rep) && curve.scale(order, rep)?.is_none(),
            format!("representative {col} invalid"),
        )?;
        ensure(
            curve.scale(lambda, rep)? == curve.frob(rep),
            format!("representative {col} Frobenius mismatch"),
        )?;
        let (mut current, mut coefficient) = (rep, 1u64);
        for k in 0..n as usize {
            let i = 2 * (col * n as usize + k);
            ensure(
                points[i] == current && points[i + 1] == curve.neg(current),
                format!("orbit {col}:{k} point mismatch"),
            )?;
            ensure(
                labels[i] == (col, coefficient)
                    && labels[i + 1] == (col, mod_sub(0, coefficient, order)),
                format!("orbit {col}:{k} label mismatch"),
            )?;
            ensure(
                curve.on_curve(current),
                format!("orbit {col}:{k} off curve"),
            )?;
            current = curve.frob(current);
            coefficient = mod_mul(coefficient, lambda, order);
        }
        ensure(
            current == rep && coefficient == 1,
            format!("orbit {col} does not close"),
        )?;
    }
    Ok(Base {
        curve,
        order,
        generator,
        points,
        labels,
        reps,
        hash: base["base_hash"]
            .as_str()
            .ok_or("missing base hash")?
            .to_string(),
    })
}

fn check_pair_roots(
    curve: BinaryCurve,
    chosen: &[Point],
    codes: &[u64],
    roots: &[u64],
) -> Check<()> {
    ensure(
        chosen.len() == 4 && codes.len() == 4 && roots.len() == 2,
        "S3 pair shape mismatch",
    )?;
    for (i, &p) in chosen.iter().enumerate() {
        ensure(
            p.ok_or("identity summand")?.0 == codes[i],
            format!("summand {i} x-code mismatch"),
        )?;
    }
    for (i, &root) in roots.iter().enumerate() {
        let offset = 2 * i;
        let x1 = curve.add(chosen[offset], chosen[offset + 1])?.map(|p| p.0);
        let x2 = curve
            .add(chosen[offset], curve.neg(chosen[offset + 1]))?
            .map(|p| p.0);
        ensure(
            x1 == Some(root) || x2 == Some(root),
            format!("pair {i} S3 root mismatch"),
        )?;
    }
    Ok(())
}

fn insert_row(
    mut row: Vec<u64>,
    pivots: &mut BTreeMap<usize, Vec<u64>>,
    order: u64,
    columns: usize,
) -> Check<bool> {
    ensure(row.len() == columns + 1, "rank row width mismatch")?;
    for col in 0..columns {
        let value = row[col];
        if value == 0 {
            continue;
        }
        if let Some(pivot) = pivots.get(&col) {
            for j in col..=columns {
                row[j] = mod_sub(row[j], mod_mul(value, pivot[j], order), order);
            }
        } else {
            let inverse = mod_pow(value, order - 2, order);
            ensure(mod_mul(value, inverse, order) == 1, "noninvertible pivot")?;
            for item in row.iter_mut().take(columns + 1).skip(col) {
                *item = mod_mul(*item, inverse, order);
            }
            pivots.insert(col, row);
            return Ok(true);
        }
    }
    ensure(row[columns] == 0, "inconsistent relation system")?;
    Ok(false)
}

fn solve(pivots: &BTreeMap<usize, Vec<u64>>, order: u64, columns: usize) -> Check<Vec<u64>> {
    ensure(pivots.len() == columns, "rank not full")?;
    let mut logs = vec![0u64; columns];
    for col in (0..columns).rev() {
        let row = pivots.get(&col).ok_or("missing pivot")?;
        let mut rhs = row[columns];
        for j in col + 1..columns {
            rhs = mod_sub(rhs, mod_mul(row[j], logs[j], order), order);
        }
        logs[col] = rhs;
    }
    Ok(logs)
}

fn replay_rank(path: &Path, base: &Base, mut rank_seed: u64) -> Check<(Value, Vec<u64>)> {
    let input = File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut lines = BufReader::new(input).lines();
    let first = lines
        .next()
        .ok_or("empty rank trace")?
        .map_err(|e| e.to_string())?;
    let header: Value = serde_json::from_str(&first).map_err(|e| e.to_string())?;
    ensure(
        header["kind"] == "compact_orbit_rank_header",
        "missing rank header",
    )?;
    let columns = base.reps.len();
    let mut pivots = BTreeMap::new();
    let mut equations = Vec::new();
    let (mut attempts, mut relations, mut failures, mut without_gain) = (0u64, 0u64, 0u64, 0u64);
    let mut solution: Option<Value> = None;
    for (offset, line) in lines.enumerate() {
        let line = line.map_err(|e| format!("{}:{}: {e}", path.display(), offset + 2))?;
        if line.trim().is_empty() {
            continue;
        }
        let parsed: Value = serde_json::from_str(&line)
            .map_err(|e| format!("{}:{}: {e}", path.display(), offset + 2))?;
        let row = &parsed;
        if row["kind"] == "compact_orbit_rank_solution" {
            ensure(solution.is_none(), "duplicate rank solution")?;
            solution = Some(parsed);
            continue;
        }
        ensure(
            solution.is_none() && row["kind"] == "compact_orbit_rank_attempt",
            "rank trace kind/order mismatch",
        )?;
        ensure(
            integer(row, "attempt_index")? == attempts,
            "rank attempt index mismatch",
        )?;
        rank_seed = rank_seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let scalar = ((rank_seed >> 11) % (base.order - 1)) + 1;
        ensure(
            integer(row, "scalar")? == scalar,
            format!("attempt {attempts} scalar mismatch"),
        )?;
        let column = (0..columns)
            .find(|c| !pivots.contains_key(c))
            .ok_or("attempt after full rank")?;
        ensure(
            integer(row, "pivotless_column")? == column as u64,
            "pivotless column mismatch",
        )?;
        ensure(
            integer(row, "rank_before")? == pivots.len() as u64,
            "rank_before mismatch",
        )?;
        if row["found"] == true {
            let indices = numbers(row, "point_indices")?;
            ensure(
                indices.len() == 4 && indices.iter().all(|&i| (i as usize) < base.points.len()),
                "invalid rank summand indices",
            )?;
            let chosen: Vec<_> = indices.iter().map(|&i| base.points[i as usize]).collect();
            let target = point(&row["target"])?;
            ensure(
                base.curve.add(target, base.reps[column])?
                    == base.curve.scale(scalar, base.generator)?,
                "rank target group mismatch",
            )?;
            ensure(
                base.curve.sum(chosen.iter().copied())? == target,
                "rank summands do not add to target",
            )?;
            check_pair_roots(
                base.curve,
                &chosen,
                &numbers(row, "x_codes")?,
                &numbers(row, "pinned_intermediates")?,
            )?;
            let mut expected = vec![0u64; columns + 1];
            for i in indices {
                let (col, coeff) = base.labels[i as usize];
                expected[col] = mod_add(expected[col], coeff, base.order);
            }
            expected[column] = mod_add(expected[column], 1, base.order);
            expected[columns] = scalar;
            ensure(
                numbers(row, "row")? == expected,
                format!("attempt {attempts} row mismatch"),
            )?;
            let gained = insert_row(expected.clone(), &mut pivots, base.order, columns)?;
            ensure(
                row["gained"].as_bool() == Some(gained),
                "rank gain mismatch",
            )?;
            relations += 1;
            without_gain += u64::from(!gained);
            equations.push(expected);
        } else {
            ensure(row["found"] == false, "missing found status")?;
            failures += 1;
        }
        attempts += 1;
        ensure(
            integer(row, "rank_after")? == pivots.len() as u64,
            "rank_after mismatch",
        )?;
    }
    let solution = solution.as_ref().ok_or("missing rank solution")?;
    ensure(
        pivots.len() == columns && integer(solution, "rank")? == columns as u64,
        "incomplete rank",
    )?;
    for (field, value) in [
        ("attempts", attempts),
        ("relations", relations),
        ("failures", failures),
        ("rows_without_gain", without_gain),
    ] {
        ensure(
            integer(solution, field)? == value,
            format!("solution {field} mismatch"),
        )?;
    }
    let logs = solve(&pivots, base.order, columns)?;
    ensure(numbers(solution, "logs")? == logs, "solution logs mismatch")?;
    for row in &equations {
        let lhs = row
            .iter()
            .take(columns)
            .zip(logs.iter())
            .fold(0, |acc, (&a, &b)| {
                mod_add(acc, mod_mul(a, b, base.order), base.order)
            });
        ensure(lhs == row[columns], "solved equation mismatch")?;
    }
    for (col, rep) in base.reps.iter().enumerate() {
        ensure(
            base.curve.scale(logs[col], base.generator)? == *rep,
            format!("representative {col} log mismatch"),
        )?;
    }
    Ok((
        json!({"attempts":attempts,"relations":relations,"failures":failures,"rows_without_gain":without_gain,"rank":columns}),
        logs,
    ))
}

fn replay_target(
    path: &Path,
    base: &Base,
    logs: &[u64],
    workload: Option<&Value>,
) -> Check<(Point, u64)> {
    let row = one(path, "compact_orbit_dlp_target")?;
    let target = point(&row["published_q"])?;
    ensure(
        target == point(&row["target"])?,
        "published and internal target differ",
    )?;
    if let Some(w) = workload {
        ensure(
            target == point(&w["primary_target"])?,
            "workload target mismatch",
        )?;
    }
    let indices = numbers(&row, "point_indices")?;
    ensure(
        indices.len() == 4 && indices.iter().all(|&i| (i as usize) < base.points.len()),
        "invalid target summand indices",
    )?;
    let chosen: Vec<_> = indices.iter().map(|&i| base.points[i as usize]).collect();
    ensure(
        base.curve.sum(chosen.iter().copied())? == target,
        "target summands do not add to Q",
    )?;
    check_pair_roots(
        base.curve,
        &chosen,
        &numbers(&row, "x_codes")?,
        &numbers(&row, "pinned_intermediates")?,
    )?;
    let scalar = indices.iter().fold(0, |acc, &i| {
        let (col, coeff) = base.labels[i as usize];
        mod_add(acc, mod_mul(coeff, logs[col], base.order), base.order)
    });
    ensure(
        scalar == integer(&row, "recovered_scalar")?,
        "target recovered scalar mismatch",
    )?;
    ensure(
        base.curve.scale(scalar, base.generator)? == target,
        "target scalar replay failed",
    )?;
    ensure(
        row["group_verified"] == true,
        "producer target verification absent",
    )?;
    if let Some(w) = workload {
        ensure(
            scalar == integer(w, "verification_scalar")?,
            "workload scalar mismatch",
        )?;
    }
    Ok((target, scalar))
}

fn replay(run_dir: &Path, workload_path: Option<&Path>, require_receipt: bool) -> Check<Value> {
    let base_record = one(&run_dir.join("base.jsonl"), "point_defined_factor_base")?;
    let rank_header = one(&run_dir.join("rank.jsonl"), "compact_orbit_rank_header")?;
    let base = check_base(&base_record, &rank_header)?;
    let workload: Option<Value> = workload_path
        .map(|p| {
            fs::read(p)
                .map_err(|e| e.to_string())
                .and_then(|b| serde_json::from_slice(&b).map_err(|e| e.to_string()))
        })
        .transpose()?;
    if let Some(w) = &workload {
        for name in [
            "n",
            "a",
            "subgroup_order",
            "generator",
            "field_modulus_low_terms",
        ] {
            let rhs = if name == "generator" {
                rank_header[name].clone()
            } else {
                base_record[name].clone()
            };
            same(&w[name], &rhs, &format!("workload {name}"))?;
        }
    }
    let seed = workload
        .as_ref()
        .map(|w| integer(w, "rank_seed"))
        .transpose()?
        .unwrap_or(20261009);
    let (rank, logs) = replay_rank(&run_dir.join("rank.jsonl"), &base, seed)?;
    let (target, scalar) = replay_target(
        &run_dir.join("ic-target.jsonl"),
        &base,
        &logs,
        workload.as_ref(),
    )?;
    let summary = one(
        &run_dir.join("ic.stdout.jsonl"),
        "compact_orbit_dlp_summary",
    )?;
    for (name, value) in [
        ("rank", rank["rank"].clone()),
        ("rank_attempts", rank["attempts"].clone()),
        ("rank_relations", rank["relations"].clone()),
        ("rank_failures", rank["failures"].clone()),
        ("factor_base_points", json!(base.points.len())),
        ("base_hash", json!(base.hash)),
    ] {
        same(&summary[name], &value, &format!("summary {name}"))?;
    }
    if run_dir.join("rho.stdout.jsonl").exists() {
        let rho = one(&run_dir.join("rho.stdout.jsonl"), "rho_public_fixture")?;
        ensure(point(&rho["published_q"])? == target, "rho point mismatch")?;
        ensure(
            integer(&rho, "recovered_fixture_scalar")? == scalar && rho["verified"] == true,
            "rho scalar mismatch",
        )?;
        ensure(
            base.curve.scale(scalar, base.generator)? == target,
            "rho scalar replay failed",
        )?;
    }
    if require_receipt {
        let receipt: Value = serde_json::from_slice(
            &fs::read(run_dir.join("receipt.json")).map_err(|e| e.to_string())?,
        )
        .map_err(|e| e.to_string())?;
        ensure(
            receipt["ic"]["status"] == "success" && receipt["rho"]["status"] == "success",
            "arm receipt status mismatch",
        )?;
        ensure(
            receipt["analysis"]["status"] == "pending_independent_relation_replay",
            "analysis receipt status mismatch",
        )?;
        let files = receipt["files_sha256"]
            .as_object()
            .ok_or("missing file hashes")?;
        for (name, digest) in files {
            ensure(
                !name.contains('/') && !name.contains('\\'),
                "unsafe receipt filename",
            )?;
            ensure(
                sha256_file(&run_dir.join(name))? == digest.as_str().ok_or("invalid file hash")?,
                format!("{name} hash mismatch"),
            )?;
        }
    }
    Ok(
        json!({"status":"PASS","curve_n":base.curve.n,"subgroup_order":base.order,
        "factor_base_points":base.points.len(),"orbit_columns":base.reps.len(),"rank":rank,
        "target":target,"recovered_scalar":scalar,
        "base_sha256":sha256_file(&run_dir.join("base.jsonl"))?,
        "rank_trace_sha256":sha256_file(&run_dir.join("rank.jsonl"))?,
        "target_sha256":sha256_file(&run_dir.join("ic-target.jsonl"))?}),
    )
}

fn main() {
    if let Err(e) = entry() {
        eprintln!("native replay failed: {e}");
        std::process::exit(1);
    }
}

fn entry() -> Check<()> {
    let mut run_dir: Option<PathBuf> = None;
    let mut workload: Option<PathBuf> = None;
    let mut out: Option<PathBuf> = None;
    let mut require_receipt = false;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--run-dir" => {
                run_dir = Some(PathBuf::from(args.next().ok_or("missing --run-dir value")?))
            }
            "--workload" => {
                workload = Some(PathBuf::from(
                    args.next().ok_or("missing --workload value")?,
                ))
            }
            "--out" => out = Some(PathBuf::from(args.next().ok_or("missing --out value")?)),
            "--require-receipt" => require_receipt = true,
            _ => return Err(format!("unknown argument {arg}")),
        }
    }
    let result = replay(
        &run_dir.ok_or("--run-dir is required")?,
        workload.as_deref(),
        require_receipt,
    )?;
    let output = format!(
        "{}\n",
        serde_json::to_string_pretty(&result).map_err(|e| e.to_string())?
    );
    if let Some(path) = out {
        let mut file = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .map_err(|e| format!("{}: {e}", path.display()))?;
        use std::io::Write;
        file.write_all(output.as_bytes())
            .map_err(|e| e.to_string())?;
    } else {
        print!("{output}");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT: AtomicU64 = AtomicU64::new(0);

    struct Fixture(PathBuf);

    impl Fixture {
        fn copy() -> Self {
            let archived = Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("experiments/koblitz-n53-compact-cold-20261009/validation/n13-smoke");
            let dir = std::env::temp_dir().join(format!(
                "n13-native-replay-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            fs::create_dir(&dir).unwrap();
            for name in [
                "base.jsonl",
                "rank.jsonl",
                "ic-target.jsonl",
                "ic.stdout.jsonl",
                "rho.stdout.jsonl",
            ] {
                fs::copy(archived.join(name), dir.join(name)).unwrap();
            }
            Self(dir)
        }

        fn edit(&self, name: &str, index: usize, mutate: impl FnOnce(&mut Value)) {
            let path = self.0.join(name);
            let mut records = reference::rows(&path).unwrap();
            mutate(&mut records[index]);
            let mut content = records
                .iter()
                .map(|v| serde_json::to_string(v).unwrap())
                .collect::<Vec<_>>()
                .join("\n");
            content.push('\n');
            fs::write(path, content).unwrap();
        }
    }

    impl Drop for Fixture {
        fn drop(&mut self) {
            fs::remove_dir_all(&self.0).unwrap();
        }
    }

    #[test]
    fn archived_n13_certificate_matches() {
        let fixture = Fixture::copy();
        let actual = replay(&fixture.0, None, false).unwrap();
        let archived = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("experiments/koblitz-n53-compact-cold-20261009/validation/n13-smoke/replay.json");
        let expected: Value = serde_json::from_slice(&fs::read(archived).unwrap()).unwrap();
        assert_eq!(actual, expected);
    }

    #[test]
    fn changed_modular_row_is_rejected() {
        let fixture = Fixture::copy();
        fixture.edit("rank.jsonl", 1, |row| row["row"][0] = json!(1736));
        assert!(replay(&fixture.0, None, false).is_err());
    }

    #[test]
    fn changed_s3_root_is_rejected() {
        let fixture = Fixture::copy();
        fixture.edit("rank.jsonl", 1, |row| {
            row["pinned_intermediates"][0] = json!(8192)
        });
        assert!(replay(&fixture.0, None, false).is_err());
    }

    #[test]
    fn changed_target_scalar_is_rejected() {
        let fixture = Fixture::copy();
        fixture.edit("ic-target.jsonl", 0, |row| {
            row["recovered_scalar"] = json!(8)
        });
        assert!(replay(&fixture.0, None, false).is_err());
    }
}
