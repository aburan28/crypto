//! Native replay of the disclosed n17 SAT source control, without a solver.
//! Source pins are historical evidence, never a new registration or yield sample.
#[cfg(unix)]
use std::os::unix::fs::OpenOptionsExt;
use std::{
    fs,
    io::{Read, Write},
    path::Path,
};

use serde_json::{json, Value};

use super::{identity, json as records, oracle::Curve};

const PREPARATION_SHA: &str = "91856ab78550436d3f668367f9aebd9e2c0604bd64b1472d9d19ec318e2b144e";
const SOURCE_DIR: &str = "research/ic_candidate_tournament_20260915/goal_20260924/prepared-report-contract-v1/sat-source-check";
const PREPARATION: &str = "research/ic_candidate_tournament_20260915/goal_20260924/prepared-ic-state-v1/sat-preparation.json";
const PINS: &[(&str, usize, &str)] = &[
    (
        "manifest.json",
        4034,
        "9f516b3b6e1600a63ae42bd6e06342c3a8539c8f84d5bccc86b62d11c500822b",
    ),
    (
        "instance.anf",
        28877,
        "e5d49feddef960a0ea8a74836f517563a4f63f987abfc45defbe36a337c83057",
    ),
    (
        "instance.xor.cnf",
        39608,
        "d46c9d1edea83cb1ea2d09ce6cab4f42ff15bce99263b676430e61f4930ab3a0",
    ),
    (
        "instance.magma",
        46815,
        "30c623c4366c36bea854ab4ac08ff0deed0da4c37347395a3c442f525f1e4803",
    ),
];

fn require(condition: bool, message: &str) -> Result<(), String> {
    if condition {
        Ok(())
    } else {
        Err(message.into())
    }
}

fn check_pin(bytes: &[u8], length: usize, digest: &str) -> Result<(), String> {
    require(
        bytes.len() == length && identity::sha256_hex(bytes) == digest,
        "retained source bytes differ from accepted pin",
    )
}

fn read_regular(path: &Path, max_bytes: u64) -> Result<Vec<u8>, String> {
    let meta = fs::symlink_metadata(path).map_err(|e| format!("{}: {e}", path.display()))?;
    require(
        meta.is_file() && !meta.file_type().is_symlink() && meta.len() <= max_bytes,
        "input must be a bounded regular file, not a symlink",
    )?;
    let mut options = fs::OpenOptions::new();
    options.read(true);
    #[cfg(unix)]
    options.custom_flags(libc::O_NOFOLLOW);
    let file = options.open(path).map_err(|e| e.to_string())?;
    require(
        file.metadata().map_err(|e| e.to_string())?.is_file(),
        "input is not a regular file",
    )?;
    let mut bytes = Vec::new();
    file.take(max_bytes + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| e.to_string())?;
    require(
        bytes.len() as u64 <= max_bytes,
        "input exceeded its byte limit",
    )?;
    Ok(bytes)
}

fn variable(token: &str, variables: usize) -> Result<usize, String> {
    let index = token.parse::<usize>().map_err(|_| "invalid variable")?;
    require(index > 0 && index <= variables, "variable outside model")?;
    Ok(index - 1)
}

fn header(line: &str) -> Result<(usize, usize), String> {
    let tokens = line.split_whitespace().collect::<Vec<_>>();
    require(
        tokens.len() == 4 && tokens[..2] == ["p", "cnf"],
        "invalid source header",
    )?;
    let n = tokens[2]
        .parse::<usize>()
        .map_err(|_| "invalid variable count")?;
    let m = tokens[3]
        .parse::<usize>()
        .map_err(|_| "invalid row count")?;
    require(
        n > 0 && n <= 4096 && m > 0 && m <= 100_000,
        "source size outside replay envelope",
    )?;
    Ok((n, m))
}

/// WDSAT convention: a row has constant one exactly when T is absent.
fn verify_anf(text: &str, model: &[bool]) -> Result<usize, String> {
    let mut lines = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('c'));
    let (n, expected) = header(lines.next().ok_or("missing ANF header")?)?;
    require(model.len() == n, "ANF model width differs")?;
    let mut rows = 0;
    for line in lines {
        let tokens = line.split_whitespace().collect::<Vec<_>>();
        require(
            tokens.len() >= 2 && tokens[0] == "x" && tokens.last() == Some(&"0"),
            "malformed ANF row",
        )?;
        let terms = &tokens[1..tokens.len() - 1];
        let (mut at, mut parity, mut has_t) = (0, false, false);
        while at < terms.len() {
            match terms[at] {
                "T" => {
                    require(!has_t, "duplicate ANF T marker")?;
                    has_t = true;
                    at += 1;
                }
                token if token.starts_with('.') => {
                    let degree = token[1..]
                        .parse::<usize>()
                        .map_err(|_| "bad monomial degree")?;
                    require(
                        degree > 0 && degree <= n && degree < terms.len() - at,
                        "truncated or oversized monomial",
                    )?;
                    let mut value = true;
                    for term in &terms[at + 1..at + degree + 1] {
                        // Do not short circuit: even a false product validates every variable.
                        value &= model[variable(term, n)?];
                    }
                    parity ^= value;
                    at += degree + 1;
                }
                token => {
                    parity ^= model[variable(token, n)?];
                    at += 1;
                }
            }
        }
        require(!(parity ^ !has_t), "ANF equation is not satisfied")?;
        rows += 1;
    }
    require(rows == expected, "ANF row count differs")?;
    Ok(rows)
}

#[derive(Clone)]
struct Row {
    xor: bool,
    literals: Vec<i32>,
}

fn parse_cnf(text: &str) -> Result<(usize, Vec<Row>), String> {
    let mut lines = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('c'));
    let (variables, expected) = header(lines.next().ok_or("missing CNF header")?)?;
    let mut rows = Vec::new();
    for line in lines {
        let tokens = line.split_whitespace().collect::<Vec<_>>();
        let xor = tokens.first() == Some(&"x");
        let terms = &tokens[usize::from(xor)..];
        require(
            !terms.is_empty() && terms.last() == Some(&"0"),
            "unterminated CNF row",
        )?;
        let literals = terms[..terms.len() - 1]
            .iter()
            .map(|token| {
                let literal = token.parse::<i32>().map_err(|_| "bad CNF literal")?;
                require(
                    literal != 0 && literal.unsigned_abs() as usize <= variables,
                    "CNF literal outside model",
                )?;
                Ok(literal)
            })
            .collect::<Result<Vec<_>, String>>()?;
        rows.push(Row { xor, literals });
        require(rows.len() <= expected, "too many CNF rows")?;
    }
    require(rows.len() == expected, "CNF row count differs")?;
    Ok((variables, rows))
}

fn expanded_model(source: &[bool], variables: usize, rows: &[Row]) -> Result<Vec<bool>, String> {
    require(
        variables >= source.len(),
        "source wider than expanded model",
    )?;
    let mut model = source.to_vec();
    model.resize(variables, false);
    let mut defined = vec![false; variables - source.len()];
    for row in rows.iter().filter(|r| !r.xor) {
        let Some((&output, arguments)) = row.literals.split_last() else {
            continue;
        };
        if arguments.len() >= 2
            && output > source.len() as i32
            && arguments
                .iter()
                .all(|&a| a < 0 && a.unsigned_abs() as usize <= source.len())
        {
            let output = output as usize - 1;
            require(
                output < variables && !defined[output - source.len()],
                "duplicate or out-of-range auxiliary AND definition",
            )?;
            model[output] = arguments
                .iter()
                .all(|a| source[a.unsigned_abs() as usize - 1]);
            defined[output - source.len()] = true;
        }
    }
    require(
        defined.iter().all(|&d| d),
        "missing auxiliary AND definition",
    )?;
    Ok(model)
}

fn verify_cnf(rows: &[Row], model: &[bool]) -> Result<(usize, usize), String> {
    let (mut cnf, mut xor) = (0, 0);
    for row in rows {
        let mut parity = false;
        let mut any = false;
        for &literal in &row.literals {
            let value = *model
                .get(literal.unsigned_abs() as usize - 1)
                .ok_or("model too short")?
                ^ (literal < 0);
            parity ^= value;
            any |= value;
        }
        require(
            if row.xor { parity } else { any },
            "CNF/XOR row is not satisfied",
        )?;
        if row.xor {
            xor += 1
        } else {
            cnf += 1
        }
    }
    Ok((cnf, xor))
}

/// Independently validate an expanded native result; no shared producer parser.
pub fn verify_native_model(anf: &str, cnf: &str, stdout: &str) -> Result<Vec<bool>, String> {
    let (variables, rows) = parse_cnf(cnf)?;
    require(
        variables >= 51,
        "expanded source is narrower than n17 model",
    )?;
    let mut model = vec![None; variables];
    let mut terminated = false;
    for line in stdout.lines().filter_map(|line| line.strip_prefix("v ")) {
        for token in line.split_whitespace() {
            require(!terminated, "model continues after terminator")?;
            let number = token
                .parse::<i32>()
                .map_err(|_| "bad native model literal")?;
            if number == 0 {
                terminated = true;
                continue;
            }
            let index = number.unsigned_abs() as usize;
            require(
                index > 0 && index <= variables,
                "model literal outside source",
            )?;
            require(model[index - 1].is_none(), "model repeats a variable")?;
            model[index - 1] = Some(number > 0);
        }
    }
    require(terminated, "native model lacks terminator")?;
    let model = model
        .into_iter()
        .map(|v| v.ok_or("native model is incomplete".into()))
        .collect::<Result<Vec<_>, String>>()?;
    require(
        verify_anf(anf, &model[..51])? == 50,
        "source ANF equation count differs",
    )?;
    let (_, xor) = verify_cnf(&rows, &model)?;
    require(xor == 50, "source XOR equation count differs")?;
    Ok(model)
}

/// Substitute the accepted witness; never consume a registration or invoke a solver.
pub fn replay(root: &Path) -> Result<Value, String> {
    // Bind before arithmetic, and reject a rebuild/replacement during replay.
    // This is a local executable gate, not full scientific-runtime admission.
    let executable = std::env::current_exe().map_err(|e| e.to_string())?;
    let binary_sha256 = identity::sha256_hex(&fs::read(&executable).map_err(|e| e.to_string())?);
    let preparation_bytes = read_regular(&root.join(PREPARATION), 2_000_000)?;
    let prep_text = std::str::from_utf8(&preparation_bytes).map_err(|e| e.to_string())?;
    let prep = records::parse(prep_text)?;
    require(
        identity::sha256(&prep)? == PREPARATION_SHA,
        "preparation seal differs",
    )?;
    let preparation_admission = super::prepared_sat::verify(&prep)?;
    let inputs = prep.at("certificate")?.at("inputs")?;
    let curve = Curve::new(inputs.at("fixture")?)?;
    require(
        curve.n == 17 && curve.a == 1 && curve.r == 65587 && curve.modulus == 131081,
        "only the retained toy curve is supported",
    )?;
    let base = inputs.at("base")?.as_arr().ok_or("base must be an array")?;
    require(base.len() == 63, "geometric base count differs")?;
    let points = [29, 51, 2].map(|i| curve.decode(&base[i]));
    let points = points.into_iter().collect::<Result<Vec<_>, _>>()?;
    let target = Some((62577, 27783));
    require(
        curve.add(curve.add(points[0], points[1]), points[2]) == target,
        "retained witness does not readd to query",
    )?;
    let xs = points
        .iter()
        .map(|p| p.map(|p| p.0).ok_or("identity witness point"))
        .collect::<Result<Vec<_>, _>>()?;
    let symmetric = [
        xs[0] ^ xs[1] ^ xs[2],
        curve.fm(xs[0], xs[1]) ^ curve.fm(xs[0], xs[2]) ^ curve.fm(xs[1], xs[2]),
        curve.fm(curve.fm(xs[0], xs[1]), xs[2]),
    ];
    let mut source = Vec::new();
    for (value, width) in xs
        .iter()
        .copied()
        .chain(symmetric)
        .zip([6, 6, 6, 6, 11, 16])
    {
        require(value < 1 << width, "source coefficient overflow")?;
        source.extend((0..width).map(|bit| value & (1 << bit) != 0));
    }
    let dir = root.join(SOURCE_DIR);
    let mut files = Vec::new();
    let mut inventory = serde_json::Map::new();
    for &(name, length, digest) in PINS {
        let bytes = read_regular(&dir.join(name), length as u64)?;
        check_pin(&bytes, length, digest)?;
        inventory.insert(name.into(), json!({"bytes":length,"sha256":digest}));
        files.push(bytes);
    }
    let manifest: Value = serde_json::from_slice(&files[0]).map_err(|e| e.to_string())?;
    require(
        manifest["n"] == 17
            && manifest["curve_a"] == 1
            && manifest["ell"] == 6
            && manifest["m"] == 3
            && manifest["representation"] == "symmetrised_s4"
            && manifest["source_variables"] == 51
            && manifest["source_equations"] == 50
            && manifest["factor_base_basis_bitmasks"] == json!(["1", "2", "4", "8", "16", "32"])
            && manifest["target"] == json!({"x":"62577","y":"27783"}),
        "retained source layout differs",
    )?;
    for (name, index) in [
        ("wdsat_anf", 1),
        ("cryptominisat_xor_dimacs", 2),
        ("magma_boolean_f4", 3),
    ] {
        let descriptor = &manifest["exports"][name];
        require(
            descriptor["path"] == PINS[index].0
                && descriptor["bytes"] == files[index].len()
                && descriptor["blake3"] == blake3::hash(&files[index]).to_hex().as_str(),
            "export descriptor differs from native source hash",
        )?;
    }
    let anf = std::str::from_utf8(&files[1]).map_err(|e| e.to_string())?;
    let cnf = std::str::from_utf8(&files[2]).map_err(|e| e.to_string())?;
    let equations = verify_anf(anf, &source)?;
    let (variables, rows) = parse_cnf(cnf)?;
    require(variables == 767, "expanded variable count differs")?;
    let model = expanded_model(&source, variables, &rows)?;
    let (clauses, xors) = verify_cnf(&rows, &model)?;
    require(
        (equations, clauses, xors) == (50, 2364, 50),
        "retained source counts differ",
    )?;
    let digest_bits = |bits: &[bool]| {
        identity::sha256_hex(&bits.iter().map(|&b| u8::from(b)).collect::<Vec<_>>())
    };
    require(
        identity::sha256_hex(&fs::read(&executable).map_err(|e| e.to_string())?) == binary_sha256,
        "checker executable changed during replay",
    )?;
    Ok(json!({
        "schema_version":1,"status":"PASS_NATIVE_RETAINED_SAT_SOURCE_REPLAY",
        "implementation":"Rust icprog; independent shift-and-add curve checker",
        "checker_source_sha256":identity::sha256_hex(include_bytes!("sat_source.rs")),
        "oracle_source_sha256":identity::sha256_hex(include_bytes!("oracle.rs")),
        "checker_binary_sha256":binary_sha256,
        "checker_binary_unchanged_before_after":true,
        "host_architecture":std::env::consts::ARCH,"host_os":std::env::consts::OS,
        "preparation_sha256":PREPARATION_SHA,
        "preparation_admission":preparation_admission,
        "historical_archive_reopened":false,
        "source_custody":"committed byte-exact copies checked against accepted pins; original archive is not reopened",
        "historical_archive_sha256":"c94aba5c67afbe60b5109169d2d57c1c44ce89b218c63d7b023ed023bdb64b24",
        "historical_execution_sha256":"43539f7d440289dae1ad4867ba1ca951bd4664bf01147b9070eaa68f3bad07cb",
        "source_files":inventory,"query_point":[62577,27783],"witness_indices":[29,51,2],
        "witness_points":points,"group_replay":true,
        "source_assignment_sha256":digest_bits(&source),"expanded_assignment_sha256":digest_bits(&model),
        "source_variables":source.len(),"auxiliary_and_definitions":variables-source.len(),
        "anf_equations":equations,"cnf_clauses":clauses,"xor_rows":xors,
        "native_solvers_executed":0,"fresh_targets_generated":0,
        "original_native_outcome":"CONFLICT_BUDGET_INCONCLUSIVE",
        "source_bound_scientific_runtime_admitted":false,"original_ic_target_complete":false,
        "candidate_id":null,"online_wall_ns":null,"online_speedup":null,"promotion_eligible":false,
        "scope":"native postexecution substitution of one disclosed witness; no solver search or natural-yield estimate"
    }))
}

pub fn run(root: &Path, out: Option<&Path>) -> Result<String, String> {
    let report = replay(root)?;
    let text = serde_json::to_string_pretty(&report).map_err(|e| e.to_string())? + "\n";
    if let Some(path) = out {
        let mut file = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)
            .map_err(|e| format!("create immutable {}: {e}", path.display()))?;
        file.write_all(text.as_bytes()).map_err(|e| e.to_string())?;
        file.sync_all().map_err(|e| e.to_string())?;
    }
    Ok(text)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_length_mutation_and_truncated_sources_fail_the_pin() {
        let bytes=include_bytes!("../../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-report-contract-v1/sat-source-check/instance.anf");
        let (_, length, digest) = PINS[1];
        assert!(check_pin(bytes, length, digest).is_ok());
        let mut changed = bytes.to_vec();
        changed[100] ^= 1;
        assert!(check_pin(&changed, length, digest).is_err());
        assert!(check_pin(&bytes[..bytes.len() - 1], length, digest).is_err());
    }

    #[test]
    fn retained_model_corruption_is_rejected_by_source_rows() {
        let historical: Value=serde_json::from_str(include_str!("../../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-report-contract-v1/sat-source-check/result.json")).unwrap();
        let mut model = historical["cnf_assignment"]
            .as_array()
            .unwrap()
            .iter()
            .map(|b| b.as_bool().unwrap())
            .collect::<Vec<_>>();
        let (_,rows)=parse_cnf(include_str!("../../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-report-contract-v1/sat-source-check/instance.xor.cnf")).unwrap();
        assert!(verify_cnf(&rows, &model).is_ok());
        for index in 0..model.len() {
            model[index] = !model[index];
            assert!(
                verify_cnf(&rows, &model).is_err(),
                "accepted corrupted bit {index}"
            );
            model[index] = !model[index];
        }
    }

    #[test]
    fn retained_native_source_replay_matches_independent_historical_assignment() {
        let report = replay(Path::new(env!("CARGO_MANIFEST_DIR"))).unwrap();
        assert_eq!(
            report["source_assignment_sha256"],
            "054d401221bdfc267c17817d5faad7a1cebc9c6d9ea26214695ec57c3a98e650"
        );
        assert_eq!(
            report["expanded_assignment_sha256"],
            "bac1e84db816e7547c555eaa2ac60a79d55769d8bfe8661149fc67f122409dc2"
        );
        assert_eq!(report["native_solvers_executed"], 0);
        assert_eq!(report["promotion_eligible"], false);
    }

    #[test]
    fn anf_constant_and_monomial_semantics_and_truncation() {
        assert!(verify_anf("p cnf 2 1\nx .2 1 2 0\n", &[true, true]).is_ok());
        assert!(verify_anf("p cnf 2 1\nx .2 1 2 T 0\n", &[true, false]).is_ok());
        for text in [
            "p cnf 2 1\nx .2 1 0",
            "p cnf 2 1\nx 1 T T 0",
            "p cnf 2 1\nx .2 1 3 T 0",
            "p cnf 2 1\nx 0 T 0",
            "p cnf 2 2\nx 1 T 0",
            "p cnf 2 1\nx 1 T",
        ] {
            assert!(
                verify_anf(text, &[false, false]).is_err(),
                "accepted {text}"
            );
        }
    }

    #[test]
    fn cnf_rejects_invalid_headers_literals_counts_and_terminators() {
        for text in [
            "1 0",
            "p cnf 2 1\n3 0",
            "p cnf 2 1\n-2147483648 0",
            "p cnf 2 1\n1 0 2 0",
            "p cnf 2 1\n1",
            "p cnf 2 2\n1 0",
            "p cnf 2 1\n1 0\n2 0",
            "p cnf 2 1\np cnf 2 1",
        ] {
            assert!(parse_cnf(text).is_err(), "accepted {text}");
        }
    }

    #[test]
    fn signed_xors_and_empty_rows_are_checked() {
        let (_, rows) = parse_cnf("p cnf 2 2\nx -1 2 0\n1 0").unwrap();
        assert!(verify_cnf(&rows, &[true, true]).is_ok());
        assert!(verify_cnf(&rows, &[true, false]).is_err());
        for line in ["0", "x 0"] {
            let (_, rows) = parse_cnf(&format!("p cnf 2 1\n{line}")).unwrap();
            assert!(verify_cnf(&rows, &[true, false]).is_err());
        }
    }

    #[test]
    fn auxiliary_definitions_must_be_complete_unique_and_satisfy_reverse_clauses() {
        let (_, rows) = parse_cnf("p cnf 3 3\n-1 -2 3 0\n1 -3 0\n2 -3 0").unwrap();
        let model = expanded_model(&[true, false], 3, &rows).unwrap();
        assert_eq!(model, vec![true, false, false]);
        assert!(verify_cnf(&rows, &model).is_ok());
        assert!(verify_cnf(&rows, &[true, false, true]).is_err());
        assert!(expanded_model(&[true, false], 4, &rows).is_err());
        let mut duplicate = rows.clone();
        duplicate.push(rows[0].clone());
        assert!(expanded_model(&[true, false], 3, &duplicate).is_err());
    }
}
