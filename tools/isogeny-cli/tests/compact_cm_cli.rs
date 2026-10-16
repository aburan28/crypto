use serde_json::{json, Value};
use std::{fs, path::Path, process::Command};
fn call(args: &[&str], code: i32) -> Value {
    let out = Command::new(env!("CARGO_BIN_EXE_isogeny"))
        .args(args)
        .current_dir(std::env::temp_dir())
        .output()
        .unwrap();
    assert_eq!(
        out.status.code(),
        Some(code),
        "{}\n{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    serde_json::from_slice(&out.stdout).unwrap()
}
fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}
#[test]
fn installed_compact_command_retains_controls_and_rejects_forged_certificates() {
    let retained = std::env::var_os("ISOGENY_COMPACT_EVIDENCE");
    let root = retained
        .as_ref()
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| {
            std::env::temp_dir().join(format!(
                "compact-cm-{}-{}",
                std::process::id(),
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .as_nanos()
            ))
        });
    fs::create_dir(&root).unwrap();
    let root = fs::canonicalize(root).unwrap();
    let construct = root.join("positive");
    assert_eq!(
        call(
            &[
                "compact-cm",
                "construct",
                "--exponent",
                "3",
                "--out",
                s(&construct)
            ],
            0
        )["status"],
        "CONSTRUCTED"
    );
    let path = construct.join("certificate.json");
    let before = fs::read(&path).unwrap();
    let c: Value = serde_json::from_slice(&before).unwrap();
    let verified = call(
        &[
            "compact-cm",
            "verify",
            "--input",
            s(&path),
            "--out",
            s(&root.join("replay")),
        ],
        0,
    );
    assert_eq!(verified["geometric_prime_degree"], "7");
    assert_eq!(verified["resources"]["status"], "PASS");
    assert!(
        verified["resources"]["sampled_peak_resident_bytes"]
            .as_str()
            .unwrap()
            .parse::<u64>()
            .unwrap()
            > 0
    );
    assert_eq!(
        verified["resources"]["resident_memory_limit_bytes"],
        "8589934592"
    );

    assert_eq!(verified["public_point_checks"].as_array().unwrap().len(), 8);
    assert!(verified["expanded_map_challenge_score"].is_null());
    call(
        &[
            "compact-cm",
            "construct",
            "--exponent",
            "3",
            "--out",
            s(&construct),
        ],
        2,
    );
    assert_eq!(fs::read(&path).unwrap(), before);
    for exponent in ["11", "9", "2"] {
        let out = root.join(format!("negative-{exponent}"));
        assert_eq!(
            call(
                &[
                    "compact-cm",
                    "construct",
                    "--exponent",
                    exponent,
                    "--out",
                    s(&out)
                ],
                1
            )["status"],
            "REJECTED"
        );
    }
    let mutations = [
        ("degree", json!("31")),
        ("a", json!("4")),
        ("b", json!("3")),
        ("beta", json!("1")),
        ("differential", json!("0")),
        ("subgroup_scalar", json!("0")),
        ("subgroup_inverse", json!("0")),
    ];
    for (key, value) in mutations {
        let mut bad = c.clone();
        bad[key] = value;
        let input = root.join(format!("mutated-{key}.json"));
        fs::write(&input, serde_json::to_vec(&bad).unwrap()).unwrap();
        let out = root.join(format!("judged-{key}"));
        assert_eq!(
            call(
                &[
                    "compact-cm",
                    "verify",
                    "--input",
                    s(&input),
                    "--out",
                    s(&out)
                ],
                1
            )["status"],
            "REJECTED"
        );
    }
    for key in ["terminal", "transitions", "states_sha256"] {
        let mut bad = c.clone();
        bad["primality"][key] = if key == "transitions" {
            json!(2)
        } else {
            json!("1")
        };
        let input = root.join(format!("ll-{key}.json"));
        fs::write(&input, serde_json::to_vec(&bad).unwrap()).unwrap();
        call(
            &[
                "compact-cm",
                "verify",
                "--input",
                s(&input),
                "--out",
                s(&root.join(format!("judged-ll-{key}"))),
            ],
            1,
        );
    }
    let gx = c["domain"]["generator"][0].as_str().unwrap();
    let gy = c["domain"]["generator"][1].as_str().unwrap();
    let direct = call(
        &[
            "compact-cm",
            "evaluate",
            "--input",
            s(&path),
            "--x",
            gx,
            "--y",
            gy,
            "--mode",
            "geometric",
        ],
        0,
    );
    let reduced = call(
        &[
            "compact-cm",
            "evaluate",
            "--input",
            s(&path),
            "--x",
            gx,
            "--y",
            gy,
            "--mode",
            "subgroup",
        ],
        0,
    );
    assert_eq!(direct["point"], reduced["point"]);
    call(
        &[
            "compact-cm",
            "evaluate",
            "--input",
            s(&path),
            "--x",
            "0",
            "--y",
            "0",
            "--mode",
            "geometric",
        ],
        2,
    );
    if retained.is_none() {
        fs::remove_dir_all(root).unwrap();
    }
}
