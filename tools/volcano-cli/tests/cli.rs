//! Exercise the actual binary, including an interrupted SQLite writer.
use rusqlite::Connection;
use serde_json::Value;
use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Stdio},
    time::{Duration, Instant},
};

struct Scratch(PathBuf);
impl Scratch {
    fn new(name: &str) -> Self {
        let path = std::env::temp_dir().join(format!("volcano-cli-{name}-{}", std::process::id()));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        if std::thread::panicking() {
            eprintln!("retained failed control: {}", self.0.display());
        } else {
            fs::remove_dir_all(&self.0).unwrap();
        }
    }
}
fn invoke(args: &[&str], code: i32) -> Value {
    let output = Command::new(env!("CARGO_BIN_EXE_isogeny-class"))
        .args(args)
        .output()
        .unwrap();
    assert_eq!(
        output.status.code(),
        Some(code),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).unwrap()
}
fn content(
    path: &Path,
) -> (
    Vec<(u64, u64, u64, u64, u64)>,
    Vec<(u64, u64, String, u64, u64)>,
) {
    let db = Connection::open(path).unwrap();
    let nodes = db
        .prepare("SELECT id,a,b,j,expanded FROM nodes ORDER BY id")
        .unwrap()
        .query_map([], |r| {
            Ok((r.get(0)?, r.get(1)?, r.get(2)?, r.get(3)?, r.get(4)?))
        })
        .unwrap()
        .collect::<std::result::Result<Vec<_>, _>>()
        .unwrap();
    let edges = db
        .prepare("SELECT source,ell,kernel,target,u FROM edges ORDER BY source,ell,kernel")
        .unwrap()
        .query_map([], |r| {
            Ok((r.get(0)?, r.get(1)?, r.get(2)?, r.get(3)?, r.get(4)?))
        })
        .unwrap()
        .collect::<std::result::Result<Vec<_>, _>>()
        .unwrap();
    (nodes, edges)
}
#[test]
fn killed_writer_recovers_and_resumes_to_the_uninterrupted_graph() {
    let dir = Scratch::new("kill");
    let input = dir.0.join("curve.json");
    fs::write(
        &input,
        include_bytes!("../../../research/full_volcano_classes_20261016/inputs/A.json"),
    )
    .unwrap();
    let interrupted = dir.0.join("interrupted.sqlite");
    let reference = dir.0.join("reference.sqlite");
    let stderr = fs::File::create(dir.0.join("writer.stderr")).unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_isogeny-class"))
        .args([
            "build",
            "--curve",
            input.to_str().unwrap(),
            "--out",
            interrupted.to_str().unwrap(),
            "--steps",
            "all",
        ])
        .stdout(Stdio::null())
        .stderr(stderr)
        .spawn()
        .unwrap();
    let deadline = Instant::now() + Duration::from_secs(20);
    loop {
        assert!(
            child.try_wait().unwrap().is_none(),
            "large writer exited before interruption"
        );
        let committed = if interrupted.is_file() {
            Connection::open_with_flags(&interrupted, rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY)
                .ok()
                .and_then(|db| {
                    db.busy_timeout(Duration::from_millis(1)).ok()?;
                    db.query_row::<u64, _, _>(
                        "SELECT count(*) FROM nodes WHERE expanded=1",
                        [],
                        |r| r.get(0),
                    )
                    .ok()
                })
                .unwrap_or(0)
        } else {
            0
        };
        if committed > 0 {
            break;
        }
        if Instant::now() > deadline {
            child.kill().unwrap();
            child.wait().unwrap();
            panic!("writer did not commit within correctness-control deadline");
        }
        std::thread::sleep(Duration::from_millis(5));
    }
    child.kill().unwrap();
    child.wait().unwrap();
    let state = invoke(&["status", "--input", interrupted.to_str().unwrap()], 0);
    assert_eq!(state["status"], "PARTIAL");
    let expanded = state["expanded_vertices"].as_u64().unwrap();
    assert!(expanded > 0);
    invoke(
        &[
            "resume",
            "--input",
            interrupted.to_str().unwrap(),
            "--steps",
            "8",
        ],
        3,
    );
    invoke(&["verify", "--input", interrupted.to_str().unwrap()], 3);
    invoke(
        &[
            "build",
            "--curve",
            input.to_str().unwrap(),
            "--out",
            reference.to_str().unwrap(),
            "--steps",
            &(expanded + 8).to_string(),
        ],
        3,
    );
    invoke(&["verify", "--input", reference.to_str().unwrap()], 3);
    assert_eq!(content(&interrupted), content(&reference));
}
#[test]
fn invalid_budget_and_existing_output_never_overwrite_files() {
    let dir = Scratch::new("usage");
    let input = dir.0.join("curve.json");
    fs::write(&input, b"{\"p\":113,\"a\":19,\"b\":2}").unwrap();
    let out = dir.0.join("graph.sqlite");
    invoke(
        &[
            "build",
            "--curve",
            input.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
            "--steps",
            "bad",
        ],
        2,
    );
    assert!(!out.exists());
    fs::write(&out, b"preserved existing bytes").unwrap();
    invoke(
        &[
            "build",
            "--curve",
            input.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
            "--steps",
            "all",
        ],
        2,
    );
    assert_eq!(fs::read(out).unwrap(), b"preserved existing bytes");
}
