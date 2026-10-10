//! Native cap for the existing CM order-support preflight; no class label on failure.
use crypto_lib::hash::sha256::sha256;
use std::{
    fs,
    path::Path,
    process::{Command, Stdio},
    thread,
    time::{Duration, Instant},
};
fn hash(p: &Path) -> String {
    hex::encode(sha256(&fs::read(p).unwrap()))
}
fn main() {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    assert_eq!(args.len(), 1);
    let root = Path::new(&args[0]);
    let raw = root.join("cm_preflight.stdout");
    let err = root.join("cm_preflight.stderr");
    assert!(!raw.exists() && !err.exists());
    let now = Instant::now();
    let mut child = Command::new("/opt/homebrew/bin/gp")
        .args(["-q", "-f", "-s", "256M"])
        .arg(root.join("cm_preflight.gp"))
        .stdout(Stdio::from(fs::File::create(&raw).unwrap()))
        .stderr(Stdio::from(fs::File::create(&err).unwrap()))
        .spawn()
        .unwrap();
    let (code, cap) = loop {
        if let Some(c) = child.try_wait().unwrap() {
            break (c.code().unwrap_or(-1), false);
        }
        if now.elapsed() >= Duration::from_secs(60) {
            child.kill().unwrap();
            break (child.wait().unwrap().code().unwrap_or(-1), true);
        }
        thread::sleep(Duration::from_millis(100));
    };
    let status = if cap {
        "TIMEOUT_UNRESOLVED"
    } else if fs::metadata(&err).unwrap().len() > 0 || code != 0 {
        "IMPLEMENTATION_ERROR_UNRESOLVED"
    } else {
        "FIRST_ORDER_COMPLETE_CLASS_UNRESOLVED"
    };
    fs::write(root.join("cm_preflight_receipt.json"),serde_json::to_string_pretty(&serde_json::json!({"status":status,"cap_seconds":60,"exit_code":code,"wall_ms":now.elapsed().as_millis(),"script_sha256":hash(&root.join("cm_preflight.gp")),"driver_sha256":hash(&root.join("probe.rs")),"binary_sha256":hash(&std::env::current_exe().unwrap()),"stdout_sha256":hash(&raw),"stderr_sha256":hash(&err),"class_label":"UNRESOLVED"})).unwrap()+"\n").unwrap();
    println!("{status}");
}
