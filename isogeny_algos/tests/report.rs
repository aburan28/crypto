/// The `report` binary regenerates the committed Markdown summaries byte for byte (they were
/// produced by the earlier Python script, which it replaces).
#[test]
fn report_regenerates_committed_tables() {
    let bin = env!("CARGO_BIN_EXE_report");
    for name in ["baseline", "baseline-v1-all", "baseline-v2"] {
        let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/results/");
        let out = std::process::Command::new(bin)
            .arg(format!("{dir}{name}.jsonl"))
            .output()
            .unwrap();
        assert!(out.status.success());
        let want = std::fs::read_to_string(format!("{dir}{name}.md")).unwrap();
        assert_eq!(String::from_utf8(out.stdout).unwrap(), want, "{name}.md");
    }
}
