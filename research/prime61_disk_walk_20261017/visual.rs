//! Render recorded capacity diagnostics only; no reruns or timing inference.
use serde_json::{json, Value};
use std::{fs, path::Path};
fn last(path: &Path) -> Value {
    let t = fs::read_to_string(path).unwrap();
    serde_json::from_str(t.lines().last().unwrap()).unwrap()
}
fn main() {
    let a: Vec<String> = std::env::args().collect();
    assert_eq!(a.len(), 3, "PILOT_DIRECTORY OUTPUT_DIRECTORY");
    let input = Path::new(&a[1]);
    let out = Path::new(&a[2]);
    fs::create_dir_all(out).unwrap();
    let mut stages = Vec::new();
    let mut peak = 0;
    for n in [1, 8, 32] {
        let v = last(&input.join(format!("expand-{n}.stdout")));
        assert_eq!(v["expanded"], n);
        assert_eq!(v["failures"], 0);
        peak = peak.max(v["process_peak_rss_kib"].as_u64().unwrap());
        stages.push(v);
    }
    for phase in [
        "init",
        "measure-1",
        "measure-8",
        "measure-32",
        "verify-1",
        "verify-8",
        "verify-32",
        "snapshot-1",
        "snapshot-8",
        "snapshot-32",
        "final-status",
    ] {
        let v = last(&input.join(format!("{phase}.stdout")));
        peak = peak.max(v["process_peak_rss_kib"].as_u64().unwrap());
    }
    let status = last(&input.join("final-status.stdout"));
    assert_eq!(status["discovered"], 250);
    assert_eq!(status["edges"], 480);
    assert_eq!(status["runs"], 768);
    assert_eq!(status["target_trials"], 76800000);
    assert_eq!(status["relations"], 0);
    assert_eq!(status["verified_solves"], 0);
    let summary = json!({"schema":"prime61.disk-walk-summary/v1","status":status,"expansion_checkpoints":stages,"maximum_recorded_worker_peak_rss_kib":peak,"full_class_exhausted":false,"old_panel_workloads_matched":16,"old_panel_run_rows_matched":192,"cost_scope":"native relation-stage counters and worker capacity diagnostics; no admitted wall performance or full cost"});
    fs::write(
        out.join("summary.json"),
        serde_json::to_string_pretty(&summary).unwrap() + "\n",
    )
    .unwrap();
    let mut svg = String::from(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="650" viewBox="0 0 1000 650"><rect width="1000" height="650" fill="#f8fafc"/><g font-family="DejaVu Sans,sans-serif" fill="#16304a"><text x="28" y="35" font-size="25" font-weight="bold">Disk-backed prime61 walk: checkpointed capacity pilot</text><text x="28" y="62" font-size="15">2026-10-17 • single native worker • three cumulative checkpoints • full class still open</text><rect x="28" y="82" width="944" height="82" rx="7" fill="#e1edf5"/><text x="45" y="112" font-size="17" font-weight="bold">One source in RAM → SQLite frontier / witnesses / runs → compressed snapshot → Drive</text><text x="45" y="143" font-size="15">Atomic expansion and run commits; exact j deduplication; hash-checked restore and resume.</text><text x="28" y="200" font-size="18" font-weight="bold">Discovered curves (count)</text><text x="528" y="200" font-size="18" font-weight="bold">Expansion-process peak RSS (KiB)</text>"##,
    );
    for (i, v) in stages.iter().enumerate() {
        let y = 225 + i * 76;
        let expanded = v["expanded"].as_u64().unwrap();
        let nodes = v["discovered"].as_u64().unwrap();
        let rss = v["process_peak_rss_kib"].as_u64().unwrap();
        svg.push_str(&format!(r##"<text x="28" y="{}" font-size="14">{expanded} expanded</text><rect x="145" y="{y}" width="{}" height="24" fill="#367ba5"/><text x="{}" y="{}" font-size="14">{nodes}</text><text x="528" y="{}" font-size="14">{expanded} expanded</text><rect x="650" y="{y}" width="{}" height="24" fill="#638c72"/><text x="905" y="{}" font-size="14">{rss}</text>"##,y+17,nodes*280/250,153+nodes*280/250,y+17,y+17,rss*230/20000,y+17));
    }
    svg.push_str(&format!(r##"<text x="28" y="469" font-size="14">32 expanded • 250 discovered • 480 verified edges • 218 pending • 0 failures</text><rect x="28" y="494" width="944" height="86" rx="7" fill="#fff0d1"/><text x="44" y="522" font-size="17" font-weight="bold">64 measured curves × 12 streams = 768 capped runs</text><text x="44" y="549" font-size="15">76,800,000 target trials • 0 relations • 0 solves • peak worker RSS across phases: {peak} KiB</text><text x="28" y="608" font-size="13">Capacity diagnostics, not a speedup comparison. SQLite retains the graph on disk; disk use still grows.</text><text x="28" y="633" font-size="12">Sources: results/pilot/expand-*.stdout and final-status.stdout; restore control and native baseline audit.</text></g></svg>"##));
    fs::write(out.join("capacity.svg"), svg).unwrap();
}
