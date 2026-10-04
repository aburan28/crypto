//! Frozen target-blind n37 folded-table rank producer. The matching
//! independent general-law replay is `n37_shared_rank_replay`.

use crypto_lib::cryptanalysis::ic_boundary::{koblitz_instance, Calibration};
use crypto_lib::cryptanalysis::ic_framework::shared_rank::{run_shared_rank, SharedRankSpec};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args()
        .nth(1)
        .ok_or("usage: n37_shared_rank_folded OUTPUT.json")?;
    let inst = koblitz_instance(0, 37).ok_or("registered n37 curve unavailable")?;
    let spec = SharedRankSpec {
        columns: 42,
        raw_x_cap: 1_000_000,
        rank_seed: 202610031137,
        max_trials: 1_000_000,
    };
    let report = run_shared_rank(&inst, &spec, &Calibration::default())?;
    let file = std::fs::File::create(path)?;
    serde_json::to_writer_pretty(file, &report)?;
    if !report.verified {
        return Err(format!(
            "target-blind rank gate failed at rank {}/{} after {} trials",
            report.rank, report.columns, report.trials
        )
        .into());
    }
    Ok(())
}
