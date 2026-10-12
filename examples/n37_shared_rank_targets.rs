//! Point-only target descent over the verified, reusable n37 rank table.
//! This executable does not open or receive the scalar fixture.

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::ic_boundary::{koblitz_instance, Calibration};
use crypto_lib::cryptanalysis::ic_framework::shared_rank::{
    run_shared_rank_targets, SharedRankSpec, SharedTargetSpec,
};
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use std::fs;
use std::path::Path;

const POINTS_SHA: &str = "81e49c561bd208fa1f51f7bc38f62d430e1758ad20230be989f151c0168d3775";
const COUNT: usize = 16;

fn run(points_path: &Path, output: &Path) -> Result<(), String> {
    let points_bytes = fs::read(points_path).map_err(|error| error.to_string())?;
    let digest = sha256_hex(&points_bytes);
    if digest != POINTS_SHA {
        return Err("point-only target file differs from frozen input".into());
    }
    let mut points = Vec::new();
    for line in points_bytes
        .split(|&b| b == b'\n')
        .filter(|line| !line.is_empty())
    {
        let [x, y]: [u64; 2] =
            serde_json::from_slice(line).map_err(|error| format!("target coordinates: {error}"))?;
        points.push(FastPoint::affine(x, y));
    }
    if points.len() != COUNT {
        return Err(format!("expected {COUNT} targets, got {}", points.len()));
    }
    let inst = koblitz_instance(0, 37).ok_or("registered n37 curve unavailable")?;
    let rank_spec = SharedRankSpec {
        columns: 42,
        raw_x_cap: 1_000_000,
        rank_seed: 202610031137,
        max_trials: 1_000_000,
    };
    let target_spec = SharedTargetSpec {
        residual_seed: 202610031649,
        max_attempts: 64,
    };
    let report = run_shared_rank_targets(
        &inst,
        &rank_spec,
        &target_spec,
        &points,
        digest,
        &Calibration::default(),
    )?;
    serde_json::to_writer_pretty(
        fs::File::create(output).map_err(|error| error.to_string())?,
        &report,
    )
    .map_err(|error| error.to_string())?;
    if !report.verified {
        return Err("not all point-only target logs were recovered".into());
    }
    Ok(())
}

fn main() -> Result<(), String> {
    let mut args = std::env::args().skip(1);
    let points = args
        .next()
        .ok_or("usage: n37_shared_rank_targets POINTS RAW")?;
    let output = args.next().ok_or("missing RAW output")?;
    if args.next().is_some() {
        return Err("too many arguments".into());
    }
    run(Path::new(&points), Path::new(&output))
}
