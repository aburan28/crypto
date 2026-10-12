//! Deterministic controls for the independent Python query-law checker.
//! Uses the real pinned rand crates; no timed candidate or DLP is run here.
use crypto_lib::cryptanalysis::koblitz_index_calculus::probe_scalar;
use rand::{rngs::StdRng, Rng, RngCore, SeedableRng};
use serde_json::json;

fn main() {
    let seeds = [0, 1, 2026092556, u64::MAX];
    let orders = [2, 31, 127, 65587, 1439393, u64::MAX];
    let trials = [0, 1, 7, 63, 64, 65, 127, 128, 65535, u64::MAX];
    let mut streams = Vec::new();
    let mut probes = Vec::new();
    for seed in seeds {
        let mut rng = StdRng::seed_from_u64(seed);
        streams.push(json!({"seed": seed, "order": null,
            "values": (0..257).map(|_| rng.next_u64()).collect::<Vec<_>>()}));
        for order in orders {
            let mut rng = StdRng::seed_from_u64(seed);
            streams.push(json!({"seed": seed, "order": order,
                "values": (0..257).map(|_| rng.gen_range(1..order)).collect::<Vec<_>>()}));
            probes.push(json!({"seed": seed, "order": order, "trials": trials,
                "values": trials.map(|t| probe_scalar(seed, t, order))}));
        }
    }
    println!(
        "{}",
        json!({"schema_version": 1, "streams": streams, "probes": probes})
    );
}
