//! Summarise `research/pdp_shape_replay_20261015`: for each ladder (the
//! pilot's frozen one and the three registered ones), every cell's word
//! operations under the default and `KIC_F4_ROOT_REPLAY=1` arms, their
//! ratio, the splits, and whether the verdicts agree; ladder totals.
//! Writes `summary.json` in the study directory.
//!
//! ```sh
//! cargo run --release --example pdp_shape_replay_summary -- research/pdp_shape_replay_20261015
//! ```

use std::fs;
use std::path::Path;

use serde_json::{json, Value};

fn load(p: &Path) -> Value {
    serde_json::from_str(&fs::read_to_string(p).expect("stage.json")).expect("JSON")
}

fn main() {
    let dir = std::env::args().nth(1).expect("the study directory");
    let dir = Path::new(&dir);
    let arms = [
        (
            "frozen (pilot)",
            dir.join("pilot/frozen-default/stage.json"),
            dir.join("pilot/frozen-replay/stage.json"),
        ),
        (
            "chain",
            dir.join("raw/chain-default/stage.json"),
            dir.join("raw/chain-replay/stage.json"),
        ),
        (
            "chain-holdout",
            dir.join("raw/chain-holdout-default/stage.json"),
            dir.join("raw/chain-holdout-replay/stage.json"),
        ),
        (
            "holdout",
            dir.join("raw/holdout-default/stage.json"),
            dir.join("raw/holdout-replay/stage.json"),
        ),
    ];
    let mut ladders = Vec::new();
    for (name, a, b) in arms {
        let (a, b) = (load(&a), load(&b));
        let (mut ta, mut tb) = (0u64, 0u64);
        let mut cells = Vec::new();
        for (ra, rb) in a["rows"]
            .as_array()
            .unwrap()
            .iter()
            .zip(b["rows"].as_array().unwrap())
        {
            assert_eq!(ra["curve"], rb["curve"]);
            assert_eq!(ra["m"], rb["m"]);
            let (wa, wb) = (
                ra["word_ops"].as_u64().unwrap(),
                rb["word_ops"].as_u64().unwrap(),
            );
            ta += wa;
            tb += wb;
            let same = ra["verdict_digest"] == rb["verdict_digest"]
                && ra["decomposed"] == rb["decomposed"]
                && ra["exhausted"] == rb["exhausted"];
            cells.push(json!({
                "curve": ra["curve"], "m": ra["m"], "factor_index": ra["factor_index"],
                "first_target": ra["first_target"], "targets": ra["targets"],
                "word_ops_default": wa, "word_ops_replay": wb,
                "default_over_replay": wa as f64 / wb as f64,
                "splits_default": ra["splits"], "splits_replay": rb["splits"],
                "f4_calls_default": ra["f4_calls"], "f4_calls_replay": rb["f4_calls"],
                "verdicts_identical": same,
            }));
        }
        ladders.push(json!({
            "ladder": name, "cells": cells,
            "total_default": ta, "total_replay": tb,
            "default_over_replay": ta as f64 / tb as f64,
        }));
    }
    let out = json!({"schema": "pdp_shape_replay.summary/v1", "ladders": ladders});
    let path = dir.join("summary.json");
    fs::write(&path, serde_json::to_string_pretty(&out).unwrap() + "\n").unwrap();
    eprintln!("wrote {}", path.display());
}
