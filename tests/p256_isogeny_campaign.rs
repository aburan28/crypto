use crypto_lib::cryptanalysis::p256_isogeny_campaign::{
    CairnShardReceipt, ObjectRef, RunManifest, ShardComplete, COMPLETE_SCHEMA, P256_ICV1,
    RECEIPT_SCHEMA, RUN_SCHEMA,
};

const A: &str = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
const B: &str = "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb";
const C: &str = "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc";
const D: &str = "dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd";
const E: &str = "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee";

fn run() -> RunManifest {
    RunManifest {
        schema: RUN_SCHEMA.into(),
        run_id: "p256-2p40-v1".into(),
        curve: P256_ICV1.into(),
        shard_bits: 20,
        records_per_shard_bits: 20,
        config_sha256: A.into(),
        result_class: "stage-diagnostic".into(),
    }
}

fn marker(run: &RunManifest, shard_id: u64) -> ShardComplete {
    let attempt_id = "batch-42-attempt-1";
    let prefix = run.attempt_prefix(shard_id, attempt_id).unwrap();
    ShardComplete {
        schema: COMPLETE_SCHEMA.into(),
        run_id: run.run_id.clone(),
        config_sha256: run.config_sha256.clone(),
        shard_id,
        attempt_id: attempt_id.into(),
        records: 1 << 20,
        distinct_j_within_shard: 1 << 20,
        raw: ObjectRef {
            key: format!("{prefix}/raw.jsonl.zst"),
            sha256: B.into(),
            bytes: 4096,
        },
        stage1: ObjectRef {
            key: format!("{prefix}/stage1.jsonl.zst"),
            sha256: C.into(),
            bytes: 512,
        },
        merkle_root_sha256: D.into(),
    }
}

#[test]
fn production_layout_is_exactly_2p40_emitted_records() {
    let run = run();
    run.validate_2p40().unwrap();
    assert_eq!(run.shard_count().unwrap(), 1 << 20);
    assert_eq!(run.records_per_shard().unwrap(), 1 << 20);
    assert_eq!(run.total_records().unwrap(), 1u128 << 40);
    assert_eq!(
        run.complete_key(0xabc).unwrap(),
        "p256-isogeny-search/runs/p256-2p40-v1/shards/00000abc/complete.json"
    );
}

#[test]
fn strided_assignments_cover_each_shard_once() {
    let run = run();
    let task_count = 10_000;
    for shard in [0, 1, 9_999, 10_000, (1 << 20) - 1] {
        let owner = shard % task_count;
        let assignment = run.assignment(owner, task_count).unwrap();
        assert!(assignment.contains(shard));
        let position = shard / task_count;
        assert_eq!(assignment.shard_at(position), Some(shard));
    }
    assert_eq!(run.assignment(0, task_count).unwrap().len(), 105);
    assert_eq!(run.assignment(9_999, task_count).unwrap().len(), 104);
}

#[test]
fn a_commit_marker_can_only_name_its_immutable_attempt_objects() {
    let run = run();
    let mut changed_key = marker(&run, 7);
    changed_key.validate(&run).unwrap();

    changed_key.raw.key = changed_key.raw.key.replace("00000007", "00000008");
    assert!(changed_key.validate(&run).unwrap_err().contains("raw.key"));

    let mut marker = marker(&run, 7);
    marker.distinct_j_within_shard -= 1;
    assert!(marker
        .validate(&run)
        .unwrap_err()
        .contains("repeated j-invariant"));
}

#[test]
fn cairn_receipt_is_content_addressed_but_explicitly_not_a_work_proof() {
    let run = run();
    let receipt = marker(&run, 7).cairn_receipt(&run, E).unwrap();
    receipt.validate(&run).unwrap();
    assert_eq!(receipt.schema, RECEIPT_SCHEMA);
    assert_eq!(receipt.claim_class, "coordination-receipt");
    assert_eq!(receipt.verification, "s3-content-address-only");
    assert_eq!(receipt.claim_key(), "p256-2p40-v1:00000007");

    let artifact = receipt.artifact();
    assert_eq!(artifact["shard_id"], 7);
    assert_eq!(artifact["complete_sha256"], E);
    let hash = receipt
        .commitment_hash(
            "sha256:objective",
            "alice",
            "00112233445566778899aabbccddeeff",
        )
        .unwrap();
    assert!(hash.starts_with("sha256:"));
    assert_eq!(hash.len(), 71);
}

#[test]
fn receipt_cannot_be_relabelled_as_verified_piecework() {
    let run = run();
    let receipt = marker(&run, 7).cairn_receipt(&run, E).unwrap();
    let encoded = serde_json::to_string(&receipt).unwrap();
    let mut changed: CairnShardReceipt = serde_json::from_str(&encoded).unwrap();
    changed.verification = "verified-isogeny-work".into();
    assert!(changed.validate(&run).unwrap_err().contains("overstates"));
}
