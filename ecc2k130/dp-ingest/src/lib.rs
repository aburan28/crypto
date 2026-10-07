//! Pure core of the ECC2K-130 distinguished-point ingest.
//!
//! This is the offline half of `ecc2k130/aws/dp_ingest.py`: record decoding,
//! key shapes, envelope checks, coverage counting, and the weight-32 cutoff
//! verdict. The live daemon is still Python; Postgres and S3 land in a later
//! phase under the same crate. See
//! `research/ecc2k130_dp_ingest_rust_20261007/PROTOCOL.md`.

pub mod coverage;
pub mod cutoff;
pub mod decode;
pub mod envelope;
pub mod keys;
pub mod status;

pub use coverage::{coverage_rows, covered_records, CoverageRow};
pub use cutoff::{
    dp_weight_verdict, estimate_dp_weight, theoretical_iter_per_dp_log2, CAMPAIGN_DP_WEIGHT,
    CAMPAIGN_ITER_PER_DP_LOG2, DP_RATIO_MIN_RECORDS, DP_RATIO_TOLERANCE_LOG2, FIELD_BITS,
};
pub use decode::{decode, seed_of, Decoded, COEFF_BYTES, RECORD_BYTES};
pub use envelope::{
    check_envelope, is_witness_v2, table_record_body, DP_MAGIC_TABLE3, DP_MAGIC_V2, ENVELOPE_SUFFIX,
};
pub use keys::{
    classify_key, found_at, slot_of_key, worker_id, Classified, KeyKind, LEGACY_KEY_RE,
    ORBIT_KEY_RE,
};
pub use status::{
    campaign_state, public_status, walk_rate_between, window_sum, CampaignState, WalkRate,
    MIN_WALK_RATE_SPAN_S,
};
