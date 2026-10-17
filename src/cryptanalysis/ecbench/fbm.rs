//! FBM1: factor-base *measurement* identity.
//!
//! FB1 names a point set on a curve. FBM1 names the measurement
//! configuration that produced yield / DoR / cost numbers against that
//! base: curve × FB1 × oracle × solver × linalg × instrument version.
//! Results (S, yield, solving degree, …) attach to an FBM1 id; they are
//! never part of the hash (same rule as ICV1: the identity is the object,
//! not a run's outcome).
//!
//! Spec: `docs/ecbench/FACTOR_BASE_PROTOCOL.md` §FBM1.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};

use crate::cryptanalysis::ecbench::canonical::short_id;
use crate::cryptanalysis::ecbench::methods::{FactorBaseFacts, ResolvedMethod, SolverStats};

/// Schema string frozen into every FBM1 preimage.
pub const FBM_SCHEMA: &str = "ecbench.factor_base_measurement/v1";

/// Instrument tag for the native `ic.pipeline` path.
pub const INSTRUMENT_IC_PIPELINE: &str = "ecbench.ic_pipeline/v1";

/// Instrument tag for Sage DoR cells joined by FB1.
pub const INSTRUMENT_SAGE_DOR: &str = "sage.factor_base_dor/v1";

/// The configuration half of an FBM1 identity (no measured numbers).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct FactorBaseMeasurementConfig {
    pub schema: String,
    pub curve: String,
    pub fb_id: String,
    pub fb_sha256: String,
    pub family: String,
    pub fb_params: BTreeMap<String, String>,
    pub oracle: String,
    pub solver: String,
    pub linalg: String,
    pub instrument: String,
}

/// A sealed FBM1 id plus the config that produced it.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct FactorBaseMeasurementId {
    /// `FBM1h` + 12 hex.
    pub fbm_id: String,
    pub fbm_sha256: String,
    pub config: FactorBaseMeasurementConfig,
}

impl FactorBaseMeasurementConfig {
    /// Preimage object (sorted keys via [`short_id`]).
    pub fn preimage(&self) -> Value {
        json!({
            "schema": self.schema,
            "curve": self.curve,
            "fb_id": self.fb_id,
            "fb_sha256": self.fb_sha256,
            "family": self.family,
            "fb_params": self.fb_params,
            "oracle": self.oracle,
            "solver": self.solver,
            "linalg": self.linalg,
            "instrument": self.instrument,
        })
    }

    pub fn seal(&self) -> Result<FactorBaseMeasurementId, String> {
        let (fbm_id, fbm_sha256) = short_id("FBM1", &self.preimage())?;
        Ok(FactorBaseMeasurementId {
            fbm_id,
            fbm_sha256,
            config: self.clone(),
        })
    }
}

/// Build an FBM1 config from an IC method and the FB1 facts it used.
pub fn config_from_ic_method(
    curve_slug: &str,
    fb: &FactorBaseFacts,
    method: &ResolvedMethod,
    instrument: &str,
) -> FactorBaseMeasurementConfig {
    let p = &method.params;
    FactorBaseMeasurementConfig {
        schema: FBM_SCHEMA.into(),
        curve: curve_slug.into(),
        fb_id: fb.fb_id.clone(),
        fb_sha256: fb.fb_sha256.clone(),
        family: fb.family.clone(),
        fb_params: fb.params.clone(),
        oracle: p.get("oracle").cloned().unwrap_or_default(),
        solver: p.get("solver").cloned().unwrap_or_default(),
        linalg: p
            .get("linalg")
            .cloned()
            .unwrap_or_else(|| "incremental-gauss".into()),
        instrument: instrument.into(),
    }
}

/// One measured cell under an FBM1 id (results, not identity).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FactorBaseMeasurementCell {
    pub schema: String,
    pub fbm_id: String,
    pub fbm_sha256: String,
    pub config: FactorBaseMeasurementConfig,
    /// Source record id when the cell came from an ecbench session.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub record_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub session_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub s: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub total_gae: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub trials: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub relations: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub yield_per_trial: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub solving_degree_mean: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub solving_degree_max: Option<u32>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub semi_regular_degree: Option<u32>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub macaulay_rows: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub macaulay_columns: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub status: Option<String>,
    /// Free-form extras (Sage D* fields, etc.).
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub extra: BTreeMap<String, Value>,
}

impl FactorBaseMeasurementCell {
    pub fn from_sealed_config(
        id: FactorBaseMeasurementId,
        record_id: Option<String>,
        session_id: Option<String>,
        s: Option<f64>,
        total_gae: Option<f64>,
        trials: Option<u64>,
        relations: Option<u64>,
        solver: Option<&SolverStats>,
        status: Option<String>,
    ) -> Self {
        let yield_per_trial = match (relations, trials) {
            (Some(r), Some(t)) if t > 0 => Some(r as f64 / t as f64),
            _ => None,
        };
        Self {
            schema: "ecbench.factor_base_measurement_cell/v1".into(),
            fbm_id: id.fbm_id,
            fbm_sha256: id.fbm_sha256,
            config: id.config,
            record_id,
            session_id,
            s,
            total_gae,
            trials,
            relations,
            yield_per_trial,
            solving_degree_mean: solver.and_then(|s| s.solving_degree_mean),
            solving_degree_max: solver.and_then(|s| s.solving_degree_max),
            semi_regular_degree: solver.and_then(|s| s.semi_regular_degree),
            macaulay_rows: solver.and_then(|s| s.macaulay_rows),
            macaulay_columns: solver.and_then(|s| s.macaulay_columns),
            status,
            extra: BTreeMap::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fbm1_is_stable_under_param_reorder() {
        let mut a = BTreeMap::new();
        a.insert("w".into(), "1".into());
        a.insert("columns".into(), "4".into());
        let mut b = BTreeMap::new();
        b.insert("columns".into(), "4".into());
        b.insert("w".into(), "1".into());
        let ca = FactorBaseMeasurementConfig {
            schema: FBM_SCHEMA.into(),
            curve: "icv1-f2m23-tm5197-1f85e9e1".into(),
            fb_id: "FB1h089756e11c3b".into(),
            fb_sha256: "aa".repeat(32),
            family: "cyclotomic".into(),
            fb_params: a,
            oracle: "mitm-frobenius:m=2".into(),
            solver: String::new(),
            linalg: "incremental-gauss".into(),
            instrument: INSTRUMENT_IC_PIPELINE.into(),
        };
        let mut cb = ca.clone();
        cb.fb_params = b;
        assert_eq!(ca.seal().unwrap().fbm_id, cb.seal().unwrap().fbm_id);
    }

    #[test]
    fn different_oracle_changes_fbm1() {
        let base = FactorBaseMeasurementConfig {
            schema: FBM_SCHEMA.into(),
            curve: "icv1-f2m17-tm101-00378d4e".into(),
            fb_id: "FB1hdd4118e0e739".into(),
            fb_sha256: "bb".repeat(32),
            family: "koblitz-orbit".into(),
            fb_params: BTreeMap::from([("divisor".into(), "0;1".into())]),
            oracle: "mitm-frobenius:m=2".into(),
            solver: String::new(),
            linalg: "incremental-gauss".into(),
            instrument: INSTRUMENT_IC_PIPELINE.into(),
        };
        let mut other = base.clone();
        other.oracle = "mitm:m=2".into();
        assert_ne!(base.seal().unwrap().fbm_id, other.seal().unwrap().fbm_id);
    }
}
