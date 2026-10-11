//! Status-document helpers that need no database or S3.

use std::collections::BTreeMap;

pub const MIN_WALK_RATE_SPAN_S: f64 = 600.0;

/// The page's `state` field, from collisions, recent points and ingest lag.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CampaignState {
    CollisionRecorded,
    Collecting,
    IngestBehind,
    IdleOrStale,
    Empty,
}

impl CampaignState {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::CollisionRecorded => "COLLISION_RECORDED",
            Self::Collecting => "COLLECTING",
            Self::IngestBehind => "INGEST_BEHIND",
            Self::IdleOrStale => "IDLE_OR_STALE",
            Self::Empty => "EMPTY",
        }
    }
}

/// Priority: collision → collecting → ingest behind → idle → empty.
pub fn campaign_state(
    collisions: u64,
    dps_last_hour: u64,
    outstanding: u64,
    unrecognised: u64,
    dps: u64,
) -> CampaignState {
    let behind = outstanding > 0 || unrecognised > 0;
    if collisions > 0 {
        CampaignState::CollisionRecorded
    } else if dps_last_hour > 0 {
        CampaignState::Collecting
    } else if behind {
        CampaignState::IngestBehind
    } else if dps > 0 {
        CampaignState::IdleOrStale
    } else {
        CampaignState::Empty
    }
}

/// Points in the last `hours`, from whole-hour rollup buckets.
///
/// `hourly` is `(unix_seconds_of_hour_start, count)`. `now` is injected so
/// tests do not depend on the wall clock.
pub fn window_sum(hourly: &[(i64, u64)], hours: i64, now: i64) -> u64 {
    let floor = now - hours * 3600;
    let floor = floor - floor.rem_euclid(3600);
    hourly
        .iter()
        .filter(|(hour, _)| *hour >= floor)
        .map(|(_, n)| *n)
        .sum()
}

/// Drop `work.per_slot` from a status payload (the public document).
pub fn public_status(
    mut work: BTreeMap<String, String>,
    has_per_slot: bool,
) -> BTreeMap<String, String> {
    if has_per_slot {
        work.remove("per_slot");
    }
    work
}

#[derive(Clone, Debug, PartialEq)]
pub struct WalkRate {
    pub iterations_per_second: f64,
    pub window_seconds: i64,
    pub measured_from: String,
    pub measured_to: String,
    pub iterations_from: u128,
    pub iterations_to: u128,
}

/// Iterations per second between two ingest publishes, or `None`.
pub fn walk_rate_between(
    cur_iter: u128,
    cur_at_unix: i64,
    cur_generated_at: &str,
    prev_iter: u128,
    prev_at_unix: i64,
    prev_generated_at: &str,
) -> Option<WalkRate> {
    if cur_iter == 0 || prev_iter == 0 || cur_iter < prev_iter {
        return None;
    }
    let span = (cur_at_unix - prev_at_unix) as f64;
    if span < MIN_WALK_RATE_SPAN_S {
        return None;
    }
    Some(WalkRate {
        iterations_per_second: (cur_iter - prev_iter) as f64 / span,
        window_seconds: span.round() as i64,
        measured_from: prev_generated_at.into(),
        measured_to: cur_generated_at.into(),
        iterations_from: prev_iter,
        iterations_to: cur_iter,
    })
}

/// Refuse a public blob that would leak private fields.
pub fn refuse_banned_fields(blob: &str) -> Result<(), String> {
    for banned in ["point_key", "walk_seed", "password", "secret"] {
        if blob.contains(banned) {
            return Err(format!("refusing to publish field {banned}"));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn state_priority() {
        assert_eq!(campaign_state(1, 0, 0, 0, 0).as_str(), "COLLISION_RECORDED");
        assert_eq!(campaign_state(0, 5, 9, 0, 0).as_str(), "COLLECTING");
        assert_eq!(campaign_state(0, 0, 3, 0, 0).as_str(), "INGEST_BEHIND");
        assert_eq!(campaign_state(0, 0, 0, 1, 10).as_str(), "INGEST_BEHIND");
        assert_eq!(campaign_state(0, 0, 0, 0, 10).as_str(), "IDLE_OR_STALE");
        assert_eq!(campaign_state(0, 0, 0, 0, 0).as_str(), "EMPTY");
    }

    #[test]
    fn window_sum_aligns_to_the_hour() {
        // now = 2026-09-18T17:30:00Z → floor for 1h is 16:00.
        let now = 1_726_678_200i64; // approx; use explicit buckets relative to now
        let hour = |h: i64| now - now.rem_euclid(3600) - h * 3600;
        let hourly = [(hour(0), 10), (hour(1), 20), (hour(2), 40)];
        assert_eq!(window_sum(&hourly, 1, now), 10 + 20);
        assert_eq!(window_sum(&hourly, 2, now), 10 + 20 + 40);
    }

    #[test]
    fn walk_rate_needs_the_minimum_span() {
        assert!(walk_rate_between(200, 1_000, "b", 100, 500, "a").is_none());
        let rate = walk_rate_between(200, 1_700, "b", 100, 1_000, "a").unwrap();
        assert!((rate.iterations_per_second - 100.0 / 700.0).abs() < 1e-12);
        assert_eq!(rate.window_seconds, 700);
    }

    #[test]
    fn banned_fields_are_refused() {
        assert!(refuse_banned_fields(r#"{"dps":1}"#).is_ok());
        assert!(refuse_banned_fields(r#"{"point_key":"x"}"#).is_err());
    }
}
