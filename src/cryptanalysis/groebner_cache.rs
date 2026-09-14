//! Optional Redis/ElastiCache storage for reusable Boolean F4 results.
//!
//! The cache is deliberately fail-open: a missing endpoint, an unreachable
//! ElastiCache node, or a malformed cached value never changes the local
//! solver's answer.  Set `IC_GROEBNER_CACHE_URL` to enable it.  For an
//! ElastiCache replication group with in-transit encryption, use a `rediss://`
//! URL and provide the auth token in the URL through the deployment secret.

use crate::cryptanalysis::pq_groebner_f2::F2BoolPoly;
use redis::Commands;
use serde::{Deserialize, Serialize};
use std::sync::OnceLock;

const CACHE_VERSION: u8 = 1;
const DEFAULT_TTL_SECS: u64 = 7 * 24 * 60 * 60;
const DEFAULT_MAX_BYTES: usize = 16 * 1024 * 1024;

#[derive(Debug, Serialize, Deserialize)]
struct CachedRows {
    version: u8,
    n_vars: usize,
    rows: Vec<Vec<u64>>,
}

/// A lightweight Redis client.  Connections are short-lived so a dead node or
/// a failover does not poison a long-lived solver process.
pub(crate) struct GroebnerCache {
    client: redis::Client,
    namespace: String,
    ttl_secs: u64,
    max_bytes: usize,
}

static CACHE: OnceLock<Option<GroebnerCache>> = OnceLock::new();
static CONFIG_WARNING_EMITTED: OnceLock<()> = OnceLock::new();

impl GroebnerCache {
    fn from_env() -> Option<Self> {
        let url = std::env::var("IC_GROEBNER_CACHE_URL")
            .ok()
            .filter(|value| !value.trim().is_empty())?;
        let client = match redis::Client::open(url) {
            Ok(client) => client,
            Err(error) => {
                warn_once(format!("invalid IC_GROEBNER_CACHE_URL; cache disabled: {error}"));
                return None;
            }
        };
        let ttl_secs = parse_env_u64("IC_GROEBNER_CACHE_TTL_SECS", DEFAULT_TTL_SECS).max(1);
        let max_bytes =
            parse_env_usize("IC_GROEBNER_CACHE_MAX_BYTES", DEFAULT_MAX_BYTES).max(1);
        let namespace = std::env::var("IC_GROEBNER_CACHE_NAMESPACE")
            .ok()
            .filter(|value| !value.trim().is_empty())
            .unwrap_or_else(|| "crypto:ic:groebner:v1".to_owned());
        Some(Self {
            client,
            namespace,
            ttl_secs,
            max_bytes,
        })
    }

    fn key(&self, equations: &[F2BoolPoly], n_vars: usize, degree: u32) -> String {
        // Serialize the exact Boolean system, including equation and term
        // order.  The caller's deterministic construction order is part of
        // the solver input and therefore belongs in the cache key.
        let input = (
            CACHE_VERSION,
            n_vars,
            degree,
            equations
                .iter()
                .map(|poly| {
                    (
                        poly.n_vars,
                        poly.terms.iter().map(|term| term.mask).collect::<Vec<_>>(),
                    )
                })
                .collect::<Vec<_>>(),
        );
        let encoded = serde_json::to_vec(&input).expect("cache key input is serializable");
        let digest = blake3::hash(&encoded);
        format!("{}:f4:{}", self.namespace, hex::encode(digest.as_bytes()))
    }

    fn get(
        &self,
        equations: &[F2BoolPoly],
        n_vars: usize,
        degree: u32,
    ) -> Option<Vec<F2BoolPoly>> {
        let key = self.key(equations, n_vars, degree);
        let mut connection = match self.client.get_connection() {
            Ok(connection) => connection,
            Err(error) => {
                warn_once(format!("ElastiCache read connection failed; bypassing cache: {error}"));
                return None;
            }
        };
        let bytes: Option<Vec<u8>> = match connection.get(&key) {
            Ok(bytes) => bytes,
            Err(error) => {
                warn_once(format!("ElastiCache read failed; bypassing cache: {error}"));
                return None;
            }
        };
        let Some(bytes) = bytes else {
            return None;
        };
        match serde_json::from_slice::<CachedRows>(&bytes) {
            Ok(payload) if payload.version == CACHE_VERSION && payload.n_vars == n_vars => {
                let rows = payload
                    .rows
                    .into_iter()
                    .map(|terms| {
                        F2BoolPoly::from_monos(
                            terms
                                .into_iter()
                                .map(crate::cryptanalysis::pq_groebner_f2::F2BoolMono::from_mask)
                                .collect(),
                            n_vars,
                        )
                    })
                    .collect();
                Some(rows)
            }
            Ok(_) | Err(_) => {
                // A rolling deployment can encounter an old or truncated
                // value.  Treat it as a miss and let the local solver repair it.
                None
            }
        }
    }

    fn put(
        &self,
        equations: &[F2BoolPoly],
        n_vars: usize,
        degree: u32,
        rows: &[F2BoolPoly],
    ) {
        let payload = CachedRows {
            version: CACHE_VERSION,
            n_vars,
            rows: rows
                .iter()
                .map(|poly| poly.terms.iter().map(|term| term.mask).collect())
                .collect(),
        };
        let Ok(bytes) = serde_json::to_vec(&payload) else {
            return;
        };
        if bytes.len() > self.max_bytes {
            return;
        }
        let key = self.key(equations, n_vars, degree);
        let Ok(mut connection) = self.client.get_connection() else {
            return;
        };
        if let Err(error) = connection.set_ex::<_, _, ()>(&key, bytes, self.ttl_secs) {
            warn_once(format!("ElastiCache write failed; continuing without cache: {error}"));
        }
    }
}

/// Compute a matrix-F4 result, consulting ElastiCache when configured.
pub(crate) fn get_or_compute(
    equations: &[F2BoolPoly],
    n_vars: usize,
    degree: u32,
    compute: impl FnOnce() -> Option<Vec<F2BoolPoly>>,
) -> Option<Vec<F2BoolPoly>> {
    let cache = CACHE.get_or_init(GroebnerCache::from_env);
    let Some(cache) = cache.as_ref() else {
        return compute();
    };
    if let Some(rows) = cache.get(equations, n_vars, degree) {
        return Some(rows);
    }
    let rows = compute()?;
    cache.put(equations, n_vars, degree, &rows);
    Some(rows)
}

fn parse_env_u64(name: &str, default: u64) -> u64 {
    std::env::var(name)
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(default)
}

fn parse_env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(default)
}

fn warn_once(message: String) {
    if CONFIG_WARNING_EMITTED.set(()).is_ok() {
        eprintln!("warning: {message}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};

    #[test]
    fn cached_rows_round_trip_without_redis() {
        let rows = vec![F2BoolPoly::from_monos(
            vec![F2BoolMono::from_mask(0), F2BoolMono::from_mask(5)],
            3,
        )];
        let payload = CachedRows {
            version: CACHE_VERSION,
            n_vars: 3,
            rows: rows
                .iter()
                .map(|poly| poly.terms.iter().map(|term| term.mask).collect())
                .collect(),
        };
        let bytes = serde_json::to_vec(&payload).unwrap();
        let decoded: CachedRows = serde_json::from_slice(&bytes).unwrap();
        assert_eq!(decoded.version, CACHE_VERSION);
        assert_eq!(decoded.n_vars, 3);
        assert_eq!(decoded.rows, vec![vec![5, 0]]);
    }
}
