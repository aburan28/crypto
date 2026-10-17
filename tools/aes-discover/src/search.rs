//! Search is exploratory; validation consumes a frozen, hashed event family.
use crate::cipher::{xor, Aes, Block, Config};
use crate::{block, digest, hex, stats, stream};
use rand::RngCore;
use rand_chacha::ChaCha20Rng;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::{HashMap, HashSet};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, clap::ValueEnum)]
#[serde(rename_all = "snake_case")]
pub enum OracleKind {
    Aes,
    Random,
}

/// Lazy permutation, including consistent inverse queries. Outputs are sampled
/// without replacement. The public fixture PRNG is explicitly ChaCha20.
pub struct RandomPermutation {
    rng: ChaCha20Rng,
    forward: HashMap<Block, Block>,
    backward: HashMap<Block, Block>,
}
impl RandomPermutation {
    pub fn new(seed: u64, domain: &str, index: u64) -> Self {
        Self {
            rng: stream(seed, domain, index),
            forward: HashMap::new(),
            backward: HashMap::new(),
        }
    }
    pub fn encrypt(&mut self, x: Block) -> Block {
        if let Some(y) = self.forward.get(&x) {
            return *y;
        }
        let y = loop {
            let mut y = [0; 16];
            self.rng.fill_bytes(&mut y);
            if !self.backward.contains_key(&y) {
                break y;
            }
        };
        self.forward.insert(x, y);
        self.backward.insert(y, x);
        y
    }
    pub fn decrypt(&mut self, y: Block) -> Block {
        if let Some(x) = self.backward.get(&y) {
            return *x;
        }
        let x = loop {
            let mut x = [0; 16];
            self.rng.fill_bytes(&mut x);
            if !self.forward.contains_key(&x) {
                break x;
            }
        };
        self.forward.insert(x, y);
        self.backward.insert(y, x);
        x
    }
}

pub enum Oracle {
    Aes(Aes),
    Random(Box<RandomPermutation>),
}
impl Oracle {
    pub fn encrypt(&mut self, x: Block) -> Block {
        match self {
            Self::Aes(a) => a.encrypt(x),
            Self::Random(p) => p.encrypt(x),
        }
    }
}

pub fn public_key(seed: u64, domain: &str, index: u64, bits: usize) -> Vec<u8> {
    let mut key = vec![0; bits / 8];
    stream(seed, domain, index).fill_bytes(&mut key);
    key
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Projection {
    /// XOR these distinct output-difference bytes, then apply mask.
    pub bytes: Vec<usize>,
    pub mask: u8,
    pub value: u8,
}
impl Projection {
    pub fn check(&self) -> Result<(), String> {
        if self.bytes.is_empty()
            || self.bytes.len() > 16
            || self.mask == 0
            || self.value & !self.mask != 0
        {
            return Err("invalid projection".into());
        }
        let mut seen = HashSet::new();
        for &b in &self.bytes {
            if b >= 16 || !seen.insert(b) {
                return Err("projection bytes must be distinct and in 0..16".into());
            }
        }
        Ok(())
    }
    pub fn evaluate(&self, d: &Block) -> bool {
        self.bytes.iter().fold(0, |acc, &i| acc ^ d[i]) & self.mask == self.value
    }
    pub fn null_probability(&self) -> f64 {
        stats::null_probability(self.mask.count_ones(), self.value == 0)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Candidate {
    pub id: String,
    pub input_difference: String,
    pub event: Projection,
    pub direction: String,
    pub discovery_hits: u64,
    pub discovery_samples: u64,
    pub discovery_kl_score: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct FrozenSpec {
    pub schema: String,
    pub cipher: Config,
    pub discovery_oracle: OracleKind,
    pub discovery_seed: u64,
    pub discovery_key_sha256: Vec<String>,
    pub searched_input_differences: usize,
    pub searched_events_per_difference: usize,
    pub discovery_transcript_sha256: String,
    pub candidates: Vec<Candidate>,
    pub validation: ValidationPlan,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct ValidationPlan {
    pub seed: u64,
    pub samples_per_key_per_event: u64,
    pub keys: usize,
    pub alpha: f64,
    pub primary_oracle: OracleKind,
    pub stopping_rule: String,
}

impl ValidationPlan {
    fn check(&self, discovery_seed: u64, candidates: usize, max_pairs: u64) -> Result<(), String> {
        if self.seed == discovery_seed {
            return Err("validation seed must differ from discovery".into());
        }
        if !self.alpha.is_finite() || self.alpha <= 0.0 || self.alpha >= 1.0 {
            return Err("alpha must be between zero and one".into());
        }
        if self.stopping_rule != "fixed_samples_no_optional_stopping" {
            return Err("unsupported stopping rule".into());
        }
        budget(
            self.samples_per_key_per_event,
            self.keys,
            candidates,
            max_pairs,
        )
    }
}

impl FrozenSpec {
    pub fn check(&self, max_pairs: u64) -> Result<(), String> {
        self.cipher.check()?;
        if self.schema != "aes-discover/spec/v1" {
            return Err("unknown spec schema".into());
        }
        if self.discovery_key_sha256.is_empty()
            || self.candidates.is_empty()
            || self.candidates.len() > 4096
        {
            return Err("spec must have discovery keys and 1..4096 candidates".into());
        }
        self.validation
            .check(self.discovery_seed, self.candidates.len(), max_pairs)?;
        let mut ids = HashSet::new();
        for c in &self.candidates {
            c.event.check()?;
            if block(&c.input_difference)? == [0; 16] {
                return Err("zero input difference is not admissible".into());
            }
            if c.direction != "increased" && c.direction != "decreased" {
                return Err("invalid direction".into());
            }
            if c.id != candidate_id(self.cipher, &c.input_difference, &c.event)?
                || !ids.insert(&c.id)
            {
                return Err("candidate ID mismatch or duplicate".into());
            }
            if !c.discovery_kl_score.is_finite()
                || c.discovery_samples == 0
                || c.discovery_hits > c.discovery_samples
            {
                return Err("invalid discovery diagnostics".into());
            }
        }
        Ok(())
    }
}

pub fn budget(samples: u64, keys: usize, candidates: usize, max_pairs: u64) -> Result<(), String> {
    if samples == 0 || keys == 0 || keys > 1024 || candidates == 0 || candidates > 4096 {
        return Err("samples must be positive; keys in 1..1024; candidates in 1..4096".into());
    }
    let total = samples
        .checked_mul(keys as u64)
        .and_then(|x| x.checked_mul(candidates as u64))
        .ok_or("pair count overflow")?;
    total.checked_mul(2).ok_or("oracle query count overflow")?;
    if total > max_pairs {
        return Err(format!("{total} pairs exceed --max-pairs={max_pairs}; raise it explicitly for a larger process"));
    }
    Ok(())
}

fn candidate_id(
    config: Config,
    difference: &str,
    projection: &Projection,
) -> Result<String, String> {
    let bytes = serde_json::to_vec(&(config, difference, projection)).map_err(|e| e.to_string())?;
    Ok(format!("AESD1-{}", &digest(&bytes)[..24]))
}

fn groups() -> Vec<Vec<usize>> {
    let mut g: Vec<Vec<usize>> = (0..16).map(|i| vec![i]).collect();
    g.extend((0..4).map(|c| (0..4).map(|r| 4 * c + r).collect()));
    g.extend((0..4).map(|r| (0..4).map(|c| 4 * c + r).collect()));
    g.extend((0..4).map(|c| (0..4).map(|r| 4 * ((c + 4 - r) % 4) + r).collect()));
    g.push((0..16).collect());
    g
}

fn difference(index: usize, rng: &mut ChaCha20Rng) -> Block {
    let mut d = [0; 16];
    if index < 48 {
        d[(index / 3) % 16] = [1, 0x80, 0x1f][index % 3];
    } else {
        let n = 2 + ((index - 48) / 32).min(14);
        while d.iter().filter(|&&b| b != 0).count() < n {
            let p = rng.next_u32() as usize % 16;
            d[p] = (rng.next_u32() % 255 + 1) as u8;
        }
    }
    d
}

fn next_pair(rng: &mut ChaCha20Rng, used: &mut HashSet<Block>, d: Block) -> (Block, Block) {
    loop {
        let mut p = [0; 16];
        rng.fill_bytes(&mut p);
        let q = xor(p, d);
        if !used.contains(&p) && !used.contains(&q) {
            used.insert(p);
            used.insert(q);
            return (p, q);
        }
    }
}

pub struct DiscoveryOptions {
    pub config: Config,
    pub oracle: OracleKind,
    pub seed: u64,
    pub samples: u64,
    pub keys: usize,
    pub candidates: usize,
    pub top: usize,
    pub max_pairs: u64,
    pub validation: ValidationPlan,
}

pub fn discover(o: DiscoveryOptions) -> Result<FrozenSpec, String> {
    o.config.check()?;
    budget(o.samples, o.keys, o.candidates, o.max_pairs)?;
    if o.top == 0 || o.top > o.candidates {
        return Err("top must be in 1..=candidates".into());
    }
    o.validation.check(o.seed, o.top, o.max_pairs)?;
    let gs = groups();
    let mut diff_rng = stream(o.seed, "discovery/differences", 0);
    let keys: Vec<_> = (0..o.keys)
        .map(|k| public_key(o.seed, "discovery/keys", k as u64, o.config.key_bits))
        .collect();
    let mut oracles: Vec<_> = keys
        .iter()
        .enumerate()
        .map(|(k, key)| {
            Ok(match o.oracle {
                OracleKind::Aes => Oracle::Aes(Aes::new(key, o.config)?),
                OracleKind::Random => Oracle::Random(Box::new(RandomPermutation::new(
                    o.seed,
                    "discovery/primary-prp",
                    k as u64,
                ))),
            })
        })
        .collect::<Result<_, String>>()?;
    let mut used: Vec<HashSet<Block>> = (0..o.keys).map(|_| HashSet::new()).collect();
    let mut transcript = Sha256::new();
    let mut candidates = Vec::new();
    for i in 0..o.candidates {
        let d = difference(i, &mut diff_rng);
        let mut hist = vec![[0u64; 256]; gs.len()];
        for k in 0..o.keys {
            let mut rng = stream(o.seed, "discovery/plaintexts", (i * o.keys + k) as u64);
            for _ in 0..o.samples {
                let (p, q) = next_pair(&mut rng, &mut used[k], d);
                let a = oracles[k].encrypt(p);
                let b = oracles[k].encrypt(q);
                transcript.update(p);
                transcript.update(q);
                transcript.update(a);
                transcript.update(b);
                let out = xor(a, b);
                for (g, counts) in gs.iter().zip(&mut hist) {
                    counts[g.iter().fold(0u8, |acc, &j| acc ^ out[j]) as usize] += 1;
                }
            }
        }
        let n = o.samples * o.keys as u64;
        let mut best = (
            f64::NEG_INFINITY,
            Projection {
                bytes: vec![0],
                mask: 255,
                value: 0,
            },
            0u64,
        );
        for (g, counts) in gs.iter().zip(&hist) {
            for (value, &count) in counts.iter().enumerate() {
                let p = stats::null_probability(8, value == 0);
                let score = n as f64 * stats::kl(count as f64 / n as f64, p);
                if score > best.0 {
                    best = (
                        score,
                        Projection {
                            bytes: g.clone(),
                            mask: 255,
                            value: value as u8,
                        },
                        count,
                    );
                }
            }
            for bit in 0..8 {
                let mask = 1 << bit;
                let count: u64 = counts
                    .iter()
                    .enumerate()
                    .filter(|(v, _)| v & mask == 0)
                    .map(|(_, n)| *n)
                    .sum();
                let score = n as f64 * stats::kl(count as f64 / n as f64, 0.5);
                if score > best.0 {
                    best = (
                        score,
                        Projection {
                            bytes: g.clone(),
                            mask: mask as u8,
                            value: 0,
                        },
                        count,
                    );
                }
            }
        }
        let input_difference = hex(&d);
        let direction = if best.2 as f64 / n as f64 >= best.1.null_probability() {
            "increased"
        } else {
            "decreased"
        }
        .into();
        candidates.push(Candidate {
            id: candidate_id(o.config, &input_difference, &best.1)?,
            input_difference,
            event: best.1,
            direction,
            discovery_hits: best.2,
            discovery_samples: n,
            discovery_kl_score: best.0,
        });
    }
    candidates.sort_by(|a, b| {
        b.discovery_kl_score
            .total_cmp(&a.discovery_kl_score)
            .then(a.id.cmp(&b.id))
    });
    candidates.truncate(o.top);
    let spec = FrozenSpec {
        schema: "aes-discover/spec/v1".into(),
        cipher: o.config,
        discovery_oracle: o.oracle,
        discovery_seed: o.seed,
        discovery_key_sha256: keys.iter().map(|k| digest(k)).collect(),
        searched_input_differences: o.candidates,
        searched_events_per_difference: gs.len() * (256 + 8),
        discovery_transcript_sha256: hex(&transcript.finalize()),
        candidates,
        validation: o.validation,
    };
    spec.check(o.max_pairs)?;
    Ok(spec)
}

#[derive(Debug, Serialize)]
pub struct ValidationRow {
    pub candidate_id: String,
    pub key_index: usize,
    pub public_key_sha256: String,
    pub samples: u64,
    pub primary_hits: u64,
    pub control_hits: u64,
    pub primary_rate: f64,
    pub control_rate: f64,
    pub ideal_pair_probability: f64,
    pub effect_from_ideal: f64,
    pub hoeffding_family_deviation_threshold: f64,
    pub sufficient_effect_for_80pct_power_iid_model: f64,
    pub primary_family_bound: f64,
    pub control_family_bound: f64,
    pub direction_matches: bool,
    pub candidate_pass: bool,
    pub transcript_sha256: String,
}

#[derive(Debug, Serialize)]
pub struct ValidationReport {
    pub schema: String,
    pub spec_sha256: String,
    pub validation_plan: ValidationPlan,
    pub cipher: Config,
    pub family_tests: usize,
    pub oracle_queries_per_arm: u64,
    pub prp_coupling_slack_per_test: f64,
    pub statistical_method: String,
    pub random_control_passed: bool,
    pub validated_candidates: Vec<String>,
    pub rows: Vec<ValidationRow>,
    pub scope: String,
}

pub fn validate(
    spec: &FrozenSpec,
    exact_spec_bytes: &[u8],
    max_pairs: u64,
) -> Result<ValidationReport, String> {
    let bound_spec: FrozenSpec =
        serde_json::from_slice(exact_spec_bytes).map_err(|e| e.to_string())?;
    if &bound_spec != spec {
        return Err("specification object does not match frozen bytes".into());
    }
    spec.check(max_pairs)?;
    let plan = &spec.validation;
    let m = spec.candidates.len() * plan.keys * 2; // Both primary and null, all keys/events.
    let n = plan.samples_per_key_per_event;
    let q_per_key = 2 * n * spec.candidates.len() as u64;
    let mut rows = Vec::new();
    let mut key_ids = HashSet::new();
    for k in 0..plan.keys {
        let key = public_key(plan.seed, "validation/keys", k as u64, spec.cipher.key_bits);
        let key_id = digest(&key);
        if spec.discovery_key_sha256.contains(&key_id) || !key_ids.insert(key_id.clone()) {
            return Err("validation key overlap detected".into());
        }
        let mut primary = match plan.primary_oracle {
            OracleKind::Aes => Oracle::Aes(Aes::new(&key, spec.cipher)?),
            OracleKind::Random => Oracle::Random(Box::new(RandomPermutation::new(
                plan.seed,
                "validation/primary-prp",
                k as u64,
            ))),
        };
        let mut control = RandomPermutation::new(plan.seed, "validation/control-prp", k as u64);
        let mut used = HashSet::new();
        for (i, c) in spec.candidates.iter().enumerate() {
            let d = block(&c.input_difference)?;
            let mut rng = stream(
                plan.seed,
                "validation/plaintexts",
                (k * spec.candidates.len() + i) as u64,
            );
            let mut hits = 0;
            let mut null_hits = 0;
            let mut transcript = Sha256::new();
            for _ in 0..n {
                let (p, q) = next_pair(&mut rng, &mut used, d);
                let a = primary.encrypt(p);
                let b = primary.encrypt(q);
                let na = control.encrypt(p);
                let nb = control.encrypt(q);
                hits += u64::from(c.event.evaluate(&xor(a, b)));
                null_hits += u64::from(c.event.evaluate(&xor(na, nb)));
                for v in [p, q, a, b, na, nb] {
                    transcript.update(v);
                }
            }
            let p0 = c.event.null_probability();
            let primary_family_bound = stats::corrected_bound(hits, n, p0, m, q_per_key);
            let control_family_bound = stats::corrected_bound(null_hits, n, p0, m, q_per_key);
            let rate = hits as f64 / n as f64;
            let direction_matches = if c.direction == "increased" {
                rate > p0
            } else {
                rate < p0
            };
            rows.push(ValidationRow {
                candidate_id: c.id.clone(),
                key_index: k,
                public_key_sha256: key_id.clone(),
                samples: n,
                primary_hits: hits,
                control_hits: null_hits,
                primary_rate: rate,
                control_rate: null_hits as f64 / n as f64,
                ideal_pair_probability: p0,
                effect_from_ideal: rate - p0,
                hoeffding_family_deviation_threshold: ((2.0 * m as f64 / plan.alpha).ln()
                    / (2.0 * n as f64))
                    .sqrt(),
                sufficient_effect_for_80pct_power_iid_model: ((2.0 * m as f64 / plan.alpha).ln()
                    / (2.0 * n as f64))
                    .sqrt()
                    + (5.0_f64.ln() / (2.0 * n as f64)).sqrt(),
                primary_family_bound,
                control_family_bound,
                direction_matches,
                candidate_pass: primary_family_bound <= plan.alpha && direction_matches,
                transcript_sha256: hex(&transcript.finalize()),
            });
        }
    }
    let random_control_passed = rows.iter().all(|r| r.control_family_bound > plan.alpha);
    let validated_candidates = if random_control_passed {
        spec.candidates
            .iter()
            .filter(|c| {
                rows.iter()
                    .filter(|r| r.candidate_id == c.id)
                    .all(|r| r.candidate_pass)
            })
            .map(|c| c.id.clone())
            .collect()
    } else {
        Vec::new()
    };
    Ok(ValidationReport { schema:"aes-discover/validation/v1".into(),spec_sha256:digest(exact_spec_bytes),
        validation_plan:plan.clone(),cipher:spec.cipher,family_tests:m,
        oracle_queries_per_arm:q_per_key * plan.keys as u64,
        prp_coupling_slack_per_test:stats::coupling_slack(q_per_key),
        statistical_method:"two-sided Bernoulli Chernoff upper bound + finite-query PRP coupling slack; Bonferroni over both arms, every frozen event and key; fixed sample counts".into(),
        random_control_passed,validated_candidates,rows,
        scope:"fresh public synthetic fixed keys tested individually; no universal-key, novelty, optimality, impossibility, key-recovery or full-AES break claim".into() })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn options(rounds: usize, oracle: OracleKind) -> DiscoveryOptions {
        DiscoveryOptions {
            config: Config {
                key_bits: 128,
                rounds,
                final_mix_columns: false,
            },
            oracle,
            seed: 7,
            samples: 64,
            keys: 1,
            candidates: 4,
            top: 2,
            max_pairs: 100_000,
            validation: ValidationPlan {
                seed: 9,
                samples_per_key_per_event: 1024,
                keys: 3,
                alpha: 0.01,
                primary_oracle: oracle,
                stopping_rule: "fixed_samples_no_optional_stopping".into(),
            },
        }
    }
    #[test]
    fn permutation_is_bijective_and_query_order_consistent() {
        let mut p = RandomPermutation::new(4, "test", 0);
        let mut outputs = HashSet::new();
        for i in 0..2048u64 {
            let mut x = [0; 16];
            x[..8].copy_from_slice(&i.to_le_bytes());
            let y = p.encrypt(x);
            assert_eq!(p.encrypt(x), y);
            assert_eq!(p.decrypt(y), x);
            assert!(outputs.insert(y));
        }
        let x = p.decrypt([0; 16]);
        assert_eq!(p.encrypt(x), [0; 16]);
    }
    #[test]
    fn known_signal_survives_fresh_keys_and_random_winner_does_not() {
        let spec = discover(options(1, OracleKind::Aes)).unwrap();
        let bytes = serde_json::to_vec(&spec).unwrap();
        let a = validate(&spec, &bytes, 100_000).unwrap();
        assert_eq!(a.validated_candidates.len(), 2);
        assert!(a.random_control_passed);
        let spec = discover(options(10, OracleKind::Random)).unwrap();
        let bytes = serde_json::to_vec(&spec).unwrap();
        let a = validate(&spec, &bytes, 100_000).unwrap();
        assert!(a.validated_candidates.is_empty());
        assert!(a.random_control_passed);
    }
    #[test]
    fn duplicate_bytes_zero_differences_and_reused_seed_rejected() {
        assert!(Projection {
            bytes: vec![0, 0],
            mask: 255,
            value: 0
        }
        .check()
        .is_err());
        let mut spec = discover(options(1, OracleKind::Aes)).unwrap();
        spec.validation.seed = spec.discovery_seed;
        assert!(spec.check(100_000).is_err());
        spec.validation.seed = 99;
        spec.candidates[0].input_difference = "00".repeat(16);
        assert!(spec.check(100_000).is_err());
        assert!(budget(u64::MAX, 1024, 4096, u64::MAX).is_err());
        assert!(budget(u64::MAX, 1, 1, u64::MAX).is_err());
    }
    #[test]
    fn library_rejects_a_forged_frozen_byte_binding() {
        let mut spec = discover(options(1, OracleKind::Aes)).unwrap();
        let bytes = serde_json::to_vec(&spec).unwrap();
        spec.validation.seed += 1;
        assert!(validate(&spec, &bytes, 100_000)
            .unwrap_err()
            .contains("frozen bytes"));
    }
}
