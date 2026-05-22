use crate::buchberger::BuchbergerState;
use rand::Rng;

// The hook a learned policy plugs into. `select` returns an index into state.pairs.
pub trait Strategy {
    fn select(&mut self, state: &BuchbergerState) -> usize;
    fn name(&self) -> &str;
}

pub struct NormalStrategy;
impl Strategy for NormalStrategy {
    fn select(&mut self, state: &BuchbergerState) -> usize {
        (0..state.pairs.len())
            .min_by_key(|&i| state.pairs[i].lcm.degree())
            .expect("no pairs to select from")
    }
    fn name(&self) -> &str { "normal" }
}

pub struct SugarStrategy;
impl Strategy for SugarStrategy {
    fn select(&mut self, state: &BuchbergerState) -> usize {
        (0..state.pairs.len())
            .min_by_key(|&i| (state.pairs[i].sugar, state.pairs[i].lcm.degree()))
            .expect("no pairs to select from")
    }
    fn name(&self) -> &str { "sugar" }
}

pub struct RandomStrategy<R: Rng>(pub R);
impl<R: Rng> Strategy for RandomStrategy<R> {
    fn select(&mut self, state: &BuchbergerState) -> usize {
        self.0.gen_range(0..state.pairs.len())
    }
    fn name(&self) -> &str { "random" }
}

// Stub for a learned policy. In practice the RL trainer sits on the Python side
// and chooses pairs via the JSON-RPC env; this struct exists so the trait surface
// is unified inside Rust (useful for offline rollouts of an exported policy).
pub struct LearnedPolicy {
    // TODO: load weights (ONNX / torch traced) and run inference here.
}
impl Strategy for LearnedPolicy {
    fn select(&mut self, _state: &BuchbergerState) -> usize {
        unimplemented!("wire up a Rust-side inference runtime (tract / candle / ort)")
    }
    fn name(&self) -> &str { "learned" }
}
