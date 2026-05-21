use crate::buchberger::BuchbergerState;
use crate::poly::Poly;
use serde::{Deserialize, Serialize};

// Per-polynomial features. Cheap scalars first; a real GNN feature pipeline would
// also expose support vectors and leading-monomial signatures.
#[derive(Serialize, Deserialize, Debug)]
pub struct PolyFeat {
    pub idx: usize,
    pub total_degree: u32,
    pub n_terms: usize,
    pub lm_degree: u32,
    pub sugar: u32,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct PairFeat {
    pub idx: usize,
    pub i: usize,
    pub j: usize,
    pub lcm_degree: u32,
    pub sugar: u32,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Observation {
    pub nvars: usize,
    pub basis_size: usize,
    pub polys: Vec<PolyFeat>,
    pub pairs: Vec<PairFeat>,
    pub step_count: usize,
    pub arith_ops: usize,
    pub done: bool,
}

pub struct Env {
    pub state: BuchbergerState,
    pub max_steps: usize,
}

impl Env {
    pub fn new(initial: Vec<Poly>, max_steps: usize) -> Self {
        Env { state: BuchbergerState::new(initial), max_steps }
    }

    pub fn done(&self) -> bool {
        self.state.is_done() || self.state.step_count >= self.max_steps
    }

    pub fn observe(&self) -> Observation {
        let polys = self.state.basis.iter().enumerate().map(|(idx, p)| PolyFeat {
            idx,
            total_degree: p.total_degree(),
            n_terms: p.nterms(),
            lm_degree: p.lm().map(|m| m.degree()).unwrap_or(0),
            sugar: self.state.sugars[idx],
        }).collect();

        let pairs = self.state.pairs.iter().enumerate().map(|(idx, p)| PairFeat {
            idx,
            i: p.i,
            j: p.j,
            lcm_degree: p.lcm.degree(),
            sugar: p.sugar,
        }).collect();

        Observation {
            nvars: self.state.nvars,
            basis_size: self.state.basis.len(),
            polys,
            pairs,
            step_count: self.state.step_count,
            arith_ops: self.state.arith_ops,
            done: self.done(),
        }
    }

    // Returns (reward, done). Reward shaping is intentionally minimal here —
    // -1 per step incentivizes shorter rollouts. Terminal bonuses / arith-op
    // penalties belong in the Python trainer where they're easy to ablate.
    pub fn step(&mut self, pair_idx: usize) -> (f64, bool) {
        if self.done() { return (0.0, true); }
        self.state.step(pair_idx);
        (-1.0, self.done())
    }
}
