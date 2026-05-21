use crate::monomial::Monomial;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CriticalPair {
    pub i: usize,
    pub j: usize,
    pub lcm: Monomial,
    pub sugar: u32,
}
