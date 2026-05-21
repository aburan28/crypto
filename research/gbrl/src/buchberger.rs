use std::collections::HashSet;

use crate::monomial::Monomial;
use crate::pair::CriticalPair;
use crate::poly::Poly;
use crate::reduce::reduce;
use crate::spoly::spoly;

// Buchberger state. The RL action at each step is "which critical pair to process next" —
// the index into `pairs`.
pub struct BuchbergerState {
    pub basis: Vec<Poly>,
    pub sugars: Vec<u32>,
    pub pairs: Vec<CriticalPair>,
    pub nvars: usize,
    pub step_count: usize,
    pub arith_ops: usize,
}

pub struct StepResult {
    pub added: bool,
    pub new_poly_idx: Option<usize>,
    pub reduce_steps: usize,
}

impl BuchbergerState {
    pub fn new(initial: Vec<Poly>) -> Self {
        let nvars = initial.first().map(|p| p.nvars).unwrap_or(0);
        let sugars: Vec<u32> = initial.iter().map(|p| p.total_degree()).collect();
        let mut state = BuchbergerState {
            basis: initial,
            sugars,
            pairs: vec![],
            nvars,
            step_count: 0,
            arith_ops: 0,
        };
        // Install initial pairs in order, applying GM at each insertion. This
        // way the F0 input gets the same pruning the algorithm uses later.
        let n_initial = state.basis.len();
        if n_initial > 0 {
            // First pair: nothing to prune.
            for i in 0..n_initial {
                for j in (i + 1)..n_initial {
                    state.install_pair_raw(i, j);
                }
            }
            // Run GM "from the perspective of" each later basis element so the
            // initial queue is consistent with what step() would have produced
            // if we'd added them one-by-one.
            // (Simpler: just call install_pair_raw for all initial pairs; GM
            // will prune as new polys get added during step(). We accept a
            // slightly larger initial queue for cleaner semantics.)
        }
        state
    }

    // Install one critical pair. Only the Buchberger first criterion (coprime LM)
    // is applied here. Full GM pruning happens in `apply_gebauer_moller` after a
    // new poly is added — it needs the whole set of candidate pairs to reason
    // about M/F dominance, so per-pair pruning isn't enough.
    fn install_pair_raw(&mut self, i: usize, j: usize) {
        let fi = &self.basis[i];
        let fj = &self.basis[j];
        if fi.is_zero() || fj.is_zero() { return; }
        let lmi = fi.lm().unwrap();
        let lmj = fj.lm().unwrap();
        if lmi.gcd_is_one(lmj) { return; }
        let lcm = lmi.lcm(lmj);
        let dfi = lcm.degree() - lmi.degree();
        let dfj = lcm.degree() - lmj.degree();
        let sugar = (self.sugars[i] + dfi).max(self.sugars[j] + dfj);
        self.pairs.push(CriticalPair { i, j, lcm, sugar });
    }

    // Gebauer-Möller installation: given a freshly-added basis element at index
    // `idx`, generate new pairs (k, idx) and prune both new and old pairs using
    // the M, F, and B criteria.
    //
    // Reference: Gebauer & Möller (1988), "On an installation of Buchberger's
    // algorithm". The three criteria together cut S-polynomial work by ~order
    // of magnitude on cyclic-n and similar growth-bound families.
    fn apply_gebauer_moller(&mut self, idx: usize) {
        let lm_h = self.basis[idx].lm().expect("nonzero new poly").clone();
        let sugar_h = self.sugars[idx];

        // 1. Enumerate candidate new pairs (k, idx) with their lcms and sugars.
        //    Coprime pairs are dropped immediately (first criterion).
        let mut cands: Vec<(usize, Monomial, u32)> = Vec::new();
        for k in 0..idx {
            let fk = &self.basis[k];
            if fk.is_zero() { continue; }
            let lm_k = fk.lm().unwrap();
            if lm_k.gcd_is_one(&lm_h) { continue; }
            let lcm = lm_k.lcm(&lm_h);
            let dfk = lcm.degree() - lm_k.degree();
            let dfh = lcm.degree() - lm_h.degree();
            let sugar = (self.sugars[k] + dfk).max(sugar_h + dfh);
            cands.push((k, lcm, sugar));
        }

        // 2. M-criterion: drop pair (k, idx) if some other candidate (k', idx)
        //    has lcm strictly dividing this one. "Strictly" means divides AND
        //    not equal — otherwise we'd cascade-drop pairs with identical lcms.
        let n = cands.len();
        let mut keep = vec![true; n];
        for i in 0..n {
            if !keep[i] { continue; }
            for j in 0..n {
                if i == j || !keep[j] { continue; }
                if cands[j].1.divides(&cands[i].1) && cands[j].1 != cands[i].1 {
                    keep[i] = false;
                    break;
                }
            }
        }
        let cands: Vec<_> = cands.into_iter().zip(keep).filter_map(|(c, k)| if k { Some(c) } else { None }).collect();

        // 3. F-criterion: among surviving candidates sharing the same lcm, keep
        //    only one (we keep whichever comes first in basis-index order).
        let mut seen: HashSet<Monomial> = HashSet::new();
        let cands: Vec<_> = cands.into_iter().filter(|(_, lcm, _)| seen.insert(lcm.clone())).collect();

        // 4. B-criterion: prune existing queue. Discard old pair (i, j) if
        //    LM(h) divides lcm(LM(i), LM(j)) AND the "split" lcms with h are
        //    both strictly different from the original lcm. The lcm-inequality
        //    guards ensure we don't drop pairs where h gives no new info.
        let basis = &self.basis;
        let lm_h_ref = &lm_h;
        self.pairs.retain(|p| {
            if !lm_h_ref.divides(&p.lcm) { return true; }
            let lm_i = basis[p.i].lm().unwrap();
            let lm_j = basis[p.j].lm().unwrap();
            let lcm_ih = lm_i.lcm(lm_h_ref);
            let lcm_jh = lm_j.lcm(lm_h_ref);
            if lcm_ih == p.lcm || lcm_jh == p.lcm { return true; }
            false
        });

        // 5. Install the surviving new pairs.
        for (k, lcm, sugar) in cands {
            self.pairs.push(CriticalPair { i: k, j: idx, lcm, sugar });
        }
    }

    pub fn is_done(&self) -> bool { self.pairs.is_empty() }

    // Process pair at `pair_idx`. Returns whether the reduction added a new basis element.
    pub fn step(&mut self, pair_idx: usize) -> StepResult {
        assert!(pair_idx < self.pairs.len(), "pair_idx out of range");
        let pair = self.pairs.swap_remove(pair_idx);
        self.step_count += 1;

        let s = spoly(&self.basis[pair.i], &self.basis[pair.j]);
        let (r, ops) = reduce(&s, &self.basis);
        self.arith_ops += ops;

        if r.is_zero() {
            return StepResult { added: false, new_poly_idx: None, reduce_steps: ops };
        }

        let idx = self.basis.len();
        self.sugars.push(pair.sugar);
        self.basis.push(r);
        self.apply_gebauer_moller(idx);
        StepResult { added: true, new_poly_idx: Some(idx), reduce_steps: ops }
    }
}
