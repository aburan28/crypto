//! # The factor bases, oracles and matrices that ship with the framework.
//!
//! Each of these wraps a component the repository already measures, so
//! a framework run reproduces a ledger row rather than approximating
//! one.  They are also the worked examples: a new plug-in looks like
//! these, and `docs/ic/FRAMEWORK.md` walks through writing one.
//!
//! ## Factor bases
//!
//! | name | regime | what it is |
//! |:--|:--|:--|
//! | `prime-abscissa` | prime | the smallest abscissae carrying a point, one column each |
//! | `binary-subspace` | binary, Koblitz | points whose abscissa lies in an `F_2`-subspace |
//! | `koblitz-orbit` | Koblitz | the same, each signed Frobenius orbit folded onto one column |
//! | `compact-orbit-scan` | Koblitz | deterministic raw-abscissa scan, subgroup projection and signed Frobenius closure |
//! | `cyclotomic` | Koblitz, `2n+1` prime | the γ-weight-`w` abscissae `Σ_{k∈S} γ^k` (lead A2), projected and folded like `compact-orbit-scan` |
//! | `glv-orbit` | prime | the closure of the smallest abscissae under the curve's automorphism group, each orbit folded onto one column (`j = 0`: 6 points a column, `j = 1728`: 4, generic: 2) |
//! | `gls-line` | GLS over `F_{p²}` | the `ψ`-stable line `x ∈ u·s·F_p`, each `⟨−1, ψ⟩`-orbit folded onto one column (4 points a column) |
//!
//! The column map is the interesting part.  `koblitz-orbit` gives a
//! whole orbit of `2n` points one unknown, which cuts the relations
//! needed by `2n`; setting its `no_fold` parameter gives the control
//! that shows what the fold buys.  That pair is the one place this
//! ledger records an *advance* in the `AGENTS.md` §3 sense rather than
//! engineering, and it is a two-word change in a sweep file.
//!
//! ## Decomposition oracles
//!
//! | name | summands | how |
//! |:--|:--|:--|
//! | `subtract` | 2 | test `R − P ∈ F` for every signed base point; no table |
//! | `mitm` | 2 or 3 | a pair table probed once per target, or once per `R − P` |
//! | `mitm-prime-batch` | 2 | the same folded prime pair table, built with one batched inversion per abscissa row and x-only sums |
//! | `mitm-frobenius` | 2 or 3 | the same on a Frobenius-folded table (Koblitz) |
//! | `mitm-frobenius-counted` | 2 or 3 | the same table and probe, with native table setup counters |
//! | `rootfind-subtract` | 2 or 3 | prime curves: `prime_fast`'s direct-subtraction root finding, every `x(R ∓ F_i)` of a target from one Montgomery batched inversion, then a lookup (`m = 2`) or a probe of the negation-folded pair table (`m = 3`) |
//!
//! The pair table is where memory buys time, and `mitm`'s
//! `negation_folded` parameter halves the additions that build it.

use std::collections::{BTreeMap, HashSet};
use std::time::{Duration, Instant};

use num_traits::ToPrimitive;

use super::large_prime::{LargePrimeCollector, LargePrimeConfig};
use super::stages::{
    BooleanSystem, DecompositionOracle, FactorBaseBuilder, InstanceCtx, Params, SolverCost,
    SolverTotals, SolverVerdict, SystemShape, SystemSolver,
};
use crate::cryptanalysis::gls_fp2::{gls_line_base, Fp2Curve, Fp2Point, GlsInstance};
use crate::cryptanalysis::glv_invariant_base::{glv_orbit_base, AutomorphismGroup};
use crate::cryptanalysis::ic_boundary::{
    binary_subspace_factor_base, decompose_mitm, decompose_mitm_frobenius, koblitz_factor_base,
    prime_factor_base, prime_factor_base_strided, BinaryGroup, BinaryInstance, ColumnFold, CountedGroup, FactorBase,
    FrobeniusPairTable, GroupOps, OracleCounters, PairTable, PrimeCurve, PrimeInstance, PrimePoint,
    Relation,
};
use crate::cryptanalysis::koblitz_factor_base_search::subspace_is_trace_zero;
use crate::cryptanalysis::koblitz_fast::FastPoint;
use crate::cryptanalysis::koblitz_fast_arith::FastBinaryCurve;
use crate::cryptanalysis::koblitz_groebner::{f4_word_ops_thread, FieldStructure, SolverEngine};
use crate::cryptanalysis::koblitz_index_calculus::build_frobenius_factor_base_from_divisor;
use crate::cryptanalysis::koblitz_symmetrised::{
    build_symmetrised_factor_base, build_symmetrised_system, frobenius_view_of_symmetrised,
    symmetrised_groebner_decompose, SymmetrisedFactorBase,
};
use crate::cryptanalysis::pq_descent_symbolic::{descend, max_n_prime, SymbolicDescent};
use crate::cryptanalysis::prime_fast::{Fp64, MAX_BITS as PRIME_FAST_MAX_BITS};

// ── Factor bases ───────────────────────────────────────────────────

/// The smallest abscissae of a prime-field curve that carry a point,
/// one column each.
pub struct PrimeAbscissaBase<'i> {
    pub instance: &'i PrimeInstance,
}

impl FactorBaseBuilder<PrimeCurve> for PrimeAbscissaBase<'_> {
    fn name(&self) -> &str {
        "prime-abscissa"
    }

    fn describe(&self, params: &Params) -> String {
        match params.get("stride") {
            None | Some("1") => format!(
                "the {} smallest abscissae carrying a point, one column each",
                params.get("size").unwrap_or("?")
            ),
            Some(stride) => format!(
                "{} abscissae carrying a point, scanned at stride {stride} modulo p, one column each",
                params.get("size").unwrap_or("?")
            ),
        }
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            ("size", "abscissae in the base; the family optimum is about #E^(1/3)"),
            ("stride", "fixed nonzero abscissa stride modulo p; default 1"),
        ]
    }

    fn build(
        &self,
        _ctx: &InstanceCtx<PrimeCurve>,
        params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<FactorBase<PrimePoint>, String> {
        let size = params.u64("size")? as usize;
        if size == 0 {
            return Err("factor base size must be positive".into());
        }
        let stride = params.get("stride").map_or(Ok(1u64), |s| {
            s.parse::<u64>().map_err(|_| "stride must be an unsigned integer".to_string())
        })?;
        if stride == 0 || stride >= self.instance.curve.p {
            return Err("stride must lie in 1..p".into());
        }
        if stride == 1 {
            Ok(prime_factor_base(self.instance, size))
        } else {
            Ok(prime_factor_base_strided(self.instance, size, stride))
        }
    }
}

/// Points of a binary curve whose abscissa lies in an `F_2`-subspace
/// of the field, one column per abscissa.
pub struct BinarySubspaceBase<'i> {
    pub instance: &'i BinaryInstance,
}

impl<'a> FactorBaseBuilder<BinaryGroup<'a>> for BinarySubspaceBase<'_> {
    fn name(&self) -> &str {
        "binary-subspace"
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "points whose abscissa lies in the {}-dimensional subspace spanned by {{1, z, …}}",
            params.get("dimension").unwrap_or("?")
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[(
            "dimension",
            "F_2-dimension of the abscissa subspace; about 2^dimension abscissae",
        )]
    }

    fn build(
        &self,
        _ctx: &InstanceCtx<BinaryGroup<'a>>,
        params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<FactorBase<FastPoint>, String> {
        let dim = params.u64("dimension")? as u32;
        if dim == 0 || dim > 30 {
            return Err(format!("subspace dimension {dim} outside 1..=30"));
        }
        let basis: Vec<u64> = (0..dim).map(|k| 1u64 << k).collect();
        Ok(binary_subspace_factor_base(self.instance, &basis))
    }
}

/// **The Frobenius closure of a prefix subspace** (lead L5, the ordered
/// τ-slot relation space): every point whose abscissa lies in
/// `σ^t(V)` for some `t`, `V = ⟨1, z, …, z^{d−1}⟩`, each signed Frobenius
/// orbit folded onto one column with weights `±λ^t`.  `V` itself is not
/// Frobenius-stable, so the closure has about `n·|F_V|` points on about
/// `|F_V|/2` columns.  The base records `V`'s basis as its
/// `subspace_basis`, so the algebraic oracle can descend with summand
/// `i` over `σ^{j_i}(V)` (`descent-algebraic:slots=all`) and the lifted
/// abscissae are abscissae of base points; the pair-table oracles run on
/// it like on any closed base.
pub struct KoblitzClosureBase<'i> {
    pub instance: &'i BinaryInstance,
}

impl<'a> FactorBaseBuilder<BinaryGroup<'a>> for KoblitzClosureBase<'_> {
    fn name(&self) -> &str {
        "koblitz-closure"
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "the signed Frobenius closure of the {}-dimensional prefix subspace, one column per orbit",
            params.get("dimension").unwrap_or("?")
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[(
            "dimension",
            "F_2-dimension of the seed subspace V = <1, z, ...>; the base is the union of the Frobenius images of V's points",
        )]
    }

    fn build(
        &self,
        ctx: &InstanceCtx<BinaryGroup<'a>>,
        params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<FactorBase<FastPoint>, String> {
        let inst = self.instance;
        let kc = inst
            .koblitz
            .as_ref()
            .ok_or("koblitz-closure requires a Koblitz curve")?;
        let dim = params.u64("dimension")? as u32;
        if dim == 0 || dim > 24 {
            return Err(format!("subspace dimension {dim} outside 1..=24"));
        }
        let n = inst.n as usize;
        let lambda = kc
            .lambda
            .to_u64()
            .ok_or("Frobenius eigenvalue does not fit in u64")?;
        let r = ctx.r;
        let basis: Vec<u64> = (0..dim).map(|k| 1u64 << k).collect();
        let g = ctx.group;
        let mut points = Vec::new();
        let mut col_of = Vec::new();
        let mut coef_of = Vec::new();
        let mut seen: HashSet<u64> = HashSet::new();
        let mut columns = 0usize;
        let (mut solves, mut frob) = (0u64, 0u64);
        for idx in 1..(1u64 << dim) {
            let mut x = 0u64;
            for (j, &e) in basis.iter().enumerate() {
                if (idx >> j) & 1 == 1 {
                    x ^= e;
                }
            }
            solves += 1;
            for seed in inst.points_with_x(x) {
                if seen.contains(&seed.pack()) {
                    continue;
                }
                // One column for the signed orbit of `seed`.
                let column = columns;
                let mut member = seed;
                let mut coefficient = 1u64;
                let mut fresh = false;
                for _ in 0..n {
                    if seen.insert(member.pack()) {
                        fresh = true;
                        points.push(member);
                        col_of.push(column);
                        coef_of.push(coefficient);
                        let neg = g.neg(member);
                        seen.insert(neg.pack());
                        points.push(neg);
                        col_of.push(column);
                        coef_of.push((r - coefficient) % r);
                    }
                    member = inst.fast.frobenius_k(member, 1);
                    frob += 1;
                    coefficient = ((coefficient as u128 * lambda as u128) % r as u128) as u64;
                }
                if fresh {
                    columns += 1;
                }
            }
        }
        if columns == 0 {
            return Err("the closure base is empty".into());
        }
        let mut fb = FactorBase::from_column_map(
            format!("Frobenius closure of the prefix subspace of dimension {dim}"),
            points,
            col_of,
            coef_of,
            columns,
            |p| p.pack(),
            |p| inst.fast.neg(*p).pack(),
            |p| p.x,
        )?;
        fb.subspace_basis = Some(basis);
        fb.dimension = Some(dim);
        fb.cost.count("abscissae_scanned", solves);
        fb.cost.count("as_solves", solves);
        fb.cost.count("frobenius_maps", frob);
        Ok(fb)
    }
}

/// A Koblitz base with each signed Frobenius orbit folded onto one
/// column: `2n` points, one unknown.
pub struct KoblitzOrbitBase<'i> {
    pub instance: &'i BinaryInstance,
}

/// A signed-Frobenius-orbit base whose abscissa subspace is contained in
/// `ker Tr_{F_{2^n}/F_2}`.  This is an executable factor-base policy rather
/// than merely a score attached to a search candidate: a divisor that is not
/// trace-zero is rejected before relation collection starts.
pub struct KoblitzTraceZeroBase<'i> {
    pub instance: &'i BinaryInstance,
}

impl<'a> FactorBaseBuilder<BinaryGroup<'a>> for KoblitzTraceZeroBase<'_> {
    fn name(&self) -> &str {
        "koblitz-trace-zero"
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "the trace-zero Frobenius-invariant abscissa subspace selected by divisor {}, each signed orbit folded onto one column",
            params.get("divisor").unwrap_or("?")
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            (
                "divisor",
                "`;`-separated factors of x^n - 1; the resulting invariant subspace must lie in the absolute-trace kernel",
            ),
            (
                "projected",
                "1 to fold by the signed Frobenius orbits of the cofactor projections [h]P instead of \
                 of the raw points: orbits whose projections coincide or are Frobenius translates \
                 share one column (on a base closed under +T, P and P + T do), which is the column \
                 count the relation matrix actually has to determine",
            ),
        ]
    }

    fn build(
        &self,
        _ctx: &InstanceCtx<BinaryGroup<'a>>,
        params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<FactorBase<FastPoint>, String> {
        let kc = self
            .instance
            .koblitz
            .as_ref()
            .ok_or("koblitz-trace-zero requires a Koblitz curve")?;
        let idx = divisor_indices(params)?;
        let frob = build_frobenius_factor_base_from_divisor(kc, &idx)
            .ok_or_else(|| format!("no invariant subspace for divisor {idx:?}"))?;
        if !subspace_is_trace_zero(kc, &frob) {
            return Err(format!(
                "divisor {idx:?} is not trace-zero: at least one abscissa has absolute trace one"
            ));
        }
        let fold = if params.flag("projected") {
            ColumnFold::ProjectedSignedFrobeniusOrbit
        } else {
            ColumnFold::SignedFrobeniusOrbit
        };
        koblitz_factor_base(
            self.instance,
            &frob,
            fold,
            format!(
                "trace-zero invariant subspace, divisor {idx:?}, dimension {}{}",
                frob.ell,
                if params.flag("projected") {
                    ", columns merged by cofactor projection"
                } else {
                    ""
                }
            ),
        )
        .ok_or_else(|| "the trace-zero subspace produced no usable factor base".into())
    }
}

impl<'a> FactorBaseBuilder<BinaryGroup<'a>> for KoblitzOrbitBase<'_> {
    fn name(&self) -> &str {
        "koblitz-orbit"
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "the Frobenius-invariant subspace spanned by divisor factors {}, {}",
            params.get("divisor").unwrap_or("?"),
            if params.flag("no_fold") {
                "one column per abscissa (the unfolded control)"
            } else {
                "each signed Frobenius orbit folded onto one column"
            }
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            (
                "divisor",
                "`;`-separated indices into the factorisation of the Frobenius characteristic \
                 polynomial, e.g. `divisor=1;2` (a `,` would end the parameter on the command \
                 line); the invariant subspace they span becomes the abscissa set. This is \
                 the knob that chooses the base, so sweeping it sweeps factor bases",
            ),
            (
                "no_fold",
                "1 to give every abscissa its own column: the control that shows what the fold buys",
            ),
            (
                "projected",
                "1 to fold by the signed Frobenius orbits of the cofactor projections [h]P instead of \
                 of the raw points: orbits whose projections coincide or are Frobenius translates \
                 share one column (on a base closed under +T, P and P + T do), which is the column \
                 count the relation matrix actually has to determine",
            ),
        ]
    }

    fn build(
        &self,
        _ctx: &InstanceCtx<BinaryGroup<'a>>,
        params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<FactorBase<FastPoint>, String> {
        let kc = self
            .instance
            .koblitz
            .as_ref()
            .ok_or("this instance is not a Koblitz curve, so it has no Frobenius to fold by")?;
        let idx = divisor_indices(params)?;
        let frob = build_frobenius_factor_base_from_divisor(kc, &idx)
            .ok_or_else(|| format!("no invariant subspace for divisor {idx:?} on this curve"))?;
        let fold = if params.flag("no_fold") {
            ColumnFold::Abscissa
        } else if params.flag("projected") {
            ColumnFold::ProjectedSignedFrobeniusOrbit
        } else {
            ColumnFold::SignedFrobeniusOrbit
        };
        koblitz_factor_base(
            self.instance,
            &frob,
            fold,
            format!(
                "invariant subspace, divisor {idx:?}, dimension {}",
                frob.ell
            ),
        )
        .ok_or_else(|| "the invariant subspace produced no usable factor base".into())
    }
}

/// The archived compact-orbit producer's public, deterministic raw-x scan.
/// All group work is charged to `ops`; native field and container work is
/// retained on the base as priced counters or explicit lower-bound gaps.
pub struct CompactOrbitScanBase<'i> {
    pub instance: &'i BinaryInstance,
}

impl<'a> FactorBaseBuilder<BinaryGroup<'a>> for CompactOrbitScanBase<'_> {
    fn name(&self) -> &str {
        "compact-orbit-scan"
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "the first {} distinct projected signed Frobenius orbits from raw x <= {}",
            params.get("columns").unwrap_or("?"),
            params.get("raw_x_cap").unwrap_or("?")
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            ("columns", "number of distinct signed Frobenius orbits"),
            (
                "raw_x_cap",
                "inclusive maximum raw abscissa; exhaustion is an error",
            ),
            (
                "norm",
                "Z[tau]-multiplier closure (lead A1): also take alpha(P) for every alpha = a + b*tau of norm a^2 - ab + 2b^2 <= norm, with weight (a + b*lambda) on P's column; 0 or absent for the plain {+-tau^k} orbits",
            ),
        ]
    }

    fn build(
        &self,
        ctx: &InstanceCtx<BinaryGroup<'a>>,
        params: &Params,
        ops: &mut GroupOps,
    ) -> Result<FactorBase<FastPoint>, String> {
        let inst = self.instance;
        let kc = inst
            .koblitz
            .as_ref()
            .ok_or("compact-orbit-scan requires a Koblitz curve")?;
        if ctx.r != inst.r || ctx.cofactor != inst.cofactor || ctx.field_degree != Some(inst.n) {
            return Err("compact-orbit-scan context does not match its binary instance".into());
        }
        let columns =
            usize::try_from(params.u64("columns")?).map_err(|_| "columns does not fit in usize")?;
        let raw_x_cap = params.u64("raw_x_cap")?;
        let norm_cap = match params.get("norm") {
            None => 0i64,
            Some(_) => params.u64("norm")? as i64,
        };
        // Multipliers a + b*tau of norm a^2 - ab + 2b^2 <= norm_cap, one per
        // associate class up to sign and Frobenius (the orbit loop below
        // supplies +-tau^k alpha(P); duplicates are dropped by key).
        let mut multipliers: Vec<(i64, i64)> = Vec::new();
        if norm_cap > 0 {
            let bound = (norm_cap as f64).sqrt() as i64 + 2;
            for a in -bound..=bound {
                for b in -bound..=bound {
                    if (a, b) == (0, 0) || a * a - a * b + 2 * b * b > norm_cap {
                        continue;
                    }
                    multipliers.push((a, b));
                }
            }
        }
        let mut multiplier_points = 0u64;
        if columns == 0 || raw_x_cap == 0 || inst.n == 0 || inst.n >= 64 || inst.n.is_multiple_of(2)
        {
            return Err(
                "compact-orbit-scan requires positive columns and cap and odd degree < 64".into(),
            );
        }
        // The archived source selects the half-trace root first. The
        // framework's linear Artin-Schreier solver can choose its other
        // valid root, changing representatives, point order and base hash.
        let source_curve = FastBinaryCurve::new(&inst.irreducible, inst.a)
            .ok_or("compact-orbit source curve is outside its field limit")?;
        let orbit_size = (inst.n as usize)
            .checked_mul(2)
            .ok_or("orbit size overflows usize")?;
        let capacity = columns
            .checked_mul(orbit_size)
            .ok_or("factor base size overflows usize")?;
        let lambda = kc
            .lambda
            .to_u64()
            .ok_or("Frobenius eigenvalue does not fit in u64")?;
        let g = ctx.group;
        let mut seen = HashSet::new();
        let mut points = Vec::with_capacity(capacity);
        let mut col_of = Vec::with_capacity(capacity);
        let mut coef_of = Vec::with_capacity(capacity);
        let mut scanned = 0u64;
        let mut lifts = 0u64;
        let mut projected = 0u64;
        let mut key_squares = 0u64;
        let mut frobenius_maps = 0u64;
        let mut equation_checks = 0u64;
        let max_x = raw_x_cap.min((1u64 << inst.n) - 1);
        for raw_x in 1..=max_x {
            scanned += 1;
            for source_point in source_curve.points_with_x(inst.b, raw_x) {
                let (x, y) = source_point.ok_or("raw-x lift unexpectedly returned infinity")?;
                let point = FastPoint::affine(x, y);
                lifts += 1;
                let seed = g.mul(ops, point, ctx.cofactor);
                if g.is_identity(&seed) || seed.x <= 1 {
                    continue;
                }
                projected += 1;
                equation_checks += 1;
                if !inst.fast.is_on_curve(seed) {
                    return Err(format!("projected point from raw x {raw_x} is off curve"));
                }
                let mut x = seed.x;
                let mut key = x;
                for _ in 1..inst.n {
                    x = inst.gf.sqr(x);
                    key = key.min(x);
                    key_squares += 1;
                }
                if !seen.insert(key) {
                    continue;
                }
                if !g.is_identity(&g.mul(ops, seed, ctx.r)) {
                    return Err(format!(
                        "projected point from raw x {raw_x} has wrong order"
                    ));
                }
                if g.mul(ops, seed, lambda) != inst.fast.frobenius_k(seed, 1) {
                    return Err(format!(
                        "Frobenius eigenvalue check failed for raw x {raw_x}"
                    ));
                }
                frobenius_maps += 1;
                let column = seen.len() - 1;
                let mut member = seed;
                let mut coefficient = 1u64;
                for _ in 0..inst.n {
                    equation_checks += 1;
                    if !inst.fast.is_on_curve(member) {
                        return Err(format!("Frobenius image from raw x {raw_x} is off curve"));
                    }
                    points.push(member);
                    col_of.push(column);
                    coef_of.push(coefficient);
                    points.push(g.neg(member));
                    col_of.push(column);
                    coef_of.push((ctx.r - coefficient) % ctx.r);
                    member = inst.fast.frobenius_k(member, 1);
                    frobenius_maps += 1;
                    coefficient = ((coefficient as u128 * lambda as u128) % ctx.r as u128) as u64;
                }
                if member != seed || coefficient != 1 {
                    return Err(format!("Frobenius orbit from raw x {raw_x} is not closed"));
                }
                // Lead A1: the orbits of alpha(P) = [a]P + [b]tau(P) for the
                // small-norm multipliers, weight (a + b*lambda)*lambda^t.
                if !multipliers.is_empty() {
                    let tau_seed = inst.fast.frobenius_k(seed, 1);
                    let r128 = ctx.r as u128;
                    let to_scalar = |v: i64| -> u64 { v.rem_euclid(ctx.r as i64) as u64 };
                    let mut local_seen: HashSet<u64> = HashSet::new();
                    for &(a, b) in &multipliers {
                        let sa = to_scalar(a);
                        let sb = to_scalar(b);
                        let pa = if sa == 0 {
                            g.identity()
                        } else {
                            g.mul(ops, seed, sa)
                        };
                        let pb = if sb == 0 {
                            g.identity()
                        } else {
                            g.mul(ops, tau_seed, sb)
                        };
                        let image = g.add(ops, pa, pb);
                        if g.is_identity(&image) {
                            continue;
                        }
                        // Skip images already in the base (an associate of a
                        // Frobenius power, or a repeat), by orbit key.
                        let mut xk = image.x;
                        let mut key = xk;
                        for _ in 1..inst.n {
                            xk = inst.gf.sqr(xk);
                            key = key.min(xk);
                        }
                        if !local_seen.insert(key) || seen.contains(&key) {
                            continue;
                        }
                        let weight =
                            ((sa as u128 + (sb as u128 * lambda as u128) % r128) % r128) as u64;
                        let mut member = image;
                        let mut coefficient = weight;
                        for _ in 0..inst.n {
                            points.push(member);
                            col_of.push(column);
                            coef_of.push(coefficient);
                            points.push(g.neg(member));
                            col_of.push(column);
                            coef_of.push((ctx.r - coefficient) % ctx.r);
                            member = inst.fast.frobenius_k(member, 1);
                            frobenius_maps += 1;
                            coefficient = ((coefficient as u128 * lambda as u128) % r128) as u64;
                            multiplier_points += 2;
                        }
                        if member != image {
                            return Err(format!(
                                "multiplier orbit from raw x {raw_x} is not closed"
                            ));
                        }
                    }
                }
                if seen.len() == columns {
                    break;
                }
            }
            if seen.len() == columns {
                break;
            }
        }
        if seen.len() != columns {
            return Err(format!(
                "raw-x cap {raw_x_cap} exhausted after {scanned} abscissae with {} of {columns} orbits",
                seen.len()
            ));
        }
        let mut fb = FactorBase::from_column_map(
            format!("first {columns} projected signed Frobenius orbits, raw x <= {raw_x_cap}"),
            points,
            col_of,
            coef_of,
            columns,
            |p| p.pack(),
            |p| inst.fast.neg(*p).pack(),
            |p| p.x,
        )?;
        fb.cost.count("abscissae_scanned_uncharged", scanned);
        fb.cost.count("inversions", scanned);
        fb.cost.count("legacy_as_solves_uncharged", scanned);
        fb.cost.count("point_lifts_uncharged", lifts);
        fb.cost.count(
            "lift_field_squares_uncharged",
            scanned * inst.n as u64 + lifts / 2 * (inst.n as u64 - 1),
        );
        fb.cost
            .count("lift_field_multiplies_uncharged", scanned + lifts / 2);
        fb.cost.count("cofactor_projected_points", projected);
        fb.cost
            .count("curve_equation_checks_uncharged", equation_checks);
        fb.cost
            .count("min_key_field_squares_uncharged", key_squares);
        fb.cost.count("frobenius_maps", frobenius_maps);
        // Only with multipliers: a plain scan's counters (and the frozen
        // evidence built from them) stay exactly what they were.
        if !multipliers.is_empty() {
            fb.cost.count("multiplier_points", multiplier_points);
        }
        fb.cost.count(
            "hash_probes_uncharged",
            projected + fb.points.len() as u64 * 3,
        );
        fb.cost.count(
            "point_negations_uncharged",
            lifts / 2 + columns as u64 * inst.n as u64 + fb.points.len() as u64,
        );
        fb.cost
            .count("factor_point_inserts_uncharged", fb.points.len() as u64);
        fb.cost.count(
            "factor_base_vec_bytes_floor_uncharged",
            (fb.points.len() * std::mem::size_of::<FastPoint>()
                + fb.col_of.len() * std::mem::size_of::<usize>()
                + fb.coef_of.len() * std::mem::size_of::<u64>()
                + fb.neg_index.len() * std::mem::size_of::<usize>()) as u64,
        );
        Ok(fb)
    }
}

/// **The cyclotomic-sparse base** (lead A2): abscissae `x = Σ_{k∈S} γ^k`
/// with `|S| = w`, `γ` a primitive `(2n+1)`-th root of unity — the second
/// Frobenius-stable family at prime `n`, with a polynomial membership
/// description (`cyclotomic_factor_base`).  Each abscissa's rational
/// points are projected into the subgroup by the cofactor and the
/// projected signed Frobenius orbits are folded onto one column each
/// with weights `±λ^t`, exactly as `compact-orbit-scan` does for raw
/// abscissae.  Needs `2n + 1` prime with `ord_{2n+1}(2) = n`
/// (`n ∈ {11, 23, 35, 39, …}` in the single-word range).
pub struct CyclotomicBase<'i> {
    pub instance: &'i BinaryInstance,
}

impl<'a> FactorBaseBuilder<BinaryGroup<'a>> for CyclotomicBase<'_> {
    fn name(&self) -> &str {
        "cyclotomic"
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "projected signed Frobenius orbits of the gamma-weight-{} abscissae (gamma a primitive (2n+1)-th root of unity)",
            params.get("w").unwrap_or("?")
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            (
                "w",
                "gamma-weight of the abscissae: x = sum of w distinct powers of gamma; the set has C(2n+1, w) elements",
            ),
            (
                "columns",
                "optional cap on the number of distinct projected signed orbits, in canonical abscissa order; 0 or absent takes them all",
            ),
        ]
    }

    fn build(
        &self,
        ctx: &InstanceCtx<BinaryGroup<'a>>,
        params: &Params,
        ops: &mut GroupOps,
    ) -> Result<FactorBase<FastPoint>, String> {
        use crate::cryptanalysis::cyclotomic_factor_base::CyclotomicField;
        let inst = self.instance;
        let kc = inst
            .koblitz
            .as_ref()
            .ok_or("cyclotomic requires a Koblitz curve")?;
        if ctx.r != inst.r || ctx.cofactor != inst.cofactor || ctx.field_degree != Some(inst.n) {
            return Err("cyclotomic context does not match its binary instance".into());
        }
        let w = usize::try_from(params.u64("w")?).map_err(|_| "w does not fit in usize")?;
        let column_cap = match params.get("columns") {
            None => 0usize,
            Some(_) => usize::try_from(params.u64("columns")?)
                .map_err(|_| "columns does not fit in usize")?,
        };
        if inst.n == 0 || inst.n >= 64 || inst.n.is_multiple_of(2) {
            return Err("cyclotomic requires an odd degree < 64".into());
        }
        let field = CyclotomicField::new(&inst.irreducible, inst.n)?;
        let abscissae = field.gamma_weight_set(w)?;
        let lambda = kc
            .lambda
            .to_u64()
            .ok_or("Frobenius eigenvalue does not fit in u64")?;
        let g = ctx.group;
        let mut seen = HashSet::new();
        let mut points = Vec::new();
        let mut col_of = Vec::new();
        let mut coef_of = Vec::new();
        let (mut lifts, mut projected, mut frobenius_maps) = (0u64, 0u64, 0u64);
        'outer: for x in &abscissae {
            let raw_x = x.raw_bits().first().copied().unwrap_or(0);
            for point in inst.points_with_x(raw_x) {
                lifts += 1;
                let seed = g.mul(ops, point, ctx.cofactor);
                if g.is_identity(&seed) {
                    continue;
                }
                projected += 1;
                let mut xk = seed.x;
                let mut key = xk;
                for _ in 1..inst.n {
                    xk = inst.gf.sqr(xk);
                    key = key.min(xk);
                }
                if !seen.insert(key) {
                    continue;
                }
                if !g.is_identity(&g.mul(ops, seed, ctx.r)) {
                    return Err(format!("projected point from x {raw_x} has wrong order"));
                }
                let column = seen.len() - 1;
                let mut member = seed;
                let mut coefficient = 1u64;
                for _ in 0..inst.n {
                    if !inst.fast.is_on_curve(member) {
                        return Err(format!("Frobenius image from x {raw_x} is off curve"));
                    }
                    points.push(member);
                    col_of.push(column);
                    coef_of.push(coefficient);
                    points.push(g.neg(member));
                    col_of.push(column);
                    coef_of.push((ctx.r - coefficient) % ctx.r);
                    member = inst.fast.frobenius_k(member, 1);
                    frobenius_maps += 1;
                    coefficient = ((coefficient as u128 * lambda as u128) % ctx.r as u128) as u64;
                }
                if member != seed || coefficient != 1 {
                    return Err(format!("Frobenius orbit from x {raw_x} is not closed"));
                }
                if column_cap > 0 && seen.len() == column_cap {
                    break 'outer;
                }
            }
        }
        let columns = seen.len();
        if columns == 0 {
            return Err(format!(
                "no gamma-weight-{w} abscissa carries a subgroup point on {}",
                inst.name
            ));
        }
        let mut fb = FactorBase::from_column_map(
            format!(
                "{columns} projected signed Frobenius orbits of the {} gamma-weight-{w} abscissae",
                abscissae.len()
            ),
            points,
            col_of,
            coef_of,
            columns,
            |p| p.pack(),
            |p| inst.fast.neg(*p).pack(),
            |p| p.x,
        )?;
        fb.cost.count("abscissae_enumerated_uncharged", abscissae.len() as u64);
        fb.cost.count("point_lifts_uncharged", lifts);
        fb.cost.count("cofactor_projected_points", projected);
        fb.cost.count("frobenius_maps", frobenius_maps);
        Ok(fb)
    }
}

/// The prime-field analogue of `koblitz-orbit`: the `size` smallest
/// abscissae, closed under the curve's automorphism group and folded
/// one orbit to a column.  `group` picks the group (`auto` takes what
/// the curve has); `no_fold` keeps the closed point set and folds by
/// negation only, which is the control that shows what the fold buys
/// on the same points.  See `glv_invariant_base`.
pub struct GlvOrbitBase<'i> {
    pub instance: &'i PrimeInstance,
}

impl FactorBaseBuilder<PrimeCurve> for GlvOrbitBase<'_> {
    fn name(&self) -> &str {
        "glv-orbit"
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "the closure of the {} smallest abscissae under the {} automorphism group, {}",
            params.get("size").unwrap_or("?"),
            params.get("group").unwrap_or("auto"),
            if params.flag("no_fold") {
                "folded by negation only (the control)"
            } else {
                "each orbit folded onto one column"
            }
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            (
                "size",
                "seed abscissae; the base is their closure under the group (3× as many on j = 0, 2× on j = 1728)",
            ),
            (
                "group",
                "auto (what the curve has), negation, j0 or j1728",
            ),
            (
                "no_fold",
                "1 to fold the same points by negation only: the control that shows what the fold buys",
            ),
        ]
    }

    fn build(
        &self,
        _ctx: &InstanceCtx<PrimeCurve>,
        params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<FactorBase<PrimePoint>, String> {
        let size = params.u64("size")? as usize;
        if size == 0 {
            return Err("factor base size must be positive".into());
        }
        let group = AutomorphismGroup::parse(params.get("group").unwrap_or("auto"))?;
        let (fb, _) = glv_orbit_base(self.instance, size, group, !params.flag("no_fold"))?;
        Ok(fb)
    }
}

/// The `ψ`-stable line of a GLS twist over `F_{p²}`, folded by
/// `⟨−1, ψ⟩`; `no_fold` folds the same points by negation only.  See
/// `gls_fp2`.
pub struct GlsLineBase<'i> {
    pub instance: &'i GlsInstance,
}

impl FactorBaseBuilder<Fp2Curve> for GlsLineBase<'_> {
    fn name(&self) -> &str {
        "gls-line"
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "the ψ-stable line x ∈ u·s·F_p, {}",
            if params.flag("no_fold") {
                "folded by negation only (the control)"
            } else {
                "each ⟨−1, ψ⟩-orbit folded onto one column"
            }
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[(
            "no_fold",
            "1 to fold the same points by negation only: the control that shows what the fold buys",
        )]
    }

    fn build(
        &self,
        _ctx: &InstanceCtx<Fp2Curve>,
        params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<FactorBase<Fp2Point>, String> {
        let (fb, _) = gls_line_base(self.instance, !params.flag("no_fold"))?;
        Ok(fb)
    }
}

// ── Decomposition oracles ──────────────────────────────────────────

/// Test `R − P ∈ F` for every signed base point: no table, `|F|` group
/// operations a target, and the honest baseline every table oracle is
/// an improvement on.
pub struct SubtractOracle;

impl<G: CountedGroup> DecompositionOracle<G> for SubtractOracle {
    fn name(&self) -> &str {
        "subtract"
    }

    fn summands(&self) -> u32 {
        2
    }

    fn describe(&self, _params: &Params) -> String {
        "for each signed base point P, test whether R − P is in the base".into()
    }

    fn decompose(
        &mut self,
        ctx: &InstanceCtx<G>,
        fb: &FactorBase<G::Elt>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: G::Elt,
    ) -> Option<Vec<usize>> {
        for (i, p) in fb.points.iter().enumerate() {
            let rest = ctx.group.add(ops, point, ctx.group.neg(*p));
            counters.lookups += 1;
            let key = ctx.group.key(&rest);
            if let Some(j) = fb.index_of_key(key) {
                return Some(vec![i, j]);
            }
        }
        None
    }
}

/// Meet in the middle over a pair table: memory for time, the lever
/// that makes a two-summand search `O(1)` a target instead of `O(|F|)`.
pub struct MitmOracle {
    m: u32,
    table: Option<PairTable>,
}

impl MitmOracle {
    pub fn new(m: u32) -> Self {
        Self { m, table: None }
    }
}

impl<G: CountedGroup> DecompositionOracle<G> for MitmOracle {
    fn name(&self) -> &str {
        "mitm"
    }

    fn summands(&self) -> u32 {
        self.m
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "a table of pair sums probed once per target (m = 2) or once per R − P (m = 3), {}",
            if params.flag("negation_folded") {
                "built with the negation fold"
            } else {
                "built in full"
            }
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[(
            "negation_folded",
            "1 to build only one of each ±pair, halving the table's additions",
        )]
    }

    fn prepare(
        &mut self,
        ctx: &InstanceCtx<G>,
        fb: &FactorBase<G::Elt>,
        params: &Params,
        ops: &mut GroupOps,
    ) -> Result<(), String> {
        let table = if params.flag("negation_folded") {
            PairTable::build_negation_folded(ctx.group, fb)
        } else {
            PairTable::build(ctx.group, fb)
        };
        // The table's construction is part of the pipeline's cost.
        ops.merge(table.build_ops);
        self.table = Some(table);
        Ok(())
    }

    fn decompose(
        &mut self,
        ctx: &InstanceCtx<G>,
        fb: &FactorBase<G::Elt>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: G::Elt,
    ) -> Option<Vec<usize>> {
        let table = self.table.as_ref()?;
        decompose_mitm(ctx.group, fb, table, self.m, ops, counters, point)
    }
}

/// The two-summand prime MITM search with the same probe and table keys as
/// `mitm`, but a batched-inversion, x-only table construction.
pub struct PrimeBatchMitmOracle {
    table: Option<PairTable>,
}

impl PrimeBatchMitmOracle {
    pub fn new() -> Self {
        Self { table: None }
    }
}

impl DecompositionOracle<PrimeCurve> for PrimeBatchMitmOracle {
    fn name(&self) -> &str {
        "mitm-prime-batch"
    }

    fn summands(&self) -> u32 {
        2
    }

    fn describe(&self, _params: &Params) -> String {
        "negation-folded pair table built with batched prime-field inversions; identical MITM probe"
            .into()
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[("m", "must be 2")]
    }

    fn prepare(
        &mut self,
        ctx: &InstanceCtx<PrimeCurve>,
        fb: &FactorBase<PrimePoint>,
        params: &Params,
        ops: &mut GroupOps,
    ) -> Result<(), String> {
        if params.u64_or("m", 2)? != 2 || params.0.keys().any(|key| key != "m") {
            return Err("mitm-prime-batch accepts only m=2".into());
        }
        let table = PairTable::build_negation_folded_prime_batch(ctx.group, fb)?;
        ops.merge(table.build_ops);
        self.table = Some(table);
        Ok(())
    }

    fn decompose(
        &mut self,
        ctx: &InstanceCtx<PrimeCurve>,
        fb: &FactorBase<PrimePoint>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: PrimePoint,
    ) -> Option<Vec<usize>> {
        decompose_mitm(ctx.group, fb, self.table.as_ref()?, 2, ops, counters, point)
    }
}

/// One abscissa of the base as the root-finding sweep sees it: both
/// coordinates of its positive representative in Montgomery form, and the
/// indices of `P` and `−P` in the run's base.
#[derive(Clone, Copy, Debug)]
struct SweepRep {
    x: u64,
    y: u64,
    plus: usize,
    minus: usize,
}

/// Direct-subtraction root finding on a prime curve: the decomposition
/// oracle of `prime_fast::collect_relations`
/// (`RESEARCH_SEMAEV_ROOTFIND_VS_GROEBNER.md` §2), inside `ic.pipeline`.
///
/// The Semaev `S₃` root in the third variable is literally `x(R ∓ F_i)`,
/// so a two-summand sweep over the base needs no summation polynomial,
/// no Gröbner basis and no square root: the `x_R − x_i` of every
/// abscissa are inverted together by Montgomery's trick (one field
/// inversion a target), `λ` for `R − F_i` and `R + F_i` share that
/// inverse, and each resulting point is looked up in the base (`m = 2`)
/// or probed against the negation-folded pair table (`m = 3`).
///
/// It is the `subtract` oracle's search — the same `2·|F_x|` candidate
/// points a target, found in the same order — on `prime_fast`'s
/// arithmetic.  It is charged exactly as `subtract` is: one group addition
/// per candidate point, so the two arms read the same `S` and differ in
/// the constant the unit does not see (one batched inversion against
/// `2·|F_x|`), which the record keeps as `inversions` (counted, priced
/// only where the curve has a pinned ratio) and as wall time.  Every hit
/// is re-added in the group on the ledger before it is returned.
pub struct RootfindSubtractOracle<'i> {
    m: u32,
    instance: &'i PrimeInstance,
    field: Option<Fp64>,
    reps: Vec<SweepRep>,
    table: Option<PairTable>,
    diffs: Vec<u64>,
    degenerate: Vec<bool>,
    scratch: Vec<u64>,
    /// Hits whose summands did not re-add to the target in the group.
    /// Must stay zero.
    pub mismatches: u64,
}

impl<'i> RootfindSubtractOracle<'i> {
    pub fn new(m: u32, instance: &'i PrimeInstance) -> Self {
        Self {
            m,
            instance,
            field: None,
            reps: Vec::new(),
            table: None,
            diffs: Vec::new(),
            degenerate: Vec::new(),
            scratch: Vec::new(),
            mismatches: 0,
        }
    }
}

impl DecompositionOracle<PrimeCurve> for RootfindSubtractOracle<'_> {
    fn name(&self) -> &str {
        "rootfind-subtract"
    }

    fn summands(&self) -> u32 {
        self.m
    }

    fn describe(&self, _params: &Params) -> String {
        match self.m {
            2 => "every x(R ∓ F_i) from one Montgomery batched inversion, each looked up in the \
                  base (prime_fast's direct-subtraction root finding)"
                .into(),
            _ => "every R ∓ F_i from one Montgomery batched inversion, each probed against the \
                  negation-folded pair table (prime_fast's direct-subtraction root finding)"
                .into(),
        }
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[("m", "summands, 2 or 3; 3 builds the negation-folded pair table")]
    }

    fn prepare(
        &mut self,
        ctx: &InstanceCtx<PrimeCurve>,
        fb: &FactorBase<PrimePoint>,
        params: &Params,
        ops: &mut GroupOps,
    ) -> Result<(), String> {
        if !(2..=3).contains(&self.m) {
            return Err(format!("rootfind-subtract takes m = 2 or 3, got {}", self.m));
        }
        for key in params.0.keys() {
            if key != "m" {
                return Err(format!(
                    "rootfind-subtract reads only `m`; `{key}` would change the method's identity \
                     without changing the search"
                ));
            }
        }
        let p = self.instance.curve.p;
        if p != ctx.group.p {
            return Err("rootfind-subtract was built for another curve".into());
        }
        if p < 5 || p & 1 == 0 || p >= (1u64 << PRIME_FAST_MAX_BITS) {
            return Err(format!(
                "rootfind-subtract needs an odd prime 3 < p < 2^{PRIME_FAST_MAX_BITS}"
            ));
        }
        let f = Fp64::new(p);
        let mut reps = Vec::with_capacity(fb.abscissa_list.len());
        for x in &fb.abscissa_list {
            let idx = fb
                .x_index
                .get(x)
                .ok_or("every listed abscissa is indexed")?;
            let plus = *idx.first().ok_or("an abscissa with no point")?;
            let minus = idx.get(1).copied().unwrap_or(plus);
            let pt = fb.points[plus];
            if pt.infinity || pt.x != *x {
                return Err("the base's abscissa index does not match its points".into());
            }
            reps.push(SweepRep {
                x: f.to_mont(pt.x),
                y: f.to_mont(pt.y),
                plus,
                minus,
            });
        }
        self.diffs = vec![0; reps.len()];
        self.degenerate = vec![false; reps.len()];
        self.reps = reps;
        self.field = Some(f);
        if self.m == 3 {
            let table = PairTable::build_negation_folded(ctx.group, fb);
            ops.merge(table.build_ops);
            self.table = Some(table);
        }
        Ok(())
    }

    fn decompose(
        &mut self,
        ctx: &InstanceCtx<PrimeCurve>,
        fb: &FactorBase<PrimePoint>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: PrimePoint,
    ) -> Option<Vec<usize>> {
        let f = self.field?;
        if point.infinity {
            return None;
        }
        let g = ctx.group;
        let xr = f.to_mont(point.x);
        let yr = f.to_mont(point.y);
        ops.field_ops.muls += 2;
        // One inversion for the whole sweep: `(x_R − x_i)⁻¹` for every
        // abscissa.  An abscissa equal to `x_R` means `R = ±F_i`; that
        // row is taken on the ledger below, exactly as `subtract` takes
        // it, so the two oracles examine the same candidates in the same
        // order whatever the target.
        let k = self.reps.len();
        for (i, rep) in self.reps.iter().enumerate() {
            let d = f.sub(xr, rep.x);
            self.degenerate[i] = d == 0;
            self.diffs[i] = if d == 0 { f.one } else { d };
        }
        f.batch_inv(&mut self.diffs, &mut self.scratch);
        counters.inversions += 1;
        ops.field_ops.invs += 1;
        ops.field_ops.muls += 3 * k as u64;
        for i in 0..k {
            let rep = self.reps[i];
            if self.degenerate[i] {
                for summand in [rep.plus, rep.minus] {
                    let t = g.add(ops, point, g.neg(fb.points[summand]));
                    if let Some(found) = self.lookup(ctx, fb, ops, counters, point, summand, t) {
                        return Some(found);
                    }
                }
                continue;
            }
            let inv = self.diffs[i];
            // `T = R − ε·F_i`: ε = +1 adds `−F_i = (x_i, −y_i)`, so
            // `λ = (−y_i − y_R)/(x_i − x_R) = (y_R + y_i)/(x_R − x_i)`;
            // ε = −1 adds `F_i`, `λ = (y_R − y_i)/(x_R − x_i)`.  Each `T`
            // is one group addition's result and is charged as one; its
            // field work (two multiplications, one squaring, two
            // conversions out of Montgomery form counted as
            // multiplications) is tallied beside the ledger's.
            let lam_minus = f.mul(f.add(yr, rep.y), inv);
            let lam_plus = f.mul(f.sub(yr, rep.y), inv);
            for (summand, lam) in [(rep.plus, lam_minus), (rep.minus, lam_plus)] {
                let x3 = f.sub(f.sub(f.sqr(lam), xr), rep.x);
                let y3 = f.sub(f.mul(lam, f.sub(xr, x3)), yr);
                ops.adds += 1;
                ops.field_ops.muls += 4;
                ops.field_ops.sqrs += 1;
                let t = PrimePoint::affine(f.from_mont(x3), f.from_mont(y3));
                if let Some(found) = self.lookup(ctx, fb, ops, counters, point, summand, t) {
                    return Some(found);
                }
            }
        }
        None
    }
}

impl RootfindSubtractOracle<'_> {
    /// Is `t = R − F_summand` a base point (`m = 2`), or a pair sum
    /// (`m = 3`)?  A hit is re-added in the group before it is returned.
    fn lookup(
        &mut self,
        ctx: &InstanceCtx<PrimeCurve>,
        fb: &FactorBase<PrimePoint>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: PrimePoint,
        summand: usize,
        t: PrimePoint,
    ) -> Option<Vec<usize>> {
        let g = ctx.group;
        if g.is_identity(&t) {
            counters.lookups += 1;
            return None;
        }
        match self.m {
            2 => {
                counters.lookups += 1;
                let j = fb.index_of_key(g.key(&t))?;
                // Fail closed: `F_summand + T` must be `R` in the group,
                // on the ledger.
                if g.add(ops, fb.points[summand], fb.points[j]) == point {
                    return Some(vec![summand, j]);
                }
                self.mismatches += 1;
                counters.frobfold_mismatches += 1;
                None
            }
            _ => {
                let table = self.table.as_ref()?;
                let (j, k) = table.probe(g, fb, ops, &mut counters.lookups, t)?;
                Some(vec![summand, j, k])
            }
        }
    }
}

/// Meet in the middle on a table folded by the signed Frobenius: the
/// same probe over a table `2n` times smaller.
///
/// With `large_primes=1` or `2` (and `m=3`) the oracle runs the same
/// probe loop through the large-prime collector
/// ([`LargePrimeCollector`]), which keeps partial relations whose
/// residual lies in a bounded domain `D` of `lp_domain_factor · |F|`
/// signed points and returns the complete relations they merge into;
/// see `ic_framework::large_prime`.
pub struct FrobeniusMitmOracle<'i> {
    m: u32,
    instance: &'i BinaryInstance,
    table: Option<FrobeniusPairTable>,
    large: Option<LargePrimeCollector<'i>>,
    count_setup: bool,
}

impl<'i> FrobeniusMitmOracle<'i> {
    pub fn new(m: u32, instance: &'i BinaryInstance) -> Self {
        Self {
            m,
            instance,
            table: None,
            large: None,
            count_setup: false,
        }
    }

    pub fn new_counted(m: u32, instance: &'i BinaryInstance) -> Self {
        Self {
            m,
            instance,
            table: None,
            large: None,
            count_setup: true,
        }
    }

    /// The large-prime collector of this run, when `large_primes` is set.
    pub fn large_prime_collector(&self) -> Option<&LargePrimeCollector<'i>> {
        self.large.as_ref()
    }

    /// The first error the collector's exact eliminator returned, if any.
    pub fn large_prime_error(&self) -> Option<&str> {
        self.large.as_ref().and_then(|lp| lp.last_error())
    }
}

impl<'a> DecompositionOracle<BinaryGroup<'a>> for FrobeniusMitmOracle<'_> {
    fn name(&self) -> &str {
        if self.count_setup {
            "mitm-frobenius-counted"
        } else {
            "mitm-frobenius"
        }
    }

    fn summands(&self) -> u32 {
        self.m
    }

    fn describe(&self, params: &Params) -> String {
        match LargePrimeConfig::from_params(params) {
            Ok(Some(config)) => format!(
                "a pair table folded by the signed Frobenius, one entry per orbit, plus a large-prime collector with up to {} residuals in a domain of {} × |F| signed points",
                config.large_primes, config.domain_factor
            ),
            _ => "a pair table folded by the signed Frobenius, one entry per orbit".into(),
        }
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            ("m", "summands, 2 or 3"),
            (
                "large_primes",
                "0 (default), 1 or 2 residual points outside the base a partial relation may carry; needs m=3",
            ),
            (
                "lp_domain_factor",
                "size of the large-prime domain D as a multiple of |F| in signed points (default 10); D is the next projected signed Frobenius orbits of the raw-x scan after the base's",
            ),
            (
                "lp_raw_x_cap",
                "inclusive raw-abscissa cap of the scan that materialises D (default 2^40)",
            ),
        ]
    }

    fn prepare(
        &mut self,
        ctx: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        params: &Params,
        ops: &mut GroupOps,
    ) -> Result<(), String> {
        if let Some(config) = LargePrimeConfig::from_params(params)? {
            if self.m != 3 {
                return Err("large_primes needs the three-summand probe loop (m=3)".into());
            }
            // The collector builds the classic table itself (its direct
            // path is the plain probe) and the large table beside it.
            self.large = Some(LargePrimeCollector::new(
                self.instance,
                ctx,
                fb,
                config,
                ops,
            )?);
            self.table = None;
            return Ok(());
        }
        let table = FrobeniusPairTable::build(self.instance, fb)
            .ok_or("this instance has no Koblitz structure to fold by")?;
        ops.merge(table.build_ops);
        self.table = Some(table);
        Ok(())
    }

    fn setup_native(&self) -> BTreeMap<String, u64> {
        let mut counts = BTreeMap::new();
        if let Some(lp) = &self.large {
            counts.extend(lp.setup_counters());
        }
        if self.count_setup {
            if let Some(table) = &self.table {
                counts.insert("canonicalisations".into(), table.build_canonicalisations);
                counts.insert("frobenius_maps".into(), table.build_frobenius_maps);
                counts.insert("lookups".into(), table.build_lookups);
                counts.insert("pair_table_entries".into(), table.entries);
                counts.insert("pair_table_representatives".into(), table.representatives);
            }
        }
        counts
    }

    fn decompose(
        &mut self,
        ctx: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: FastPoint,
    ) -> Option<Vec<usize>> {
        let table = self.table.as_ref()?;
        decompose_mitm_frobenius(ctx.group, fb, table, self.m, ops, counters, point)
    }

    fn decompose_relations(
        &mut self,
        ctx: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: FastPoint,
        a: u64,
        b: u64,
    ) -> Vec<Relation> {
        match self.large.as_mut() {
            Some(lp) => lp.decompose(ctx, fb, ops, counters, point, a, b),
            None => match self.decompose(ctx, fb, ops, counters, point) {
                Some(summands) => vec![Relation::Summands(summands)],
                None => Vec::new(),
            },
        }
    }

    fn relation_counters(&self) -> BTreeMap<String, u64> {
        self.large
            .as_ref()
            .map(|lp| lp.relation_counters())
            .unwrap_or_default()
    }
}

// ── The bridge to the solver plug point ────────────────────────────

/// Holds a [`SystemSolver`] and the record of what it cost, so an
/// algebraic oracle can be written without repeating the bookkeeping.
///
/// This is what makes the solver plug point reachable from a whole
/// pipeline: swap `buchberger-f2` for `sat-cdcl` and the change
/// propagates through the relation phase to the recovered key, with
/// every stage still counted.
pub struct SolverHarness {
    pub solver: Box<dyn SystemSolver>,
    pub params: Params,
    pub budget: Option<Duration>,
    last_shape: Option<SystemShape>,
    last_cost: Option<SolverCost>,
    /// Targets whose solver hit its budget.  Counted, never silently
    /// folded into "did not decompose": an oracle that did that would
    /// quietly lower its own hit rate and make the relation phase it
    /// feeds unmeasurable.
    pub budget_exceeded: u64,
}

impl SolverHarness {
    pub fn new(solver: Box<dyn SystemSolver>, params: Params, budget: Option<Duration>) -> Self {
        Self {
            solver,
            params,
            budget,
            last_shape: None,
            last_cost: None,
            budget_exceeded: 0,
        }
    }

    /// Solve one descended system, recording what it cost.
    pub fn solve(&mut self, system: &BooleanSystem) -> SolverVerdict {
        let shape = system.shape();
        self.last_shape = Some(shape.clone());
        if !self.solver.accepts(&shape) {
            self.last_cost = None;
            self.budget_exceeded += 1;
            return SolverVerdict::BudgetExceeded;
        }
        let (verdict, cost) = self.solver.solve(system, &self.params, self.budget);
        self.last_cost = Some(cost);
        if matches!(verdict, SolverVerdict::BudgetExceeded) {
            self.budget_exceeded += 1;
        }
        verdict
    }

    pub fn shape(&self) -> Option<&SystemShape> {
        self.last_shape.as_ref()
    }

    pub fn cost(&self) -> Option<&SolverCost> {
        self.last_cost.as_ref()
    }
}

// ── The algebraic oracle ───────────────────────────────────────────

/// Descend the Semaev polynomial over the base's abscissa subspace and
/// hand the boolean system to the chosen [`SystemSolver`].
///
/// This is the oracle that makes the solver plug point reach the `S`
/// column: swap `buchberger-f2` for `sat-cdcl` and the change
/// propagates through the relation phase to the recovered logarithm,
/// with every stage still counted.  It requires a factor base whose
/// abscissae form an `F_2`-subspace — `binary-subspace` and
/// `koblitz-orbit` both record theirs in `FactorBase::subspace_basis`.
///
/// ## What it charges
///
/// The descent and the solver are engine work, not group operations,
/// and are reported through [`DecompositionOracle::solver_totals`] in
/// the solver's own unit; the runner prices them into `S`.  Lifting a
/// solution back to signed base points costs real group additions,
/// which are charged to `ops` like any other oracle's.
///
/// ## The cap
///
/// The descent is symbolic — the summation polynomial expanded in
/// `F_{2^n}[v]/(v² − v)` term by term — so it builds no table and is
/// capped only by the monomial mask: `m·n' ≤ 64` boolean variables,
/// `n' ≤ 32` at two summands and `n' ≤ 21` at three.  What limits a
/// row above that is the engine, which says so through `accepts`.
pub struct DescentAlgebraicOracle<'i> {
    summands: u32,
    instance: &'i BinaryInstance,
    v_basis: Vec<u64>,
    /// Ordered τ-slots (lead L5): `None` is the plain descent (every
    /// summand over `V`); `Some(n)` tries, per target, the systems with
    /// summand 2 over `σ^{j}(V)` for `j = 0, …, n−1` (and summand 3 over
    /// `V`), stopping at the first that lifts.
    slot_rotations: Option<u32>,
    harness: SolverHarness,
    totals: SolverTotals,
    /// Systems the solver decided satisfiable none of whose solutions
    /// lifted to base points summing to the target.  Expected, not a
    /// defect: a summation polynomial vanishes over the algebraic
    /// closure, so a solution may name abscissae whose points live on
    /// the quadratic twist (`y ∈ F_{2^{2n}}`, never in the base), and
    /// the pair table never sees those.  Reported as
    /// `unliftable_systems` beside `lift_failures` so the hit rate can
    /// be read against what the solver actually found.
    pub unliftable: u64,
}

impl<'i> DescentAlgebraicOracle<'i> {
    pub fn new(
        summands: u32,
        instance: &'i BinaryInstance,
        solver: Box<dyn SystemSolver>,
        solver_params: Params,
        budget: Option<Duration>,
    ) -> Self {
        let totals = SolverTotals {
            solver: solver.name().to_string(),
            ..Default::default()
        };
        Self {
            summands,
            instance,
            v_basis: Vec::new(),
            slot_rotations: None,
            harness: SolverHarness::new(solver, solver_params, budget),
            totals,
            unliftable: 0,
        }
    }

    /// The descent with summand 2 over `σ^{j}(V)`: the τ-slot system.
    fn descend_slot(&self, x_r: u64, j: u32) -> SymbolicDescent {
        let gf = &self.instance.gf;
        let shifted: Vec<u64> = self
            .v_basis
            .iter()
            .map(|&e| (0..j).fold(e, |acc, _| gf.sqr(acc)))
            .collect();
        let mut bases = vec![self.v_basis.clone(), shifted];
        if self.summands == 3 {
            bases.push(self.v_basis.clone());
        }
        crate::cryptanalysis::pq_descent_symbolic::descend_with_bases(
            gf,
            self.instance.b,
            x_r,
            &bases,
        )
        .expect("prepare checked the dimension")
    }

    pub fn solver_name(&self) -> &str {
        self.harness.solver.name()
    }

    /// The descended system for a target abscissa; its `lift` maps a
    /// solution back to abscissae.  `prepare` must have run.
    fn descend(&self, x_r: u64) -> SymbolicDescent {
        descend(
            &self.instance.gf,
            self.instance.b,
            x_r,
            &self.v_basis,
            self.summands,
        )
        .expect("prepare checked the dimension")
    }
}

fn system_of(sys: &SymbolicDescent) -> BooleanSystem {
    BooleanSystem {
        equations: sys.equations.clone(),
        n_vars: sys.n_vars,
    }
}

impl<'a> DecompositionOracle<BinaryGroup<'a>> for DescentAlgebraicOracle<'_> {
    fn name(&self) -> &str {
        "descent-algebraic"
    }

    fn summands(&self) -> u32 {
        self.summands
    }

    fn describe(&self, _params: &Params) -> String {
        format!(
            "Weil-descend S_{} over the base's abscissa subspace and solve with {}",
            self.summands + 1,
            self.harness.solver.name()
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            ("m", "summands, 2 (descends S3) or 3 (descends S4)"),
            (
                "slots",
                "ordered tau-slots (Koblitz): `all` tries summand 2 over every Frobenius image of V per target; absent for the plain descent",
            ),
        ]
    }

    fn prepare(
        &mut self,
        _ctx: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<(), String> {
        self.slot_rotations = match params.get("slots") {
            None | Some("0") | Some("none") => None,
            Some("all") => {
                if self.instance.koblitz.is_none() {
                    return Err(
                        "slots need a Koblitz curve (the Frobenius must be an endomorphism)".into(),
                    );
                }
                Some(self.instance.n)
            }
            Some(other) => return Err(format!("unknown slots value `{other}`; use all")),
        };
        if !(2..=3).contains(&self.summands) {
            return Err(format!(
                "descent-algebraic takes m = 2 or 3, got {}",
                self.summands
            ));
        }
        let basis = fb.subspace_basis.as_ref().ok_or(
            "descent-algebraic needs a factor base whose abscissae form an F_2-subspace; \
             use binary-subspace or koblitz-orbit",
        )?;
        let cap = max_n_prime(self.summands);
        if basis.len() as u32 > cap {
            return Err(format!(
                "the descent's monomial mask holds n' ≤ {cap} at m = {}; this base has n' = {}",
                self.summands,
                basis.len()
            ));
        }
        self.v_basis = basis.clone();
        // Every system on one base has the same shape, so descend one
        // abscissa now and let the engine decline the shape here, once,
        // with a reason a sweep can report — rather than target by
        // target, which would run the whole relation phase for nothing.
        let shape = system_of(&self.descend(self.instance.generator.x)).shape();
        if !self.harness.solver.accepts(&shape) {
            return Err(format!(
                "solver `{}` declines a system of {} unknowns and {} equations; \
                 use a smaller base or another engine",
                self.harness.solver.name(),
                shape.n_vars,
                shape.n_equations
            ));
        }
        Ok(())
    }

    fn decompose(
        &mut self,
        ctx: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: FastPoint,
    ) -> Option<Vec<usize>> {
        if self.v_basis.is_empty() {
            return None;
        }
        // Descend, solve, and lift every solution until one sums to R;
        // with slots, one system per rotation of the second summand.
        let rotations = self.slot_rotations.unwrap_or(1);
        for j in 0..rotations {
            let sys = if self.slot_rotations.is_some() {
                self.descend_slot(point.x, j)
            } else {
                self.descend(point.x)
            };
            let system = system_of(&sys);
            let shape = system.shape();
            let verdict = self.harness.solve(&system);
            let exceeded = matches!(verdict, SolverVerdict::BudgetExceeded);
            self.totals.absorb(&shape, self.harness.cost(), exceeded);
            *self.totals.extra.entry("slot_systems".into()).or_insert(0) += 1;
            let SolverVerdict::Solved(solutions) = verdict else {
                continue;
            };
            for v in solutions {
                let xs = sys.lift(v);
                if let Some(indices) = crate::cryptanalysis::ic_boundary::lift_abscissae(
                    ctx.group, fb, ops, &xs, point,
                ) {
                    *self.totals.extra.entry("slot_hits".into()).or_insert(0) += 1;
                    *self
                        .totals
                        .extra
                        .entry(format!("slot_hit_rotation_{j}"))
                        .or_insert(0) += 1;
                    return Some(indices);
                }
                counters.lift_failures += 1;
            }
            self.unliftable += 1;
            counters.unliftable_systems += 1;
        }
        None
    }

    fn last_system(&self) -> Option<SystemShape> {
        self.totals.shape.clone()
    }

    fn solver_totals(&self) -> Option<SolverTotals> {
        Some(self.totals.clone())
    }
}

// ── The symmetrised Koblitz base and oracle ────────────────────────

/// Parse a `divisor` parameter: indices into
/// [`crate::cryptanalysis::koblitz_index_calculus::all_factors_of_x_n_minus_1`],
/// separated by `;` (a `,` would split the plug-in spec itself; a `,` is
/// still accepted here for a [`Params`] built in code, but
/// [`Params::parse_spec`] never lets one through from a command line).
fn divisor_indices(params: &Params) -> Result<Vec<usize>, String> {
    let raw = params
        .get("divisor")
        .ok_or("missing parameter `divisor`; try 0;1 (index 0 is x + 1)")?;
    let idx: Result<Vec<usize>, String> = raw
        .split([',', ';'])
        .filter(|s| !s.is_empty())
        .map(|s| {
            s.trim()
                .parse::<usize>()
                .map_err(|_| format!("divisor index `{s}` is not a number"))
        })
        .collect();
    let idx = idx?;
    if idx.is_empty() {
        return Err("parameter `divisor` selected no factors".into());
    }
    Ok(idx)
}

/// `F_u = { P : u(P) ∈ V }` for the frame `u = 1/(x + 1)` and a
/// Frobenius-invariant subspace `V ∋ 1`, each signed Frobenius orbit
/// folded onto one column.
///
/// This is the base of `research/notes/index-calculus/RESEARCH_EXOTIC_COORDINATES.md`
/// §8: the rational 2-torsion point `T = (0, 1)` acts on the `x`-line
/// as `x ↦ 1/x`, the frame turns that into `u ↦ u + 1`, and a subspace
/// containing `1` is closed under it for free.  `T` itself (`u = 1`) is
/// a base point.  The abscissae of `F_u` are *not* an `F_2`-subspace,
/// so the `descent-algebraic` oracle cannot run on this base; the
/// `symmetrised` oracle is the one that can, and `mitm-frobenius` runs
/// on it like on any other.  The same `divisor` handed to
/// `koblitz-orbit` gives the `x`-frame base `F_x` for the same `V`,
/// which is the matched control.
pub struct KoblitzSymmetrisedBase<'i> {
    pub instance: &'i BinaryInstance,
}

impl<'a> FactorBaseBuilder<BinaryGroup<'a>> for KoblitzSymmetrisedBase<'_> {
    fn name(&self) -> &str {
        "koblitz-symmetrised"
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "F_u = {{P : 1/(x + 1) ∈ V}} for the Frobenius-invariant V ∋ 1 of divisor {}, {}",
            params.get("divisor").unwrap_or("?"),
            if params.flag("no_fold") {
                "one column per abscissa (the unfolded control)"
            } else {
                "each signed Frobenius orbit folded onto one column"
            }
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            (
                "divisor",
                "`;`-separated indices into the factorisation of xⁿ − 1, e.g. `divisor=0;1` (a `,` \
                 would end the parameter on the command line); index 0 is x + 1, which must be \
                 included so that 1 ∈ V; the same indices on koblitz-orbit give the x-frame base \
                 for the same V",
            ),
            (
                "no_fold",
                "1 to give every abscissa its own column: the control that shows what the fold buys",
            ),
            (
                "projected",
                "1 to fold by the signed Frobenius orbits of the cofactor projections [h]P instead of \
                 of the raw points: orbits whose projections coincide or are Frobenius translates \
                 share one column (on a base closed under +T, P and P + T do), which is the column \
                 count the relation matrix actually has to determine",
            ),
        ]
    }

    fn build(
        &self,
        _ctx: &InstanceCtx<BinaryGroup<'a>>,
        params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<FactorBase<FastPoint>, String> {
        let kc = self
            .instance
            .koblitz
            .as_ref()
            .ok_or("this instance is not a Koblitz curve, so it has no 1/(x + 1) frame")?;
        let idx = divisor_indices(params)?;
        let fb_u = build_symmetrised_factor_base(kc, &idx).ok_or_else(|| {
            format!(
                "no symmetrised base for divisor {idx:?}: the subspace must contain 1 \
                 (include index 0, the factor x + 1) and be a proper invariant subspace"
            )
        })?;
        let view = frobenius_view_of_symmetrised(kc, &fb_u)
            .ok_or("the u-frame base has no Frobenius-orbit view; this is a bug")?;
        let fold = if params.flag("no_fold") {
            ColumnFold::Abscissa
        } else if params.flag("projected") {
            ColumnFold::ProjectedSignedFrobeniusOrbit
        } else {
            ColumnFold::SignedFrobeniusOrbit
        };
        let mut fb = koblitz_factor_base(
            self.instance,
            &view,
            fold,
            format!(
                "u-frame invariant subspace, divisor {idx:?}, dimension {}, {} points including T",
                fb_u.ell,
                fb_u.points.len()
            ),
        )
        .ok_or_else(|| "the u-frame subspace produced no usable factor base".to_string())?;
        // The base is the subspace V of the u-frame, so a frozen copy
        // carries its dimension and basis (single words, as the
        // x-frame bases do) rather than only a description.
        fb.dimension = Some(fb_u.ell as u32);
        fb.subspace_basis = Some(
            fb_u.v_basis
                .iter()
                .map(|v| u64::try_from(v.to_biguint()).unwrap_or(0))
                .collect(),
        );
        Ok(fb)
    }
}

/// Decompose over `F_u` by solving the symmetrised summation polynomial
/// in the invariants `w = u² + u`, `s = Σu` with the repository's
/// boolean matrix-F4 and splitting (`koblitz_symmetrised`).
///
/// The system has `m(ℓ − 1) + 1` unknowns of boolean degree 2 at
/// `m = 2` and 4 at `m = 3`, against `m·ℓ` unknowns of degree 2 and 6
/// for the direct `x`-frame descent on the same `V`.  A root lifts to
/// points summing to `R` or to `R + T`; both are relations over `F_u`,
/// and the second carries `T` as an extra summand.
///
/// ## What it charges
///
/// The solve is engine work, reported through `solver_totals` in the
/// same unit and by the same rule as the `inherited-f4` adapter (word
/// XORs counted, priced by the runner); lifting a root back to signed
/// base points goes through the framework's own `lift_abscissae`, so
/// its group additions are charged to `ops` exactly as the
/// `descent-algebraic` oracle's are.
///
/// ## Parameters
///
/// `m` (2 or 3), `divisor` (the base's; copied from the factor base
/// when omitted on the command line), `engine`
/// (`inherited-f4` default, `matrix-f4`, `matrix-f5`), `max_degree`
/// (Macaulay cap, default 3: the cap is absolute, so the degree-4
/// `m = 3` system still builds its own degree at the root) and
/// `node_budget` (algebraic reduction calls before a solve gives up,
/// default 4096; splitting decisions are counted separately).
pub struct SymmetrisedOracle<'i> {
    summands: u32,
    instance: &'i BinaryInstance,
    fb_u: Option<SymmetrisedFactorBase>,
    st: Option<FieldStructure>,
    engine: SolverEngine,
    node_budget: usize,
    shape: Option<SystemShape>,
    totals: SolverTotals,
    /// Solves whose root lifted in the `u`-frame but whose abscissae
    /// did not re-lift through the framework's sign resolution; must
    /// stay zero, and is reported so that it can be seen to.
    pub unliftable: u64,
    /// Targets with `u(R) ∈ {0, ∞}` or an oversize layout, which the
    /// system cannot be built for.
    pub unbuildable: u64,
}

impl<'i> SymmetrisedOracle<'i> {
    pub fn new(summands: u32, instance: &'i BinaryInstance) -> Self {
        Self {
            summands,
            instance,
            fb_u: None,
            st: None,
            engine: SolverEngine::InheritedF4 { max_degree: 3 },
            node_budget: 4096,
            shape: None,
            totals: SolverTotals::default(),
            unliftable: 0,
            unbuildable: 0,
        }
    }
}

impl<'a> DecompositionOracle<BinaryGroup<'a>> for SymmetrisedOracle<'_> {
    fn name(&self) -> &str {
        "symmetrised"
    }

    fn summands(&self) -> u32 {
        self.summands
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "solve the symmetrised S_{} in w = u² + u, s = Σu over F_u with {} (Macaulay cap {}, {} reduction calls)",
            self.summands + 1,
            params.get("engine").unwrap_or("inherited-f4"),
            params.get("max_degree").unwrap_or("3"),
            params.get("node_budget").unwrap_or("4096"),
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            ("m", "summands, 2 (symmetrised S3) or 3 (symmetrised S4)"),
            (
                "divisor",
                "the koblitz-symmetrised base's divisor; copied from the base when omitted",
            ),
            ("engine", "inherited-f4 (default), matrix-f4 or matrix-f5"),
            (
                "max_degree",
                "highest Macaulay degree built before splitting (default 3)",
            ),
            (
                "node_budget",
                "algebraic reduction calls before a solve gives up (default 4096)",
            ),
        ]
    }

    fn prepare(
        &mut self,
        _ctx: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<(), String> {
        if !(2..=3).contains(&self.summands) {
            return Err(format!(
                "symmetrised takes m = 2 or 3, got {}",
                self.summands
            ));
        }
        let kc = self
            .instance
            .koblitz
            .as_ref()
            .ok_or("symmetrised needs a Koblitz instance")?;
        let idx = divisor_indices(params)?;
        let fb_u = build_symmetrised_factor_base(kc, &idx)
            .ok_or_else(|| format!("no symmetrised base for divisor {idx:?}"))?;
        // The base the runner built must be this F_u, point for point:
        // a mismatched divisor would make every relation a lie.
        if fb_u.points.len() != fb.points.len() {
            return Err(format!(
                "symmetrised needs the koblitz-symmetrised base with the same divisor: \
                 F_u has {} points, the base has {}",
                fb_u.points.len(),
                fb.points.len()
            ));
        }
        for p in &fb_u.points {
            let key = self.instance.fast.lift(p).pack();
            if fb.index_of_key(key).is_none() {
                return Err(
                    "symmetrised needs the koblitz-symmetrised base with the same divisor: \
                     a point of F_u is not in the base"
                        .into(),
                );
            }
        }
        let max_degree = params.u64_or("max_degree", 3)? as u32;
        self.engine = match params.get("engine").unwrap_or("inherited-f4") {
            "inherited-f4" => SolverEngine::InheritedF4 { max_degree },
            "matrix-f4" => SolverEngine::MatrixF4 { max_degree },
            "matrix-f5" => SolverEngine::MatrixF5 { max_degree },
            other => {
                return Err(format!(
                    "unknown engine `{other}`; try inherited-f4, matrix-f4 or matrix-f5"
                ))
            }
        };
        self.node_budget = params.u64_or("node_budget", 4096)? as usize;
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        // Every system on one base has the same shape; build one for the
        // generator so the report's degree columns come from a real
        // system, and so an oversize layout is refused here, once.
        let g = self.instance.fast.lower(self.instance.generator);
        let sys = build_symmetrised_system(kc, &fb_u, &g, self.summands as usize, &st).ok_or_else(
            || {
                format!(
                    "the symmetrised system for dimension {} at m = {} exceeds the 64-variable \
                     layout or the generator has u ∈ {{0, ∞}}",
                    fb_u.ell, self.summands
                )
            },
        )?;
        self.shape = Some(
            BooleanSystem {
                equations: sys.equations,
                n_vars: sys.n_vars,
            }
            .shape(),
        );
        self.totals = SolverTotals {
            solver: format!(
                "symmetrised {} (cap {}, budget {})",
                match self.engine {
                    SolverEngine::InheritedF4 { .. } => "inherited-f4",
                    SolverEngine::MatrixF4 { .. } => "matrix-f4",
                    SolverEngine::MatrixF5 { .. } => "matrix-f5",
                    _ => "engine",
                },
                max_degree,
                self.node_budget
            ),
            ..Default::default()
        };
        self.fb_u = Some(fb_u);
        self.st = Some(st);
        Ok(())
    }

    fn decompose(
        &mut self,
        ctx: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: FastPoint,
    ) -> Option<Vec<usize>> {
        let (Some(fb_u), Some(st), Some(shape)) = (&self.fb_u, &self.st, &self.shape) else {
            return None;
        };
        let kc = self.instance.koblitz.as_ref()?;
        let target = self.instance.fast.lower(point);
        let before = f4_word_ops_thread();
        let started = Instant::now();
        let outcome = symmetrised_groebner_decompose(
            kc,
            fb_u,
            st,
            &target,
            self.summands as usize,
            self.engine,
            self.node_budget,
        );
        let word_xors = f4_word_ops_thread().wrapping_sub(before);
        let wall_ns = started.elapsed().as_nanos() as u64;
        let Some(out) = outcome else {
            self.unbuildable += 1;
            return None;
        };
        let mut extra = BTreeMap::new();
        extra.insert("splits".to_string(), out.effort);
        extra.insert("oversize".to_string(), out.oversize as u64);
        extra.insert("used_t".to_string(), u64::from(out.used_t));
        let cost = SolverCost {
            ops: word_xors,
            // The same unit string as the `inherited-f4` adapter, so the
            // runner prices both arms by the same rule.
            op_unit: "word XORs (elimination, specialisation and linear elimination only)".into(),
            wall_ns,
            peak_bytes: 0,
            degree_reached: Some(out.built_degree),
            solving_degree: None,
            timed_out: !out.complete,
            extra,
        };
        self.totals.absorb(shape, Some(&cost), !out.complete);
        let indices = out.relation?;
        // Re-lift through the framework so the sign resolution is charged
        // exactly as it is for the x-frame oracle.  `T` (x = 0) is a base
        // point like any other, so a relation through `R + T` lifts too.
        let xs: Vec<u64> = indices
            .iter()
            .map(|&i| self.instance.fast.lift(&fb_u.points[i]).x)
            .collect();
        match crate::cryptanalysis::ic_boundary::lift_abscissae(ctx.group, fb, ops, &xs, point) {
            Some(found) => Some(found),
            None => {
                counters.lift_failures += 1;
                self.unliftable += 1;
                None
            }
        }
    }

    fn last_system(&self) -> Option<SystemShape> {
        self.shape.clone()
    }

    fn solver_totals(&self) -> Option<SolverTotals> {
        Some(self.totals.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ic_boundary::koblitz_instance;

    /// The command-line form of a multi-factor divisor, end to end from
    /// the plug-in string to the indices both Koblitz bases read.
    #[test]
    fn divisor_indices_parse_the_semicolon_form_from_a_plugin_string() {
        for spec in [
            "koblitz-orbit:divisor=1;2",
            "koblitz-symmetrised:divisor=1;2",
        ] {
            let (_, params) = Params::parse_spec(spec).unwrap();
            assert_eq!(divisor_indices(&params).unwrap(), vec![1, 2], "{spec}");
        }
        let (_, params) = Params::parse_spec("koblitz-orbit:divisor=0,no_fold=1").unwrap();
        assert_eq!(divisor_indices(&params).unwrap(), vec![0]);
    }

    /// The comma form never reaches `divisor_indices` from a command
    /// line: the spec parser refuses it first, and says why.
    #[test]
    fn the_comma_form_is_refused_before_the_plugin_sees_it() {
        let err = Params::parse_spec("koblitz-orbit:divisor=1,2").unwrap_err();
        assert!(err.contains("`divisor=1;2`"), "{err}");
    }

    #[test]
    fn divisor_indices_reject_a_missing_empty_or_non_numeric_value() {
        let mut p = Params::default();
        assert!(divisor_indices(&p)
            .unwrap_err()
            .contains("missing parameter `divisor`"));
        p.set("divisor", ";");
        assert!(divisor_indices(&p)
            .unwrap_err()
            .contains("selected no factors"));
        p.set("divisor", "1;x");
        assert!(divisor_indices(&p)
            .unwrap_err()
            .contains("`x` is not a number"));
    }

    #[test]
    fn trace_zero_factor_base_is_an_executable_policy() {
        let inst = koblitz_instance(1, 9).expect("registered n9 Koblitz instance");
        let group = BinaryGroup(&inst.fast);
        let ctx = InstanceCtx {
            group: &group,
            generator: inst.generator,
            target: inst.generator,
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: Some(inst.n),
        };
        let mut ops = GroupOps::default();
        let (_, trace_zero) = Params::parse_spec("koblitz-trace-zero:divisor=1").unwrap();
        let base = KoblitzTraceZeroBase { instance: &inst }
            .build(&ctx, &trace_zero, &mut ops)
            .expect("a divisor omitting x + 1 is trace-zero");
        assert!(base.columns > 0);

        let (_, trace_one) = Params::parse_spec("koblitz-trace-zero:divisor=0").unwrap();
        let error =
            match (KoblitzTraceZeroBase { instance: &inst }).build(&ctx, &trace_one, &mut ops) {
                Err(error) => error,
                Ok(_) => panic!("a divisor containing x + 1 must be rejected"),
            };
        assert!(error.contains("not trace-zero"), "{error}");
    }

    /// Freeze the entire producer support, not merely its size or first
    /// representative. This is the preregistered gate before fresh Q.
    #[test]
    fn compact_orbit_scan_matches_archived_n37_source42() {
        let inst = koblitz_instance(0, 37).expect("registered n37 Koblitz instance");
        let g = BinaryGroup(&inst.fast);
        let ctx = InstanceCtx {
            group: &g,
            generator: inst.generator,
            target: inst.generator,
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: Some(inst.n),
        };
        let (_, params) =
            Params::parse_spec("compact-orbit-scan:columns=42,raw_x_cap=1000000").unwrap();
        let mut ops = GroupOps::default();
        let fb = CompactOrbitScanBase { instance: &inst }
            .build(&ctx, &params, &mut ops)
            .unwrap();
        let source = std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/research/notes/ecc2k130/n37_four_policy_support_20261003/SOURCE42.jsonl"
        ))
        .unwrap();
        let frozen: serde_json::Value =
            serde_json::from_str(source.lines().next().unwrap()).unwrap();
        let coordinates: Vec<[u64; 2]> =
            serde_json::from_value(frozen["factor_base_point_coordinates"].clone()).unwrap();
        let labels: Vec<[u64; 2]> =
            serde_json::from_value(frozen["factor_base_point_labels"].clone()).unwrap();
        let representatives: Vec<[u64; 2]> =
            serde_json::from_value(frozen["factor_base_representatives"].clone()).unwrap();
        assert_eq!(fb.columns, 42);
        assert_eq!(fb.points.len(), 3108);
        assert_eq!(coordinates.len(), fb.points.len());
        assert_eq!(labels.len(), fb.points.len());
        assert_eq!(representatives.len(), fb.columns);
        for (i, p) in fb.points.iter().enumerate() {
            assert_eq!([p.x, p.y], coordinates[i], "point {i}");
            assert_eq!([fb.col_of[i] as u64, fb.coef_of[i]], labels[i], "label {i}");
        }
        for (column, rep) in representatives.iter().enumerate() {
            let p = fb.points[column * 2 * inst.n as usize];
            assert_eq!([p.x, p.y], *rep, "representative {column}");
        }
        let mut base_hash = blake3::Hasher::new();
        base_hash.update(b"compact-orbit-constructed-base-v1");
        for (i, p) in fb.points.iter().enumerate() {
            for word in [p.x, p.y, fb.col_of[i] as u64, fb.coef_of[i]] {
                base_hash.update(&word.to_le_bytes());
            }
        }
        assert_eq!(
            base_hash.finalize().to_hex().as_str(),
            frozen["base_hash"].as_str().unwrap()
        );
        assert_eq!(
            frozen["base_hash"].as_str().unwrap(),
            "8423b135df3515b0e284126d9eeb95f71e3901ae12908b3cad8eb03d1f31d4bc"
        );
        assert!(ops.scalar_mults > 0 && fb.cost.get("inversions") > 0);
        assert!(fb.cost.get("min_key_field_squares_uncharged") > 0);
    }
}
