//! # `ecbench`: one harness for every ECDLP method.
//!
//! Pollard rho in its variants, baby-step giant-step, the kangaroo and
//! index calculus run here behind one interface, on the same single-target
//! workloads, charged in the same unit, under the same isolation, with
//! the same records.  `docs/ecbench/README.md` is the standard this
//! implements; `src/bin/ecbench.rs` is the command line.
//!
//! | module | what it holds |
//! |---|---|
//! | [`canonical`] | canonical JSON, SHA-256 identities, derived seeds |
//! | [`workload`] | curve constructions → ICV1-named single-target workloads |
//! | [`methods`] | the method registry and the one `solve` every method goes through |
//! | [`generic`] | counted BSGS (three forms) and the kangaroo on any [`CountedGroup`] |
//! | [`spec`] | the experiment spec and its deterministic expansion |
//! | [`host`] | the host capsule and its environment class |
//! | [`isolation`] | CPU reservation, eviction, pinning, NUMA binding, kernel counters |
//! | [`record`] | the record, the child protocol, isolation grading |
//! | [`runner`] | sessions: interleaved measured children, sealed append-only records |
//! | [`signals`] | interruption: the host is restored and no child outlives the runner |
//! | [`compare`] | paired ratios with bootstrap intervals; wall time gated by level |
//! | [`audit`] | re-derive a session from its files; replay runs exactly |
//! | [`db`] | SQL that loads sessions into the schema in `docs/ecbench/schema.sql` |
//! | [`isolab`] | an `isolab.job/v1` that runs a spec on an independent lab worker |
//!
//! [`CountedGroup`]: crate::cryptanalysis::ic_boundary::CountedGroup

pub mod audit;
pub mod canonical;
pub mod compare;
pub mod db;
pub mod generic;
pub mod host;
pub mod isolab;
pub mod isolation;
pub mod methods;
pub mod record;
pub mod runner;
pub mod signals;
pub mod spec;
pub mod stats;
pub mod workload;
