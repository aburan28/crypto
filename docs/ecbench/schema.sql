-- ecbench database schema, version 1 (SQLite 3.24+; PostgreSQL with minor
-- type changes).
--
-- This database is an INDEX, rebuilt from committed session directories by
--   ecbench db sql <session dirs...> | sqlite3 ecbench.db
-- It is never the record of a measurement: the session directory (sealed
-- records, host capsule, hashes) is, and `ecbench verify` audits it.  Every
-- statement `ecbench db sql` emits is idempotent, so loading a session twice
-- changes nothing.
--
-- Conventions
--   * Curves are keyed by ICV1 slug (AGENTS.md §11, docs/curves/ICV1.md).
--     EC1 aliases and curve UIDs are joined per representation, never inferred.
--   * 64-bit unsigned values that may exceed 2^63 (seeds, planted scalars,
--     group orders) are TEXT in decimal; counts that cannot are INTEGER.
--   * S = total group-addition equivalents / sqrt(r), over the whole method,
--     cold, one target.  A row whose `lower_bound` is 1 left work unpriced.
--   * Only runs with status = 'verified' and warmup = 0 are results.

PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS schema_info (
  key   TEXT PRIMARY KEY,
  value TEXT NOT NULL
);

-- ── What is measured on ──────────────────────────────────────────────

CREATE TABLE IF NOT EXISTS curves (
  slug                    TEXT PRIMARY KEY,           -- icv1-...
  icv1                    TEXT NOT NULL UNIQUE,       -- full ICV1 string
  family                  TEXT NOT NULL CHECK (family IN ('prime', 'koblitz', 'binary')),
  field_degree            INTEGER,                    -- m for GF(2^m)
  field_bits              INTEGER,                    -- bits of p for GF(p)
  group_order             TEXT NOT NULL,
  r                       TEXT NOT NULL,              -- prime subgroup order
  cofactor                TEXT NOT NULL,
  log2_r                  REAL NOT NULL,
  automorphisms_available INTEGER NOT NULL,           -- 2, or 2m on a Koblitz curve
  floor_s                 REAL NOT NULL,              -- sqrt(pi / 2A)
  registered              INTEGER NOT NULL CHECK (registered IN (0, 1))
);

-- One generator of one subgroup: the exact representation measured.
CREATE TABLE IF NOT EXISTS curve_representations (
  slug        TEXT NOT NULL REFERENCES curves(slug),
  generator_x TEXT NOT NULL,
  generator_y TEXT NOT NULL,
  ec1         TEXT,                                   -- registry EC1 alias, when it has this one
  curve_uid   TEXT,
  PRIMARY KEY (slug, generator_x, generator_y)
);

-- The constructor calls that produced a curve (a search, the roster, explicit).
CREATE TABLE IF NOT EXISTS curve_constructions (
  slug         TEXT NOT NULL REFERENCES curves(slug),
  construction TEXT NOT NULL,
  PRIMARY KEY (slug, construction)
);

-- One single-target workload: curve, generator, target point.
CREATE TABLE IF NOT EXISTS workloads (
  workload_id     TEXT PRIMARY KEY,                   -- W + 12 hex
  workload_sha256 TEXT NOT NULL UNIQUE,
  curve_slug      TEXT NOT NULL REFERENCES curves(slug),
  generator_x     TEXT NOT NULL,
  generator_y     TEXT NOT NULL,
  target_x        TEXT NOT NULL,
  target_y        TEXT NOT NULL,
  target_law      TEXT NOT NULL,
  target_seed     TEXT NOT NULL,
  target_index    INTEGER NOT NULL,
  planted         TEXT                                -- the known answer; NULL for a public target
);

-- ── What is measured ─────────────────────────────────────────────────

-- An algorithm configuration: a registry method with its parameters,
-- defaults written out.  Rho, BSGS, kangaroo and index calculus share it.
CREATE TABLE IF NOT EXISTS algorithms (
  method_id     TEXT PRIMARY KEY,                     -- ECM1h + 12 hex
  method_sha256 TEXT NOT NULL UNIQUE,
  method        TEXT NOT NULL,                        -- e.g. rho.negation, ic.pipeline
  family        TEXT NOT NULL CHECK (family IN ('rho', 'bsgs', 'kangaroo', 'ic')),
  params_json   TEXT NOT NULL,
  entry         TEXT NOT NULL                         -- the code path it runs
);

-- A factor base, identified by its points.
CREATE TABLE IF NOT EXISTS factor_bases (
  fb_id         TEXT PRIMARY KEY,                     -- FB1h + 12 hex
  fb_sha256     TEXT NOT NULL UNIQUE,
  curve_slug    TEXT NOT NULL REFERENCES curves(slug),
  family        TEXT NOT NULL,                        -- prime-abscissa, koblitz-orbit, ...
  params_json   TEXT NOT NULL,
  description   TEXT,
  signed_points INTEGER NOT NULL,
  abscissae     INTEGER NOT NULL,
  columns       INTEGER NOT NULL,
  dimension     INTEGER,
  points_sha256 TEXT NOT NULL                         -- sha256 of the sorted point keys
);

-- The points themselves, when a factor-base dump (`ecbench fb`) is loaded.
CREATE TABLE IF NOT EXISTS factor_base_points (
  fb_id  TEXT NOT NULL REFERENCES factor_bases(fb_id),
  idx    INTEGER NOT NULL,
  x      TEXT NOT NULL,
  y      TEXT NOT NULL,
  col    INTEGER NOT NULL,                            -- relation-matrix column
  coef   TEXT NOT NULL,                               -- coefficient of the point in its column
  PRIMARY KEY (fb_id, idx)
);

-- ── Where and how ────────────────────────────────────────────────────

CREATE TABLE IF NOT EXISTS hosts (
  env_class_id     TEXT PRIMARY KEY,                  -- ECBENV2h (or v1 ECBENV1h) + 12 hex of the stable facts
  env_class_sha256 TEXT NOT NULL UNIQUE,
  os               TEXT NOT NULL,
  arch             TEXT NOT NULL,
  kernel_release   TEXT,
  cpu_model        TEXT,
  features         TEXT NOT NULL,                     -- space-separated
  logical_cpus     INTEGER NOT NULL,
  physical_cores   INTEGER,
  numa_nodes       INTEGER NOT NULL,
  virtual          INTEGER,                           -- 1, 0, or NULL when unknown
  rustc            TEXT,
  capsule_json     TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS sessions (
  session_id       TEXT PRIMARY KEY,                  -- ECBS1h + 12 hex
  label            TEXT NOT NULL,
  spec_id          TEXT NOT NULL,                     -- ECS1h + 12 hex
  spec_sha256      TEXT NOT NULL,
  env_class_id     TEXT NOT NULL REFERENCES hosts(env_class_id),
  git_commit       TEXT,
  git_dirty        INTEGER,
  binary_sha256    TEXT,
  cpu_request      TEXT NOT NULL,
  run_cpu          INTEGER,
  reserved_cpus    TEXT,
  numa_node        INTEGER,
  preflight_quiet  INTEGER,
  started_unix_ms  INTEGER NOT NULL,
  finished_unix_ms INTEGER,
  status           TEXT NOT NULL,                     -- complete | interrupted | running
  records_sha256   TEXT,
  path             TEXT NOT NULL                      -- where the session directory was loaded from
);

CREATE TABLE IF NOT EXISTS arms (
  session_id TEXT NOT NULL REFERENCES sessions(session_id),
  arm        TEXT NOT NULL,
  role       TEXT NOT NULL CHECK (role IN ('reference', 'baseline', 'candidate', 'control')),
  method_id  TEXT NOT NULL REFERENCES algorithms(method_id),
  PRIMARY KEY (session_id, arm)
);

-- ── What happened ────────────────────────────────────────────────────

CREATE TABLE IF NOT EXISTS runs (
  record_id            TEXT PRIMARY KEY,              -- ECR1h + 12 hex, seals the record
  session_id           TEXT NOT NULL REFERENCES sessions(session_id),
  run_id               TEXT NOT NULL,                 -- <method_id><workload_id>R<n>
  seq                  INTEGER NOT NULL,
  round                INTEGER NOT NULL,
  warmup               INTEGER NOT NULL CHECK (warmup IN (0, 1)),
  arm                  TEXT NOT NULL,
  method_id            TEXT NOT NULL REFERENCES algorithms(method_id),
  workload_id          TEXT NOT NULL REFERENCES workloads(workload_id),
  fb_id                TEXT REFERENCES factor_bases(fb_id),
  algorithm_seed       TEXT NOT NULL,
  status               TEXT NOT NULL CHECK (status IN ('verified', 'wrong_answer', 'exhausted', 'error', 'timeout', 'crashed')),
  recovered            TEXT,
  total_gae            REAL,
  s                    REAL,
  ratio_to_floor       REAL,
  automorphisms_used   INTEGER,
  lower_bound          INTEGER NOT NULL,
  unpriced             TEXT NOT NULL,                 -- space-separated counter names
  deterministic        INTEGER NOT NULL,
  solve_wall_ns        INTEGER,
  process_wall_ns      INTEGER NOT NULL,
  user_ns              INTEGER NOT NULL,
  sys_ns               INTEGER NOT NULL,
  max_rss_kib          INTEGER NOT NULL,
  involuntary_switches INTEGER NOT NULL,
  run_delay_ns         INTEGER,
  isolation_level      TEXT NOT NULL CHECK (isolation_level IN ('L0', 'L1', 'L2', 'L3')),
  foreign_busy_ticks   INTEGER,
  steal_ticks          INTEGER,
  error                TEXT,
  UNIQUE (session_id, seq),
  UNIQUE (session_id, run_id)
);

CREATE TABLE IF NOT EXISTS phases (
  record_id     TEXT NOT NULL REFERENCES runs(record_id),
  ord           INTEGER NOT NULL,
  phase         TEXT NOT NULL,
  adds          INTEGER NOT NULL,
  doubles       INTEGER NOT NULL,
  scalar_mults  INTEGER NOT NULL,
  gae           REAL NOT NULL,
  wall_ns       INTEGER,
  native_json   TEXT NOT NULL,
  PRIMARY KEY (record_id, ord)
);

CREATE TABLE IF NOT EXISTS run_counters (
  record_id TEXT NOT NULL REFERENCES runs(record_id),
  name      TEXT NOT NULL,
  value     INTEGER NOT NULL,
  PRIMARY KEY (record_id, name)
);

CREATE TABLE IF NOT EXISTS isolation_blockers (
  record_id TEXT NOT NULL REFERENCES runs(record_id),
  ord       INTEGER NOT NULL,
  blocker   TEXT NOT NULL,                            -- "L2: 3 ticks of hypervisor steal ..."
  PRIMARY KEY (record_id, ord)
);

CREATE TABLE IF NOT EXISTS comparisons (
  comparison_id     TEXT PRIMARY KEY,                 -- ECC1h + 12 hex
  session_a         TEXT NOT NULL,
  arm_a             TEXT NOT NULL,
  session_b         TEXT NOT NULL,
  arm_b             TEXT NOT NULL,
  ops_status        TEXT NOT NULL,                    -- ok | incomplete | empty
  ops_ratio         REAL,
  ops_ci_low        REAL,
  ops_ci_high       REAL,
  ops_bounded       INTEGER NOT NULL,
  wall_status       TEXT NOT NULL,                    -- admitted | descriptive | refused
  wall_ratio        REAL,
  wall_ci_low       REAL,
  wall_ci_high      REAL,
  verdict           TEXT NOT NULL,
  json              TEXT NOT NULL
);

-- vs_rho claims (ecbench claim build): one IC run against one strong-rho
-- run on one public target, keyed by the repository's IC1 run id.  A
-- claim rebuilt with an independent receipt replaces its earlier form;
-- the identity trigger keeps the run id on one candidate, workload and
-- session.
CREATE TABLE IF NOT EXISTS claims (
  run_id                     TEXT PRIMARY KEY,                -- <candidate_id>W<workload_id>R<n>
  candidate_id               TEXT NOT NULL,                   -- IC1 label
  candidate_manifest_sha256  TEXT NOT NULL,
  claim_workload_id          TEXT NOT NULL,                   -- IC1 workload id (12 hex), not ecbench's W…
  workload_manifest_sha256   TEXT NOT NULL,
  session_id                 TEXT NOT NULL,
  ic_record_id               TEXT NOT NULL,
  rho_record_id              TEXT NOT NULL,
  n_or_bits                  INTEGER NOT NULL,
  ic_online_wall_ms          REAL NOT NULL,
  rho_online_wall_ms         REAL NOT NULL,
  online_speedup             REAL NOT NULL,                   -- rho / IC
  independent_validation     INTEGER NOT NULL,
  replay_certificate_sha256  TEXT,                            -- NULL until a receipt is cited
  check_status               TEXT NOT NULL,                   -- PASS | FAIL, the vs_rho checker at load time
  verdict                    TEXT NOT NULL,
  json                       TEXT NOT NULL
);

-- ── Identity guards ──────────────────────────────────────────────────
--
-- An insert tolerates a conflict on its primary key only when the row is
-- the same row again.  These triggers turn a key that arrives with a
-- different identity (a slug on another ICV1 string, a workload id on
-- another target, a method id on other parameters) into an error instead
-- of a silently kept first row.  They fire before conflict resolution, so
-- a true duplicate still passes to `ON CONFLICT ... DO NOTHING`.

CREATE TRIGGER IF NOT EXISTS curves_identity BEFORE INSERT ON curves
WHEN EXISTS (SELECT 1 FROM curves WHERE slug = NEW.slug AND icv1 <> NEW.icv1)
BEGIN SELECT RAISE(ABORT, 'curve slug already loaded with a different ICV1 string'); END;

CREATE TRIGGER IF NOT EXISTS workloads_identity BEFORE INSERT ON workloads
WHEN EXISTS (SELECT 1 FROM workloads WHERE workload_id = NEW.workload_id AND workload_sha256 <> NEW.workload_sha256)
BEGIN SELECT RAISE(ABORT, 'workload id already loaded with a different workload hash'); END;

CREATE TRIGGER IF NOT EXISTS algorithms_identity BEFORE INSERT ON algorithms
WHEN EXISTS (SELECT 1 FROM algorithms WHERE method_id = NEW.method_id AND method_sha256 <> NEW.method_sha256)
BEGIN SELECT RAISE(ABORT, 'method id already loaded with a different method hash'); END;

CREATE TRIGGER IF NOT EXISTS factor_bases_identity BEFORE INSERT ON factor_bases
WHEN EXISTS (SELECT 1 FROM factor_bases WHERE fb_id = NEW.fb_id AND fb_sha256 <> NEW.fb_sha256)
BEGIN SELECT RAISE(ABORT, 'factor base id already loaded with a different factor-base hash'); END;

CREATE TRIGGER IF NOT EXISTS hosts_identity BEFORE INSERT ON hosts
WHEN EXISTS (SELECT 1 FROM hosts WHERE env_class_id = NEW.env_class_id AND env_class_sha256 <> NEW.env_class_sha256)
BEGIN SELECT RAISE(ABORT, 'host class id already loaded with a different class hash'); END;

CREATE TRIGGER IF NOT EXISTS runs_identity BEFORE INSERT ON runs
WHEN EXISTS (SELECT 1 FROM runs WHERE record_id = NEW.record_id AND (session_id <> NEW.session_id OR seq <> NEW.seq))
BEGIN SELECT RAISE(ABORT, 'record id already loaded for a different session or sequence number'); END;

CREATE TRIGGER IF NOT EXISTS claims_identity BEFORE INSERT ON claims
WHEN EXISTS (SELECT 1 FROM claims WHERE run_id = NEW.run_id
             AND (candidate_manifest_sha256 <> NEW.candidate_manifest_sha256
                  OR workload_manifest_sha256 <> NEW.workload_manifest_sha256
                  OR session_id <> NEW.session_id))
BEGIN SELECT RAISE(ABORT, 'claim run id already loaded for a different candidate, workload or session'); END;

-- ── Views ────────────────────────────────────────────────────────────

-- Every arm on every workload: the table AGENTS.md §2 asks for, one unit.
CREATE VIEW IF NOT EXISTS arm_workload_summary AS
SELECT r.session_id, r.arm, a.method, a.family, r.method_id,
       w.curve_slug, c.log2_r, r.workload_id,
       COUNT(*)                                          AS measured,
       SUM(r.status = 'verified')                        AS verified,
       AVG(CASE WHEN r.status = 'verified' THEN r.s END) AS mean_s,
       MIN(CASE WHEN r.status = 'verified' THEN r.s END) AS min_s,
       MAX(CASE WHEN r.status = 'verified' THEN r.s END) AS max_s,
       c.floor_s,
       AVG(CASE WHEN r.status = 'verified' THEN r.s END) / c.floor_s AS mean_ratio_to_floor,
       MAX(r.lower_bound)                                AS any_lower_bound,
       MIN(r.isolation_level)                            AS worst_level
FROM runs r
JOIN algorithms a ON a.method_id = r.method_id
JOIN workloads w ON w.workload_id = r.workload_id
JOIN curves c ON c.slug = w.curve_slug
WHERE r.warmup = 0
GROUP BY r.session_id, r.arm, r.workload_id;

-- Each method's mean S per curve across sessions, against the curve's floor.
CREATE VIEW IF NOT EXISTS method_by_curve AS
SELECT a.method, a.family, r.method_id, w.curve_slug, c.log2_r, c.floor_s,
       COUNT(*)                     AS verified_runs,
       AVG(r.s)                     AS mean_s,
       AVG(r.s) / c.floor_s         AS mean_ratio_to_floor,
       MAX(r.lower_bound)           AS any_lower_bound
FROM runs r
JOIN algorithms a ON a.method_id = r.method_id
JOIN workloads w ON w.workload_id = r.workload_id
JOIN curves c ON c.slug = w.curve_slug
WHERE r.warmup = 0 AND r.status = 'verified'
GROUP BY r.method_id, w.curve_slug;

-- Index-calculus runs with their factor base and phase split.
CREATE VIEW IF NOT EXISTS ic_phase_split AS
SELECT r.record_id, r.session_id, r.arm, w.curve_slug, f.fb_id, f.family AS fb_family,
       f.signed_points, f.columns, p.phase, p.gae, r.total_gae,
       p.gae / r.total_gae AS share
FROM runs r
JOIN phases p ON p.record_id = r.record_id
JOIN workloads w ON w.workload_id = r.workload_id
LEFT JOIN factor_bases f ON f.fb_id = r.fb_id
WHERE r.fb_id IS NOT NULL;

INSERT OR REPLACE INTO schema_info (key, value) VALUES ('ecbench_schema_version', '1');
