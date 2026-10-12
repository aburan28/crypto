PRAGMA foreign_keys = ON;

CREATE TABLE input_file (
  id TEXT PRIMARY KEY,
  repository_path TEXT NOT NULL UNIQUE,
  sha256 TEXT NOT NULL CHECK(length(sha256) = 64),
  bytes INTEGER NOT NULL CHECK(bytes > 0),
  evidence_grade TEXT NOT NULL,
  scope TEXT NOT NULL
) STRICT;

CREATE TABLE source_ref (
  id TEXT PRIMARY KEY,
  title TEXT NOT NULL,
  url TEXT,
  repository_path TEXT,
  CHECK ((url IS NULL) != (repository_path IS NULL))
) STRICT;

CREATE TABLE family (
  id TEXT PRIMARY KEY,
  name TEXT NOT NULL,
  category TEXT NOT NULL,
  scope TEXT NOT NULL,
  evidence_status TEXT NOT NULL,
  local_status TEXT NOT NULL,
  domain TEXT NOT NULL,
  recognition TEXT NOT NULL,
  construction TEXT NOT NULL,
  exhaustive_scope TEXT NOT NULL,
  subgroup_gate TEXT NOT NULL,
  cost_gate TEXT NOT NULL,
  source_ids_json TEXT NOT NULL CHECK(json_valid(source_ids_json)),
  record_json TEXT NOT NULL CHECK(json_valid(record_json)),
  input_id TEXT NOT NULL REFERENCES input_file(id)
) STRICT;
CREATE INDEX family_category_idx ON family(category, evidence_status);

CREATE TABLE curve_model (
  icv1 TEXT PRIMARY KEY,
  model_sha256 TEXT NOT NULL UNIQUE CHECK(length(model_sha256) = 64),
  slug TEXT NOT NULL UNIQUE,
  characteristic TEXT NOT NULL,
  extension_degree INTEGER NOT NULL CHECK(extension_degree > 0),
  field_bits INTEGER NOT NULL CHECK(field_bits > 0),
  field_cardinality TEXT NOT NULL,
  trace TEXT NOT NULL,
  curve_order TEXT NOT NULL,
  subgroup_order TEXT,
  inventory_law TEXT NOT NULL,
  class_label TEXT NOT NULL,
  cubic_model_label TEXT NOT NULL,
  quadratic_model_label TEXT NOT NULL,
  subgroup_advantage TEXT NOT NULL,
  model_json TEXT NOT NULL CHECK(json_valid(model_json)),
  record_json TEXT NOT NULL CHECK(json_valid(record_json)),
  input_id TEXT NOT NULL REFERENCES input_file(id)
) STRICT;
CREATE INDEX curve_model_field_idx ON curve_model(field_bits, extension_degree, characteristic);
CREATE INDEX curve_model_labels_idx ON curve_model(class_label, cubic_model_label, quadratic_model_label);

CREATE TABLE curve_property (
  icv1 TEXT PRIMARY KEY REFERENCES curve_model(icv1),
  properties_json TEXT NOT NULL CHECK(json_valid(properties_json)),
  family_checks_json TEXT NOT NULL CHECK(json_valid(family_checks_json)),
  witnesses_json TEXT NOT NULL CHECK(json_valid(witnesses_json)),
  record_json TEXT NOT NULL CHECK(json_valid(record_json)),
  input_id TEXT NOT NULL REFERENCES input_file(id)
) STRICT;

CREATE TABLE class_observation (
  characteristic INTEGER NOT NULL CHECK(characteristic > 3),
  trace INTEGER NOT NULL,
  mechanism_scope TEXT NOT NULL CHECK(mechanism_scope IN ('two_branch_union', 'cubic_norm_one')),
  trace_status TEXT NOT NULL CHECK(trace_status IN ('ordinary', 'nonordinary_or_unrealized')),
  label_grade TEXT NOT NULL CHECK(label_grade IN ('complete_exact', 'probabilistic_diagnostic')),
  cubic_support INTEGER CHECK(cubic_support IN (0, 1)),
  quadratic_support INTEGER CHECK(quadratic_support IN (0, 1)),
  combined_support INTEGER CHECK(combined_support IN (0, 1)),
  representative_count INTEGER CHECK(representative_count >= 0),
  frobenius_2_depth INTEGER,
  fundamental_discriminant TEXT,
  two_split INTEGER,
  v2_delta_quarter INTEGER,
  class_number_odd INTEGER,
  source_line INTEGER NOT NULL CHECK(source_line >= 2),
  input_id TEXT NOT NULL REFERENCES input_file(id),
  PRIMARY KEY (characteristic, trace, mechanism_scope),
  CHECK (
    (mechanism_scope = 'two_branch_union' AND label_grade = 'complete_exact'
      AND trace_status = 'ordinary' AND cubic_support IS NOT NULL
      AND quadratic_support IS NOT NULL AND combined_support IS NOT NULL
      AND representative_count IS NULL)
    OR
    (mechanism_scope = 'cubic_norm_one' AND label_grade = 'probabilistic_diagnostic'
      AND cubic_support IS NULL AND quadratic_support IS NULL
      AND combined_support IS NULL AND representative_count IS NOT NULL)
  )
) STRICT;
CREATE INDEX class_observation_scope_idx ON class_observation(mechanism_scope, label_grade, characteristic, combined_support);

CREATE TABLE conjecture (
  id TEXT PRIMARY KEY,
  title TEXT NOT NULL,
  status TEXT NOT NULL,
  kind TEXT NOT NULL,
  domain TEXT NOT NULL,
  falsifier TEXT NOT NULL,
  lean_status TEXT NOT NULL,
  record_json TEXT NOT NULL CHECK(json_valid(record_json)),
  input_id TEXT NOT NULL REFERENCES input_file(id)
) STRICT;

CREATE TABLE heuristic (
  id TEXT PRIMARY KEY,
  statement TEXT NOT NULL,
  scope TEXT NOT NULL,
  status TEXT NOT NULL,
  calibration TEXT NOT NULL,
  falsifier TEXT NOT NULL,
  linked_conjecture TEXT REFERENCES conjecture(id),
  source_ids_json TEXT NOT NULL CHECK(json_valid(source_ids_json)),
  record_json TEXT NOT NULL CHECK(json_valid(record_json)),
  input_id TEXT NOT NULL REFERENCES input_file(id)
) STRICT;

CREATE TABLE theorem (
  id TEXT PRIMARY KEY,
  name TEXT NOT NULL,
  exact_statement TEXT NOT NULL,
  hypotheses TEXT NOT NULL,
  result_status TEXT NOT NULL,
  lean_declaration TEXT,
  lean_status TEXT NOT NULL,
  open_bridge TEXT,
  source_ids_json TEXT NOT NULL CHECK(json_valid(source_ids_json)),
  record_json TEXT NOT NULL CHECK(json_valid(record_json)),
  input_id TEXT NOT NULL REFERENCES input_file(id)
) STRICT;

CREATE TABLE mechanism (
  id TEXT PRIMARY KEY,
  family_id TEXT REFERENCES family(id),
  name TEXT NOT NULL,
  construction TEXT NOT NULL,
  recognition TEXT NOT NULL,
  witness_requirements TEXT NOT NULL,
  selected_subgroup_gate TEXT NOT NULL,
  measured_cost_status TEXT NOT NULL,
  validation_status TEXT NOT NULL,
  source_ids_json TEXT NOT NULL CHECK(json_valid(source_ids_json)),
  record_json TEXT NOT NULL CHECK(json_valid(record_json)),
  input_id TEXT NOT NULL REFERENCES input_file(id)
) STRICT;

CREATE TABLE family_proof_obligation (
  family_id TEXT PRIMARY KEY REFERENCES family(id),
  theorem_id TEXT REFERENCES theorem(id),
  heuristic_id TEXT REFERENCES heuristic(id),
  principal_claim TEXT NOT NULL,
  proof_state TEXT NOT NULL,
  next_bridge TEXT NOT NULL,
  record_json TEXT NOT NULL CHECK(json_valid(record_json)),
  input_id TEXT NOT NULL REFERENCES input_file(id)
) STRICT;

CREATE VIEW class_summary AS
SELECT characteristic, mechanism_scope, label_grade,
       COUNT(*) AS rows,
       SUM(CASE WHEN trace_status='ordinary' THEN 1 ELSE 0 END) AS ordinary_rows,
       SUM(CASE WHEN mechanism_scope='two_branch_union' AND combined_support=1 THEN 1
                WHEN mechanism_scope='cubic_norm_one' AND trace_status='ordinary' AND representative_count>0 THEN 1
                ELSE 0 END) AS positive_ordinary_rows,
       SUM(CASE WHEN mechanism_scope='two_branch_union' AND combined_support=0 THEN 1
                WHEN mechanism_scope='cubic_norm_one' AND trace_status='ordinary' AND representative_count=0 THEN 1
                ELSE 0 END) AS zero_ordinary_rows
FROM class_observation
GROUP BY characteristic, mechanism_scope, label_grade;
