"""Tests for ICMS.  Standard library only (PyYAML for the YAML specs).

    python3 -m unittest discover -s tools/icms/tests -v
"""
from __future__ import annotations

import copy
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

TOOLS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REPO = os.path.dirname(TOOLS)
sys.path.insert(0, TOOLS)

from icms import adapters  # noqa: E402
from icms.canonical import CanonicalError, canonical_json, candidate_label, sha256_hex, spec_id  # noqa: E402
from icms.compare import bootstrap_median_ratio, compare, spec_differences  # noqa: E402
from icms.environment import capsule, env_class, parse_cpu_list, stable_capsule  # noqa: E402
from icms.execute import build_env, run_measured  # noqa: E402
from icms.gates import DEFAULT_THRESHOLDS, effective_thresholds, evaluate  # noqa: E402
from icms.registry import Registry  # noqa: E402
from icms.schemacheck import Validator, audit_schema  # noqa: E402
from icms.audit import audit_session, validators  # noqa: E402
from icms.session import SessionError, auto_cpus, plan, run_session  # noqa: E402
from icms.spec import SpecError, identity_view, load, validate  # noqa: E402

def _read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def _dump(obj, path):
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh)


FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
SPECS = os.path.join(REPO, "docs", "ic", "measurement", "specs")


def base_spec() -> dict:
    return {
        "icms": "icms.spec/v1",
        "label": "test",
        "instance": {"curve": {"regime": "koblitz", "degree": 23, "koblitz_a": 1},
                     "workload": {"targets": 1, "law": "known_answer", "seeds": [7]}},
        "factor_base": {"family": "binary_subspace", "params": {"dimension": 8, "basis_law": "prefix"}, "quotient": ["negation"]},
        "decomposition": {"arity": 2, "method": "mitm"},
        "relations": {"collector": "walk", "stop": "first_log"},
        "linear_algebra": {"method": "incremental_gauss"},
        "reference": {"rho": "rho.measured_matched", "rho_runs": 2},
        "accounting": {"unit": "crypto.S.gae_pinned"},
        "measurement": {"window": "cold_end_to_end", "repetitions": 5},
        "execution": {"adapter": "crypto.ic_bench", "threads": 1},
    }


class CanonicalTests(unittest.TestCase):
    def test_sorted_compact_ascii(self):
        self.assertEqual(canonical_json({"b": 1, "a": [2, "é"]}), '{"a":[2,"\\u00e9"],"b":1}')

    def test_floats_refused_in_identity(self):
        with self.assertRaises(CanonicalError):
            sha256_hex({"x": 0.5})

    def test_candidate_label_needs_an_actual_point_count(self):
        self.assertIsNone(candidate_label(23, "kb1", None, 2, "mitm", "walk", "gauss", "pdp", 0, {"a": 1}))
        lab = candidate_label(23, "kb1", 236, 2, "mitm", "walk", "gauss", "pdp", 0, {"a": 1})
        self.assertRegex(lab, r"^IC1N23Ckb1fb236PDP2mitmRCwalkLAgaussTDpdpISO0h[0-9a-f]{12}$")


class SchemaTests(unittest.TestCase):
    def setUp(self):
        self.schema = json.loads(_read(os.path.join(REPO, "docs", "ic", "measurement", "schema", "spec.v1.json")))

    def test_schema_uses_only_supported_keywords(self):
        self.assertEqual(audit_schema(self.schema), [])

    def test_unknown_keyword_is_refused_at_load(self):
        with self.assertRaises(ValueError):
            Validator({"type": "object", "propertyNames": {"pattern": "x"}})

    def test_base_spec_is_valid(self):
        self.assertEqual(Validator(self.schema).errors(base_spec()), [])

    def test_unknown_key_and_missing_required(self):
        s = base_spec()
        s["bogus"] = 1
        del s["relations"]
        errs = Validator(self.schema).errors(s)
        self.assertTrue(any("unknown key 'bogus'" in e for e in errs))
        self.assertTrue(any("missing required 'relations'" in e for e in errs))

    def test_algebraic_needs_encoding_and_solver(self):
        s = base_spec()
        s["decomposition"] = {"arity": 2, "method": "algebraic"}
        errs = Validator(self.schema).errors(s)
        self.assertTrue(any("encoding" in e for e in errs) and any("solver" in e for e in errs))


class SpecTests(unittest.TestCase):
    def test_defaults_do_not_change_identity(self):
        a = validate(base_spec())
        explicit = base_spec()
        explicit["measurement"].update({"warmup": 1, "interleave": "abab", "isolation_required": "L2"})
        b = validate(explicit)
        self.assertEqual(load_id(a), load_id(b))

    def test_label_is_not_identity(self):
        a, b = base_spec(), base_spec()
        b["label"] = "something else"
        self.assertEqual(load_id(validate(a)), load_id(validate(b)))

    def test_dimension_is_identity(self):
        b = base_spec()
        b["factor_base"]["params"]["dimension"] = 7
        self.assertNotEqual(load_id(validate(base_spec())), load_id(validate(b)))

    def test_unknown_family_refused(self):
        s = base_spec()
        s["factor_base"]["family"] = "made_up"
        with self.assertRaises(SpecError):
            validate(s)

    def test_family_params_required(self):
        s = base_spec()
        s["factor_base"]["params"] = {}
        with self.assertRaises(SpecError) as cm:
            validate(s)
        self.assertTrue(any("dimension" in p for p in cm.exception.problems))

    def test_primary_window_is_one_target(self):
        s = base_spec()
        s["measurement"]["window"] = "online_one_target"
        s["instance"]["workload"]["targets"] = 4
        with self.assertRaises(SpecError):
            validate(s)

    def test_sat_solver_must_pin_search_options(self):
        s = base_spec()
        s["decomposition"] = {"arity": 2, "method": "algebraic", "encoding": "expanded_semaev",
                              "solver": {"name": "sat-cdcl", "options": {"sat_conflict_budget": 1000}}}
        with self.assertRaises(SpecError) as cm:
            validate(s)
        self.assertTrue(any("sat_xor_encoding" in p for p in cm.exception.problems))

    def test_threshold_floats_hash_exactly(self):
        s = base_spec()
        s["measurement"]["thresholds"] = {"max_run_delay_fraction": 0.001}
        self.assertIn("'0.001'", repr(identity_view(validate(s))["measurement"]["thresholds"]))

    def test_shipped_specs_validate_and_adapters_accept_them(self):
        reg = Registry.load()
        files = sorted(f for f in os.listdir(SPECS) if f.endswith((".yaml", ".yml", ".json")))
        self.assertTrue(files)
        for f in files:
            sp = load(os.path.join(SPECS, f), reg)
            self.assertEqual(adapters.get(sp["spec"]["execution"]["adapter"]).check(sp["spec"]), [], f)

    def test_workloads_one_per_seed(self):
        s = base_spec()
        s["instance"]["workload"]["seeds"] = [1, 2, 3]
        from icms.spec import workloads
        w = workloads(validate(s))
        self.assertEqual(len({x["workload_id"] for x in w}), 3)


def load_id(norm):
    from icms.spec import spec_id as sid
    return sid(norm)


class GateTests(unittest.TestCase):
    STABLE = {"topology": {"isolated": "", "nohz_full": None, "smt_active": "0", "smt_control": "notsupported",
                           "intel_pstate_no_turbo": None, "cpufreq_boost": None, "cpus": {"3": {"governor": None}}},
              "virtualization": {"detect_virt": "kvm"}}
    SESSION = {"reservation": {"evicted": True, "threads_moved": 3, "left_on_reserved": {}},
               "preflight": {"quiet": True, "other_cpu_seconds": 0.0, "psi_some_avg10_max": 0.0}}

    def execution(self, **over):
        ex = {"pinned_cpus": [3], "child_affinity_observed": [3], "wall_ns": 1_000_000_000,
              "schedstat": {"run_ns": 990_000_000, "wait_ns": 1_000_000, "threads_seen": 1},
              "rusage": {"nivcsw": 3},
              "contention": {"samples": 20, "samples_with_foreign_runnable": 0,
                             "pinned_cpu_jiffies_delta": {"user": 99, "steal": 0, "idle": 1},
                             "psi_total_delta_us": {"memory": {"some": 0.0, "full": 0.0}}},
              "max_threads_observed": 2}
        ex.update(over)
        return ex

    def test_quiet_vm_run_earns_l2_not_l3(self):
        g = evaluate(self.execution(), self.SESSION, self.STABLE, 1)
        self.assertEqual(g["earned_level"], "L2")
        self.assertIn("bare_metal", g["blocking_next_level"])

    def test_cpu0_never_earns_l1(self):
        g = evaluate(self.execution(pinned_cpus=[0], child_affinity_observed=[0]), self.SESSION, self.STABLE, 1)
        self.assertEqual(g["earned_level"], "L0")

    def test_steal_blocks_l2(self):
        ex = self.execution()
        ex["contention"]["pinned_cpu_jiffies_delta"]["steal"] = 2
        self.assertEqual(evaluate(ex, self.SESSION, self.STABLE, 1)["earned_level"], "L1")

    def test_oversubscribed_threads_fail_run_delay(self):
        ex = self.execution(schedstat={"run_ns": 600_000_000, "wait_ns": 1_800_000_000, "threads_seen": 5})
        g = evaluate(ex, self.SESSION, self.STABLE, 1)
        self.assertEqual(g["earned_level"], "L1")
        self.assertIn("run_delay", g["blocking_next_level"])

    def test_unknown_never_passes(self):
        ex = self.execution(schedstat={})
        self.assertEqual(evaluate(ex, self.SESSION, self.STABLE, 1)["earned_level"], "L0")

    def test_busy_preflight_caps_at_l1(self):
        sess = dict(self.SESSION, preflight={"quiet": False})
        self.assertEqual(evaluate(self.execution(), sess, self.STABLE, 1)["earned_level"], "L1")

    def test_isolated_bare_metal_earns_l3(self):
        st = copy.deepcopy(self.STABLE)
        st["topology"].update({"isolated": "2-3", "intel_pstate_no_turbo": "1", "cpus": {"3": {"governor": "performance"}}})
        st["virtualization"]["detect_virt"] = "none"
        self.assertEqual(evaluate(self.execution(), self.SESSION, st, 1)["earned_level"], "L3")

    def test_thresholds_can_only_tighten(self):
        self.assertEqual(effective_thresholds({"max_run_delay_fraction": 0.001})["max_run_delay_fraction"], 0.001)
        with self.assertRaises(ValueError):
            effective_thresholds({"max_run_delay_fraction": 0.5})
        with self.assertRaises(ValueError):
            effective_thresholds({"not_a_threshold": 1})
        self.assertEqual(set(DEFAULT_THRESHOLDS), set(effective_thresholds(None)))


class EnvironmentTests(unittest.TestCase):
    def test_capsule_has_every_section_and_a_stable_class(self):
        a, b = capsule(), capsule()
        for k in ("cpu", "topology", "kernel", "virtualization", "cgroup", "memory", "os", "toolchain"):
            self.assertIn(k, a["stable"])
        self.assertEqual(a["env_class_id"], b["env_class_id"])
        self.assertRegex(a["env_class_id"], r"^ENV1h[0-9a-f]{12}$")
        self.assertIn("vm/swappiness", a["stable"]["kernel"]["sysctl"])

    def test_cpu_list_parsing(self):
        self.assertEqual(parse_cpu_list("0-2,5"), [0, 1, 2, 5])
        self.assertEqual(parse_cpu_list(""), [])
        self.assertEqual(parse_cpu_list(None), [])

    def test_env_class_ignores_volatile_facts(self):
        st = stable_capsule()
        st2 = copy.deepcopy(st)
        st2["captured_unix"] = 0
        self.assertEqual(env_class(st), env_class(st2))


class ExecuteTests(unittest.TestCase):
    def test_pinned_child_is_measured_and_its_environment_is_built(self):
        cpu = max(os.sched_getaffinity(0))
        os.environ["KIC_TEST_KNOB"] = "1"
        try:
            env, rec = build_env({"DECLARED": "yes"})
        finally:
            del os.environ["KIC_TEST_KNOB"]
        self.assertNotIn("KIC_TEST_KNOB", env)
        self.assertIn("KIC_TEST_KNOB", rec["dropped_engine_knobs"])
        self.assertEqual(env["RAYON_NUM_THREADS"], "1")
        with tempfile.TemporaryDirectory() as d:
            code = ("import os,sys\nassert os.environ['DECLARED']=='yes'\nassert 'KIC_TEST_KNOB' not in os.environ\n"
                    "s=0\nfor i in range(300000): s+=i\nprint(sorted(os.sched_getaffinity(0)))")
            ex = run_measured([sys.executable, "-c", code], {cpu}, env, d, os.path.join(d, "o"), os.path.join(d, "e"),
                              timeout_s=60, interval=0.01)
            self.assertEqual(ex["exit"]["returncode"], 0, _read(os.path.join(d, "e")))
            self.assertEqual(ex["child_affinity_observed"], [cpu])
            self.assertEqual(_read(os.path.join(d, "o")).strip(), str([cpu]))
            self.assertGreater(ex["schedstat"]["run_ns"], 0)
            self.assertGreaterEqual(ex["schedstat"]["threads_seen"], 1)
            self.assertIsNotNone(ex["outputs"]["stdout"]["sha256"])

    def test_timeout_kills_the_process_group(self):
        cpu = max(os.sched_getaffinity(0))
        with tempfile.TemporaryDirectory() as d:
            env, _ = build_env({})
            ex = run_measured([sys.executable, "-c", "import time; time.sleep(30)"], {cpu}, env, d,
                              os.path.join(d, "o"), os.path.join(d, "e"), timeout_s=0.5, interval=0.05)
            self.assertTrue(ex["exit"]["timed_out"])
            self.assertLess(ex["wall_ns"], 10_000_000_000)


class CryptoIcBenchAdapterTests(unittest.TestCase):
    def setUp(self):
        self.ad = adapters.get("crypto.ic_bench")
        self.ctx = adapters.Context(repo_root=REPO, exec_dir="/nonexistent")

    def test_command_renders_the_plugin_grammar(self):
        spec = validate(base_spec())
        cmd = self.ad.command(spec, {"record": {"seed": 7}}, self.ctx)
        argv = cmd["argv"]
        self.assertEqual(argv[1], "bench")
        self.assertIn("binary-subspace:dimension=8", argv)
        self.assertIn("mitm:m=2", argv)
        self.assertEqual(argv[argv.index("--seed") + 1], "7")
        self.assertIn("--json", argv)

    def test_divisor_lists_join_with_semicolons(self):
        s = base_spec()
        s["factor_base"] = {"family": "frobenius_stable_subspace", "params": {"divisor": [0, 2]}, "quotient": ["negation"]}
        cmd = self.ad.command(validate(s), {"record": {"seed": 1}}, self.ctx)
        self.assertIn("koblitz-orbit:divisor=0;2,no_fold=1", cmd["argv"])

    def test_refusals(self):
        s = base_spec()
        s["relations"]["stop"] = "full_rank"
        s["factor_base"]["large_primes"] = {"mode": "double"}
        probs = self.ad.check(validate(s))
        self.assertTrue(any("first_log" in p for p in probs))
        self.assertTrue(any("large-prime" in p for p in probs))

    def test_parse_pair_table_row(self):
        out = self.ad.parse(validate(base_spec()), {"record": {"seed": 7}}, self.ctx,
                            os.path.join(FIX, "ic-bench-k23-mitm2.json"))
        self.assertEqual(out["outcome"]["status"], "complete")
        fb = out["metrics"]["factor_base"]
        self.assertEqual((fb["signed_points"], fb["abscissae_with_points"], fb["columns"]), (236, 118, 118))
        self.assertTrue(out["units"]["crypto.S.gae_pinned"]["deterministic"])
        self.assertEqual(out["windows"]["online_one_target"]["conformance"], "derived")
        self.assertEqual(out["consistency"], [])
        cold = out["windows"]["cold_end_to_end"]["ops"]
        parts = sum(p["ops"] for p in out["phases"].values())
        self.assertAlmostEqual(cold, parts, places=6)

    def test_parse_sat_row_is_host_dependent(self):
        s = base_spec()
        s["instance"]["curve"]["degree"] = 17
        s["factor_base"]["params"]["dimension"] = 6
        s["decomposition"] = {"arity": 2, "method": "algebraic", "encoding": "expanded_semaev",
                              "solver": {"name": "sat-cdcl", "options": {"sat_xor_encoding": "native", "sat_conflict_budget": 200000}}}
        out = self.ad.parse(validate(s), {"record": {"seed": 1}}, self.ctx, os.path.join(FIX, "ic-bench-k17-sat.json"))
        u = out["units"]["crypto.S.gae_pinned"]
        self.assertFalse(u["deterministic"])
        self.assertTrue(any("wall time" in r for r in u["host_dependent_because"]))
        self.assertEqual(out["metrics"]["system"]["n_vars"], 12)
        self.assertEqual(out["metrics"]["solver"]["op_unit"], "conflicts")
        self.assertIsNone(out["metrics"]["solver"]["peak_bytes"], "an unmeasured 0 must become null")

    def test_declared_quotient_mismatch_is_reported(self):
        s = base_spec()
        s["factor_base"]["quotient"] = []
        out = self.ad.parse(validate(s), {"record": {"seed": 7}}, self.ctx, os.path.join(FIX, "ic-bench-k23-mitm2.json"))
        self.assertTrue(out["consistency"])


class ReviewRegressionTests(unittest.TestCase):
    """Defects an adversarial review reproduced; each must stay fixed."""

    def toy(self, **edit):
        s = toy_spec(8)
        for path, value in edit.items():
            cur = s
            keys = path.split("__")
            for k in keys[:-1]:
                cur = cur.setdefault(k, {})
            cur[keys[-1]] = value
        return s

    def test_koblitz_curve_must_name_a(self):
        s = self.toy()
        del s["instance"]["curve"]["koblitz_a"]
        s["instance"]["curve"]["ref"] = "icv1-f2m11-x"
        with self.assertRaises(SpecError):
            validate(s)

    def test_boolean_is_not_an_integer_enum_value(self):
        with self.assertRaises(SpecError):
            validate(self.toy(instance__curve__koblitz_a=True))
        self.assertTrue(Validator({"enum": [0, 1]}).errors(False))
        self.assertEqual(Validator({"enum": [1]}).errors(1.0), [], "1 and 1.0 are one JSON number")

    def test_thresholds_checked_before_anything_runs(self):
        for th in ({"max_steal_fraction": 0.5}, {"max_steal_fractoin": 0.0}, {"max_steal_fraction": "0"}):
            with self.assertRaises(SpecError, msg=th):
                validate(self.toy(measurement__thresholds=th))

    def test_written_out_defaults_share_an_id(self):
        same = [(self.toy(factor_base__large_primes={"mode": "double"}),
                 self.toy(factor_base__large_primes={"mode": "double", "max_large_primes": 2})),
                (self.toy(), self.toy(measurement__thresholds={"max_steal_fraction": 0.0})),
                (self.toy(measurement__thresholds={"max_involuntary_switches_per_s": 40}),
                 self.toy(measurement__thresholds={"max_involuntary_switches_per_s": 40.0})),
                (self.toy(), self.toy(reference={}))]
        for a, b in same:
            self.assertEqual(spec_id(identity_view(validate(a))), spec_id(identity_view(validate(b))))

    def test_fixed_file_law_needs_its_file(self):
        with self.assertRaises(SpecError):
            validate(self.toy(instance__workload__law="fixed_file"))
        with self.assertRaises(SpecError):
            validate(self.toy(instance__workload__file={"path": "t.dat", "sha256": "0" * 64}))

    def test_null_family_parameter_is_missing(self):
        s = self.toy(factor_base={"family": "binary_subspace", "params": {"dimension": 8, "basis_law": None}})
        with self.assertRaises(SpecError):
            validate(s)

    def test_symmetrised_solver_becomes_oracle_parameters(self):
        ad = adapters.get("crypto.ic_bench")
        s = base_spec()
        s["factor_base"] = {"family": "frobenius_symmetrised", "params": {"divisor": [1]}, "quotient": ["negation", "frobenius"]}
        s["decomposition"] = {"arity": 2, "method": "symmetrised", "encoding": "symmetrized", "limits": {"per_call_seconds": 120},
                              "solver": {"name": "matrix-f4", "options": {"max_degree": 5, "split": "auto", "node_budget": 100}}}
        s["relations"]["max_trials"] = 1000
        norm = validate(s)
        self.assertEqual(ad.check(norm), [])
        argv = ad.command(norm, {"record": {"seed": 1}}, adapters.Context(repo_root=REPO))["argv"]
        self.assertNotIn("--solver", argv, "ic bench ignores --solver for the symmetrised oracle")
        self.assertIn("symmetrised:engine=matrix-f4,m=2,max_degree=5,node_budget=100", argv)
        norm["decomposition"]["solver"]["options"]["split"] = "mom"
        self.assertTrue(any("split" in x for x in ad.check(norm)))

    def test_adapters_refuse_curves_they_do_not_resolve(self):
        for name, spec_file in (("crypto.ic_bench", "k23-subspace-d7-mitm2.yaml"),
                                ("cryptanalysis.ic_bench", "ca-k0n13-prefix-l3-m3.yaml"),
                                ("autoresearcher.index_calculus", "ar-p16-smallx-m2-enum.yaml")):
            s = copy.deepcopy(load(os.path.join(SPECS, spec_file))["spec"])
            s["instance"]["curve"]["ref"] = "icv1-anything"
            self.assertTrue(any("ref" in x for x in adapters.get(name).check(s)), name)

    def test_validator_keyword_semantics(self):
        v = Validator({"type": "object", "properties": {"a.b": {"type": "integer"}},
                       "patternProperties": {"^a\\.": {"minimum": 10}}})
        self.assertTrue(v.errors({"a.b": 1}), "patternProperties also applies to a key in properties")
        with self.assertRaises(ValueError):
            Validator({"type": "array", "items": [{"type": "string"}]})
        self.assertTrue(Validator({"type": "array", "uniqueItems": True}).errors([{"a": 1, "b": 2}, {"b": 2, "a": 1}]))


class RunnerRegressionTests(unittest.TestCase):
    """Execution defects an adversarial review reproduced; each must stay fixed."""

    def setUp(self):
        try:
            self.cpu = min(auto_cpus())
        except SessionError as exc:
            self.skipTest(str(exc))
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def run_cmd(self, argv, timeout=20):
        d = self.tmp.name
        env, _ = build_env({})
        return run_measured(argv, {self.cpu}, env, d, os.path.join(d, "out"), os.path.join(d, "err"), timeout_s=timeout)

    def test_a_producer_that_repins_itself_fails_affinity(self):
        code = ("import os,time; os.sched_setaffinity(0, set(range(os.cpu_count()))); t=time.time()\n"
                "while time.time()-t<0.3: pass")
        ex = self.run_cmd([sys.executable, "-c", code])
        self.assertGreater(ex["tree_affinity"]["tasks_outside"], 0)
        gate = evaluate(ex, {}, {"topology": {}}, 1)
        self.assertEqual([c["status"] for c in gate["checks"] if c["id"] == "affinity_observed"], ["fail"])

    def test_nothing_the_child_starts_outlives_the_record(self):
        for script in ("sleep 30 & echo started", "setsid sleep 30 & echo started"):
            ex = self.run_cmd(["sh", "-c", script])
            self.assertGreaterEqual(ex["descendants_killed"], 1, script)
            out = subprocess.run(["pgrep", "-f", "^sleep 30$"], capture_output=True, text=True).stdout.split()
            self.assertEqual(out, [], f"{script}: a sleep survived")

    def test_launch_failure_has_the_full_record_shape(self):
        ex = self.run_cmd(["/nonexistent/producer"])
        ex["env"], ex["argv0"] = build_env({})[1], {"name": "/nonexistent/producer", "resolved": None, "sha256": None}
        schema = json.loads(_read(os.path.join(REPO, "docs", "ic", "measurement", "schema", "record.v1.json")))
        self.assertEqual(Validator(schema).errors(ex, schema["$defs"]["execution"]), [])
        self.assertIn("launch_error", ex)

    def test_own_children_are_not_foreign(self):
        pids = os.path.join(self.tmp.name, "pids")
        script = f"for i in $(seq 60); do sh -c 'echo $$ >> {pids}; i=0; while [ $i -lt 2000 ]; do i=$((i+1)); done'; done"
        ex = self.run_cmd(["sh", "-c", script])
        ours = set(_read(pids).split())
        self.assertFalse([e for e in ex["contention"]["foreign_examples"] if str(e["pid"]) in ours])

    def test_undersampled_schedstat_cannot_pass_parallelism_or_run_delay(self):
        ex = {"pinned_cpus": [3], "child_affinity_observed": [3], "wall_ns": 351_000_000,
              "schedstat": {"run_ns": 19_400_000, "wait_ns": 0, "threads_seen": 3},
              "rusage": {"user_s": 0.40, "sys_s": 0.055, "nivcsw": 0}, "contention": {}}
        checks = {c["id"]: c for c in evaluate(ex, {}, {"topology": {}}, 1)["checks"]}
        self.assertEqual(checks["parallelism"]["status"], "fail")
        self.assertEqual(checks["run_delay"]["status"], "unknown")

    def test_a_user_thread_left_on_the_reservation_fails_l1(self):
        session = {"reservation": {"evicted": True, "threads_moved": 3,
                                   "left_on_reserved": {"user_threads": ["python3[2739]"], "kernel_threads": 2}}}
        checks = {c["id"]: c for c in evaluate({"pinned_cpus": [3]}, session, {"topology": {}}, 1)["checks"]}
        self.assertEqual(checks["reservation"]["status"], "fail")

    def test_a_signal_is_not_a_completed_run(self):
        spec = toy_spec(8)
        spec["execution"]["args"]["argv"] = [
            "python3", "-c", "import os,runpy,sys; sys.argv=['p', sys.argv[1], '8']; "
            "runpy.run_path(os.environ['TOY']); os.kill(os.getpid(), 9)", "{seed}"]
        spec["execution"]["env"] = {"TOY": os.path.join(FIX, "toy_producer.py")}
        spec["measurement"].update({"repetitions": 1, "warmup": 0})
        sp = os.path.join(self.tmp.name, "killed.json")
        _dump(spec, sp)
        out = os.path.join(self.tmp.name, "s")
        run_session([sp], out, auto_cpus(), lock=os.path.join(self.tmp.name, "lock"), settle=0.1, allow_busy=True,
                    log=lambda *_: None)
        rec = json.loads(_read(os.path.join(out, "records.jsonl")).splitlines()[0])
        self.assertEqual(rec["execution"]["exit"]["signal"], 9)
        self.assertEqual(rec["outcome"]["status"], "error")
        self.assertFalse(rec["outcome"]["verified"])

    def test_refusals_leave_no_directory(self):
        a, b = os.path.join(self.tmp.name, "a"), os.path.join(self.tmp.name, "b")
        os.makedirs(a), os.makedirs(b)
        for d, size in ((a, 8), (b, 9)):
            _dump(toy_spec(size), os.path.join(d, "spec.json"))
        out = os.path.join(self.tmp.name, "s")
        with self.assertRaises(SessionError):
            run_session([os.path.join(a, "spec.json"), os.path.join(b, "spec.json")], out, auto_cpus(),
                        lock=os.path.join(self.tmp.name, "lock"), allow_busy=True, log=lambda *_: None)
        self.assertFalse(os.path.exists(out))
        import fcntl
        with open(os.path.join(self.tmp.name, "lock"), "a+") as held:
            fcntl.flock(held, fcntl.LOCK_EX)
            with self.assertRaises(SessionError):
                run_session([os.path.join(a, "spec.json")], out, auto_cpus(), lock=os.path.join(self.tmp.name, "lock"),
                            allow_busy=True, log=lambda *_: None)
        self.assertFalse(os.path.exists(out))


class OtherAdapterTests(unittest.TestCase):
    """The cryptanalysis and autoresearcher adapters against frozen producer output."""

    def setUp(self):
        self.reg = Registry.load()
        self.ctx = adapters.Context(repo_root=REPO, exec_dir="/nonexistent", out_dir="/nonexistent")

    def spec(self, name):
        return load(os.path.join(SPECS, name), self.reg)

    def record_shape(self, parsed, spec):
        """The normalised parse must fit the record schema's metric and unit sections."""
        out = adapters.normalise(parsed, spec, self.reg)
        schema = json.loads(_read(os.path.join(REPO, "docs", "ic", "measurement", "schema", "record.v1.json")))
        v = Validator(schema)
        for key, ref in (("metrics", "metrics"), ("outcome", "outcome")):
            self.assertEqual(v.errors(out[key], schema["$defs"][ref]), [], key)
        for uid, u in out["units"].items():
            self.assertTrue(self.reg.has_unit(uid), uid)
            self.assertEqual(v.errors(u, schema["$defs"]["unit_value"]), [], uid)
        for name, w in out["windows"].items():
            self.assertEqual(v.errors(w, schema["$defs"]["window"]), [], name)
        for name, ph in out["phases"].items():
            self.assertIn(name, self.reg.phases)
            self.assertEqual(v.errors(ph, schema["$defs"]["phase"]), [], name)
        return out

    def test_cryptanalysis_specs_accepted_and_cell_rendered(self):
        ad = adapters.get("cryptanalysis.ic_bench")
        sp = self.spec("ca-k0n13-prefix-l3-m3.yaml")
        self.assertEqual(ad.check(sp["spec"]), [])
        cell = ad.cell(sp["spec"], sp["workloads"][0])
        self.assertEqual(cell, {"n": 13, "m": 3, "l": 3, "family": "prefix", "seed": 1, "mode": "mxl",
                                "workload_seed": 1, "targets": 1, "max_attempts": 200000, "run": 1})
        argv = ad.command(sp["spec"], sp["workloads"][0], self.ctx)["argv"]
        self.assertNotIn("--record", argv, "the runner must never touch the cryptanalysis baseline")

    def test_cryptanalysis_refuses_crypto_conventions(self):
        ad = adapters.get("cryptanalysis.ic_bench")
        s = copy.deepcopy(self.spec("ca-k0n13-prefix-l3-m3.yaml")["spec"])
        s["reference"]["rho"] = "rho.measured_matched"
        s["relations"]["stop"] = "first_log"
        s["instance"]["curve"]["koblitz_a"] = 1
        probs = " ".join(ad.check(s))
        for needle in ("rho.measured_matched", "full_rank", "K_0"):
            self.assertIn(needle, probs)

    def test_cryptanalysis_parse(self):
        ad = adapters.get("cryptanalysis.ic_bench")
        sp = self.spec("ca-k0n13-prefix-l3-m3.yaml")
        out = self.record_shape(ad.parse(sp["spec"], sp["workloads"][0], self.ctx,
                                         os.path.join(FIX, "cryptanalysis-n13-prefix-l3-m3-cell.json")), sp["spec"])
        rec = json.loads(_read(os.path.join(FIX, "cryptanalysis-n13-prefix-l3-m3-cell.json")))["receipt"]
        self.assertEqual(out["outcome"]["status"], "complete")
        self.assertEqual(out["units"]["cryptanalysis.rps"]["total"], rec["total_operations"])
        self.assertEqual(sum(p["ops"] for p in out["phases"].values()), rec["total_operations"],
                         "the eleven exclusive phases sum to the cold total")
        fb = out["metrics"]["factor_base"]
        self.assertEqual((fb["usable_points"], fb["columns"], fb["abscissae_allowed"]), (8, 4, 8))
        self.assertIsNone(out["metrics"]["system"], "no Boolean-system shape in the receipt: null, not zero")
        self.assertEqual(out["reference"]["id"], "rho.signed_frobenius")
        self.assertEqual(out["windows"]["online_one_target"]["conformance"], "derived")

    def test_autoresearcher_spec_accepted_and_seeds_shared(self):
        ad = adapters.get("autoresearcher.index_calculus")
        sp = self.spec("ar-p16-smallx-m2-enum.yaml")
        self.assertEqual(ad.check(sp["spec"]), [])
        argv = ad.command(sp["spec"], sp["workloads"][0], self.ctx)["argv"]
        for flag in ("--curve-seed", "--target-seed", "--seed"):
            self.assertEqual(argv[argv.index(flag) + 1], "0")
        self.assertIn("--no-accel", argv)
        self.assertEqual(argv[argv.index("--la-pivot") + 1], "min_fill")

    def test_autoresearcher_refuses_a_strong_rho_claim(self):
        ad = adapters.get("autoresearcher.index_calculus")
        s = copy.deepcopy(self.spec("ar-p16-smallx-m2-enum.yaml")["spec"])
        s["reference"]["rho"] = "rho.negation"
        del s["linear_algebra"]["pivot"]
        probs = " ".join(ad.check(s))
        self.assertIn("rho.plain", probs)
        self.assertIn("pivot", probs)

    def test_autoresearcher_parse_counts_signed_points(self):
        ad = adapters.get("autoresearcher.index_calculus")
        sp = self.spec("ar-p16-smallx-m2-enum.yaml")
        out = self.record_shape(ad.parse(sp["spec"], sp["workloads"][0], self.ctx,
                                         os.path.join(FIX, "autoresearcher-p16-smallx-m2-enum.json")), sp["spec"])
        doc = json.loads(_read(os.path.join(FIX, "autoresearcher-p16-smallx-m2-enum.json")))
        ic = doc["index_calculus"]
        fb = out["metrics"]["factor_base"]
        self.assertEqual(fb["abscissae_with_points"], ic["factor_base"]["size"])
        self.assertEqual(fb["usable_points"], 2 * ic["factor_base"]["size"], "signed points, as crypto and cryptanalysis count")
        self.assertEqual(out["units"]["count.s3_solves"]["total"], ic["s3_solves"] + ic["table_s3_solves"])
        self.assertEqual(out["reference"]["id"], "rho.plain")


class CompareTests(unittest.TestCase):
    def test_plan_interleaves_after_warmups(self):
        arms = [{"spec": validate(base_spec())}, {"spec": validate(base_spec())}]
        order = plan(arms, 0)
        self.assertEqual(order[:2], [(0, -1, True), (1, -1, True)])
        self.assertEqual([a for a, r, w in order[2:6]], [0, 1, 0, 1])

    def test_bootstrap(self):
        e = bootstrap_median_ratio([(1.0, 1.1)] * 3 + [(1.0, 1.2)] * 3)
        self.assertEqual(e["n_pairs"], 6)
        self.assertTrue(e["excludes_1"])
        self.assertIsNone(bootstrap_median_ratio([(1.0, 1.1)] * 2)["ci95"])

    def test_spec_differences_classify(self):
        a, b = base_spec(), base_spec()
        b["factor_base"]["params"]["dimension"] = 7
        b["instance"]["curve"]["degree"] = 29
        d = spec_differences(a, b, ["factor_base.params.dimension"])
        self.assertEqual([x["field"] for x in d["declared"]], ["factor_base.params.dimension"])
        self.assertEqual([x["field"] for x in d["forbidden"]], ["instance.curve.degree"])
        d2 = spec_differences(a, b, [])
        self.assertIn("factor_base.params.dimension", [x["field"] for x in d2["confound"]])

    def test_compare_on_a_synthetic_session(self):
        with tempfile.TemporaryDirectory() as d:
            sa, sb = base_spec(), base_spec()
            sb["factor_base"]["params"]["dimension"] = 7
            os.makedirs(os.path.join(d, "specs"))
            for name, s in (("a.json", sa), ("b.json", sb)):
                _dump(s, os.path.join(d, "specs", name))
            session = {"session_id": "S", "env_class_id": "ENV1hx", "spec_files": {"a.json": "", "b.json": ""},
                       "records_sha256": "0" * 64,
                       "arms": [{"index": 0, "spec_path": "a.json", "spec_id": "A", "workload_id": "W", "label": "a"},
                                {"index": 1, "spec_path": "b.json", "spec_id": "B", "workload_id": "W", "label": "b"},
                                {"index": 2, "spec_path": "a.json", "spec_id": "A", "workload_id": "W", "label": "a"}]}
            _dump(session, os.path.join(d, "session.json"))
            with open(os.path.join(d, "records.jsonl"), "w") as fh:
                for rnd in range(5):
                    for arm, ops, wall in ((0, 100.0, 10.0 + rnd * 0.01), (1, 80.0, 8.0), (2, 100.0, 10.0)):
                        fh.write(json.dumps({"record_id": f"r{arm}{rnd}", "arm": arm, "round": rnd, "warmup": False,
                                             "outcome": {"status": "complete", "verified": True},
                                             "units": {"crypto.S.gae_pinned": {"total": ops, "deterministic": True}},
                                             "windows": {"cold_end_to_end": {"ops": ops, "wall_ns": wall, "conformance": "exact_operations"}},
                                             "isolation": {"earned_level": "L2", "wall_admissible": True, "blocking_next_level": []},
                                             "consistency": [], "metrics": {"factor_base": {"columns": ops}}}) + "\n")
            res = compare(d, 0, 1, ["factor_base.params.dimension"])
            self.assertEqual(res["refusals"], [])
            self.assertAlmostEqual(res["ops"]["ratio_b_over_a"], 0.8)
            self.assertTrue(res["ops"]["deterministic"])
            self.assertTrue(res["wall"]["admitted"])
            self.assertLess(res["wall"]["estimate"]["ci95"][1], 1)
            self.assertEqual(res["aa_noise"][0]["arms"], [0, 2])
            refused = compare(d, 0, 1, [])
            self.assertTrue(refused["refusals"])
            self.assertFalse(refused["ops"]["admitted"])


def toy_spec(size: int) -> dict:
    return {
        "icms": "icms.spec/v1",
        "label": f"toy producer, base size {size}",
        "role": "candidate",
        "instance": {"curve": {"regime": "koblitz", "degree": 11, "koblitz_a": 0},
                     "workload": {"targets": 1, "law": "known_answer", "seeds": [5]}},
        "factor_base": {"family": "explicit", "params": {"recipe": f"toy-{size}"}},
        "decomposition": {"arity": 2, "method": "mitm"},
        "relations": {"collector": "walk", "stop": "first_log"},
        "linear_algebra": {"method": "incremental_gauss"},
        "accounting": {"unit": "count.group_additions"},
        "measurement": {"window": "cold_end_to_end", "repetitions": 5, "warmup": 1, "isolation_required": "L1",
                        "timeout_seconds": 60},
        "execution": {"adapter": "command", "threads": 1,
                      "args": {"argv": ["python3", "{repo_root}/tools/icms/tests/fixtures/toy_producer.py",
                                        "{seed}", str(size)]}},
    }


class SessionAuditTests(unittest.TestCase):
    """A real session through the command adapter, then the audit, then forgeries."""

    @classmethod
    def setUpClass(cls):
        try:
            cls.cpus = auto_cpus()
        except SessionError as exc:
            raise unittest.SkipTest(f"no reservable core: {exc}")
        cls.tmp = tempfile.TemporaryDirectory()
        d = cls.tmp.name
        for size in (8, 9):
            _dump(toy_spec(size), os.path.join(d, f"toy{size}.json"))
        cls.session_dir = os.path.join(d, "session")
        a, b = os.path.join(d, "toy8.json"), os.path.join(d, "toy9.json")
        cls.session = run_session([a, b, a], cls.session_dir, cls.cpus, lock=os.path.join(d, "lock"), settle=0.2,
                                  allow_busy=True, session_label="unit-test", log=lambda *_: None)
        cls.declared = ["factor_base.params.recipe", "execution.args.argv"]
        res = compare(cls.session_dir, 0, 1, cls.declared)
        os.makedirs(os.path.join(cls.session_dir, "comparisons"))
        with open(os.path.join(cls.session_dir, "comparisons", "toy9-vs-toy8.json"), "w") as fh:
            fh.write(json.dumps(res, indent=1, sort_keys=True) + "\n")
        cls.vals = validators()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def forged(self) -> str:
        dst = tempfile.mkdtemp(dir=self.tmp.name)
        dst = os.path.join(dst, "session")
        shutil.copytree(self.session_dir, dst)
        return dst

    def records(self, d):
        with open(os.path.join(d, "records.jsonl")) as fh:
            return [json.loads(x) for x in fh]

    def rewrite(self, d, records, reseal: bool):
        """Write records back; with reseal, also recompute every id and hash a forger could."""
        from icms.canonical import record_id, sha256_file
        if reseal:
            for r in records:
                r["record_id"] = record_id({k: v for k, v in r.items() if k != "record_id"})
        with open(os.path.join(d, "records.jsonl"), "w") as fh:
            for r in records:
                fh.write(json.dumps(r, sort_keys=True) + "\n")
        if reseal:
            sp = os.path.join(d, "session.json")
            s = json.loads(_read(sp))
            s["records_sha256"] = sha256_file(os.path.join(d, "records.jsonl"))
            _dump(s, sp)

    def test_session_records_and_comparison_audit_clean(self):
        self.assertEqual(audit_session(self.session_dir, vals=self.vals), [])
        recs = self.records(self.session_dir)
        self.assertEqual(len(recs), 3 + 3 * 5)
        self.assertEqual({r["outcome"]["status"] for r in recs}, {"complete"})
        fb = recs[0]["metrics"]["factor_base"]
        self.assertEqual(fb["usable_points"], 16)
        self.assertIsNone(fb["abscissae_allowed"], "an unreported size field is null, never absent or zero")
        self.assertEqual(fb["native"], {"producer_specific": fb["native"]["producer_specific"]})
        self.assertEqual(recs[0]["metrics"]["solver"]["conflicts"], 5, "search effort moves from system to solver")
        self.assertEqual(recs[0]["metrics"]["system"]["cnf_clauses"], 24)
        for r in recs:
            self.assertEqual(r["execution"]["pinned_cpus"], [min(self.cpus)])
            self.assertEqual(r["execution"]["child_affinity_observed"], [min(self.cpus)])
            self.assertEqual(r["execution"]["env"]["defaults"]["RAYON_NUM_THREADS"], "1")

    def test_comparison_is_exact_and_carries_the_aa_arm(self):
        res = compare(self.session_dir, 0, 1, self.declared)
        self.assertEqual(res["refusals"], [])
        self.assertTrue(res["ops"]["deterministic"])
        self.assertEqual((res["ops"]["a_median"], res["ops"]["b_median"]), (8005.0, 9005.0))
        self.assertEqual(res["aa_noise"][0]["arms"], [0, 2])
        self.assertEqual(res["structure"]["a"]["system.cnf_clauses"]["median"], 24.0)
        refused = compare(self.session_dir, 0, 1, ["factor_base.params.recipe"])
        self.assertIn("execution.args.argv", refused["refusals"][0]["fields"], "an undeclared argv change is a confound")

    def test_existing_session_is_never_overwritten(self):
        with self.assertRaises(SessionError):
            run_session([os.path.join(self.tmp.name, "toy8.json")], self.session_dir, self.cpus,
                        lock=os.path.join(self.tmp.name, "lock"), allow_busy=True, log=lambda *_: None)

    def test_edited_number_breaks_the_record_hash(self):
        d = self.forged()
        recs = self.records(d)
        recs[4]["units"]["count.group_additions"]["total"] = 1.0
        self.rewrite(d, recs, reseal=False)
        problems = " ".join(audit_session(d, vals=self.vals))
        self.assertIn("record_id", problems)
        self.assertIn("records_sha256", problems)

    def test_resealed_forged_isolation_level_is_caught_by_the_gate(self):
        d = self.forged()
        recs = self.records(d)
        for r in recs:
            r["isolation"]["earned_level"] = "L3"
            r["isolation"]["blocking_next_level"] = []
            r["isolation"]["wall_admissible"] = True
        self.rewrite(d, recs, reseal=True)
        problems = audit_session(d, vals=self.vals)
        self.assertTrue(problems)
        self.assertTrue(all("isolation differs from the gate" in p or "comparisons/" in p for p in problems), problems)

    def test_resealed_forged_observation_is_caught_by_the_raw_output_hash(self):
        d = self.forged()
        with open(os.path.join(d, "exec", "0003", "stdout"), "a") as fh:
            fh.write("tampered\n")
        self.assertIn("does not hash to the record", " ".join(audit_session(d, vals=self.vals)))

    def test_edited_comparison_is_refused(self):
        d = self.forged()
        cp = os.path.join(d, "comparisons", "toy9-vs-toy8.json")
        c = json.loads(_read(cp))
        c["ops"]["ratio_b_over_a"] = 0.5
        _dump(c, cp)
        self.assertIn("recomputing from the records", " ".join(audit_session(d, vals=self.vals)))

    def test_edited_spec_copy_is_refused(self):
        d = self.forged()
        with open(os.path.join(d, "specs", "toy8.json"), "a") as fh:
            fh.write(" ")
        self.assertIn("does not hash to session.spec_files", " ".join(audit_session(d, vals=self.vals)))

    def test_record_schema_fields_follow_the_registry(self):
        reg = Registry.load()
        schema = json.loads(_read(os.path.join(REPO, "docs", "ic", "measurement", "schema", "record.v1.json")))
        m = schema["$defs"]["metrics"]["properties"]
        self.assertEqual(set(reg.doc["factor_base_size_fields"]["fields"]),
                         set(m["factor_base"]["required"]) - {"family"})
        solver_side = {"conflicts", "decisions", "propagations", "restarts", "matrix_rows_max", "matrix_cols_max",
                       "solving_degree"}
        self.assertEqual(set(reg.doc["pdp_metrics"]["fields"]) - solver_side, set(m["system"]["required"]))
        self.assertTrue(solver_side - {"solving_degree"} <= set(m["solver"]["properties"]))
        checks = schema["$defs"]["check"]["properties"]["id"]["enum"]
        gate = evaluate({}, {}, {"topology": {}}, 1)
        self.assertEqual([c["id"] for c in gate["checks"]], checks)


class CliTests(unittest.TestCase):
    def test_validate_cli(self):
        files = [os.path.join(SPECS, f) for f in sorted(os.listdir(SPECS)) if f.endswith(".yaml")]
        out = subprocess.run([sys.executable, os.path.join(TOOLS, "icms"), "validate", *files],
                             capture_output=True, text=True, cwd=REPO)
        self.assertEqual(out.returncode, 0, out.stdout + out.stderr)


if __name__ == "__main__":
    unittest.main()
