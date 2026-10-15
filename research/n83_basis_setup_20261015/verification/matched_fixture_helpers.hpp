#pragma once
#include "fresh_fixture_helpers.hpp"

constexpr uint64_t MATCHED_SEED = 202610146301ULL;
constexpr size_t MATCHED_TARGET_COUNT = 32;
constexpr size_t MATCHED_INITIAL_ORBITS = 617;

inline std::string matched_id(size_t index) {
  require(index < MATCHED_TARGET_COUNT, "Matched fixture index out of bounds");
  std::ostringstream out;
  out << "target_" << std::setfill('0') << std::setw(2) << index + 1;
  return out.str();
}

inline void verify_and_exclude_both_prior_panels(
    const files::path& prior_dir, const files::path& speed_dir,
    Instance& instance) {
  require(instance.excluded_orbits.size() == 601,
          "Expected K600 and original public target exclusions");
  const auto prior = verify_prior_fixtures(prior_dir, instance);
  require(prior.size() == 8 && instance.excluded_orbits.size() == 609,
          "Prior optimization exclusion count mismatch");
  const auto speed = verify_fresh_fixtures(speed_dir, prior_dir, instance);
  require(speed.size() == 8, "Speed panel fixture count mismatch");
  for (const auto& fixture : speed)
    require(instance.excluded_orbits.insert(orbit_key(fixture.point)).second,
            "Speed fixture overlaps previous exclusions");
  require(instance.excluded_orbits.size() == MATCHED_INITIAL_ORBITS,
          "Expected 617 complete prior signed-orbit exclusions");
}

inline std::vector<Fixture> verify_matched_fixtures(
    const files::path& dir, const files::path& prior_dir,
    const files::path& speed_dir, const Instance& instance) {
  const auto receipt = read_json((dir / "generator_receipt.json").string());
  const auto sidecar = read_json((dir / "validation_sidecar.json").string());
  require(receipt.at("status").scalar == "PASS" &&
              num(receipt.at("seed")) == MATCHED_SEED,
          "Matched generation status/seed mismatch");
  check_modulus(receipt);
  require(num(receipt.at("n")) == 83 && num(receipt.at("a")) == 1 &&
              num(receipt.at("b")) == 1 &&
              num(receipt.at("subgroup_order")) == instance.r &&
              point(receipt.at("generator")) == instance.g,
          "Matched curve/G/r mismatch");
  require(instance.excluded_orbits.size() == MATCHED_INITIAL_ORBITS &&
              num(receipt.at("initial_excluded_signed_orbits")) ==
                  MATCHED_INITIAL_ORBITS &&
              num(receipt.at("target_count")) == MATCHED_TARGET_COUNT &&
              num(receipt.at("holdout_target_count")) == MATCHED_TARGET_COUNT &&
              num(receipt.at("training_target_count")) == 0,
          "Matched fixture/exclusion counts mismatch");
  require(digest(prior_dir / "generator_receipt.json") ==
              receipt.at("prior_optimization_generator_receipt_sha256").scalar &&
              digest(speed_dir / "generator_receipt.json") ==
                  receipt.at("prior_speed_generator_receipt_sha256").scalar,
          "Prior generator receipt binding mismatch");
  require(digest(dir / "validation_sidecar.json") ==
              receipt.at("validation_sidecar_sha256").scalar &&
              sidecar.at("kind").scalar == "validation_only_scalar_sidecar" &&
              truth(sidecar.at("must_not_be_passed_to_measured_children")),
          "Matched validation-only sidecar binding mismatch");
  const auto& targets = receipt.at("targets").array;
  const auto& validation = sidecar.at("entries").array;
  const auto& trace = sidecar.at("proposal_trace").array;
  require(targets.size() == MATCHED_TARGET_COUNT &&
              validation.size() == MATCHED_TARGET_COUNT,
          "Matched target/validation entry count mismatch");
  const uint64_t bound = instance.r - 1;
  const uint64_t threshold = (uint64_t(0) - bound) % bound;
  require(num(receipt.at("scalar_bound_r_minus_one")) == bound &&
              num(receipt.at("uniform_rejection_threshold_2pow64_mod_bound")) ==
                  threshold,
          "Matched unbiased sampling parameters mismatch");
  auto excluded = instance.excluded_orbits;
  SplitMix64 rng{MATCHED_SEED};
  uint64_t range_rejections = 0, orbit_rejections = 0;
  size_t trace_index = 0;
  std::string combined;
  std::vector<Fixture> result;
  for (size_t index = 0; index < MATCHED_TARGET_COUNT; ++index) {
    const auto& target = targets[index];
    const auto& witness = validation[index];
    const std::string id = matched_id(index);
    require(target.at("fixture_id").scalar == id &&
                target.at("split").scalar == "holdout" &&
                witness.at("fixture_id").scalar == id &&
                witness.at("split").scalar == "holdout" &&
                target.at("point_file").scalar == id + ".jsonl",
            "Matched fixture IDs/order/split/path mismatch");
    const P q = point(target.at("point"));
    require(!q.inf && on(q) && mul(instance.r, q).inf,
            "Matched curve/subgroup validity failure");
    const U key = orbit_key(q);
    require(!excluded.count(key) &&
                dec(key) == target.at("signed_frobenius_orbit_key_decimal").scalar,
            "Matched fixture overlaps excluded orbit or recorded key mismatch");
    const uint64_t d = uint64_t(num(witness.at("scalar_validation_only")));
    require(d > 0 && d < instance.r && point(witness.at("point")) == q &&
                mul(d, instance.g) == q,
            "Matched validation-only scalar replay failure");
    while (true) {
      const uint64_t raw = rng.next();
      require(trace_index < trace.size(), "Missing proposal trace row");
      const auto& event = trace[trace_index++];
      require(num(event.at("draw_number")) == rng.draws &&
                  num(event.at("raw_u64_validation_only")) == raw,
              "Matched raw proposal trace mismatch");
      if (raw < threshold) {
        ++range_rejections;
        require(event.at("decision").scalar == "REJECT_RANGE" &&
                    event.at("scalar_validation_only").scalar == "null" &&
                    event.at("orbit_key_decimal").scalar == "null" &&
                    event.at("fixture_id").scalar == "null",
                "Matched range-rejection trace mismatch");
        continue;
      }
      const uint64_t proposal = raw % bound + 1;
      const P candidate = mul(proposal, instance.g);
      const U proposed_key = orbit_key(candidate);
      require(num(event.at("scalar_validation_only")) == proposal &&
                  num(event.at("orbit_key_decimal")) == proposed_key,
              "Matched scalar/orbit proposal trace mismatch");
      if (excluded.count(proposed_key)) {
        ++orbit_rejections;
        require(event.at("decision").scalar == "REJECT_ORBIT" &&
                    event.at("fixture_id").scalar == "null",
                "Matched orbit-rejection trace mismatch");
        continue;
      }
      require(event.at("decision").scalar == "ACCEPT" &&
                  event.at("fixture_id").scalar == id && proposal == d &&
                  candidate == q &&
                  rng.draws == num(target.at("SplitMix64_draw_number")),
              "Matched accepted SplitMix/rejection stream mismatch");
      break;
    }
    excluded.insert(key);
    const auto file = dir / (id + ".jsonl");
    const auto hash = digest(file);
    const auto rows = lines(file);
    require(hash == target.at("point_file_sha256").scalar && rows.size() == 1 &&
                rows[0].type == "array" && rows[0].array.size() == 2 &&
                rows[0].at(0).type == "string" &&
                rows[0].at(1).type == "string" && point(rows[0]) == q,
            "Matched point-only file bytes/schema mismatch");
    combined += bytes(file);
    result.push_back({id, "holdout", q, d, hash});
  }
  require(trace_index == trace.size() && trace_index == rng.draws &&
              rng.draws == num(receipt.at("SplitMix64_draws")) &&
              range_rejections == num(receipt.at("uniform_range_rejections")) &&
              orbit_rejections == num(receipt.at("excluded_orbit_rejections")),
          "Matched complete stream/rejection counts mismatch");
  require(receipt.at("holdout_combined_file").scalar == "holdout_32_targets.jsonl" &&
              bytes(dir / "holdout_32_targets.jsonl") == combined &&
              digest(dir / "holdout_32_targets.jsonl") ==
                  receipt.at("holdout_combined_sha256").scalar,
          "Matched combined file order/bytes/hash mismatch");
  return result;
}
