#include "matched_fixture_helpers.hpp"

int main(int argc, char** argv) {
  try {
    require(argc == 5,
            "usage: generate_targets ARCHIVE PRIOR_OPTIMIZATION_FIXTURES PRIOR_SPEED_FIXTURES NEW_OUTPUT_DIR");
    const files::path archive = argv[1], prior_dir = argv[2], speed_dir = argv[3],
                      out = argv[4];
    require(!files::exists(out), "Refusing existing output directory");
    require(sha256("") ==
                "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855" &&
                sha256("abc") ==
                    "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
            "SHA256 known-answer failure");
    auto instance = verify_instance(archive);
    verify_and_exclude_both_prior_panels(prior_dir, speed_dir, instance);
    require(files::create_directory(out), "Could not exclusively create output directory");
    const uint64_t bound = instance.r - 1;
    const uint64_t threshold = (uint64_t(0) - bound) % bound;
    SplitMix64 rng{MATCHED_SEED};
    auto excluded = instance.excluded_orbits;
    uint64_t range_rejections = 0, orbit_rejections = 0;
    std::vector<P> points;
    std::vector<uint64_t> scalars, draws;
    std::vector<U> keys;
    std::ostringstream trace;
    while (points.size() < MATCHED_TARGET_COUNT) {
      const uint64_t raw = rng.next();
      if (rng.draws > 1) trace << ",\n";
      trace << "    {\"draw_number\":" << rng.draws
            << ",\"raw_u64_validation_only\":" << qstr(std::to_string(raw));
      if (raw < threshold) {
        ++range_rejections;
        trace << ",\"decision\":\"REJECT_RANGE\",\"scalar_validation_only\":null,"
                 "\"orbit_key_decimal\":null,\"fixture_id\":null}";
        continue;
      }
      const uint64_t d = raw % bound + 1;
      const P q = mul(d, instance.g);
      require(!q.inf && on(q) && mul(instance.r, q).inf && mul(d, instance.g) == q,
              "Generated point/scalar/subgroup validity failure");
      const U key = orbit_key(q);
      trace << ",\"scalar_validation_only\":" << qstr(std::to_string(d))
            << ",\"orbit_key_decimal\":" << qstr(dec(key));
      if (excluded.count(key)) {
        ++orbit_rejections;
        trace << ",\"decision\":\"REJECT_ORBIT\",\"fixture_id\":null}";
        continue;
      }
      trace << ",\"decision\":\"ACCEPT\",\"fixture_id\":"
            << qstr(matched_id(points.size())) << "}";
      excluded.insert(key);
      points.push_back(q);
      scalars.push_back(d);
      draws.push_back(rng.draws);
      keys.push_back(key);
    }
    std::vector<files::path> artifacts;
    std::string combined;
    std::ostringstream target_entries, side_entries;
    for (size_t index = 0; index < MATCHED_TARGET_COUNT; ++index) {
      const auto id = matched_id(index);
      const auto filename = id + ".jsonl";
      const auto text = point_json(points[index]) + "\n";
      create_new(out / filename, text);
      artifacts.push_back(out / filename);
      combined += text;
      if (index) {
        target_entries << ",\n";
        side_entries << ",\n";
      }
      side_entries << "    {\"fixture_id\":" << qstr(id)
                   << ",\"split\":\"holdout\",\"point\":" << point_json(points[index])
                   << ",\"scalar_validation_only\":" << qstr(std::to_string(scalars[index]))
                   << "}";
      target_entries << "    {\"fixture_id\":" << qstr(id)
                     << ",\"split\":\"holdout\",\"point_file\":" << qstr(filename)
                     << ",\"point\":" << point_json(points[index])
                     << ",\"point_file_sha256\":" << qstr(sha256(text))
                     << ",\"signed_frobenius_orbit_key_decimal\":" << qstr(dec(keys[index]))
                     << ",\"SplitMix64_draw_number\":" << draws[index]
                     << ",\"nonidentity_on_curve\":true,\"r_times_Q_is_infinity\":true,"
                        "\"native_scalar_point_replay\":true,"
                        "\"disjoint_orbit_from_all_prior_panels_base_and_accepted_fixtures\":true}";
    }
    const std::string side =
        "{\n  \"kind\":\"validation_only_scalar_sidecar\",\n"
        "  \"must_not_be_passed_to_measured_children\":true,\n  \"entries\":[\n" +
        side_entries.str() + "\n  ],\n  \"proposal_trace\":[\n" + trace.str() + "\n  ]\n}\n";
    create_new(out / "validation_sidecar.json", side);
    artifacts.push_back(out / "validation_sidecar.json");
    create_new(out / "holdout_32_targets.jsonl", combined);
    artifacts.push_back(out / "holdout_32_targets.jsonl");
    std::ostringstream receipt;
    receipt << "{\n  \"kind\":\"native_deterministic_n83_matched_target_generation\",\n"
            << "  \"status\":\"PASS\",\n  \"seed\":" << MATCHED_SEED << ",\n"
            << "  \"random_generator\":\"SplitMix64 with full validation-only draw trace\",\n"
            << "  \"sampling_law\":\"Unbiased rejection of 64-bit nonzero subgroup scalar proposals; accepted points conditioned on exclusion of K600, original target, all16 prior optimization targets and earlier accepted signed Frobenius orbits\",\n"
            << "  \"scalar_bound_r_minus_one\":" << bound << ",\n"
            << "  \"uniform_rejection_threshold_2pow64_mod_bound\":" << threshold << ",\n"
            << "  \"SplitMix64_draws\":" << rng.draws << ",\n"
            << "  \"uniform_range_rejections\":" << range_rejections << ",\n"
            << "  \"excluded_orbit_rejections\":" << orbit_rejections << ",\n"
            << "  \"target_count\":32,\n  \"training_target_count\":0,\n  \"holdout_target_count\":32,\n"
            << "  \"one_target_point_files\":true,\n  \"scalar_sidecar_is_validation_only\":true,\n"
            << "  \"curve_slug\":\"icv1-f2m83-t6151469093347-9eba6cd5\",\n"
            << "  \"n\":83,\n  \"a\":1,\n  \"b\":1,\n"
            << "  \"field_modulus_low_terms\":[0,2,4,7],\n"
            << "  \"subgroup_order\":\"8569786107849059\",\n"
            << "  \"generator\":" << point_json(instance.g) << ",\n"
            << "  \"prior_public_Q_excluded\":" << point_json(instance.prior) << ",\n"
            << "  \"tau_eigenvalue\":\"2998004823950457\",\n"
            << "  \"signed_frobenius_orbit_size\":166,\n"
            << "  \"K600_representatives_verified\":600,\n"
            << "  \"K600_signed_points_and_labels_verified\":99600,\n"
            << "  \"field_irreducibility_and_subgroup_primality_verified\":true,\n"
            << "  \"initial_excluded_signed_orbits\":617,\n"
            << "  \"prior_optimization_targets_excluded\":8,\n"
            << "  \"prior_speed_targets_excluded\":8,\n"
            << "  \"prior_panel_generation_reproduced\":true,\n"
            << "  \"source_revision_of_retained_inputs\":\"a8194d616ed8c04489ec46177d369f1b1670360a\",\n"
            << "  \"archived_fixture_sha256\":" << qstr(digest(archive / "frozen/fixture.json")) << ",\n"
            << "  \"retained_base_file_sha256\":" << qstr(digest(archive / "base_n83_K600.jsonl")) << ",\n"
            << "  \"retained_base_label_hash\":\"d4a4da26d6e7cb09367de922d9806631fe49989f70d378cae76f347f8d76a9f7\",\n"
            << "  \"prior_optimization_generator_receipt_sha256\":"
            << qstr(digest(prior_dir / "generator_receipt.json")) << ",\n"
            << "  \"prior_speed_generator_receipt_sha256\":"
            << qstr(digest(speed_dir / "generator_receipt.json")) << ",\n"
            << "  \"generator_source_sha256\":" << qstr(digest(__FILE__)) << ",\n"
            << "  \"generator_executable_sha256\":" << qstr(digest(argv[0])) << ",\n"
            << "  \"validation_sidecar_file\":\"validation_sidecar.json\",\n"
            << "  \"validation_sidecar_sha256\":" << qstr(sha256(side)) << ",\n"
            << "  \"holdout_combined_file\":\"holdout_32_targets.jsonl\",\n"
            << "  \"holdout_combined_sha256\":" << qstr(sha256(combined)) << ",\n"
            << "  \"fresh_exclusion_proof\":\"Minimum x among all83 Frobenius images identifies a signed orbit: sign leaves x unchanged, and a fixed nonzero x on this binary curve has exactly P and -P\",\n"
            << "  \"targets\":[\n" << target_entries.str() << "\n  ]\n}\n";
    create_new(out / "generator_receipt.json", receipt.str());
    artifacts.push_back(out / "generator_receipt.json");
    std::ostringstream hashes;
    for (const auto& path : artifacts)
      hashes << digest(path) << "  " << path.filename().string() << "\n";
    create_new(out / "artifact_hashes.sha256", hashes.str());
    std::cout << "PASS: generated32 fresh holdout point-only fixtures; initial exclusions=617; draws="
              << rng.draws << "; range_rejects=" << range_rejections
              << "; orbit_rejects=" << orbit_rejections << "\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "FAIL: " << error.what() << "\n";
    return 1;
  }
}
