#include "matched_fixture_helpers.hpp"

int main(int argc, char** argv) {
  try {
    require(argc == 6,
            "usage: check_fixtures ARCHIVE PRIOR_OPTIMIZATION_FIXTURES PRIOR_SPEED_FIXTURES MATCHED_FIXTURES NEW_RECEIPT");
    const files::path archive = argv[1], prior_dir = argv[2], speed_dir = argv[3],
                      fixtures_dir = argv[4], out = argv[5];
    require(!files::exists(out), "Refusing existing verification receipt");
    auto instance = verify_instance(archive);
    verify_and_exclude_both_prior_panels(prior_dir, speed_dir, instance);
    const auto fixtures = verify_matched_fixtures(fixtures_dir, prior_dir, speed_dir, instance);
    std::ostringstream receipt;
    receipt << "{\n  \"kind\":\"native_n83_matched_fixture_verification\",\n"
            << "  \"status\":\"PASS_FIXTURES_ONLY\",\n"
            << "  \"curve_slug\":\"icv1-f2m83-t6151469093347-9eba6cd5\",\n"
            << "  \"field_modulus_low_terms\":[0,2,4,7],\n"
            << "  \"fresh_generation_seed\":" << MATCHED_SEED << ",\n"
            << "  \"holdout_fixtures_verified\":" << fixtures.size() << ",\n"
            << "  \"old_optimization_fixtures_reproduced_and_excluded\":8,\n"
            << "  \"old_speed_fixtures_reproduced_and_excluded\":8,\n"
            << "  \"initial_excluded_signed_orbits\":617,\n"
            << "  \"retained_base_points_and_labels_verified\":99600,\n"
            << "  \"all_fixture_curve_subgroup_scalar_and_exact_draw_stream_checks\":true,\n"
            << "  \"all_new_signed_orbits_disjoint_from_original_Q_old16_targets_K600_and_one_another\":true,\n"
            << "  \"fixture_manifest_sha256\":" << qstr(digest(fixtures_dir / "generator_receipt.json")) << ",\n"
            << "  \"prior_fixture_manifest_sha256\":" << qstr(digest(prior_dir / "generator_receipt.json")) << ",\n"
            << "  \"prior_speed_fixture_manifest_sha256\":" << qstr(digest(speed_dir / "generator_receipt.json")) << ",\n"
            << "  \"validation_sidecar_sha256\":" << qstr(digest(fixtures_dir / "validation_sidecar.json")) << ",\n"
            << "  \"checker_source_sha256\":" << qstr(digest(__FILE__)) << ",\n"
            << "  \"checker_executable_sha256\":" << qstr(digest(argv[0])) << ",\n"
            << "  \"performance_claim\":false\n}\n";
    create_new(out, receipt.str());
    std::cout << "PASS: independently regenerated32 matched fixtures and all16 prior targets; verified99,600 K600 points/labels and617 initial signed-orbit exclusions\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "FAIL: " << error.what() << "\n";
    return 1;
  }
}
