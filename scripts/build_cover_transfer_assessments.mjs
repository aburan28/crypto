#!/usr/bin/env node

import { createHash } from "node:crypto";
import { readFile, writeFile } from "node:fs/promises";
import { fileURLToPath } from "node:url";
import { dirname, resolve } from "node:path";

const SCRIPT_DIR = dirname(fileURLToPath(import.meta.url));
const ROOT = resolve(SCRIPT_DIR, "..");
const OUTPUT_PATH = "docs/curves/cover-transfer-assessments.json";
const ASSESSMENT_DATE = "2026-10-13";
const CATALOG_COMMIT = "6600c7c6c0cfbb6a9e4e596bf8c5cfa292fb269e";
const MAP_HASH_DOMAIN = "cover-transfer-assessments/v2/map_base";
const SUBGROUP_HASH_DOMAIN = "cover-transfer-assessments/v2/subgroup_base";

function sha256Bytes(value) {
  return createHash("sha256").update(value).digest("hex");
}

function canonicalJson(value) {
  if (Array.isArray(value)) {
    return `[${value.map(canonicalJson).join(",")}]`;
  }
  if (value !== null && typeof value === "object") {
    return `{${Object.keys(value)
      .sort(compareText)
      .map((key) => `${JSON.stringify(key)}:${canonicalJson(value[key])}`)
      .join(",")}}`;
  }
  return JSON.stringify(value);
}

function membershipHash(domain, members) {
  return sha256Bytes(`${domain}\0${JSON.stringify(members)}`);
}

async function readJson(path) {
  return JSON.parse(await readFile(resolve(ROOT, path), "utf8"));
}

async function fileHash(path) {
  return sha256Bytes(await readFile(resolve(ROOT, path)));
}

function assert(condition, message) {
  if (!condition) {
    throw new Error(message);
  }
}

function compareText(left, right) {
  if (left < right) return -1;
  if (left > right) return 1;
  return 0;
}

function sortedUnique(values, label) {
  const sorted = [...values].sort(compareText);
  assert(new Set(sorted).size === sorted.length, `${label} contains duplicates`);
  return sorted;
}

function sortedUniqueTriples(values, label) {
  const sorted = [...values].sort((left, right) =>
    compareText(JSON.stringify(left), JSON.stringify(right)),
  );
  assert(
    new Set(sorted.map((value) => JSON.stringify(value))).size === sorted.length,
    `${label} contains duplicates`,
  );
  return sorted;
}

function sortedRecord(entries) {
  return Object.fromEntries([...entries].sort(([left], [right]) => compareText(left, right)));
}

function targetDeclaredFieldRecord(model) {
  if (model.field.startsWith("f2m-")) {
    return { label: model.field, modulus: model.modulus };
  }
  if (model.field.startsWith("fpk-")) {
    return {
      label: model.field,
      p: String(model.p),
      modulus: model.modulus.map(String),
    };
  }
  assert(model.field.startsWith("fp-"), `unsupported target field label ${model.field}`);
  return { label: model.field, p: String(model.p) };
}

function declaredFieldBinding(record) {
  const fieldUid = `urn:declared-field:1:sha256:${sha256Bytes(canonicalJson(record))}`;
  const binary = record.label.match(/^f2m-([1-9][0-9]*)-/u);
  if (binary !== null) {
    return {
      field_uid: fieldUid,
      record,
      characteristic: "2",
      extension_degree: Number(binary[1]),
      representation: "binary_polynomial_extension",
    };
  }
  const extension = record.label.match(/^fpk-([1-9][0-9]*)-([1-9][0-9]*)-/u);
  if (extension !== null) {
    assert(String(record.p) === extension[1], `${record.label}: characteristic mismatch`);
    assert(record.modulus.length === Number(extension[2]), `${record.label}: modulus degree mismatch`);
    return {
      field_uid: fieldUid,
      record,
      characteristic: String(record.p),
      extension_degree: Number(extension[2]),
      representation: "odd_polynomial_extension",
    };
  }
  assert(record.label === `fp-${record.p}`, `${record.label}: prime-field label mismatch`);
  return {
    field_uid: fieldUid,
    record,
    characteristic: String(record.p),
    extension_degree: 1,
    representation: "prime",
  };
}

const obligations = [
  {
    id: "general_binary_norm",
    description:
      "Implement and independently verify pushforward/norm on arbitrary binary Jacobian classes.",
    status: "unknown",
  },
  {
    id: "independent_map_oracle",
    description:
      "Replay pullback additivity and composition with an independent CAS or separately implemented divisor-class oracle.",
    status: "unknown",
  },
  {
    id: "prime_balanced_jacobian",
    description:
      "Implement balanced two-infinity odd-characteristic Jacobian arithmetic and induced maps for the degree-two family.",
    status: "unknown",
  },
  {
    id: "catalog_exact_subgroup_order",
    description:
      "Attach proof-grade exact-order evidence (including primality/factorization evidence where required) to every numeric EC1 subgroup-order claim.",
    status: "unknown",
  },
  {
    id: "catalog_generator_validation",
    description:
      "Replay generator membership and catalog-order annihilation for every linked odd-characteristic EC1 representation.",
    status: "unknown",
  },
  {
    id: "catalog_wide_packaged_binary",
    description:
      "Build and publish a standalone artifact that accepts every currently supported catalog input form, including odd-characteristic extension fields and normalized target models, then replay all 2,091 direct commands against that exact artifact.",
    status: "unknown",
  },
  {
    id: "coordinate_controls",
    description:
      "Add coordinate-isomorphism commutation controls for maps and divisor representations.",
    status: "unknown",
  },
  {
    id: "extension_descent",
    description:
      "Assess descent from each declared odd-characteristic extension field; same-declared-field cover geometry is not a descent certificate.",
    status: "unknown",
  },
  {
    id: "refresh_formula_inventory",
    description:
      "Regenerate jacobian-transfers.json against the current registry and preserve modeled formulas separately from transfer certificates.",
    status: "unknown",
  },
  {
    id: "full_solver",
    description:
      "If pursued separately, require genuine Jacobian relations, verified linear algebra and target decomposition, independent output checks, and complete preprocessing/failure/duplicate/memory/runtime accounting.",
    status: "unknown",
  },
];

const resultTemplates = {
  "map.geometry.catalog-supported": {
    status: "supported",
    scope:
      "Every exact map UID in the selected profile, over its declared field and subject to the catalog's explicit field-validity boundary.",
    rationale:
      "The generated HC1/CV1 row records an explicit source, target, formula, degree, and verified-over-declared-field geometry certificate.",
    evidence_refs: [
      "docs/curves/covers.json",
      "docs/curves/cover-links.yaml",
      "src/bin/curve_cover_check/checker.rs",
      "docs/curves/COVERS.md",
    ],
  },
  "map.binary-construction.supported": {
    status: "supported",
    scope:
      "Ordinary binary Weierstrass targets in the 160-row cubic-pullback profile, with exact Rabin-validated polynomial-basis fields and nonzero b.",
    rationale:
      "The uniform constructor emits H: v^2+u^2v=u^7+a*u^4+sqrt(b)*u and phi(u,v)=(u^3,u*v+sqrt(b)); all profile members match those coefficients, and the current source-level standalone constructor accepts every member of this profile.",
    evidence_refs: [
      "src/binary_ecc/cover_transfer.rs",
      "src/bin/curve_cover_check/checker.rs",
      "src/bin/hyperelliptic_cover.rs#every_current_catalog_model_reaches_the_standalone_constructor",
      "docs/curves/covers.json",
    ],
  },
  "map.binary-jacobian.supported": {
    status: "supported",
    scope:
      "Canonical reduced Mumford divisors for the one-rational-infinity binary source family.",
    rationale:
      "The checked implementation provides validity, reduction, equality, negation, addition, and scalar multiplication for this model boundary.",
    evidence_refs: [
      "src/binary_ecc/hyperelliptic.rs",
      "research/hyperelliptic_cover_infrastructure_20261006/evidence.json",
    ],
  },
  "map.binary-pullback.supported": {
    status: "supported",
    scope:
      "Rational elliptic point classes [P]-[O], including infinity and the ramified exceptional fiber, for the binary cubic family.",
    rationale:
      "The typed API constructs checked reduced Mumford representatives with an explicit generic and exceptional branch.",
    evidence_refs: ["src/binary_ecc/cover_transfer.rs#pullback_point_class"],
  },
  "map.binary-typed-pushforward.supported": {
    status: "supported",
    scope:
      "Rational source points and the implementation's own typed pullback witnesses; not arbitrary Jacobian classes.",
    rationale:
      "The implementation maps source points and evaluates its ramification-aware formal fiber witness on this bounded domain.",
    evidence_refs: ["src/binary_ecc/cover_transfer.rs"],
  },
  "map.binary-composition.supported": {
    status: "supported",
    scope:
      "phi_* phi^*=[3] on the implementation's typed rational point-class domain for the binary cubic family.",
    rationale:
      "The native composition certificate includes identity, exceptional, and generic fibers; the catalog-wide subgroup replay exercises every linked binary generator.",
    evidence_refs: [
      "src/binary_ecc/cover_transfer.rs#pull_push_composition_certificate",
      "src/bin/hyperelliptic_cover.rs#every_binary_ec1_generator_has_checked_typed_transfer_hypotheses",
    ],
  },
  "map.independent-oracle.unknown": {
    status: "unknown",
    scope:
      "Independent pullback additivity and composition replay for every exact catalog member.",
    rationale:
      "The current checks are native and implementation-internal; no independent CAS or separately implemented divisor-map oracle covers the full catalog.",
    evidence_refs: [],
    obligation_id: "independent_map_oracle",
  },
  "map.binary-arbitrary-pushforward.refuted": {
    status: "refuted",
    scope:
      "Claim that the current implementation accepts and norms an arbitrary reduced class in Jac(H).",
    rationale:
      "The public entry point deliberately returns Unsupported(ArbitraryJacobianPushforward); this refutes an implementation claim, not mathematical existence of a norm.",
    evidence_refs: [
      "src/binary_ecc/cover_transfer.rs#pushforward_arbitrary_jacobian_class",
    ],
    obligation_id: "general_binary_norm",
  },
  "map.prime-construction.supported": {
    status: "supported",
    scope:
      "Odd-characteristic quadratic-pullback rows over their declared prime field, with base-field primality retained as an explicit input assumption.",
    rationale:
      "The constructor chooses the first c in {0,1,2,3} with F(c) nonzero, then emits the monic squarefree sextic F(u^2+c), degree-two formula (u,v)->(u^2+c,v), and two rational infinity points; the current source-level standalone constructor accepts every member of this profile.",
    evidence_refs: [
      "src/bin/curve_cover_check/checker.rs",
      "src/bin/hyperelliptic_cover.rs#every_current_catalog_model_reaches_the_standalone_constructor",
      "docs/curves/covers.json",
      "docs/curves/COVERS.md",
    ],
  },
  "map.prime-normalization.supported": {
    status: "supported",
    scope:
      "The 30 odd-prime rows whose catalog certificate includes an explicit degree-one target-model map.",
    rationale:
      "Each member binds the short-Weierstrass cover target to its catalog Montgomery or Edwards representation through a verified algebraic identity.",
    evidence_refs: [
      "docs/curves/covers.json",
      "src/bin/curve_cover_check/checker.rs",
    ],
  },
  "map.direct-target.not-applicable": {
    status: "not_applicable",
    scope: "Rows whose cover target is already the catalog's short-Weierstrass model.",
    rationale: "No auxiliary degree-one target normalization is required for these rows.",
    evidence_refs: [],
  },
  "map.extension-construction.supported": {
    status: "supported",
    scope:
      "Odd-characteristic quadratic-pullback rows over declared polynomial-basis extensions, conditional on the stated base-prime input assumption.",
    rationale:
      "The catalog proves the extension modulus irreducible by Rabin and verifies the cover and map over that same declared extension field; the current source-level standalone constructor accepts every member of this profile.",
    evidence_refs: [
      "src/bin/curve_cover_check/checker.rs",
      "src/bin/hyperelliptic_cover.rs#every_current_catalog_model_reaches_the_standalone_constructor",
      "docs/curves/covers.json",
      "docs/curves/COVERS.md",
    ],
  },
  "map.prime-jacobian.refuted": {
    status: "refuted",
    scope:
      "Claim that the current odd-prime one-infinity Jacobian implementation accepts the even sextic sources in the quadratic-pullback catalog.",
    rationale:
      "The implementation requires monic degree 2g+1 and one rational infinity; every catalog source in this family has degree 2g+2 and two rational infinities.",
    evidence_refs: [
      "src/prime_hyperelliptic/curve.rs",
      "docs/curves/HYPERELLIPTIC_INFRASTRUCTURE.md",
    ],
  },
  "map.prime-transfer.unknown": {
    status: "unknown",
    scope:
      "Pullback, pushforward, and phi_*phi^*=[2] on divisor classes for the two-infinity odd-characteristic source family.",
    rationale:
      "No compatible balanced Jacobian representation or induced-map implementation is present; cover geometry alone does not certify these operations.",
    evidence_refs: [],
    obligation_id: "prime_balanced_jacobian",
  },
  "map.same-field.not-applicable": {
    status: "not_applicable",
    scope:
      "Cover, map, and elliptic target are defined over the same declared field.",
    rationale:
      "The construction introduces no additional field extension, so extension-to-source descent is not a transfer step.",
    evidence_refs: [],
  },
  "map.extension-descent.unknown": {
    status: "unknown",
    scope:
      "Descent from the declared odd-characteristic extension field to a proper subfield.",
    rationale:
      "The cover is verified over the same declared extension as its target; no proper-subfield descent or subgroup compatibility certificate was attempted.",
    evidence_refs: [],
    obligation_id: "extension_descent",
  },
  "map.advantage.unknown": {
    status: "unknown",
    scope:
      "End-to-end discrete-log cost relative to a matched generic baseline.",
    rationale:
      "No complete Jacobian relation solver, calibrated transfer cost, target decomposition, or matched total-runtime experiment exists.",
    evidence_refs: [],
    obligation_id: "full_solver",
  },
  "map.binary-pipeline.refuted": {
    status: "refuted",
    scope:
      "Claim that the binary implementation provides a general Jac(H)-to-E norm and end-to-end logarithm transfer.",
    rationale:
      "Only point-divisor and typed-fiber pushforward domains exist; arbitrary Jacobian pushforward is an explicit negative control and no solver is present.",
    evidence_refs: [
      "src/binary_ecc/cover_transfer.rs#pushforward_arbitrary_jacobian_class",
      "research/hyperelliptic_cover_infrastructure_20261006/README.md",
    ],
  },
  "map.prime-pipeline.refuted": {
    status: "refuted",
    scope:
      "Claim that the current implementation provides arithmetic and induced maps for the odd-characteristic two-infinity source family.",
    rationale:
      "The sources are outside the one-infinity prime Jacobian implementation and no balanced substitute or induced maps are present.",
    evidence_refs: [
      "src/prime_hyperelliptic/curve.rs",
      "docs/curves/HYPERELLIPTIC_INFRASTRUCTURE.md",
    ],
  },
  "subgroup.identity.supported": {
    status: "supported",
    scope:
      "Every exact [EC1, curve UID, map UID] triple in the selected subgroup membership.",
    rationale:
      "The current content-addressed graph binds the EC1 record, full curve UID, generator, catalog subgroup order, cofactor, target model, and exact map UID.",
    evidence_refs: ["docs/curves/cover-links.yaml"],
  },
  "subgroup.degree-coprime.supported": {
    status: "supported",
    scope:
      "The integer gcd between map degree and the catalog subgroup-order value for every selected triple.",
    rationale:
      "Every degree-three binary value is nonzero modulo 3 and every degree-two odd-characteristic value is odd.",
    evidence_refs: [
      "docs/curves/cover-links.yaml",
      "scripts/build_cover_transfer_assessments.mjs",
    ],
  },
  "subgroup.binary-generator.supported": {
    status: "supported",
    scope:
      "All 134 linked binary EC1 generators, using their exact catalog fields and curve coefficients.",
    rationale:
      "The content-addressed graph and independent semantic verifier establish each exact EC1/curve/map join; the native batch certificate then decodes that row's generator canonically and checks its binary Weierstrass equation.",
    evidence_refs: [
      "docs/curves/cover-links.yaml",
      "research/cover_transfer_assessment_20261013/verify_assessment.mjs",
      "src/bin/hyperelliptic_cover.rs#every_binary_ec1_generator_has_checked_typed_transfer_hypotheses",
    ],
  },
  "subgroup.binary-annihilation.supported": {
    status: "supported",
    scope:
      "All 134 linked binary generators and their numeric catalog subgroup-order values r.",
    rationale:
      "The native batch certificate checks [r]G=O. This proves ord(G) divides r, not by itself that ord(G)=r.",
    evidence_refs: [
      "src/bin/hyperelliptic_cover.rs#every_binary_ec1_generator_has_checked_typed_transfer_hypotheses",
    ],
  },
  "subgroup.numeric-order.unknown": {
    status: "unknown",
    scope:
      "Claim that each numeric catalog value r is the exact order of its linked generator.",
    rationale:
      "Catalog identity, [r]G=O, and a probable-prime screen do not constitute proof-grade exact-order evidence for every value.",
    evidence_refs: [],
    obligation_id: "catalog_exact_subgroup_order",
  },
  "subgroup.binary-pullback.supported": {
    status: "supported",
    scope:
      "The checked pullback image of every linked binary EC1 generator in the typed point-class domain.",
    rationale:
      "The batch certificate exercises the implementation-selected generic or exceptional branch for each generator and validates every resulting Mumford divisor.",
    evidence_refs: [
      "src/bin/hyperelliptic_cover.rs#every_binary_ec1_generator_has_checked_typed_transfer_hypotheses",
      "src/binary_ecc/cover_transfer.rs",
    ],
  },
  "subgroup.binary-image-annihilation.supported": {
    status: "supported",
    scope:
      "Every checked binary pullback image and its numeric catalog order value r.",
    rationale:
      "The native batch certificate checks [r]phi^*(G)=0 for all 134 exact triples.",
    evidence_refs: [
      "src/bin/hyperelliptic_cover.rs#every_binary_ec1_generator_has_checked_typed_transfer_hypotheses",
    ],
  },
  "subgroup.binary-composition.supported": {
    status: "supported",
    scope:
      "phi_*phi^*(G)=[3]G for every linked binary generator, within the typed fiber witness domain.",
    rationale:
      "The native batch certificate compares both sides for all 134 exact triples.",
    evidence_refs: [
      "src/bin/hyperelliptic_cover.rs#every_binary_ec1_generator_has_checked_typed_transfer_hypotheses",
    ],
  },
  "subgroup.binary-injective-generated.unknown": {
    status: "unknown",
    scope:
      "The actual generated subgroup <G> for each of the 134 linked binary triples.",
    rationale:
      "The checked generator identity would imply injectivity after pullback additivity and phi_*phi^*=[3] are certified on every element of <G>; the current batch test does not certify those all-element hypotheses.",
    evidence_refs: [
      "src/bin/hyperelliptic_cover.rs#every_binary_ec1_generator_has_checked_typed_transfer_hypotheses",
      "research/cover_transfer_assessment_20261013/lean/CoverTransfer.lean#pull_injective_on_subgroup",
    ],
    obligation_id: "independent_map_oracle",
  },
  "subgroup.binary-image-actual-order.unknown": {
    status: "unknown",
    scope:
      "Equality between the actual order of each linked binary generator and the actual order of its checked pullback image.",
    rationale:
      "Order preservation follows from injectivity of a homomorphic pullback on <G>, but pullback additivity and the all-element composition hypothesis remain uncertified catalog-wide.",
    evidence_refs: [
      "src/bin/hyperelliptic_cover.rs#every_binary_ec1_generator_has_checked_typed_transfer_hypotheses",
      "research/cover_transfer_assessment_20261013/lean/CoverTransfer.lean#pull_injective_on_subgroup",
    ],
    obligation_id: "independent_map_oracle",
  },
  "subgroup.image-numeric-order.unknown": {
    status: "unknown",
    scope:
      "Claim that the pullback image has exact order equal to the numeric catalog value r.",
    rationale:
      "The batch test proves that r annihilates the image, but neither equality of its actual order with ord(G) nor equality with the numeric value r is certified catalog-wide.",
    evidence_refs: [],
    obligation_id: "catalog_exact_subgroup_order",
  },
  "subgroup.binary-inverse-recovery.supported": {
    status: "supported",
    scope:
      "Recovery of G from the typed pull-push value [3]G using 3^{-1} mod r for all 134 linked binary triples.",
    rationale:
      "The batch certificate checks the modular inverse and exact point equality for each catalog value r; this is bounded typed-domain recovery, not a general Jacobian norm.",
    evidence_refs: [
      "src/bin/hyperelliptic_cover.rs#every_binary_ec1_generator_has_checked_typed_transfer_hypotheses",
    ],
  },
  "subgroup.prime-generator.unknown": {
    status: "unknown",
    scope:
      "Generator membership and [r]G=O replay for all 530 linked odd-characteristic EC1 representations.",
    rationale:
      "The current ledger binds their catalog identities but does not attach a catalog-wide native generator arithmetic replay.",
    evidence_refs: [],
    obligation_id: "catalog_generator_validation",
  },
  "subgroup.prime-transfer.unknown": {
    status: "unknown",
    scope:
      "Pullback validity, image order, kernel intersection, composition, and degree-inverse recovery for linked degree-two maps.",
    rationale:
      "Coprimality of 2 and the catalog order value would be useful after a compatible pullback/composition certificate, but neither induced map is implemented for the two-infinity sources.",
    evidence_refs: [],
    obligation_id: "prime_balanced_jacobian",
  },
  "subgroup.arbitrary-recovery.unknown": {
    status: "unknown",
    scope:
      "Operational recovery from an arbitrary class produced by a Jacobian algorithm.",
    rationale:
      "There is no general norm, Jacobian relation solver, or target decomposition producing such a result.",
    evidence_refs: [],
    obligation_id: "full_solver",
  },
  "subgroup.exact-order29.supported": {
    status: "supported",
    scope:
      "The subgroup generated by EC1N7Ce0hb6d297a2ca08's exact point over F_2^7.",
    rationale:
      "The fixture checks a nonidentity G with [29]G=O; primality of 29 proves exact source order 29.",
    evidence_refs: [
      "src/binary_ecc/cover_transfer.rs#subgroup_status_and_arbitrary_pushforward_boundary",
      "research/hyperelliptic_cover_infrastructure_20261006/evidence.json#subgroup_transfer",
    ],
  },
  "subgroup.exact-image-order29.supported": {
    status: "supported",
    scope:
      "The checked pullback image of EC1N7Ce0hb6d297a2ca08's generator.",
    rationale:
      "The image is nonidentity and is killed by the prime 29, hence its exact order is 29.",
    evidence_refs: [
      "src/binary_ecc/cover_transfer.rs#subgroup_status_and_arbitrary_pushforward_boundary",
    ],
  },
  "map.f2m7-geometry.supported": {
    status: "supported",
    scope:
      "The exact source v^2+u^2v=u^7+u and phi(u,v)=(u^3,u*v+1) over F_2[z]/(z^7+z+1).",
    rationale:
      "The historical detailed replay and current HC1/CV1 identities independently bind the displayed genus-three cover.",
    evidence_refs: [
      "research/hyperelliptic_cover_infrastructure_20261006/evidence.json#icv1-f2m7-t13-616700dd",
      "docs/curves/cover-links.yaml#icv1-f2m7-t13-616700dd",
    ],
  },
  "map.fp827-geometry.supported": {
    status: "supported",
    scope:
      "The exact source v^2=u^6+u^2+15 and phi(u,v)=(u^2,v) over F_827.",
    rationale:
      "The historical detailed replay and current HC1/CV1 identities independently bind the displayed genus-two two-infinity cover.",
    evidence_refs: [
      "research/hyperelliptic_cover_infrastructure_20261006/evidence.json#icv1-fp10-t5-192cb216",
      "docs/curves/cover-links.yaml#icv1-fp10-t5-192cb216",
    ],
  },
};

const profiles = {
  "map.binary-cubic.same-field": {
    kind: "map_base",
    predicate_id: "binary_cubic_same_field",
    hypotheses: [
      "certificate.construction=binary_cubic_pullback_v1",
      "field_relation=same_field",
      "field_check=irreducibility_proved_by_Rabin",
      "ordinary binary Weierstrass model with b nonzero",
    ],
    claims: {
      cover_existence: "map.geometry.catalog-supported",
      explicit_construction: "map.binary-construction.supported",
      curve_map_correctness: "map.geometry.catalog-supported",
      jacobian_arithmetic_for_cover: "map.binary-jacobian.supported",
      pullback_point_classes: "map.binary-pullback.supported",
      independent_pullback_homomorphism_oracle: "map.independent-oracle.unknown",
      pushforward_point_divisors: "map.binary-typed-pushforward.supported",
      arbitrary_jacobian_pushforward: "map.binary-arbitrary-pushforward.refuted",
      composition_on_typed_domain: "map.binary-composition.supported",
      alternate_target_model_map: "map.direct-target.not-applicable",
      field_extension_or_descent: "map.same-field.not-applicable",
      computational_advantage: "map.advantage.unknown",
    },
    full_transfer_pipeline: "map.binary-pipeline.refuted",
    boundary:
      "Supported as a same-field genus-three cover with checked one-infinity arithmetic and bounded rational-point-class transfer. Not a general Jacobian norm or solver.",
  },
  "map.odd-prime-quadratic.direct": {
    kind: "map_base",
    predicate_id: "odd_prime_quadratic_direct_same_field",
    hypotheses: [
      "certificate.construction=prime_quadratic_pullback_v1",
      "declared field has no extension modulus",
      "certificate.target_model_map is absent",
      "base-field primality remains the catalog's stated input assumption",
    ],
    claims: {
      cover_existence: "map.geometry.catalog-supported",
      explicit_construction: "map.prime-construction.supported",
      curve_map_correctness: "map.geometry.catalog-supported",
      jacobian_arithmetic_for_cover: "map.prime-jacobian.refuted",
      pullback_point_classes: "map.prime-transfer.unknown",
      independent_pullback_homomorphism_oracle: "map.independent-oracle.unknown",
      pushforward_point_divisors: "map.prime-transfer.unknown",
      arbitrary_jacobian_pushforward: "map.prime-transfer.unknown",
      composition_on_typed_domain: "map.prime-transfer.unknown",
      alternate_target_model_map: "map.direct-target.not-applicable",
      field_extension_or_descent: "map.same-field.not-applicable",
      computational_advantage: "map.advantage.unknown",
    },
    full_transfer_pipeline: "map.prime-pipeline.refuted",
    boundary:
      "Supported as explicit same-field genus-two cover geometry only. Balanced/two-infinity Jacobian arithmetic, induced maps, and transfer advantage remain unavailable or unknown.",
  },
  "map.odd-prime-quadratic.normalized": {
    kind: "map_base",
    predicate_id: "odd_prime_quadratic_normalized_same_field",
    hypotheses: [
      "certificate.construction=prime_quadratic_pullback_v1",
      "declared field has no extension modulus",
      "certificate.target_model_map is a verified degree-one algebraic map",
      "base-field primality remains the catalog's stated input assumption",
    ],
    claims: {
      cover_existence: "map.geometry.catalog-supported",
      explicit_construction: "map.prime-construction.supported",
      curve_map_correctness: "map.geometry.catalog-supported",
      jacobian_arithmetic_for_cover: "map.prime-jacobian.refuted",
      pullback_point_classes: "map.prime-transfer.unknown",
      independent_pullback_homomorphism_oracle: "map.independent-oracle.unknown",
      pushforward_point_divisors: "map.prime-transfer.unknown",
      arbitrary_jacobian_pushforward: "map.prime-transfer.unknown",
      composition_on_typed_domain: "map.prime-transfer.unknown",
      alternate_target_model_map: "map.prime-normalization.supported",
      field_extension_or_descent: "map.same-field.not-applicable",
      computational_advantage: "map.advantage.unknown",
    },
    full_transfer_pipeline: "map.prime-pipeline.refuted",
    boundary:
      "Supported cover geometry includes the degree-one target normalization. This does not supply balanced Jacobian arithmetic or induced transfer maps.",
  },
  "map.odd-extension-quadratic.direct": {
    kind: "map_base",
    predicate_id: "odd_extension_quadratic_direct_same_field",
    hypotheses: [
      "certificate.construction=prime_quadratic_pullback_v1",
      "declared field includes a Rabin-validated extension modulus",
      "certificate.target_model_map is absent",
      "base-field primality remains the catalog's stated input assumption",
    ],
    claims: {
      cover_existence: "map.geometry.catalog-supported",
      explicit_construction: "map.extension-construction.supported",
      curve_map_correctness: "map.geometry.catalog-supported",
      jacobian_arithmetic_for_cover: "map.prime-jacobian.refuted",
      pullback_point_classes: "map.prime-transfer.unknown",
      independent_pullback_homomorphism_oracle: "map.independent-oracle.unknown",
      pushforward_point_divisors: "map.prime-transfer.unknown",
      arbitrary_jacobian_pushforward: "map.prime-transfer.unknown",
      composition_on_typed_domain: "map.prime-transfer.unknown",
      alternate_target_model_map: "map.direct-target.not-applicable",
      field_extension_or_descent: "map.extension-descent.unknown",
      computational_advantage: "map.advantage.unknown",
    },
    full_transfer_pipeline: "map.prime-pipeline.refuted",
    boundary:
      "Supported as geometry over the exact declared extension field. Proper-subfield descent, balanced Jacobian arithmetic, and induced maps remain unproved or unavailable.",
  },
  "subgroup.binary-degree3": {
    kind: "subgroup_base",
    predicate_id: "linked_subgroup_degree3",
    hypotheses: [
      "the exact EC1/curve UID record links to a binary cubic map",
      "the catalog subgroup-order integer is coprime to 3",
      "the native batch certificate decodes and checks the linked generator",
    ],
    claims: {
      subgroup_identity: "subgroup.identity.supported",
      gcd_degree_catalog_order: "subgroup.degree-coprime.supported",
      generator_on_curve: "subgroup.binary-generator.supported",
      catalog_order_annihilates_generator: "subgroup.binary-annihilation.supported",
      numeric_exact_source_order: "subgroup.numeric-order.unknown",
      pullback_valid: "subgroup.binary-pullback.supported",
      image_annihilated_by_catalog_order:
        "subgroup.binary-image-annihilation.supported",
      composition_on_generator: "subgroup.binary-composition.supported",
      kernel_intersection_on_generated_subgroup:
        "subgroup.binary-injective-generated.unknown",
      image_order_equals_actual_source_order:
        "subgroup.binary-image-actual-order.unknown",
      image_order_equals_catalog_order: "subgroup.image-numeric-order.unknown",
      degree_inverse_recovery: "subgroup.binary-inverse-recovery.supported",
      recovery_from_arbitrary_jacobian_result:
        "subgroup.arbitrary-recovery.unknown",
      computational_advantage: "map.advantage.unknown",
    },
    boundary:
      "Catalog-wide binary checks prove the listed generator-level typed identities and degree-inverse recovery. Pullback additivity, all-element composition, generated-subgroup injectivity, and actual image-order preservation remain conditional; exact numeric-order equality is proved only where an overlay supplies proof-grade evidence.",
  },
  "subgroup.odd-degree2": {
    kind: "subgroup_base",
    predicate_id: "linked_subgroup_degree2",
    hypotheses: [
      "the exact EC1/curve UID record links to an odd-characteristic quadratic map",
      "the catalog subgroup-order integer is odd and therefore coprime to 2",
      "no balanced two-infinity induced-map implementation is inferred from coprimality",
    ],
    claims: {
      subgroup_identity: "subgroup.identity.supported",
      gcd_degree_catalog_order: "subgroup.degree-coprime.supported",
      generator_on_curve: "subgroup.prime-generator.unknown",
      catalog_order_annihilates_generator: "subgroup.prime-generator.unknown",
      numeric_exact_source_order: "subgroup.numeric-order.unknown",
      pullback_valid: "subgroup.prime-transfer.unknown",
      image_annihilated_by_catalog_order: "subgroup.prime-transfer.unknown",
      composition_on_generator: "subgroup.prime-transfer.unknown",
      kernel_intersection_on_generated_subgroup: "subgroup.prime-transfer.unknown",
      image_order_equals_actual_source_order: "subgroup.prime-transfer.unknown",
      image_order_equals_catalog_order: "subgroup.prime-transfer.unknown",
      degree_inverse_recovery: "subgroup.prime-transfer.unknown",
      recovery_from_arbitrary_jacobian_result:
        "subgroup.arbitrary-recovery.unknown",
      computational_advantage: "map.advantage.unknown",
    },
    boundary:
      "The exact subgroup identities and degree/order coprimality are catalog-wide facts. Without compatible induced maps, no kernel, image-order, or recovery result follows.",
  },
};

function classifyMap(coverRow) {
  const certificate = coverRow.certificate;
  if (certificate.construction === "binary_cubic_pullback_v1") {
    assert(coverRow.field_relation === "same_field", `${coverRow.slug}: binary field relation`);
    assert(
      coverRow.field_check === "irreducibility_proved_by_Rabin",
      `${coverRow.slug}: binary field check`,
    );
    return "map.binary-cubic.same-field";
  }
  assert(
    certificate.construction === "prime_quadratic_pullback_v1",
    `${coverRow.slug}: unknown construction`,
  );
  assert(coverRow.field_relation === "same_field", `${coverRow.slug}: odd field relation`);
  const extension = coverRow.field_check.includes("modulus over GF(p)");
  const normalized = certificate.target_model_map != null;
  assert(!(extension && normalized), `${coverRow.slug}: unmodeled normalized extension row`);
  if (extension) return "map.odd-extension-quadratic.direct";
  if (normalized) return "map.odd-prime-quadratic.normalized";
  return "map.odd-prime-quadratic.direct";
}

function freshness(status, scope, rationale, evidenceRefs, obligationId) {
  const result = {
    status,
    scope,
    rationale,
    evidence_refs: evidenceRefs,
  };
  if (obligationId) result.obligation_id = obligationId;
  return result;
}

async function sourceArtifact(path, role, assessment) {
  return {
    path,
    sha256: await fileHash(path),
    role,
    freshness: assessment,
  };
}

function integrityArtifact(path, role) {
  return sourceArtifact(
    path,
    role,
    freshness(
      "supported",
      "byte identity of this assessment artifact at the generated-ledger revision",
      "The SHA-256 digest binds references to the exact repository bytes; it does not by itself promote any mathematical or performance claim.",
      [path],
    ),
  );
}

function implementationVersionArtifact(path) {
  return sourceArtifact(
    path,
    `exact implementation-version dependency: ${path}`,
    freshness(
      "supported",
      "byte identity of a local source, dependency manifest, lockfile, or toolchain workflow used by the assessed implementation",
      "The digest fixes the implementation version for replay. It does not widen arithmetic, transfer, subgroup, or performance support.",
      [path],
    ),
  );
}

const IMPLEMENTATION_VERSION_PATHS = [
  "Cargo.toml",
  "Cargo.lock",
  ".github/workflows/curve-covers.yml",
  ".github/workflows/release.yml",
  "src/lib.rs",
  "src/binary_ecc/mod.rs",
  "src/binary_ecc/f2m.rs",
  "src/binary_ecc/poly_f2m.rs",
  "src/binary_ecc/curve.rs",
  "src/cryptanalysis/mod.rs",
  "src/cryptanalysis/ecc2k130_guard.rs",
  "src/hash/mod.rs",
  "src/hash/sha256.rs",
  "src/prime_hyperelliptic/mod.rs",
  "src/prime_hyperelliptic/fp2.rs",
  "src/prime_hyperelliptic/fp_poly.rs",
  "src/utils/mod.rs",
];

async function buildSourceArtifacts() {
  return Promise.all([
    sourceArtifact(
      "docs/curves/registry.json",
      "exact ICV1/EC1 input registry",
      freshness(
        "supported",
        "the exact input registry for all 2,091 assessed map rows",
        "The digest is computed from the current file and every membership is joined through its model hash.",
        ["docs/curves/registry.json"],
      ),
    ),
    sourceArtifact(
      "docs/curves/covers.json",
      "generated curve-level geometry catalog",
      freshness(
        "supported",
        "all 2,091 current verified-over-declared-field geometry rows",
        "The ledger derives every map profile predicate from this exact file.",
        ["docs/curves/covers.json"],
      ),
    ),
    sourceArtifact(
      "docs/curves/cover-links.yaml",
      "content-addressed ICV1/EC1/HC1/CV1 graph",
      freshness(
        "supported",
        "2,091 exact map links, 2,089 unique cover objects, and 664 EC1 triples",
        "The ledger records exact sorted UID memberships rather than count-only coverage.",
        ["docs/curves/cover-links.yaml"],
      ),
    ),
    sourceArtifact(
      "docs/curves/cover-transfer-assessment.schema.json",
      "machine-readable v2 ledger contract",
      freshness(
        "supported",
        "the schema revision used for this generated ledger",
        "The digest binds the output to its validation contract.",
        ["docs/curves/cover-transfer-assessment.schema.json"],
      ),
    ),
    sourceArtifact(
      "scripts/build_cover_transfer_assessments.mjs",
      "deterministic catalog-wide ledger generator",
      freshness(
        "supported",
        "exact profile predicates, joins, sorted memberships, and domain-separated membership hashes",
        "The script provides byte-exact --write and --check modes using Node built-ins only.",
        ["scripts/build_cover_transfer_assessments.mjs"],
      ),
    ),
    sourceArtifact(
      "src/bin/curve_cover_check/checker.rs",
      "native cover constructor and geometry verifier",
      freshness(
        "supported",
        "all four map profiles at the pinned source revision",
        "The source implements the fixed-family construction and declared-field geometry checks.",
        ["src/bin/curve_cover_check/checker.rs"],
      ),
    ),
    sourceArtifact(
      "src/bin/curve_cover_check/models.rs",
      "native catalog-model parser and normalization formulas",
      freshness(
        "supported",
        "all direct, extension-field, Montgomery, and Edwards inputs used by the four map profiles",
        "The digest pins the model conversions and coefficient formulas consumed by the geometry verifier.",
        ["src/bin/curve_cover_check/models.rs"],
      ),
    ),
    sourceArtifact(
      "src/bin/curve_cover_check/links.rs",
      "native content-addressed cover-link graph builder",
      freshness(
        "supported",
        "the exact HC1/CV1/ICV1/EC1 graph joined by this assessment",
        "The digest pins the graph construction whose exact memberships are independently checked by the verifier.",
        ["src/bin/curve_cover_check/links.rs"],
      ),
    ),
    sourceArtifact(
      "src/bin/curve_cover_check.rs",
      "native cover-check command entry point",
      freshness(
        "supported",
        "the deterministic catalog and graph generation command",
        "The digest pins the command wiring and byte-exact check mode used by the catalog control.",
        ["src/bin/curve_cover_check.rs"],
      ),
    ),
    sourceArtifact(
      "src/bin/hyperelliptic_cover.rs",
      "standalone CLI and native binary subgroup batch certificate",
      freshness(
        "supported",
        "all 134 exact linked binary EC1 triples at the pinned source revision",
        "The batch test is cited only for the checks it performs; it is not proof of every numeric order's primality.",
        ["src/bin/hyperelliptic_cover.rs#every_binary_ec1_generator_has_checked_typed_transfer_hypotheses"],
      ),
    ),
    sourceArtifact(
      "src/binary_ecc/cover_transfer.rs",
      "bounded binary pullback, pushforward, and composition implementation",
      freshness(
        "supported",
        "the binary cubic profile's typed point-class domain",
        "The digest pins the exact current transfer API and explicit arbitrary-Jacobian negative boundary.",
        ["src/binary_ecc/cover_transfer.rs"],
      ),
    ),
    sourceArtifact(
      "src/binary_ecc/hyperelliptic.rs",
      "one-infinity binary Jacobian implementation",
      freshness(
        "supported",
        "the binary cubic profile's source model",
        "The digest pins the checked divisor representation and arithmetic implementation.",
        ["src/binary_ecc/hyperelliptic.rs"],
      ),
    ),
    sourceArtifact(
      "src/prime_hyperelliptic/curve.rs",
      "odd-prime one-infinity Jacobian model boundary",
      freshness(
        "supported",
        "evidence that the catalog's even sextics are outside the current one-infinity model",
        "The model contract requires degree 2g+1 and cannot be silently applied to the two-infinity sources.",
        ["src/prime_hyperelliptic/curve.rs"],
      ),
    ),
    ...IMPLEMENTATION_VERSION_PATHS.map(implementationVersionArtifact),
    sourceArtifact(
      "research/hyperelliptic_cover_infrastructure_20261006/evidence.json",
      "immutable detailed evidence for the two retained exact overlays",
      freshness(
        "unknown",
        "full applicability of the historical replay to every current catalog row",
        "The two exact identities remain useful overlays, but the historical run is not promoted to catalog-wide execution evidence.",
        ["research/hyperelliptic_cover_infrastructure_20261006/evidence.json"],
        "independent_map_oracle",
      ),
    ),
    sourceArtifact(
      "docs/curves/jacobian-transfers.json",
      "legacy formula-screen inventory; not a transfer certificate",
      freshness(
        "refuted",
        "claim that this inventory is current for the 2,091-row registry",
        "It pins an older registry, omits current slugs, and contains modeled routes without exact HC1/CV1 subgroup certificates.",
        [
          "docs/curves/jacobian-transfers.json",
          "docs/curves/JACOBIAN_TRANSFERS.md",
        ],
        "refresh_formula_inventory",
      ),
    ),
    integrityArtifact(
      "docs/curves/JACOBIAN_TRANSFERS.md",
      "legacy formula-screen documentation pinned for the stale-inventory finding",
    ),
    sourceArtifact(
      "docs/curves/COVERS.md",
      "catalog geometry contract and proof boundary",
      freshness(
        "supported",
        "all current fixed-family geometry claims",
        "The documentation states the construction, verification, and unsupported boundaries used by the profiles.",
        ["docs/curves/COVERS.md"],
      ),
    ),
    sourceArtifact(
      "docs/curves/HYPERELLIPTIC_INFRASTRUCTURE.md",
      "standalone CLI and supported-model documentation",
      freshness(
        "supported",
        "the binary one-infinity and odd-prime two-infinity implementation boundary",
        "The documentation distinguishes explicit cover geometry from available Jacobian arithmetic.",
        ["docs/curves/HYPERELLIPTIC_INFRASTRUCTURE.md"],
      ),
    ),
    integrityArtifact(
      "research/hyperelliptic_cover_infrastructure_20261006/README.md",
      "historical exact-fixture report cited by bounded negative controls",
    ),
    integrityArtifact(
      "research/cover_transfer_assessment_20261013/lean/CoverTransfer.lean",
      "conditional group-theoretic transfer lemma source; kernel replay is tracked separately",
    ),
    integrityArtifact(
      "research/cover_transfer_assessment_20261013/lean/README.md",
      "Lean package execution boundary and non-implication record",
    ),
    integrityArtifact(
      "research/cover_transfer_assessment_20261013/verify_assessment.mjs",
      "independent catalog-wide assessment verifier",
    ),
    integrityArtifact(
      "research/cover_transfer_assessment_20261013/README.md",
      "source-linked catalog-wide assessment report",
    ),
    integrityArtifact(
      "research/cover_transfer_assessment_20261013/release_replay.txt",
      "published-artifact replay receipt and source/package boundary",
    ),
    integrityArtifact(
      "research/cover_transfer_assessment_20261013/typed_correspondence.dot",
      "assessment diagram source",
    ),
    integrityArtifact(
      "research/cover_transfer_assessment_20261013/typed_correspondence.svg",
      "rendered assessment diagram in SVG format",
    ),
    integrityArtifact(
      "research/cover_transfer_assessment_20261013/typed_correspondence.png",
      "rendered assessment diagram in PNG format",
    ),
    integrityArtifact(
      "research/cover_transfer_assessment_20261013/typed_transfer_context.dot",
      "typed transfer-context diagram source",
    ),
    integrityArtifact(
      "research/cover_transfer_assessment_20261013/typed_transfer_context.svg",
      "rendered typed transfer-context diagram in SVG format",
    ),
    integrityArtifact(
      "research/cover_transfer_assessment_20261013/typed_transfer_context.png",
      "rendered typed transfer-context diagram in PNG format",
    ),
    integrityArtifact(
      "research/cover_transfer_assessment_20261013/cover_transfer_assessment_20261013.pdf",
      "rendered catalog-wide assessment PDF",
    ),
  ]);
}

function membership(profileId, members, domain) {
  return {
    profile_id: profileId,
    count: members.length,
    members_sha256: membershipHash(`${domain}/${profileId}`, members),
    members,
  };
}

async function buildLedger() {
  const [covers, graph] = await Promise.all([
    readJson("docs/curves/covers.json"),
    readJson("docs/curves/cover-links.yaml"),
  ]);
  assert(covers.curves.length === 2091, "unexpected covers row count");
  assert(Object.keys(graph.maps).length === 2091, "unexpected graph map count");
  assert(Object.keys(graph.covers).length === 2089, "unexpected graph cover count");
  assert(graph.curves.length === 2091, "unexpected graph curve count");

  const coverRowsByModel = new Map(
    covers.curves.map((row) => [row.model_sha256, row]),
  );
  assert(coverRowsByModel.size === 2091, "cover rows must have unique model hashes");

  const mapMembers = new Map(
    Object.entries(profiles)
      .filter(([, profile]) => profile.kind === "map_base")
      .map(([id]) => [id, []]),
  );
  const subgroupMembers = new Map(
    Object.entries(profiles)
      .filter(([, profile]) => profile.kind === "subgroup_base")
      .map(([id]) => [id, []]),
  );
  const fieldBindings = new Map();
  const jacobianBindings = new Map();
  const typedMapBindings = new Map();

  for (const curve of graph.curves) {
    assert(curve.cover_links.length === 1, `${curve.icv1_identity.slug}: cover link count`);
    const link = curve.cover_links[0];
    const map = graph.maps[link.map_uid];
    assert(map != null, `${curve.icv1_identity.slug}: unresolved map UID`);
    assert(map.map_uid === link.map_uid, `${curve.icv1_identity.slug}: map key mismatch`);
    const sourceCover = graph.covers[link.cover_uid];
    assert(sourceCover != null, `${curve.icv1_identity.slug}: unresolved cover UID`);
    assert(
      map.record.source_cover_uid === link.cover_uid,
      `${curve.icv1_identity.slug}: source cover mismatch`,
    );
    assert(
      map.record.target_model_sha256 === curve.model_sha256,
      `${curve.icv1_identity.slug}: target model mismatch`,
    );
    const row = coverRowsByModel.get(curve.model_sha256);
    assert(row != null, `${curve.icv1_identity.slug}: missing covers row`);
    assert(row.slug === curve.icv1_identity.slug, `${curve.icv1_identity.slug}: slug mismatch`);
    assert(row.status === "verified_over_declared_field", `${row.slug}: geometry status`);
    const mapProfileId = classifyMap(row);
    mapMembers.get(mapProfileId).push(link.map_uid);

    const targetModel = JSON.parse(curve.model_json);
    const targetFieldRecord = targetDeclaredFieldRecord(targetModel);
    assert(
      canonicalJson(sourceCover.record.field) === canonicalJson(targetFieldRecord),
      `${curve.icv1_identity.slug}: source and target declared fields differ`,
    );
    const field = declaredFieldBinding(sourceCover.record.field);
    const priorField = fieldBindings.get(field.field_uid);
    assert(
      priorField == null || canonicalJson(priorField) === canonicalJson(field),
      `${curve.icv1_identity.slug}: declared-field UID collision`,
    );
    fieldBindings.set(field.field_uid, field);

    const infinityConfiguration =
      mapProfileId === "map.binary-cubic.same-field"
        ? "one_rational_point"
        : "two_rational_points";
    const jacobianRecord = {
      schema_version: "jacobian-context/v1",
      source_cover_uid: link.cover_uid,
      declared_field_uid: field.field_uid,
      dimension: sourceCover.record.genus,
      object: "Picard_degree_zero",
      infinity_configuration: infinityConfiguration,
    };
    const jacobianUid =
      `urn:jacobian-context:1:sha256:${sha256Bytes(canonicalJson(jacobianRecord))}`;
    const jacobian = { jacobian_uid: jacobianUid, record: jacobianRecord };
    const priorJacobian = jacobianBindings.get(jacobianUid);
    assert(
      priorJacobian == null || canonicalJson(priorJacobian) === canonicalJson(jacobian),
      `${curve.icv1_identity.slug}: Jacobian-context UID collision`,
    );
    jacobianBindings.set(jacobianUid, jacobian);

    const typedMap = {
      map_uid: link.map_uid,
      target_icv1_slug: curve.icv1_identity.slug,
      target_model_sha256: curve.model_sha256,
      source_cover_uid: link.cover_uid,
      source_jacobian_uid: jacobianUid,
      field_relationship: {
        kind: "same_declared_field_identity",
        source_field_uid: field.field_uid,
        target_field_uid: field.field_uid,
        degree: 1,
      },
    };
    assert(!typedMapBindings.has(link.map_uid), `${curve.icv1_identity.slug}: duplicate typed map`);
    typedMapBindings.set(link.map_uid, typedMap);

    for (const representation of curve.representations) {
      const triple = [representation.ec1, representation.curve_uid, link.map_uid];
      const subgroupOrder = BigInt(representation.curve.subgroup_order);
      if (map.record.degree === 3) {
        assert(subgroupOrder % 3n !== 0n, `${representation.ec1}: order divisible by 3`);
        subgroupMembers.get("subgroup.binary-degree3").push(triple);
      } else {
        assert(map.record.degree === 2, `${representation.ec1}: unexpected map degree`);
        assert(subgroupOrder % 2n === 1n, `${representation.ec1}: even order`);
        subgroupMembers.get("subgroup.odd-degree2").push(triple);
      }
    }
  }

  const mapBase = [...mapMembers.entries()].map(([profileId, values]) =>
    membership(profileId, sortedUnique(values, profileId), MAP_HASH_DOMAIN),
  );
  const subgroupBase = [...subgroupMembers.entries()].map(([profileId, values]) =>
    membership(
      profileId,
      sortedUniqueTriples(values, profileId),
      SUBGROUP_HASH_DOMAIN,
    ),
  );
  const mapCount = mapBase.reduce((sum, group) => sum + group.count, 0);
  const subgroupCount = subgroupBase.reduce((sum, group) => sum + group.count, 0);
  assert(mapCount === 2091, `map membership count is ${mapCount}`);
  assert(subgroupCount === 664, `subgroup membership count is ${subgroupCount}`);
  assert(fieldBindings.size === 311, `declared-field binding count is ${fieldBindings.size}`);
  assert(jacobianBindings.size === 2089, `Jacobian binding count is ${jacobianBindings.size}`);
  assert(typedMapBindings.size === 2091, `typed map binding count is ${typedMapBindings.size}`);
  const expectedMapCounts = new Map([
    ["map.binary-cubic.same-field", 160],
    ["map.odd-prime-quadratic.direct", 577],
    ["map.odd-prime-quadratic.normalized", 30],
    ["map.odd-extension-quadratic.direct", 1324],
  ]);
  for (const group of mapBase) {
    assert(
      group.count === expectedMapCounts.get(group.profile_id),
      `${group.profile_id}: expected ${expectedMapCounts.get(group.profile_id)}, got ${group.count}`,
    );
  }
  const expectedSubgroupCounts = new Map([
    ["subgroup.binary-degree3", 134],
    ["subgroup.odd-degree2", 530],
  ]);
  for (const group of subgroupBase) {
    assert(
      group.count === expectedSubgroupCounts.get(group.profile_id),
      `${group.profile_id}: expected ${expectedSubgroupCounts.get(group.profile_id)}, got ${group.count}`,
    );
  }

  const f2Map =
    "urn:curve-cover-map:1:sha256:0dc60016b03df7d55ee829880508251fe8261c92a241b76a78987eeaa9787a03";
  const fp827Map =
    "urn:curve-cover-map:1:sha256:5768ce11909c0f56d84f05e70d71dfb6a46fb175b6f9cc0ed727eec14aef9c8e";
  const f2Triple = [
    "EC1N7Ce0hb6d297a2ca08",
    "urn:ec-record:1:sha256:b6d297a2ca085e9d7b8a6858df7421083e12e740c89b1647660b413e00ed02ff",
    f2Map,
  ];
  const fp827Triple = [
    "EC1P10Cfphcf15d87ab724",
    "urn:ec-record:1:sha256:cf15d87ab724c8d1fc74d88af320937e4c1cc7e706d4e0e8a0d11734831342a3",
    fp827Map,
  ];
  assert(mapBase.some((group) => group.members.includes(f2Map)), "missing F2 map overlay target");
  assert(mapBase.some((group) => group.members.includes(fp827Map)), "missing F827 map overlay target");
  assert(
    subgroupBase.some((group) =>
      group.members.some((value) => JSON.stringify(value) === JSON.stringify(f2Triple)),
    ),
    "missing order-29 subgroup overlay target",
  );
  assert(
    subgroupBase.some((group) =>
      group.members.some((value) => JSON.stringify(value) === JSON.stringify(fp827Triple)),
    ),
    "missing order-823 subgroup overlay target",
  );

  return {
    schema_version: "cover-transfer-assessments/v2",
    status_vocabulary: ["unknown", "supported", "refuted", "not_applicable"],
    assessment_date: ASSESSMENT_DATE,
    catalog_commit: CATALOG_COMMIT,
    source_artifacts: await buildSourceArtifacts(),
    claim_contract: {
      map: {
        cover_existence: "An explicit source cover exists over the stated field.",
        explicit_construction:
          "The source equation, target, map formula, degree, and field are constructive and recorded.",
        curve_map_correctness:
          "Smoothness, genus, defining-equation identity, degree, separability, and infinity behavior are checked.",
        jacobian_arithmetic_for_cover:
          "The current Jacobian implementation accepts the source infinity configuration and provides checked class arithmetic.",
        pullback_point_classes:
          "The current implementation constructs pullback on its stated point-class domain.",
        independent_pullback_homomorphism_oracle:
          "An implementation-independent oracle checks induced-map behavior.",
        pushforward_point_divisors:
          "The current implementation constructs pushforward on its stated source-divisor domain.",
        arbitrary_jacobian_pushforward:
          "The current implementation norms arbitrary supported Jacobian classes.",
        composition_on_typed_domain:
          "The pull-push composition equals multiplication by map degree on the stated typed domain.",
        alternate_target_model_map:
          "Any required degree-one map between the construction target and catalog target is explicit and checked.",
        field_extension_or_descent:
          "Any field extension or proper-subfield descent used by the transfer is explicit and checked.",
        computational_advantage:
          "A matched end-to-end experiment establishes a computational advantage.",
      },
      subgroup: {
        subgroup_identity:
          "The exact EC1 record, curve UID, generator data, and map UID are joined.",
        gcd_degree_catalog_order:
          "The map degree is coprime to the numeric catalog subgroup-order value.",
        generator_on_curve:
          "The exact linked generator decodes canonically and lies on the target curve.",
        catalog_order_annihilates_generator:
          "Multiplication by the numeric catalog subgroup-order value annihilates the generator.",
        numeric_exact_source_order:
          "The generator's exact order equals the numeric catalog subgroup-order value.",
        pullback_valid:
          "The exact generator has a valid pullback image in the source Jacobian representation.",
        image_annihilated_by_catalog_order:
          "The numeric catalog subgroup-order value annihilates the pullback image.",
        composition_on_generator:
          "Pull-push on the exact generator equals multiplication by map degree.",
        kernel_intersection_on_generated_subgroup:
          "The pullback kernel intersects the actual subgroup generated by G trivially.",
        image_order_equals_actual_source_order:
          "The pullback image has the same actual order as G.",
        image_order_equals_catalog_order:
          "The pullback image's exact order equals the numeric catalog value.",
        degree_inverse_recovery:
          "The typed composition value recovers G using the degree inverse modulo the catalog order.",
        recovery_from_arbitrary_jacobian_result:
          "An arbitrary source-Jacobian result can be normed and recovered on the target subgroup.",
        computational_advantage:
          "The subgroup transfer participates in a matched end-to-end speedup.",
      },
    },
    result_templates: resultTemplates,
    profiles,
    memberships: {
      map_base: mapBase,
      subgroup_base: subgroupBase,
    },
    typed_bindings: {
      fields: sortedRecord(fieldBindings),
      jacobians: sortedRecord(jacobianBindings),
      maps: sortedRecord(typedMapBindings),
    },
    exact_overlays: [
      {
        overlay_id: "map.binary-f2m7-degree3",
        kind: "map",
        target: f2Map,
        claim_overrides: {
          curve_map_correctness: "map.f2m7-geometry.supported",
        },
        evidence_additions: [
          "research/hyperelliptic_cover_infrastructure_20261006/evidence.json#icv1-f2m7-t13-616700dd",
          "research/hyperelliptic_cover_infrastructure_20261006/README.md",
        ],
        boundary:
          "Exact overlay for ICV1 icv1-f2m7-t13-616700dd and HC1 2530bd70db89; it adds detailed replay evidence without widening the profile's typed-domain boundary.",
      },
      {
        overlay_id: "subgroup.binary-f2m7-order29",
        kind: "subgroup",
        target: f2Triple,
        claim_overrides: {
          numeric_exact_source_order: "subgroup.exact-order29.supported",
          image_order_equals_catalog_order:
            "subgroup.exact-image-order29.supported",
        },
        evidence_additions: [
          "research/hyperelliptic_cover_infrastructure_20261006/evidence.json#subgroup_transfer",
          "src/binary_ecc/cover_transfer.rs#subgroup_status_and_arbitrary_pushforward_boundary",
        ],
        boundary:
          "Exact order-29 evidence upgrades only this triple's numeric source- and image-order claims; arbitrary Jacobian recovery remains unknown.",
      },
      {
        overlay_id: "map.prime-fp827-degree2",
        kind: "map",
        target: fp827Map,
        claim_overrides: {
          curve_map_correctness: "map.fp827-geometry.supported",
        },
        evidence_additions: [
          "research/hyperelliptic_cover_infrastructure_20261006/evidence.json#icv1-fp10-t5-192cb216",
          "research/hyperelliptic_cover_infrastructure_20261006/README.md",
        ],
        boundary:
          "Exact overlay for ICV1 icv1-fp10-t5-192cb216 and its F_827 two-infinity source; induced maps remain unknown.",
      },
      {
        overlay_id: "subgroup.prime-fp827-order823",
        kind: "subgroup",
        target: fp827Triple,
        claim_overrides: {},
        evidence_additions: [
          "docs/curves/cover-links.yaml#EC1P10Cfphcf15d87ab724",
          "research/hyperelliptic_cover_infrastructure_20261006/evidence.json#icv1-fp10-t5-192cb216",
        ],
        boundary:
          "The exact order-823 catalog identity is retained, but no transfer-preservation claim follows from cover geometry alone.",
      },
    ],
    coverage: {
      catalog_maps: 2091,
      unique_cover_models: 2089,
      geometry_verified_maps: 2091,
      ec1_representations: 664,
      assessed_maps: 2091,
      unassessed_maps: 0,
      assessed_subgroups: 664,
      unassessed_subgroups: 0,
      map_membership_entries: mapCount,
      subgroup_membership_entries: subgroupCount,
      map_profiles: mapBase.length,
      subgroup_profiles: subgroupBase.length,
      scope_note:
        "Every current map UID and every exact [EC1, curve UID, map UID] triple has one base profile. Assessment means each claim receives an explicit supported, refuted, unknown, or not-applicable result; it does not mean every transfer capability is supported.",
    },
    global_cost_accounting: [
      ["cover_construction", null, null, "Executable constructors exist, but no calibrated operation count or timing is claimed."],
      ["failed_construction_searches", null, "attempts", "The fixed families do not establish general search accounting."],
      ["duplicate_candidates", null, "candidates", "No candidate-deduplication cost ledger was produced."],
      ["geometry_verification", null, null, "Verification is replayable but uncalibrated."],
      ["jacobian_arithmetic", null, null, "Correctness receipts are not performance measurements."],
      ["pullback_and_pushforward", null, null, "Bounded binary formulas are implemented, but measured conversion cost is absent."],
      ["subgroup_recovery", null, null, "Typed-domain correctness is not a calibrated recovery benchmark."],
      ["reusable_preprocessing", null, null, "No reusable-preprocessing or amortization ledger exists."],
      ["peak_memory", null, "bytes", "Peak memory was not measured."],
      ["total_runtime", null, "milliseconds", "No end-to-end transfer or solver run exists."],
      ["matched_generic_baseline", null, null, "No matched Pollard-rho or same-target baseline was run."],
      ["end_to_end_speedup", null, "ratio", "Cover existence and bounded transfer do not establish a computational advantage."],
    ].map(([phase, value, unit, note]) => ({
      phase,
      status: "unknown",
      value,
      unit,
      source_refs: [],
      note,
    })).concat([
      {
        phase: "jacobian_relations_linear_algebra_target_decomposition",
        status: "not_applicable",
        value: null,
        unit: null,
        source_refs: [],
        note:
          "These phases belong to the explicitly separate full Jacobian index-calculus solver, which was not implemented or run.",
      },
    ]),
    controls: [
      {
        control: "identity",
        assessment: freshness(
          "supported",
          "all 2,091 map UIDs and 664 exact EC1/curve/map triples",
          "The generator resolves exact content-addressed joins and rejects duplicates or missing rows.",
          [
            "scripts/build_cover_transfer_assessments.mjs",
            "docs/curves/cover-links.yaml",
          ],
        ),
        command: "node scripts/build_cover_transfer_assessments.mjs --check",
      },
      {
        control: "geometry",
        assessment: freshness(
          "supported",
          "all 2,091 fixed-family cover rows",
          "Native certificates cover defining identities, smoothness, genus, degree, separability, and infinity behavior.",
          ["src/bin/curve_cover_check/checker.rs", "docs/curves/covers.json"],
        ),
        command: "cargo test --bin curve_cover_check",
      },
      {
        control: "arithmetic_internal",
        assessment: freshness(
          "supported",
          "the implemented one-infinity prime and binary fixtures, not the even-sextic prime catalog sources",
          "Retained tests cover canonicality, closure, inverses, scalar multiplication, exceptional divisors, and shared support.",
          [
            "src/binary_ecc/hyperelliptic.rs",
            "research/hyperelliptic_cover_infrastructure_20261006/evidence.json",
          ],
        ),
        command: "cargo test --lib hyperelliptic",
      },
      {
        control: "arithmetic_independent_oracle",
        assessment: freshness(
          "unknown",
          "independent CAS or separately implemented divisor-class comparison across the catalog",
          "Current arithmetic and transfer checks are implementation-internal.",
          [],
          "independent_map_oracle",
        ),
        command: null,
      },
      {
        control: "binary_catalog_subgroups",
        assessment: freshness(
          "supported",
          "native generator-level checks for all 134 binary representations after independent exact EC1/curve/map identity joins",
          "The semantic verifier reconstructs the exact graph joins. The native batch certificate reads those representation rows and checks generator membership, catalog-order annihilation, pullback validity, image annihilation, composition on the tested generator, and degree-inverse recovery of that generator; it does not itself validate every content UID and does not certify subgroup-wide injectivity.",
          [
            "docs/curves/cover-links.yaml",
            "research/cover_transfer_assessment_20261013/verify_assessment.mjs",
            "src/bin/hyperelliptic_cover.rs#every_binary_ec1_generator_has_checked_typed_transfer_hypotheses",
          ],
        ),
        command:
          "cargo test --bin hyperelliptic-cover every_binary_ec1_generator_has_checked_typed_transfer_hypotheses",
      },
      {
        control: "odd_degree2_transfer",
        assessment: freshness(
          "refuted",
          "claim that the current one-infinity odd-prime implementation supports the 1,931 even-sextic source models",
          "The current model contract and the catalog infinity configuration are structurally incompatible.",
          ["src/prime_hyperelliptic/curve.rs", "docs/curves/covers.json"],
        ),
        command: null,
      },
      {
        control: "published_artifact",
        assessment: freshness(
          "refuted",
          "claim that the published v0.1.0-main.1371 artifact provides executable evidence for direct construction on all 2,091 current catalog models",
          "Direct artifact replay covered only its bundled 592-row registry. Static comparison of the pinned source allowlist with the current partition places 160 binary and 577 direct prime-field forms inside that syntactic domain, and 1,324 extension-field plus 30 normalized-target forms outside it; this 737/1,354 split is not a 2,091-command artifact receipt. Current source adds those paths but has not been published in that artifact.",
          [
            "research/cover_transfer_assessment_20261013/release_replay.txt",
            "src/bin/hyperelliptic_cover.rs#every_current_catalog_model_reaches_the_standalone_constructor",
          ],
          "catalog_wide_packaged_binary",
        ),
        command: null,
      },
      {
        control: "cross_repository_core_parity",
        assessment: freshness(
          "supported",
          "the seven bindings listed by the cryptanalysis hyperelliptic source-parity manifest",
          "The offline parity replay verified the five byte-exact core files, one binary-module rustdoc adapter, and one selected-domain polynomial binding. It excludes the standalone CLI, Cargo metadata, catalogs, this assessment, workflows, packages, and live-HEAD binding.",
          ["research/cover_transfer_assessment_20261013/release_replay.txt"],
        ),
        command:
          "cargo run --locked --manifest-path suite/Cargo.toml --bin hyperelliptic-source-parity -- --crypto-root /workspace/crypto-transfer-assessment",
      },
      {
        control: "coordinate_controls",
        assessment: freshness(
          "unknown",
          "commutation of induced maps with alternate elliptic coordinates and divisor representations",
          "The catalog records the 30 degree-one target normalizations, but no induced-divisor coordinate-commutation control has been executed.",
          [],
          "coordinate_controls",
        ),
        command: null,
      },
      {
        control: "paired_performance",
        assessment: freshness(
          "not_applicable",
          "this correctness and evidence-coverage assessment",
          "No performance or speedup experiment was conducted.",
          [],
        ),
        command: null,
      },
    ],
    exploration_boundary: [
      "No search for production targets, private keys, or discrete logarithms was performed.",
      "All 2,091 map rows and all 664 linked EC1 representations are classified, but unknown remains unknown and refuted implementation support is not a mathematical nonexistence claim.",
      "The 160 binary maps have bounded typed point-class operations; arbitrary source-Jacobian norm remains unsupported.",
      "The 1,931 odd-characteristic maps have verified cover geometry but no compatible balanced/two-infinity Jacobian transfer implementation.",
      "Generator checks and [r]G=O do not by themselves prove that every catalog integer r is the exact source order.",
      "Formula-screen outputs in jacobian-transfers.json were not promoted to executable or measured claims.",
      "The full Jacobian index-calculus solver, relation collection, linear algebra, and target decomposition remain a separate deliverable.",
    ],
    open_obligations: obligations,
    report: {
      markdown: "research/cover_transfer_assessment_20261013/README.md",
      pdf:
        "research/cover_transfer_assessment_20261013/cover_transfer_assessment_20261013.pdf",
      graph_source:
        "research/cover_transfer_assessment_20261013/typed_correspondence.dot",
      graph_svg:
        "research/cover_transfer_assessment_20261013/typed_correspondence.svg",
      graph_png:
        "research/cover_transfer_assessment_20261013/typed_correspondence.png",
      context_graph_source:
        "research/cover_transfer_assessment_20261013/typed_transfer_context.dot",
      context_graph_svg:
        "research/cover_transfer_assessment_20261013/typed_transfer_context.svg",
      context_graph_png:
        "research/cover_transfer_assessment_20261013/typed_transfer_context.png",
      release_receipt:
        "research/cover_transfer_assessment_20261013/release_replay.txt",
    },
    migration: {
      from_schema: "cover-transfer-assessments/v1",
      summary:
        "Replaces two duplicated exact assessment records with reusable claim profiles, exact catalog-wide memberships, and retained exact-evidence overlays.",
      compatibility_notes: [
        "Consumers must apply exactly one base profile to each membership target, then apply zero or more matching overlays.",
        "A supported geometry profile does not imply a supported Jacobian transfer profile.",
        "Unknown is an evidence state, not a failed search or nonexistence result.",
        "Map membership hashes use domain cover-transfer-assessments/v2/map_base and subgroup membership hashes use cover-transfer-assessments/v2/subgroup_base; each hashes that domain plus '/' and the profile ID, a NUL byte, and compact JSON of the sorted members array.",
      ],
    },
  };
}

async function main() {
  const mode = process.argv[2];
  if (!['--write', '--check'].includes(mode) || process.argv.length !== 3) {
    console.error(
      "usage: node scripts/build_cover_transfer_assessments.mjs (--write|--check)",
    );
    process.exitCode = 2;
    return;
  }
  const rendered = `${JSON.stringify(await buildLedger(), null, 2)}\n`;
  const output = resolve(ROOT, OUTPUT_PATH);
  if (mode === "--write") {
    await writeFile(output, rendered);
    console.log(`wrote ${OUTPUT_PATH} (${Buffer.byteLength(rendered)} bytes)`);
    return;
  }
  const current = await readFile(output, "utf8");
  if (current !== rendered) {
    console.error(
      `${OUTPUT_PATH} is stale; run node scripts/build_cover_transfer_assessments.mjs --write`,
    );
    process.exitCode = 1;
    return;
  }
  console.log(`${OUTPUT_PATH} is byte-exact and current`);
}

await main();
