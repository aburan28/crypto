#!/usr/bin/env node

import { createHash } from "node:crypto";
import { existsSync, readFileSync, realpathSync } from "node:fs";
import { dirname, isAbsolute, relative, resolve, sep } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const repo = resolve(here, "../..");
const readText = (path) => readFileSync(resolve(repo, path), "utf8");
const readJson = (path) => JSON.parse(readText(path));
const sha256 = (path) =>
  createHash("sha256").update(readFileSync(resolve(repo, path))).digest("hex");

const assessmentPath = "docs/curves/cover-transfer-assessments.json";
const expectedCatalogCommit = "6600c7c6c0cfbb6a9e4e596bf8c5cfa292fb269e";
const expectedAssessmentSha256 =
  "edd6101f4a9339d14ac02c259f7c323a718da6dee7918097f1053f016b4db9a8";
const expectedSourceArtifactPaths = [
  "docs/curves/registry.json",
  "docs/curves/covers.json",
  "docs/curves/cover-links.yaml",
  "docs/curves/cover-links.schema.json",
  "docs/curves/jacobian-transfers.json",
  "docs/curves/JACOBIAN_TRANSFERS.md",
  "research/hyperelliptic_cover_infrastructure_20261006/evidence.json",
  "research/hyperelliptic_cover_infrastructure_20261006/README.md",
  "src/binary_ecc/cover_transfer.rs",
  "src/binary_ecc/hyperelliptic.rs",
  "src/prime_hyperelliptic/curve.rs",
  "src/bin/curve_cover_check/checker.rs",
  "src/bin/hyperelliptic_cover.rs",
  "docs/curves/COVERS.md",
  "docs/curves/HYPERELLIPTIC_INFRASTRUCTURE.md",
].sort();
const expectedObligationIds = [
  "general_binary_norm",
  "independent_map_oracle",
  "prime_balanced_jacobian",
  "catalog_subgroup_certificates",
  "coordinate_controls",
  "refresh_formula_inventory",
  "catalog_wide_assessment",
  "full_solver",
].sort();
const expectedPriorImplementations = {
  "src/binary_ecc/cover_transfer.rs":
    "implementation:binary_transfer_feature_snapshot",
  "src/binary_ecc/hyperelliptic.rs":
    "implementation:binary_mumford_feature_snapshot",
  "src/prime_hyperelliptic/curve.rs":
    "implementation:prime_mumford_feature_snapshot",
};
const expectedReportPaths = {
  markdown: "research/cover_transfer_assessment_20261013/README.md",
  pdf: "research/cover_transfer_assessment_20261013/cover_transfer_assessment_20261013.pdf",
  graph_source:
    "research/cover_transfer_assessment_20261013/typed_correspondence.dot",
  graph_svg: "research/cover_transfer_assessment_20261013/typed_correspondence.svg",
  graph_png: "research/cover_transfer_assessment_20261013/typed_correspondence.png",
  release_receipt:
    "research/cover_transfer_assessment_20261013/release_replay.txt",
};
const expectedReportSha256 = {
  markdown: "f91c44b690004cd140feb84876adbf1912abdc07ecf283fbc8f5a335b68d3c5e",
  pdf: "2fa5784cef2cb12f6eead867787177fb0d8c0914460e36c362690ebec2bb2c7b",
  graph_source:
    "2f25e173860dfe8eae36a2c8e742de5319239af8d21aaba172efc17d5c2b88b1",
  graph_svg: "c2b7169ff173b61c1c1e9779239398d5de0834f91026f80eaa20a51ea93caa6b",
  graph_png: "05c418deec8ed6a972cf22cfb74a94480ca5e3f72e9e81bc3c5eb831434988d3",
  release_receipt:
    "f883e9de7c6445f3c406e435c83bba55e5ef0d334636ae40153e6fa670d12036",
};
const assessment = readJson(assessmentPath);
const registry = readJson("docs/curves/registry.json");
const covers = readJson("docs/curves/covers.json");
// This YAML file is intentionally emitted in the documented JSON-compatible subset.
const graph = readJson("docs/curves/cover-links.yaml");
const prior = readJson(
  "research/hyperelliptic_cover_infrastructure_20261006/evidence.json",
);
const transferInventory = readJson("docs/curves/jacobian-transfers.json");

const expectedMaps = {
  "map.binary_f2m7_degree3": {
    icv1Slug: "icv1-f2m7-t13-616700dd",
    modelSha256: "616700ddc6fe573b005e3eda71414a83ce40b069edaf6105b8dccf5d1692a48c",
    coverUid:
      "urn:hc-model:1:sha256:2530bd70db895c112d66c1775146eddba3abce3587853f83f82d3e10cb84d28c",
    mapUid:
      "urn:curve-cover-map:1:sha256:0dc60016b03df7d55ee829880508251fe8261c92a241b76a78987eeaa9787a03",
    degree: 3,
    construction:
      "binary_cubic_pullback_v1: E/K y^2+xy=x^3+1; H/K v^2+u^2v=u^7+u; phi(u,v)=(u^3,u*v+1)",
    typedObjects: [
      {
        id: "EC1N7Ce0hb6d297a2ca08",
        type: "elliptic_curve",
        field: "F_2[z]/(z^7+z+1)",
        dimension_or_genus: 1,
        representation: "binary Weierstrass y^2+xy=x^3+1",
      },
      {
        id: "HC1G3h2530bd70db89",
        type: "hyperelliptic_curve",
        field: "F_2[z]/(z^7+z+1)",
        dimension_or_genus: 3,
        representation: "v^2+u^2v=u^7+u with one rational infinity",
      },
      {
        id: "jac.binary_cubic_pullback.f2m7",
        type: "jacobian",
        field: "F_2[z]/(z^7+z+1)",
        dimension_or_genus: 3,
        representation: "checked reduced Mumford divisors for a one-infinity binary model",
      },
    ],
    coverRecord: {
      f: ["0x0", "0x1", "0x0", "0x0", "0x0", "0x0", "0x0", "0x1"],
      field: { label: "f2m-7-e1a6f45a", modulus: "0x83" },
      form: "v^2+h(u)*v=f(u)",
      genus: 3,
      h: ["0x0", "0x0", "0x1"],
      schema_version: "hyperelliptic-model/v1",
    },
    mapRecord: {
      degree: 3,
      schema_version: "curve-cover-map/v1",
      source_cover_uid:
        "urn:hc-model:1:sha256:2530bd70db895c112d66c1775146eddba3abce3587853f83f82d3e10cb84d28c",
      target_model_json:
        '{"a":"0x0","b":"0x1","field":"f2m-7-e1a6f45a","form":"y^2+xy=x^3+a*x^2+b","modulus":"0x83","v":"1"}',
      target_model_map: null,
      target_model_sha256:
        "616700ddc6fe573b005e3eda71414a83ce40b069edaf6105b8dccf5d1692a48c",
      x: ["0x0", "0x0", "0x0", "0x1"],
      y_0: ["0x1"],
      y_v: ["0x0", "0x1"],
    },
    claims: {
      cover_existence: "supported",
      explicit_construction: "supported",
      curve_map_correctness: "supported",
      jacobian_arithmetic: "supported",
      pullback_point_classes: "supported",
      independent_pullback_homomorphism_oracle: "unknown",
      pushforward_point_divisors: "supported",
      arbitrary_jacobian_pushforward: "refuted",
      composition_on_typed_domain: "supported",
      field_extension_or_descent: "not_applicable",
      computational_advantage: "unknown",
    },
    fullTransferPipeline: "refuted",
  },
  "map.prime_fp827_degree2": {
    icv1Slug: "icv1-fp10-t5-192cb216",
    modelSha256: "192cb216d592d1a2e502bebb1e4154f5d89ed9893de389c2dd0c2b89de94b4e4",
    coverUid:
      "urn:hc-model:1:sha256:9985829083364f51c643cdf27023049e63c9dafbb734ffa87f2c5caf187ca59f",
    mapUid:
      "urn:curve-cover-map:1:sha256:5768ce11909c0f56d84f05e70d71dfb6a46fb175b6f9cc0ed727eec14aef9c8e",
    degree: 2,
    construction:
      "prime_quadratic_pullback_v1: E/F_827 y^2=x^3+x+15; H/F_827 v^2=u^6+u^2+15; phi(u,v)=(u^2,v)",
    typedObjects: [
      {
        id: "EC1P10Cfphcf15d87ab724",
        type: "elliptic_curve",
        field: "F_827",
        dimension_or_genus: 1,
        representation: "short Weierstrass y^2=x^3+x+15",
      },
      {
        id: "HC1G2h998582908336",
        type: "hyperelliptic_curve",
        field: "F_827",
        dimension_or_genus: 2,
        representation: "v^2=u^6+u^2+15 with two rational infinity points",
      },
      {
        id: "jac.prime_quadratic_pullback.fp827",
        type: "jacobian",
        field: "F_827",
        dimension_or_genus: 2,
        representation: "balanced two-infinity divisor representation required but not implemented",
      },
    ],
    coverRecord: {
      f: ["0xf", "0x0", "0x1", "0x0", "0x0", "0x0", "0x1"],
      field: { label: "fp-827", p: "827" },
      form: "v^2+h(u)*v=f(u)",
      genus: 2,
      h: [],
      schema_version: "hyperelliptic-model/v1",
    },
    mapRecord: {
      degree: 2,
      schema_version: "curve-cover-map/v1",
      source_cover_uid:
        "urn:hc-model:1:sha256:9985829083364f51c643cdf27023049e63c9dafbb734ffa87f2c5caf187ca59f",
      target_model_json:
        '{"a":"1","b":"15","field":"fp-827","form":"y^2=x^3+a*x+b","p":"827","v":"1"}',
      target_model_map: null,
      target_model_sha256:
        "192cb216d592d1a2e502bebb1e4154f5d89ed9893de389c2dd0c2b89de94b4e4",
      x: ["0x0", "0x0", "0x1"],
      y_0: [],
      y_v: ["0x1"],
    },
    claims: {
      cover_existence: "supported",
      explicit_construction: "supported",
      curve_map_correctness: "supported",
      jacobian_arithmetic_for_cover: "refuted",
      pullback: "unknown",
      pushforward: "unknown",
      composition: "unknown",
      field_extension_or_descent: "not_applicable",
      subgroup_transfer: "unknown",
      computational_advantage: "unknown",
    },
    fullTransferPipeline: "refuted",
  },
};

const expectedSubgroups = {
  "subgroup.binary_f2m7_order29": {
    mapUid: expectedMaps["map.binary_f2m7_degree3"].mapUid,
    ec1: "EC1N7Ce0hb6d297a2ca08",
    curveUid:
      "urn:ec-record:1:sha256:b6d297a2ca085e9d7b8a6858df7421083e12e740c89b1647660b413e00ed02ff",
    subgroupOrder: "29",
    generator: ["0x26", "0x5"],
    claims: {
      exact_source_order: "supported",
      kernel_intersection: "supported",
      image_order: "supported",
      degree_inverse_recovery: "supported",
      recovery_from_arbitrary_jacobian_result: "unknown",
    },
  },
  "subgroup.prime_fp827_order823": {
    mapUid: expectedMaps["map.prime_fp827_degree2"].mapUid,
    ec1: "EC1P10Cfphcf15d87ab724",
    curveUid:
      "urn:ec-record:1:sha256:cf15d87ab724c8d1fc74d88af320937e4c1cc7e706d4e0e8a0d11734831342a3",
    subgroupOrder: "823",
    generator: ["0x2", "0x336"],
    claims: {
      subgroup_identity: "supported",
      kernel_intersection: "unknown",
      image_order: "unknown",
      degree_inverse_recovery: "unknown",
    },
  },
};

const failures = [];
const check = (condition, message) => {
  if (!condition) failures.push(message);
};
const exact = (actual, expected) => JSON.stringify(actual) === JSON.stringify(expected);
const repoPath = (path) => {
  if (
    typeof path !== "string" ||
    path.length === 0 ||
    path.includes("\\") ||
    path.includes("//") ||
    path.split("/").includes("..") ||
    isAbsolute(path)
  ) {
    return null;
  }
  const absolute = resolve(repo, path);
  const fromRepo = relative(repo, absolute);
  if (
    fromRepo.length === 0 ||
    fromRepo === ".." ||
    fromRepo.startsWith(`..${sep}`) ||
    isAbsolute(fromRepo)
  ) {
    return null;
  }
  if (existsSync(absolute)) {
    const real = realpathSync(absolute);
    const realFromRepo = relative(realpathSync(repo), real);
    if (
      realFromRepo.length === 0 ||
      realFromRepo === ".." ||
      realFromRepo.startsWith(`..${sep}`) ||
      isAbsolute(realFromRepo)
    ) {
      return null;
    }
  }
  return absolute;
};
const uniqueBy = (items, key, label) => {
  const values = items.map((item) => item[key]);
  check(new Set(values).size === values.length, `duplicate ${label}`);
};
const checkClaimMatrix = (claims, expected, label) => {
  uniqueBy(claims, "claim", `${label} claim`);
  const actual = Object.fromEntries(
    claims.map((claim) => [claim.claim, claim.assessment.status]),
  );
  check(exact(actual, expected), `${label}: claim/status matrix changed`);
};
const markdownAnchor = (line) =>
  line
    .replace(/^#+\s+/, "")
    .trim()
    .toLowerCase()
    .replace(/[`*_{}\[\]()<>]/g, "")
    .replace(/[^a-z0-9 _-]/g, "")
    .replace(/\s+/g, "-")
    .replace(/-+/g, "-");
const evidenceRefResolves = (reference) => {
  if (/^https?:\/\//.test(reference)) return false;
  const separator = reference.indexOf("#");
  const path = separator === -1 ? reference : reference.slice(0, separator);
  let fragment = null;
  try {
    fragment =
      separator === -1 ? null : decodeURIComponent(reference.slice(separator + 1));
  } catch {
    return false;
  }
  const absolute = repoPath(path);
  if (!absolute) return false;
  if (!existsSync(absolute)) return false;
  if (!fragment) return true;
  const content = readFileSync(absolute, "utf8");
  if (content.includes(fragment)) return true;
  if (path.endsWith(".md")) {
    return content
      .split(/\r?\n/)
      .filter((line) => /^#+\s+/.test(line))
      .some((line) => markdownAnchor(line) === fragment);
  }
  return false;
};

check(
  assessment.schema_version === "cover-transfer-assessments/v1",
  "wrong assessment schema version",
);
check(
  assessment.catalog_commit === expectedCatalogCommit,
  "assessment catalog commit changed",
);
check(
  sha256(assessmentPath) === expectedAssessmentSha256,
  "assessment sidecar digest changed",
);
check(
  JSON.stringify(assessment.status_vocabulary) ===
    JSON.stringify(["unknown", "supported", "refuted", "not_applicable"]),
  "status vocabulary changed",
);

uniqueBy(assessment.source_artifacts, "path", "source artifact path");
uniqueBy(assessment.map_assessments, "assessment_id", "map assessment ID");
uniqueBy(assessment.map_assessments, "map_uid", "assessed map UID");
uniqueBy(assessment.subgroup_assessments, "assessment_id", "subgroup assessment ID");
uniqueBy(assessment.open_obligations, "id", "open-obligation ID");
check(
  exact(
    assessment.source_artifacts.map((item) => item.path).sort(),
    expectedSourceArtifactPaths,
  ),
  "source-artifact path set changed",
);
check(
  exact(
    assessment.open_obligations.map((item) => item.id).sort(),
    expectedObligationIds,
  ),
  "open-obligation ID set changed",
);
check(
  exact(
    assessment.map_assessments.map((item) => item.assessment_id).sort(),
    Object.keys(expectedMaps).sort(),
  ),
  "assessed map set changed",
);
check(
  exact(
    assessment.subgroup_assessments.map((item) => item.assessment_id).sort(),
    Object.keys(expectedSubgroups).sort(),
  ),
  "assessed subgroup set changed",
);

for (const item of assessment.map_assessments) {
  const expected = expectedMaps[item.assessment_id];
  if (!expected) continue;
  check(item.icv1_slug === expected.icv1Slug, `${item.assessment_id}: ICV1 changed`);
  check(
    item.model_sha256 === expected.modelSha256,
    `${item.assessment_id}: model digest changed`,
  );
  check(item.cover_uid === expected.coverUid, `${item.assessment_id}: cover UID changed`);
  check(item.map_uid === expected.mapUid, `${item.assessment_id}: map UID changed`);
  check(item.degree === expected.degree, `${item.assessment_id}: degree changed`);
  check(
    item.construction === expected.construction,
    `${item.assessment_id}: displayed construction changed`,
  );
  check(
    exact(item.typed_objects, expected.typedObjects),
    `${item.assessment_id}: typed object/equation set changed`,
  );
  checkClaimMatrix(item.claims, expected.claims, item.assessment_id);
  check(
    item.full_transfer_pipeline.status === expected.fullTransferPipeline,
    `${item.assessment_id}: full-pipeline status changed`,
  );
}

for (const item of assessment.subgroup_assessments) {
  const expected = expectedSubgroups[item.assessment_id];
  if (!expected) continue;
  check(item.map_uid === expected.mapUid, `${item.assessment_id}: map UID changed`);
  check(item.ec1 === expected.ec1, `${item.assessment_id}: EC1 changed`);
  check(item.curve_uid === expected.curveUid, `${item.assessment_id}: curve UID changed`);
  check(
    item.subgroup_order === expected.subgroupOrder,
    `${item.assessment_id}: subgroup order changed`,
  );
  check(
    exact(item.generator, expected.generator),
    `${item.assessment_id}: generator changed`,
  );
  checkClaimMatrix(item.claims, expected.claims, item.assessment_id);
}

for (const artifact of assessment.source_artifacts) {
  const absolute = repoPath(artifact.path);
  check(Boolean(absolute), `unsafe source artifact path: ${artifact.path}`);
  check(Boolean(absolute && existsSync(absolute)), `missing ${artifact.path}`);
  if (absolute && existsSync(absolute)) {
    check(
      createHash("sha256").update(readFileSync(absolute)).digest("hex") ===
        artifact.sha256,
      `digest mismatch for ${artifact.path}`,
    );
  }
}

const priorImplementations = prior.records.flatMap(
  (record) => record.implementations ?? [],
);
for (const [path, implementationId] of Object.entries(
  expectedPriorImplementations,
)) {
  const matches = priorImplementations.filter(
    (item) => item.implementation_id === implementationId,
  );
  check(matches.length === 1, `${implementationId}: historical record is not unique`);
  if (matches.length !== 1) continue;
  const implementation = matches[0];
  check(implementation.path === path, `${implementationId}: historical path changed`);
  check(
    sha256(path) === implementation.source_sha256,
    `${implementationId}: current source no longer matches historical snapshot`,
  );
  const sourceArtifact = assessment.source_artifacts.find(
    (item) => item.path === path,
  );
  check(Boolean(sourceArtifact), `${implementationId}: current source pin is absent`);
  if (sourceArtifact) {
    check(
      sourceArtifact.sha256 === implementation.source_sha256,
      `${implementationId}: current and historical source pins disagree`,
    );
  }
}

const graphMapCount = Object.keys(graph.maps).length;
const graphCoverCount = Object.keys(graph.covers).length;
const representations = graph.curves.flatMap((curve) => curve.representations ?? []);
const geometryVerified = covers.curves.filter(
  (curve) => curve.status === "verified_over_declared_field",
).length;
const subgroupNotTested = covers.curves.filter(
  (curve) => curve.subgroup_transfer === "not_tested",
).length;
const nullAdvantages = covers.curves.filter(
  (curve) => curve.dlp_advantage === null,
).length;

check(registry.curves.length === assessment.coverage.catalog_maps, "registry count drift");
check(covers.curves.length === assessment.coverage.catalog_maps, "cover count drift");
check(graph.curves.length === assessment.coverage.catalog_maps, "graph curve count drift");
check(graphMapCount === assessment.coverage.catalog_maps, "graph map count drift");
check(
  graphCoverCount === assessment.coverage.unique_cover_models,
  "unique cover count drift",
);
check(
  geometryVerified === assessment.coverage.geometry_verified_maps,
  "verified geometry count drift",
);
check(
  representations.length === assessment.coverage.ec1_representations,
  "EC1 representation count drift",
);
check(
  subgroupNotTested === assessment.coverage.cover_rows_subgroup_not_tested,
  "not-tested cover-row subgroup count drift",
);
check(
  nullAdvantages === assessment.coverage.catalog_maps,
  "a cover row now carries a DLP-advantage claim",
);
check(
  representations.filter(
    (representation) => representation.subgroup_transport.status === "not_tested",
  ).length === assessment.coverage.ec1_subgroups_not_tested,
  "not-tested EC1 subgroup count drift",
);
check(
  assessment.map_assessments.length === assessment.coverage.assessed_maps,
  "assessed map count mismatch",
);
check(
  graphMapCount - assessment.map_assessments.length ===
    assessment.coverage.unassessed_maps,
  "unassessed map count mismatch",
);
check(
  assessment.subgroup_assessments.length === assessment.coverage.assessed_subgroups,
  "assessed subgroup count mismatch",
);
check(covers.registry_sha256 === sha256("docs/curves/registry.json"), "covers registry pin drift");
check(graph.registry_sha256 === sha256("docs/curves/registry.json"), "graph registry pin drift");
check(transferInventory.summary.curves === 321, "formula-screen curve count changed");
check(transferInventory.summary.routes === 1156, "formula-screen route count changed");
check(
  transferInventory.registry_sha256 !== sha256("docs/curves/registry.json"),
  "formula-screen freshness claim changed; reassess before treating it as stale",
);

for (const item of assessment.map_assessments) {
  const expected = expectedMaps[item.assessment_id];
  const curve = graph.curves.find(
    (candidate) => candidate.icv1_identity.slug === item.icv1_slug,
  );
  check(Boolean(curve), `${item.assessment_id}: missing ICV1 graph row`);
  if (!curve) continue;
  check(
    curve.model_sha256 === item.model_sha256,
    `${item.assessment_id}: model hash mismatch`,
  );
  check(
    curve.cover_links.some(
      (link) => link.cover_uid === item.cover_uid && link.map_uid === item.map_uid,
    ),
    `${item.assessment_id}: cover/map link mismatch`,
  );
  const cover = graph.covers[item.cover_uid];
  const map = graph.maps[item.map_uid];
  check(Boolean(cover), `${item.assessment_id}: missing cover object`);
  check(Boolean(map), `${item.assessment_id}: missing map object`);
  if (cover && expected) {
    check(
      exact(cover.record, expected.coverRecord),
      `${item.assessment_id}: exact HC1 field/equation record changed`,
    );
  }
  if (map) {
    check(
      map.record.source_cover_uid === item.cover_uid,
      `${item.assessment_id}: map source mismatch`,
    );
    check(map.record.degree === item.degree, `${item.assessment_id}: degree mismatch`);
    check(
      map.record.target_model_sha256 === item.model_sha256,
      `${item.assessment_id}: map target mismatch`,
    );
    check(
      map.verification_status === "verified_over_declared_field",
      `${item.assessment_id}: map geometry is no longer verified`,
    );
    if (expected) {
      check(
        exact(map.record, expected.mapRecord),
        `${item.assessment_id}: exact CV1 formula record changed`,
      );
    }
  }
  const catalogRow = covers.curves.find((candidate) => candidate.slug === item.icv1_slug);
  check(Boolean(catalogRow), `${item.assessment_id}: missing covers.json row`);
  if (catalogRow) {
    check(
      item.construction.startsWith(catalogRow.certificate.construction),
      `${item.assessment_id}: construction family mismatch`,
    );
    check(catalogRow.dlp_advantage === null, `${item.assessment_id}: DLP claim changed`);
    check(
      catalogRow.subgroup_transfer === "not_tested",
      `${item.assessment_id}: catalog subgroup status changed`,
    );
  }
  check(
    prior.records.some(
      (record) => record.elliptic_input.icv1.model_sha256 === item.model_sha256,
    ),
    `${item.assessment_id}: no exact prior evidence record`,
  );
}

for (const item of assessment.subgroup_assessments) {
  const curve = graph.curves.find((candidate) =>
    candidate.cover_links.some((link) => link.map_uid === item.map_uid),
  );
  const representation = curve?.representations.find(
    (candidate) => candidate.ec1 === item.ec1,
  );
  check(Boolean(representation), `${item.assessment_id}: missing EC1 representation`);
  if (!representation) continue;
  check(
    representation.curve_uid === item.curve_uid,
    `${item.assessment_id}: curve UID mismatch`,
  );
  check(
    representation.curve.subgroup_order === item.subgroup_order,
    `${item.assessment_id}: subgroup order mismatch`,
  );
  check(
    JSON.stringify(representation.curve.generator) === JSON.stringify(item.generator),
    `${item.assessment_id}: generator mismatch`,
  );
}

const claimStatuses = [
  ...assessment.source_artifacts.map((item) => item.freshness.status),
  ...assessment.map_assessments.flatMap((item) => [
    ...item.claims.map((claim) => claim.assessment.status),
    item.full_transfer_pipeline.status,
  ]),
  ...assessment.subgroup_assessments.flatMap((item) =>
    item.claims.map((claim) => claim.assessment.status),
  ),
  ...assessment.controls.map((item) => item.assessment.status),
  ...assessment.open_obligations.map((item) => item.status),
];
const allowedStatuses = new Set(assessment.status_vocabulary);
check(
  claimStatuses.every((status) => allowedStatuses.has(status)),
  "claim status outside the four-state vocabulary",
);

const controls = Object.fromEntries(
  assessment.controls.map((item) => [item.control, item.assessment.status]),
);
uniqueBy(assessment.controls, "control", "control ID");
check(
  exact(controls, {
    identity: "supported",
    geometry: "supported",
    arithmetic_internal: "supported",
    arithmetic_independent_oracle: "unknown",
    transported_data: "supported",
    coordinate_isomorphism: "unknown",
    subgroup: "supported",
    negative: "supported",
    cross_repository: "supported",
    paired_performance: "not_applicable",
  }),
  "control/status matrix changed",
);

const costs = Object.fromEntries(
  assessment.global_cost_accounting.map((item) => [item.phase, item.status]),
);
uniqueBy(assessment.global_cost_accounting, "phase", "cost phase");
check(
  exact(costs, {
    cover_construction: "unknown",
    failed_construction_searches: "unknown",
    duplicate_candidates: "unknown",
    geometry_verification: "unknown",
    jacobian_arithmetic: "unknown",
    pullback_and_pushforward: "unknown",
    subgroup_recovery: "unknown",
    reusable_preprocessing: "unknown",
    peak_memory: "unknown",
    total_runtime: "unknown",
    matched_generic_baseline: "unknown",
    jacobian_relations_linear_algebra_target_decomposition: "not_applicable",
    end_to_end_speedup: "unknown",
  }),
  "cost-accounting phase/status matrix changed",
);
check(
  assessment.global_cost_accounting.every((item) => item.value === null),
  "an unreviewed numeric cost entered the correctness-only assessment",
);

const evidenceAssessments = [
  ...assessment.source_artifacts.map((item) => item.freshness),
  ...assessment.map_assessments.flatMap((item) => [
    ...item.claims.map((claim) => claim.assessment),
    item.full_transfer_pipeline,
  ]),
  ...assessment.subgroup_assessments.flatMap((item) =>
    item.claims.map((claim) => claim.assessment),
  ),
  ...assessment.controls.map((item) => item.assessment),
];
for (const evidenceAssessment of evidenceAssessments) {
  for (const reference of evidenceAssessment.evidence_refs) {
    check(evidenceRefResolves(reference), `unresolved evidence reference: ${reference}`);
  }
}

check(exact(assessment.report, expectedReportPaths), "report artifact paths changed");
const reportFiles = Object.fromEntries(
  Object.entries(assessment.report).map(([key, path]) => [key, repoPath(path)]),
);
for (const [key, path] of Object.entries(assessment.report)) {
  check(Boolean(reportFiles[key]), `unsafe report artifact path: ${path}`);
  check(
    Boolean(reportFiles[key] && existsSync(reportFiles[key])),
    `missing report artifact ${path}`,
  );
  if (reportFiles[key] && existsSync(reportFiles[key])) {
    const bytes = readFileSync(reportFiles[key]);
    check(
      createHash("sha256").update(bytes).digest("hex") ===
        expectedReportSha256[key],
      `report artifact digest changed: ${path}`,
    );
  }
}

if (reportFiles.pdf && existsSync(reportFiles.pdf)) {
  check(
    readFileSync(reportFiles.pdf).subarray(0, 5).toString("ascii") === "%PDF-",
    "report PDF has the wrong file signature",
  );
}
if (reportFiles.graph_png && existsSync(reportFiles.graph_png)) {
  check(
    readFileSync(reportFiles.graph_png).subarray(0, 8).toString("hex") ===
      "89504e470d0a1a0a",
    "graph PNG has the wrong file signature",
  );
}
if (reportFiles.graph_svg && existsSync(reportFiles.graph_svg)) {
  check(
    /<svg\b/.test(readFileSync(reportFiles.graph_svg, "utf8")),
    "graph SVG is not recognizable SVG",
  );
}
if (reportFiles.graph_source && existsSync(reportFiles.graph_source)) {
  check(
    /\bdigraph\b/.test(readFileSync(reportFiles.graph_source, "utf8")),
    "graph source is not recognizable DOT",
  );
}

const reportText = reportFiles.markdown && existsSync(reportFiles.markdown)
  ? readFileSync(reportFiles.markdown, "utf8")
  : "";
const reportAssessmentDigest = reportText.match(
  /Its SHA-256 is:\s*```text\s*([0-9a-f]{64})\s*```/s,
);
check(Boolean(reportAssessmentDigest), "report does not expose the sidecar SHA-256");
if (reportAssessmentDigest) {
  check(
    reportAssessmentDigest[1] === sha256(assessmentPath),
    "report's embedded sidecar SHA-256 is stale",
  );
}

const releaseReceipt =
  reportFiles.release_receipt && existsSync(reportFiles.release_receipt)
  ? readFileSync(reportFiles.release_receipt, "utf8")
  : "";
for (const receiptToken of [
  "crypto-v0.1.0-main.1371-x86_64-unknown-linux-gnu.tar.gz",
  "77e69e60a68ada60c7b5fe13503e326f50fce4594b96340eca402041e60dcb2d",
  "8f76605ad6fc76e1c32e62c15829c079e24639bba1a17bd1e0935f7e2329cba5",
  "9b98aeb8467877b8947bf1e1c8949b501cf756f2400384dffcae46b2e74c2392",
  "d504eb4f9b078cfc8e51674da49c15124847e7b4a9f42ec0c6169a77189b53d7",
  "12a78b024e60293c48ae72b8c9729683d2ba7103101f371aacad6b26c26a165e",
  "eadd838a628b5b166a90c15972d7f34dc198b65f",
  "9ab7cc8dd314050c56e334af6a116a28648a53fa0281867cebf1f9c79b9b8e3e",
  "hyperelliptic mirror manifest and crypto source: verified",
]) {
  check(releaseReceipt.includes(receiptToken), `release receipt lost ${receiptToken}`);
}

if (failures.length > 0) {
  for (const failure of failures) console.error(`ERROR: ${failure}`);
  process.exit(1);
}

console.log(
  `cover transfer assessment verified: ${assessment.map_assessments.length} maps, ` +
    `${assessment.subgroup_assessments.length} subgroups, ` +
    `${assessment.coverage.unassessed_maps} map rows explicitly unassessed`,
);
