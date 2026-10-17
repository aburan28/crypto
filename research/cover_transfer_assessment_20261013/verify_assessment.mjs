#!/usr/bin/env node

import crypto from "node:crypto";
import { execFileSync } from "node:child_process";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(process.argv[2] ?? path.resolve(here, "../.."));
const curvesDir = path.join(root, "docs", "curves");
const registryPath = path.join(curvesDir, "registry.json");
const reportPath = path.join(curvesDir, "covers.json");
const graphPath = path.join(curvesDir, "cover-links.yaml");
const assessmentPath = path.join(curvesDir, "cover-transfer-assessments.json");

const readJson = (file) => JSON.parse(fs.readFileSync(file, "utf8"));
const registry = readJson(registryPath);
const report = readJson(reportPath);
const graph = readJson(graphPath);
const assessment = readJson(assessmentPath);

function canonicalJson(value) {
  if (Array.isArray(value)) {
    return `[${value.map(canonicalJson).join(",")}]`;
  }
  if (value !== null && typeof value === "object") {
    return `{${Object.keys(value)
      .sort()
      .map((key) => `${JSON.stringify(key)}:${canonicalJson(value[key])}`)
      .join(",")}}`;
  }
  return JSON.stringify(value);
}

function sha256(value) {
  return crypto.createHash("sha256").update(value).digest("hex");
}

function mod(value, modulus) {
  return ((value % modulus) + modulus) % modulus;
}

function oddFieldContext(p, degree, modulusCoefficients) {
  const zero = () => Array(degree).fill(0n);

  function decode(value) {
    if (Array.isArray(value)) {
      if (value.length !== degree) {
        throw new Error(`field vector has ${value.length} coefficients, expected ${degree}`);
      }
      return value.map(BigInt);
    }

    let encoded = BigInt(value);
    const coefficients = zero();
    for (let index = 0; index < degree; index += 1) {
      coefficients[index] = encoded % p;
      encoded /= p;
    }
    if (encoded !== 0n) {
      throw new Error("field element encoding exceeds the declared extension degree");
    }
    return coefficients;
  }

  function add(left, right) {
    return left.map((coefficient, index) => mod(coefficient + right[index], p));
  }

  function multiply(left, right) {
    const product = Array(2 * degree - 1).fill(0n);
    for (let i = 0; i < degree; i += 1) {
      for (let j = 0; j < degree; j += 1) {
        product[i + j] = mod(product[i + j] + left[i] * right[j], p);
      }
    }
    for (let exponent = 2 * degree - 2; exponent >= degree; exponent -= 1) {
      const leading = product[exponent];
      if (leading !== 0n) {
        for (let j = 0; j < degree; j += 1) {
          product[exponent - degree + j] = mod(
            product[exponent - degree + j] - leading * modulusCoefficients[j],
            p,
          );
        }
      }
    }
    return product.slice(0, degree);
  }

  function scale(value, scalar) {
    return value.map((coefficient) => mod(coefficient * scalar, p));
  }

  function equal(left, right) {
    return left.every((coefficient, index) => coefficient === right[index]);
  }

  return {
    zero,
    decode,
    add,
    multiply,
    scale,
    equal,
    one: decode(1n),
  };
}

function binaryMultiply(left, right, degree, modulusPolynomial) {
  let product = 0n;
  let multiplicand = left;
  let multiplier = right;
  while (multiplier !== 0n) {
    if ((multiplier & 1n) !== 0n) {
      product ^= multiplicand;
    }
    multiplier >>= 1n;
    multiplicand <<= 1n;
    if ((multiplicand & (1n << BigInt(degree))) !== 0n) {
      multiplicand ^= modulusPolynomial;
    }
  }
  return product;
}

function isPlainObject(value) {
  return value !== null && typeof value === "object" && !Array.isArray(value);
}

function resolveRepoPath(relativePath) {
  if (
    typeof relativePath !== "string" ||
    relativePath.length === 0 ||
    path.isAbsolute(relativePath) ||
    relativePath.includes("\\") ||
    relativePath.includes("//") ||
    relativePath.split("/").includes("..")
  ) {
    return null;
  }
  const absolute = path.resolve(root, relativePath);
  const fromRoot = path.relative(root, absolute);
  if (
    fromRoot.length === 0 ||
    fromRoot === ".." ||
    fromRoot.startsWith(".." + path.sep) ||
    path.isAbsolute(fromRoot)
  ) {
    return null;
  }
  if (fs.existsSync(absolute)) {
    const realRoot = fs.realpathSync(root);
    const realFile = fs.realpathSync(absolute);
    const realRelative = path.relative(realRoot, realFile);
    if (
      realRelative.length === 0 ||
      realRelative === ".." ||
      realRelative.startsWith(".." + path.sep) ||
      path.isAbsolute(realRelative)
    ) {
      return null;
    }
  }
  return absolute;
}

function markdownAnchor(line) {
  return line
    .replace(/^#+\s+/u, "")
    .trim()
    .toLowerCase()
    .replace(/[\x60*_{}\[\]()<>]/gu, "")
    .replace(/[^a-z0-9 _-]/gu, "")
    .replace(/\s+/gu, "-")
    .replace(/-+/gu, "-");
}

function evidenceReferenceResolves(reference) {
  if (typeof reference !== "string" || reference.length === 0 || /^https?:\/\//u.test(reference)) {
    return false;
  }
  const separator = reference.indexOf("#");
  const relativePath = separator === -1 ? reference : reference.slice(0, separator);
  let fragment = null;
  try {
    fragment = separator === -1 ? null : decodeURIComponent(reference.slice(separator + 1));
  } catch {
    return false;
  }
  const absolute = resolveRepoPath(relativePath);
  if (absolute === null || !fs.existsSync(absolute)) {
    return false;
  }
  if (!fragment) {
    return true;
  }
  const contents = fs.readFileSync(absolute, "utf8");
  if (contents.includes(fragment)) {
    return true;
  }
  if (relativePath.endsWith(".md")) {
    return contents
      .split(/\r?\n/u)
      .filter((line) => /^#+\s+/u.test(line))
      .some((line) => markdownAnchor(line) === fragment);
  }
  return false;
}

function evidenceReferencePath(reference) {
  if (typeof reference !== "string") {
    return null;
  }
  const separator = reference.indexOf("#");
  return separator === -1 ? reference : reference.slice(0, separator);
}

function pdfLiteralMetadata(bytes, key) {
  const source = bytes.toString("latin1");
  const prefix = `/${key}(`;
  const start = source.indexOf(prefix);
  if (start === -1) {
    return null;
  }
  const decoded = [];
  let depth = 1;
  for (let index = start + prefix.length; index < source.length; index += 1) {
    const character = source[index];
    if (character === "\\") {
      const escaped = source[index + 1];
      if (escaped === undefined) {
        return null;
      }
      if (/[0-7]/u.test(escaped)) {
        let octal = escaped;
        while (octal.length < 3 && /[0-7]/u.test(source[index + 1 + octal.length] ?? "")) {
          octal += source[index + 1 + octal.length];
        }
        decoded.push(Number.parseInt(octal, 8));
        index += octal.length;
        continue;
      }
      const simpleEscapes = {
        n: 0x0a,
        r: 0x0d,
        t: 0x09,
        b: 0x08,
        f: 0x0c,
        "(": 0x28,
        ")": 0x29,
        "\\": 0x5c,
      };
      decoded.push(simpleEscapes[escaped] ?? escaped.charCodeAt(0));
      index += 1;
      continue;
    }
    if (character === "(") {
      depth += 1;
    } else if (character === ")") {
      depth -= 1;
      if (depth === 0) {
        const value = Buffer.from(decoded);
        if (value.length >= 2 && value[0] === 0xfe && value[1] === 0xff) {
          const body = value.subarray(2);
          if (body.length % 2 !== 0) {
            return null;
          }
          const littleEndian = Buffer.alloc(body.length);
          for (let offset = 0; offset < body.length; offset += 2) {
            littleEndian[offset] = body[offset + 1];
            littleEndian[offset + 1] = body[offset];
          }
          return littleEndian.toString("utf16le");
        }
        return value.toString("latin1");
      }
    }
    decoded.push(character.charCodeAt(0));
  }
  return null;
}

function sortedObjectKeys(value) {
  return Object.keys(value).sort();
}

function strictStringOrder(left, right) {
  return left < right ? -1 : left > right ? 1 : 0;
}

function sortedRecord(entries) {
  return Object.fromEntries([...entries].sort(([left], [right]) => strictStringOrder(left, right)));
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
  return { label: model.field, p: String(model.p) };
}

function expectedDeclaredFieldBinding(record) {
  const fieldUid = `urn:declared-field:1:sha256:${sha256(canonicalJson(record))}`;
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
    return {
      field_uid: fieldUid,
      record,
      characteristic: String(record.p),
      extension_degree: Number(extension[2]),
      representation: "odd_polynomial_extension",
    };
  }
  return {
    field_uid: fieldUid,
    record,
    characteristic: String(record.p),
    extension_degree: 1,
    representation: "prime",
  };
}

function compareTriples(left, right) {
  for (let index = 0; index < 3; index += 1) {
    const comparison = strictStringOrder(left[index], right[index]);
    if (comparison !== 0) {
      return comparison;
    }
  }
  return 0;
}

function memberDigest(domain, members) {
  return sha256(Buffer.concat([
    Buffer.from(domain, "utf8"),
    Buffer.from([0]),
    Buffer.from(JSON.stringify(members), "utf8"),
  ]));
}

const problems = [];
function check(condition, message) {
  if (!condition) {
    problems.push(message);
  }
}

check(registry.curves.length === 2100, "unexpected registry row count");
check(report.curves.length === 2100, "unexpected cover-report row count");
check(graph.curves.length === 2100, "unexpected graph curve-node count");
check(Object.keys(graph.maps).length === 2100, "unexpected CV1 map count");
check(Object.keys(graph.covers).length === 2098, "unexpected HC1 cover count");

const registryBytes = fs.readFileSync(registryPath);
check(report.registry_sha256 === sha256(registryBytes), "covers.json registry digest mismatch");
check(graph.registry_sha256 === report.registry_sha256, "cover-links registry digest mismatch");

const registryByModelHash = new Map();
const registrySlugs = new Set();
const icv1Identities = new Set();
for (const curve of registry.curves) {
  const modelHash = sha256(curve.model_json);
  check(!registryByModelHash.has(modelHash), `duplicate model hash ${modelHash}`);
  check(!registrySlugs.has(curve.slug), `duplicate registry slug ${curve.slug}`);
  check(!icv1Identities.has(curve.icv1), `duplicate ICV1 identity ${curve.icv1}`);
  registryByModelHash.set(modelHash, curve);
  registrySlugs.add(curve.slug);
  icv1Identities.add(curve.icv1);
}

const reportByModelHash = new Map();
for (const curve of report.curves) {
  const registryCurve = registryByModelHash.get(curve.model_sha256);
  check(registryCurve !== undefined, `cover report references unknown model ${curve.model_sha256}`);
  if (registryCurve !== undefined) {
    check(registryCurve.slug === curve.slug, `cover report slug mismatch for ${curve.slug}`);
  }
  check(!reportByModelHash.has(curve.model_sha256), `duplicate cover report row ${curve.model_sha256}`);
  reportByModelHash.set(curve.model_sha256, curve);
}

const coverReferenceCounts = new Map();
const mapReferenceCounts = new Map();
const ec1Aliases = new Set();
const curveUids = new Set();
const graphModelHashes = new Set();
let representationCount = 0;

for (const graphCurve of graph.curves) {
  check(
    typeof graphCurve.model_sha256 === "string" &&
      !graphModelHashes.has(graphCurve.model_sha256),
    `duplicate or invalid graph model hash ${String(graphCurve.model_sha256)}`,
  );
  graphModelHashes.add(graphCurve.model_sha256);
  const registryCurve = registryByModelHash.get(graphCurve.model_sha256);
  check(registryCurve !== undefined, `graph references unknown model ${graphCurve.model_sha256}`);
  if (registryCurve === undefined) {
    continue;
  }

  check(graphCurve.model_json === registryCurve.model_json, `model JSON mismatch for ${registryCurve.slug}`);
  check(
    graphCurve.icv1_identity.slug === registryCurve.slug &&
      graphCurve.icv1_identity.full === registryCurve.icv1,
    `ICV1 binding mismatch for ${registryCurve.slug}`,
  );
  check(
    graphCurve.representations.length === registryCurve.representations.length,
    `representation count mismatch for ${registryCurve.slug}`,
  );

  for (let index = 0; index < graphCurve.representations.length; index += 1) {
    const graphRepresentation = graphCurve.representations[index];
    const registryRepresentation = registryCurve.representations[index];
    check(
      graphRepresentation.ec1 === registryRepresentation.ec1 &&
        graphRepresentation.curve_uid === registryRepresentation.curve_uid &&
        canonicalJson(graphRepresentation.field) === canonicalJson(registryRepresentation.field) &&
        canonicalJson(graphRepresentation.curve) === canonicalJson(registryRepresentation.curve),
      `representation payload mismatch for ${registryCurve.slug}`,
    );
    check(
      graphRepresentation.subgroup_transport.status === "not_tested" &&
        graphRepresentation.subgroup_transport.evidence_ref === null,
      `unexpected subgroup-transfer status for ${graphRepresentation.ec1}`,
    );
    check(
      !ec1Aliases.has(graphRepresentation.ec1) && !curveUids.has(graphRepresentation.curve_uid),
      `duplicate EC1 representation ${graphRepresentation.ec1}`,
    );
    ec1Aliases.add(graphRepresentation.ec1);
    curveUids.add(graphRepresentation.curve_uid);
    representationCount += 1;
  }

  check(graphCurve.cover_links.length === 1, `expected one cover link for ${registryCurve.slug}`);
  for (const link of graphCurve.cover_links) {
    const map = graph.maps[link.map_uid];
    check(link.direction === "H -> E", `wrong cover direction for ${registryCurve.slug}`);
    check(graph.covers[link.cover_uid] !== undefined, `missing HC1 ${link.cover_uid}`);
    check(map !== undefined, `missing CV1 ${link.map_uid}`);
    if (map !== undefined) {
      check(
        map.record.source_cover_uid === link.cover_uid &&
          map.record.target_model_sha256 === graphCurve.model_sha256 &&
          map.record.target_model_json === graphCurve.model_json,
        `CV1 binding mismatch for ${registryCurve.slug}`,
      );
    }
    coverReferenceCounts.set(link.cover_uid, (coverReferenceCounts.get(link.cover_uid) ?? 0) + 1);
    mapReferenceCounts.set(link.map_uid, (mapReferenceCounts.get(link.map_uid) ?? 0) + 1);
  }
}

check(
  canonicalJson([...reportByModelHash.keys()].sort(strictStringOrder)) ===
    canonicalJson([...registryByModelHash.keys()].sort(strictStringOrder)),
  "cover-report model hashes are not exactly the registry model hashes",
);
check(
  canonicalJson([...graphModelHashes].sort(strictStringOrder)) ===
    canonicalJson([...registryByModelHash.keys()].sort(strictStringOrder)),
  "graph model hashes are not a unique exact copy of the registry model hashes",
);

check(representationCount === 664, "unexpected EC1 representation count");
check(ec1Aliases.size === 664, "unexpected unique EC1 count");
check(curveUids.size === 664, "unexpected unique curve-UID count");
check(
  mapReferenceCounts.size === 2100 && [...mapReferenceCounts.values()].every((count) => count === 1),
  "CV1 references are not bijective with curve nodes",
);
check(
  coverReferenceCounts.size === 2098 &&
    [...coverReferenceCounts.values()].every((count) => count === 1 || count === 2) &&
    [...coverReferenceCounts.values()].filter((count) => count === 2).length === 2,
  "unexpected HC1 sharing pattern",
);

for (const [uid, cover] of Object.entries(graph.covers)) {
  check(uid === cover.cover_uid, `HC1 key/payload mismatch ${uid}`);
  check(uid.endsWith(sha256(canonicalJson(cover.record))), `HC1 content hash mismatch ${uid}`);
  check(coverReferenceCounts.has(uid), `unreferenced HC1 ${uid}`);
}

for (const [uid, map] of Object.entries(graph.maps)) {
  check(uid === map.map_uid, `CV1 key/payload mismatch ${uid}`);
  check(uid.endsWith(sha256(canonicalJson(map.record))), `CV1 content hash mismatch ${uid}`);
  check(mapReferenceCounts.get(uid) === 1, `unexpected CV1 reference count ${uid}`);
}

const graphCurveByModelHash = new Map(graph.curves.map((curve) => [curve.model_sha256, curve]));
const arithmeticChecks = { binary: 0, prime: 0, extension: 0, normalized_target: 0 };
const constructionCounts = { binary_cubic_pullback_v1: 0, prime_quadratic_pullback_v1: 0 };

for (const coverResult of report.curves) {
  const registryCurve = registryByModelHash.get(coverResult.model_sha256);
  if (registryCurve === undefined) {
    continue;
  }
  const targetModel = JSON.parse(registryCurve.model_json);
  const certificate = coverResult.certificate;
  constructionCounts[certificate.construction] =
    (constructionCounts[certificate.construction] ?? 0) + 1;

  check(
    coverResult.status === "verified_over_declared_field" && coverResult.exists === true,
    `unexpected cover status for ${coverResult.slug}`,
  );
  check(
    coverResult.field_relation === "same_field" && coverResult.map_direction === "H -> E",
    `unexpected field relation or direction for ${coverResult.slug}`,
  );
  check(
    coverResult.subfield_descent === "not_tested" &&
      coverResult.subgroup_transfer === "not_tested" &&
      coverResult.dlp_advantage === null &&
      coverResult.minimal_degree === null &&
      coverResult.minimal_genus === null,
    `claim-boundary mismatch for ${coverResult.slug}`,
  );

  const graphCurve = graphCurveByModelHash.get(coverResult.model_sha256);
  if (graphCurve === undefined || graphCurve.cover_links.length !== 1) {
    continue;
  }
  const link = graphCurve.cover_links[0];
  const graphMap = graph.maps[link.map_uid];
  const graphCover = graph.covers[link.cover_uid];
  if (graphMap === undefined || graphCover === undefined) {
    continue;
  }

  check(graphCover.record.field.label === targetModel.field, `HC1 field mismatch for ${coverResult.slug}`);
  check(
    graphCover.record.genus === certificate.genus &&
      canonicalJson(graphCover.record.f) === canonicalJson(certificate.f) &&
      canonicalJson(graphCover.record.h) === canonicalJson(certificate.h),
    `HC1/report certificate mismatch for ${coverResult.slug}`,
  );
  check(
    graphMap.record.degree === certificate.degree &&
      canonicalJson(graphMap.record.x) === canonicalJson(certificate.x) &&
      canonicalJson(graphMap.record.y_v) === canonicalJson(certificate.y_v) &&
      canonicalJson(graphMap.record.y_0) === canonicalJson(certificate.y_0) &&
      canonicalJson(graphMap.record.target_model_map) ===
        canonicalJson(certificate.target_model_map ?? null),
    `CV1/report certificate mismatch for ${coverResult.slug}`,
  );

  if (certificate.construction === "binary_cubic_pullback_v1") {
    const degree = Number(targetModel.field.match(/^f2m-(\d+)-/u)?.[1]);
    const modulusPolynomial = BigInt(targetModel.modulus);
    const a = BigInt(targetModel.a);
    const b = BigInt(targetModel.b);
    const d = BigInt(certificate.y_0[0]);

    check(
      targetModel.form === "y^2+xy=x^3+a*x^2+b" &&
        targetModel.field.startsWith("f2m-") &&
        b !== 0n,
      `invalid binary constructor input for ${coverResult.slug}`,
    );
    check(
      certificate.degree === 3 &&
        certificate.genus === 3 &&
        canonicalJson(certificate.h) === canonicalJson(["0x0", "0x0", "0x1"]) &&
        canonicalJson(certificate.x) === canonicalJson(["0x0", "0x0", "0x0", "0x1"]) &&
        canonicalJson(certificate.y_v) === canonicalJson(["0x0", "0x1"]) &&
        certificate.target_model_map === undefined,
      `invalid binary map shape for ${coverResult.slug}`,
    );
    check(
      certificate.f.length === 8 &&
        certificate.f[0] === "0x0" &&
        BigInt(certificate.f[1]) === d &&
        certificate.f[2] === "0x0" &&
        certificate.f[3] === "0x0" &&
        BigInt(certificate.f[4]) === a &&
        certificate.f[5] === "0x0" &&
        certificate.f[6] === "0x0" &&
        certificate.f[7] === "0x1",
      `invalid binary source polynomial for ${coverResult.slug}`,
    );
    check(
      binaryMultiply(d, d, degree, modulusPolynomial) === b,
      `binary square-root identity failed for ${coverResult.slug}`,
    );
    arithmeticChecks.binary += 1;
    continue;
  }

  if (certificate.construction !== "prime_quadratic_pullback_v1") {
    check(false, `unknown construction ${certificate.construction} for ${coverResult.slug}`);
    continue;
  }

  const shortModel = certificate.target_model_map?.source_model ?? targetModel;
  const p = BigInt(shortModel.p);
  const extensionDegree = shortModel.k === undefined ? 1 : Number(shortModel.k);
  const modulusCoefficients = shortModel.modulus?.map(BigInt) ?? [0n];
  const field = oddFieldContext(p, extensionDegree, modulusCoefficients);
  const a = field.decode(shortModel.a);
  const b = field.decode(shortModel.b);
  const c = field.decode(BigInt(certificate.x[0]));
  const cSquared = field.multiply(c, c);
  const cCubed = field.multiply(cSquared, c);
  const fAtC = field.add(field.add(cCubed, field.multiply(a, c)), b);
  const discriminantTerm = field.add(
    field.scale(field.multiply(field.multiply(a, a), a), 4n),
    field.scale(field.multiply(b, b), 27n),
  );
  const expectedF = [
    fAtC,
    field.zero(),
    field.add(field.scale(cSquared, 3n), a),
    field.zero(),
    field.scale(c, 3n),
    field.zero(),
    field.one,
  ];
  const actualF = certificate.f.map((coefficient) => field.decode(coefficient));

  check(p > 3n && (p & 1n) === 1n, `invalid odd characteristic for ${coverResult.slug}`);
  check(
    certificate.degree === 2 &&
      certificate.genus === 2 &&
      certificate.h.length === 0 &&
      canonicalJson(certificate.y_v) === canonicalJson(["0x1"]) &&
      certificate.y_0.length === 0 &&
      certificate.x.length === 3 &&
      ["0x0", "0x1", "0x2", "0x3"].includes(certificate.x[0]) &&
      certificate.x[1] === "0x0" &&
      certificate.x[2] === "0x1",
    `invalid odd-characteristic map shape for ${coverResult.slug}`,
  );
  check(!field.equal(fAtC, field.zero()), `selected shift is a root for ${coverResult.slug}`);
  check(!field.equal(discriminantTerm, field.zero()), `singular short target for ${coverResult.slug}`);
  check(
    actualF.length === expectedF.length &&
      actualF.every((coefficient, index) => field.equal(coefficient, expectedF[index])),
    `F(u^2+c) coefficient replay failed for ${coverResult.slug}`,
  );

  const selectedShift = Number(BigInt(certificate.x[0]));
  for (let priorShift = 0; priorShift < selectedShift; priorShift += 1) {
    const candidate = field.decode(BigInt(priorShift));
    const candidateSquared = field.multiply(candidate, candidate);
    const candidateCubed = field.multiply(candidateSquared, candidate);
    const value = field.add(field.add(candidateCubed, field.multiply(a, candidate)), b);
    check(field.equal(value, field.zero()), `shift is not the first nonroot for ${coverResult.slug}`);
  }

  if (extensionDegree > 1) {
    check(
      targetModel.form === "y^2=x^3+a*x+b" &&
        targetModel.field.startsWith("fpk-") &&
        targetModel.modulus.length === extensionDegree &&
        certificate.target_model_map === undefined,
      `invalid extension-field target for ${coverResult.slug}`,
    );
    arithmeticChecks.extension += 1;
  } else {
    if (targetModel.form === "y^2=x^3+a*x+b") {
      check(certificate.target_model_map === undefined, `unexpected target model map for ${coverResult.slug}`);
    } else {
      const targetMap = certificate.target_model_map;
      check(
        targetMap !== undefined &&
          targetMap.degree === 1 &&
          targetMap.status === "verified_algebraic_identity" &&
          canonicalJson(targetMap.target_model) === canonicalJson(targetModel) &&
          targetMap.source_model.form === "y^2=x^3+a*x+b" &&
          targetMap.source_model.field === targetModel.field,
        `invalid birational target-model map for ${coverResult.slug}`,
      );
      if (targetMap !== undefined) {
        const fieldPrime = BigInt(targetModel.p);
        const montgomeryA = BigInt(targetMap.montgomery_A);
        const montgomeryB = BigInt(targetMap.montgomery_B);
        const shift = BigInt(targetMap.shift_A_over_3);
        const scale = BigInt(targetMap.edwards_scale);
        const shortA = BigInt(targetMap.source_model.a);
        const shortB = BigInt(targetMap.source_model.b);
        const square = (value) => mod(value * value, fieldPrime);
        const cube = (value) => mod(value * value * value, fieldPrime);

        check(
          targetMap.schema_version === "elliptic-model-map/v1" &&
            targetModel.field === `fp-${fieldPrime}` &&
            targetMap.source_model.p === targetModel.p &&
            targetMap.source_model.field === targetModel.field &&
            targetMap.source_model.form === "y^2=x^3+a*x+b" &&
            targetMap.exceptional_points ===
              "unique extension of the birational map to smooth projective models",
          `target normalization metadata mismatch for ${coverResult.slug}`,
        );
        check(
          montgomeryA >= 0n &&
            montgomeryA < fieldPrime &&
            montgomeryB > 0n &&
            montgomeryB < fieldPrime &&
            shift >= 0n &&
            shift < fieldPrime &&
            scale > 0n &&
            scale < fieldPrime &&
            mod(3n * shift, fieldPrime) === montgomeryA,
          `target normalization parameter range/shift mismatch for ${coverResult.slug}`,
        );
        check(
          mod(square(montgomeryB) * shortA, fieldPrime) ===
              mod(1n - montgomeryA * shift, fieldPrime) &&
            mod(cube(montgomeryB) * shortB, fieldPrime) ===
              mod(2n * cube(shift) - shift, fieldPrime),
          `short/Montgomery coefficient identity failed for ${coverResult.slug}`,
        );

        if (targetModel.form === "B*y^2=x^3+A*x^2+x") {
          check(
            targetMap.kind === "montgomery_from_short" &&
              montgomeryA === BigInt(targetModel.A) &&
              montgomeryB === BigInt(targetModel.B) &&
              scale === 1n &&
              targetMap.forward.x === "B*X-A/3" &&
              targetMap.forward.y === "B*Y" &&
              targetMap.inverse.X === "(x+A/3)/B" &&
              targetMap.inverse.Y === "y/B" &&
              mod(square(montgomeryA) - 4n, fieldPrime) !== 0n,
            `Montgomery normalization replay failed for ${coverResult.slug}`,
          );
        } else {
          const twisted = targetModel.form === "a*x^2+y^2=1+d*x^2*y^2";
          const scaled =
            targetModel.form === "x^2+y^2=c^2*(1+d*x^2*y^2)";
          const targetScale = twisted ? 1n : BigInt(targetModel.c);
          const edwardsA = twisted ? BigInt(targetModel.a) : 1n;
          const rawD = BigInt(targetModel.d);
          const edwardsD = scaled
            ? mod(rawD * square(square(targetScale)), fieldPrime)
            : rawD;
          const denominator = mod(edwardsA - edwardsD, fieldPrime);
          check(
            (twisted || scaled) &&
              targetMap.kind === "edwards_from_short" &&
              scale === targetScale &&
              edwardsA > 0n &&
              edwardsA < fieldPrime &&
              rawD > 0n &&
              rawD < fieldPrime &&
              denominator !== 0n &&
              mod(denominator * montgomeryA, fieldPrime) ===
                mod(2n * (edwardsA + edwardsD), fieldPrime) &&
              mod(denominator * montgomeryB, fieldPrime) === 4n % fieldPrime &&
              targetMap.forward.x === "c*(B*X-A/3)/(B*Y)" &&
              targetMap.forward.y === "c*(B*X-A/3-1)/(B*X-A/3+1)" &&
              targetMap.inverse.X === "((c+y)/(c-y)+A/3)/B" &&
              targetMap.inverse.Y === "c*(c+y)/(B*x*(c-y))",
            `Edwards normalization replay failed for ${coverResult.slug}`,
          );
        }
        arithmeticChecks.normalized_target += 1;
      }
    }
    arithmeticChecks.prime += 1;
  }
}

check(constructionCounts.binary_cubic_pullback_v1 === 160, "unexpected binary-construction count");
check(constructionCounts.prime_quadratic_pullback_v1 === 1940, "unexpected odd-construction count");
check(arithmeticChecks.binary === 160, "unexpected binary arithmetic-replay count");
check(arithmeticChecks.prime === 616, "unexpected prime-field arithmetic-replay count");
check(arithmeticChecks.extension === 1324, "unexpected extension-field arithmetic-replay count");
check(
  arithmeticChecks.normalized_target === 30,
  "unexpected normalized-target algebraic-replay count",
);

const expectedStatusVocabulary = ["unknown", "supported", "refuted", "not_applicable"];
const expectedCatalogCommit = "ff48e2800f62b45c746b3f391540c07b16af4e75";
const expectedRootKeys = [
  "assessment_date",
  "catalog_commit",
  "claim_contract",
  "controls",
  "coverage",
  "exact_overlays",
  "exploration_boundary",
  "global_cost_accounting",
  "memberships",
  "typed_bindings",
  "migration",
  "open_obligations",
  "profiles",
  "report",
  "result_templates",
  "schema_version",
  "source_artifacts",
  "status_vocabulary",
].sort();
const expectedMapClaimIds = [
  "alternate_target_model_map",
  "arbitrary_jacobian_pushforward",
  "composition_on_typed_domain",
  "computational_advantage",
  "cover_existence",
  "curve_map_correctness",
  "explicit_construction",
  "field_extension_or_descent",
  "independent_pullback_homomorphism_oracle",
  "jacobian_arithmetic_for_cover",
  "pullback_point_classes",
  "pushforward_point_divisors",
].sort();
const expectedSubgroupClaimIds = [
  "catalog_order_annihilates_generator",
  "composition_on_generator",
  "computational_advantage",
  "degree_inverse_recovery",
  "gcd_degree_catalog_order",
  "generator_on_curve",
  "image_annihilated_by_catalog_order",
  "image_order_equals_actual_source_order",
  "image_order_equals_catalog_order",
  "kernel_intersection_on_generated_subgroup",
  "numeric_exact_source_order",
  "pullback_valid",
  "recovery_from_arbitrary_jacobian_result",
  "subgroup_identity",
].sort();
const expectedTemplateStatuses = {
  "map.advantage.unknown": "unknown",
  "map.binary-arbitrary-pushforward.refuted": "refuted",
  "map.binary-composition.supported": "supported",
  "map.binary-construction.supported": "supported",
  "map.binary-jacobian.supported": "supported",
  "map.binary-pipeline.refuted": "refuted",
  "map.binary-pullback.supported": "supported",
  "map.binary-typed-pushforward.supported": "supported",
  "map.direct-target.not-applicable": "not_applicable",
  "map.extension-construction.supported": "supported",
  "map.extension-descent.unknown": "unknown",
  "map.f2m7-geometry.supported": "supported",
  "map.fp827-geometry.supported": "supported",
  "map.geometry.catalog-supported": "supported",
  "map.independent-oracle.unknown": "unknown",
  "map.prime-construction.supported": "supported",
  "map.prime-jacobian.refuted": "refuted",
  "map.prime-normalization.supported": "supported",
  "map.prime-pipeline.refuted": "refuted",
  "map.prime-transfer.unknown": "unknown",
  "map.same-field.not-applicable": "not_applicable",
  "subgroup.arbitrary-recovery.unknown": "unknown",
  "subgroup.binary-annihilation.supported": "supported",
  "subgroup.binary-composition.supported": "supported",
  "subgroup.binary-generator.supported": "supported",
  "subgroup.binary-image-actual-order.unknown": "unknown",
  "subgroup.binary-image-annihilation.supported": "supported",
  "subgroup.binary-injective-generated.unknown": "unknown",
  "subgroup.binary-inverse-recovery.supported": "supported",
  "subgroup.binary-pullback.supported": "supported",
  "subgroup.degree-coprime.supported": "supported",
  "subgroup.exact-image-order29.supported": "supported",
  "subgroup.exact-order29.supported": "supported",
  "subgroup.identity.supported": "supported",
  "subgroup.image-numeric-order.unknown": "unknown",
  "subgroup.numeric-order.unknown": "unknown",
  "subgroup.prime-generator.unknown": "unknown",
  "subgroup.prime-transfer.unknown": "unknown",
};
const expectedProfiles = {
  "map.binary-cubic.same-field": {
    kind: "map_base",
    predicate_id: "binary_cubic_same_field",
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
  },
  "map.odd-prime-quadratic.direct": {
    kind: "map_base",
    predicate_id: "odd_prime_quadratic_direct_same_field",
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
  },
  "map.odd-prime-quadratic.normalized": {
    kind: "map_base",
    predicate_id: "odd_prime_quadratic_normalized_same_field",
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
  },
  "map.odd-extension-quadratic.direct": {
    kind: "map_base",
    predicate_id: "odd_extension_quadratic_direct_same_field",
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
  },
  "subgroup.binary-degree3": {
    kind: "subgroup_base",
    predicate_id: "linked_subgroup_degree3",
    claims: {
      subgroup_identity: "subgroup.identity.supported",
      gcd_degree_catalog_order: "subgroup.degree-coprime.supported",
      generator_on_curve: "subgroup.binary-generator.supported",
      catalog_order_annihilates_generator: "subgroup.binary-annihilation.supported",
      numeric_exact_source_order: "subgroup.numeric-order.unknown",
      pullback_valid: "subgroup.binary-pullback.supported",
      image_annihilated_by_catalog_order: "subgroup.binary-image-annihilation.supported",
      composition_on_generator: "subgroup.binary-composition.supported",
      kernel_intersection_on_generated_subgroup: "subgroup.binary-injective-generated.unknown",
      image_order_equals_actual_source_order: "subgroup.binary-image-actual-order.unknown",
      image_order_equals_catalog_order: "subgroup.image-numeric-order.unknown",
      degree_inverse_recovery: "subgroup.binary-inverse-recovery.supported",
      recovery_from_arbitrary_jacobian_result: "subgroup.arbitrary-recovery.unknown",
      computational_advantage: "map.advantage.unknown",
    },
  },
  "subgroup.odd-degree2": {
    kind: "subgroup_base",
    predicate_id: "linked_subgroup_degree2",
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
      recovery_from_arbitrary_jacobian_result: "subgroup.arbitrary-recovery.unknown",
      computational_advantage: "map.advantage.unknown",
    },
  },
};
const expectedControlStatuses = {
  identity: ["supported", null],
  geometry: ["supported", null],
  arithmetic_internal: ["supported", null],
  arithmetic_independent_oracle: ["unknown", "independent_map_oracle"],
  binary_catalog_subgroups: ["supported", null],
  odd_degree2_transfer: ["refuted", null],
  published_artifact: ["refuted", "catalog_wide_packaged_binary"],
  cross_repository_core_parity: ["supported", null],
  current_fixed_f8_cross_repository_parity: [
    "unknown",
    "current_fixed_f8_cryptanalysis_mirror",
  ],
  coordinate_controls: ["unknown", "coordinate_controls"],
  paired_performance: ["not_applicable", null],
};
const expectedSemanticDigests = {
  claim_contract: "36a2bb9524810994c02633462e71ae08ac88da6fbb520df72e66e05150628634",
  result_templates: "74d48723d80213cef7032b499e2d36a051ebed3c23e03bfce5819ad80f5edb7c",
  profiles: "807800e78066e8c3b0caa6634072a585cf4ae924a40b13817d9ca7a3c0aa3710",
  controls: "6ef9bbf2eed4dd04ff0025f06c6d6b73e62519c54519e356c0dc7e3c22eb1850",
  open_obligations: "142a997795d925f10923d7fc4adb80bfa886ca7d5c46f5a911529d8538189629",
  global_cost_accounting:
    "b3b77a401d4dd76cabfc9a65d5a738fcfb977fd1a3abc2db3733c41820304619",
};
const expectedCostStatuses = {
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
  end_to_end_speedup: "unknown",
  jacobian_relations_linear_algebra_target_decomposition: "not_applicable",
};
const requiredSourceArtifactPaths = [
  "Cargo.toml",
  "Cargo.lock",
  ".github/workflows/curve-covers.yml",
  ".github/workflows/release.yml",
  "docs/curves/registry.json",
  "docs/curves/covers.json",
  "docs/curves/cover-links.yaml",
  "docs/curves/cover-transfer-assessment.schema.json",
  "docs/curves/fixed-f8-genus3-order-certificate.json",
  "scripts/build_cover_transfer_assessments.mjs",
  "src/bin/curve_cover_check/checker.rs",
  "src/bin/curve_cover_check/models.rs",
  "src/bin/curve_cover_check/links.rs",
  "src/bin/curve_cover_check.rs",
  "src/bin/hyperelliptic_cover.rs",
  "src/bin/ecdlp_audit/fixed_f8_order.rs",
  "src/lib.rs",
  "src/binary_ecc/mod.rs",
  "src/binary_ecc/f2m.rs",
  "src/binary_ecc/poly_f2m.rs",
  "src/binary_ecc/curve.rs",
  "src/binary_ecc/cover_transfer.rs",
  "src/binary_ecc/hyperelliptic.rs",
  "src/cryptanalysis/mod.rs",
  "src/cryptanalysis/ecc2k130_guard.rs",
  "src/hash/mod.rs",
  "src/hash/sha256.rs",
  "src/prime_hyperelliptic/mod.rs",
  "src/prime_hyperelliptic/fp2.rs",
  "src/prime_hyperelliptic/fp_poly.rs",
  "src/prime_hyperelliptic/curve.rs",
  "src/utils/mod.rs",
  "research/cover_transfer_assessment_20261013/verify_assessment.mjs",
  "research/cover_transfer_assessment_20261013/typed_transfer_context.dot",
  "research/cover_transfer_assessment_20261013/typed_transfer_context.svg",
  "research/cover_transfer_assessment_20261013/typed_transfer_context.png",
];
const closedMapPredicates = [
  "binary_cubic_same_field",
  "odd_prime_quadratic_direct_same_field",
  "odd_prime_quadratic_normalized_same_field",
  "odd_extension_quadratic_direct_same_field",
].sort();
const closedSubgroupPredicates = [
  "linked_subgroup_degree3",
  "linked_subgroup_degree2",
].sort();
const mapMembershipDomain = "cover-transfer-assessments/v2/map_base";
const subgroupMembershipDomain = "cover-transfer-assessments/v2/subgroup_base";

check(
  assessment.schema_version === "cover-transfer-assessments/v2",
  "assessment schema version is not cover-transfer-assessments/v2",
);
check(
  canonicalJson(sortedObjectKeys(assessment)) === canonicalJson(expectedRootKeys),
  "assessment root key set changed",
);
check(
  assessment.catalog_commit === expectedCatalogCommit,
  "assessment catalog commit changed",
);
for (const relativePath of [
  "docs/curves/registry.json",
  "docs/curves/covers.json",
  "docs/curves/cover-links.yaml",
]) {
  try {
    const committedBytes = execFileSync(
      "git",
      ["show", `${assessment.catalog_commit}:${relativePath}`],
      { cwd: root, encoding: "buffer", maxBuffer: 32 * 1024 * 1024 },
    );
    check(
      sha256(committedBytes) === sha256(fs.readFileSync(path.join(root, relativePath))),
      `catalog file ${relativePath} does not match catalog_commit ${assessment.catalog_commit}`,
    );
  } catch (error) {
    check(
      false,
      `cannot resolve catalog_commit ${assessment.catalog_commit} for ${relativePath}: ${error.message}`,
    );
  }
}
check(
  canonicalJson(assessment.status_vocabulary) === canonicalJson(expectedStatusVocabulary),
  "assessment status vocabulary changed",
);
for (const [section, expectedDigest] of Object.entries(expectedSemanticDigests)) {
  check(
    sha256(canonicalJson(assessment[section])) === expectedDigest,
    `authoritative semantic digest mismatch for ${section}`,
  );
}

const allowedStatuses = new Set(expectedStatusVocabulary);
const sourceArtifacts = Array.isArray(assessment.source_artifacts)
  ? assessment.source_artifacts
  : [];
check(Array.isArray(assessment.source_artifacts), "source_artifacts is not an array");
const sourceArtifactPaths = new Set();
for (const artifact of sourceArtifacts) {
  if (!isPlainObject(artifact)) {
    check(false, "source_artifacts contains a non-object");
    continue;
  }
  check(
    typeof artifact.path === "string" && !sourceArtifactPaths.has(artifact.path),
    "duplicate or invalid source-artifact path",
  );
  sourceArtifactPaths.add(artifact.path);
  const absolute = resolveRepoPath(artifact.path);
  check(absolute !== null, "unsafe source-artifact path " + String(artifact.path));
  check(absolute !== null && fs.existsSync(absolute), "missing source artifact " + String(artifact.path));
  if (absolute !== null && fs.existsSync(absolute)) {
    check(
      /^[0-9a-f]{64}$/u.test(artifact.sha256) &&
        sha256(fs.readFileSync(absolute)) === artifact.sha256,
      "source-artifact digest mismatch for " + artifact.path,
    );
  }
}
for (const requiredPath of requiredSourceArtifactPaths) {
  check(sourceArtifactPaths.has(requiredPath), "missing required source pin " + requiredPath);
}

const checkerSourceDigest = sha256(
  Buffer.concat(
    [
      "src/bin/curve_cover_check/checker.rs",
      "src/bin/curve_cover_check/models.rs",
      "src/bin/curve_cover_check/links.rs",
      "src/bin/curve_cover_check.rs",
    ].map((relativePath) => fs.readFileSync(path.join(root, relativePath))),
  ),
);
check(
  /^[0-9a-f]{64}$/u.test(report.checker_source_sha256) &&
    report.checker_source_sha256 === checkerSourceDigest,
  "cover report checker_source_sha256 does not match the complete native checker source",
);

const catalogProducerSourceDigest = sha256(
  Buffer.concat(
    [
      "src/bin/curve_cover_check/checker.rs",
      "src/bin/curve_cover_check/models.rs",
      "src/bin/hyperelliptic_cover.rs",
    ].map((relativePath) => fs.readFileSync(path.join(root, relativePath))),
  ),
);
const infrastructureDocumentation = fs.readFileSync(
  path.join(curvesDir, "HYPERELLIPTIC_INFRASTRUCTURE.md"),
  "utf8",
);
check(
  /checker source,\s+model\s+normalization source/u.test(infrastructureDocumentation) &&
    infrastructureDocumentation.includes(`\`${catalogProducerSourceDigest}\``),
  "catalog-export producer-source definition or digest is stale",
);

const resultTemplates = isPlainObject(assessment.result_templates)
  ? assessment.result_templates
  : {};
const profiles = isPlainObject(assessment.profiles) ? assessment.profiles : {};
const claimContract = isPlainObject(assessment.claim_contract)
  ? assessment.claim_contract
  : {};
const mapClaimContract = isPlainObject(claimContract.map) ? claimContract.map : {};
const subgroupClaimContract = isPlainObject(claimContract.subgroup)
  ? claimContract.subgroup
  : {};
const mapClaimIds = sortedObjectKeys(mapClaimContract);
const subgroupClaimIds = sortedObjectKeys(subgroupClaimContract);
check(isPlainObject(assessment.result_templates), "result_templates is not an object");
check(isPlainObject(assessment.profiles), "profiles is not an object");
check(isPlainObject(assessment.claim_contract), "claim_contract is not an object");
check(mapClaimIds.length > 0, "map claim contract is empty");
check(subgroupClaimIds.length > 0, "subgroup claim contract is empty");
check(
  canonicalJson(mapClaimIds) === canonicalJson(expectedMapClaimIds),
  "map claim contract is not the closed authoritative claim set",
);
check(
  canonicalJson(subgroupClaimIds) === canonicalJson(expectedSubgroupClaimIds),
  "subgroup claim contract is not the closed authoritative claim set",
);
for (const [kind, contract] of [
  ["map", mapClaimContract],
  ["subgroup", subgroupClaimContract],
]) {
  for (const [claimId, description] of Object.entries(contract)) {
    check(
      typeof claimId === "string" &&
        claimId.length > 0 &&
        typeof description === "string" &&
        description.trim().length > 0,
      kind + " claim contract has an invalid entry",
    );
  }
}

const obligations = Array.isArray(assessment.open_obligations)
  ? assessment.open_obligations
  : [];
check(Array.isArray(assessment.open_obligations), "open_obligations is not an array");
const obligationIds = new Set();
for (const obligation of obligations) {
  check(isPlainObject(obligation), "open_obligations contains a non-object");
  if (!isPlainObject(obligation)) {
    continue;
  }
  check(
    typeof obligation.id === "string" &&
      obligation.id.length > 0 &&
      !obligationIds.has(obligation.id),
    "duplicate or invalid open-obligation ID",
  );
  obligationIds.add(obligation.id);
  check(obligation.status === "unknown", "open obligation is not unknown: " + String(obligation.id));
}

const usedTemplateIds = new Set();
const referencedObligationIds = new Set();
check(
  canonicalJson(sortedObjectKeys(resultTemplates)) ===
    canonicalJson(sortedObjectKeys(expectedTemplateStatuses)),
  "result-template set changed from the authoritative verifier matrix",
);
for (const [templateId, expectedStatus] of Object.entries(expectedTemplateStatuses)) {
  check(
    resultTemplates[templateId]?.status === expectedStatus,
    `authoritative status mismatch for result template ${templateId}`,
  );
}
for (const [templateId, template] of Object.entries(resultTemplates)) {
  check(templateId.length > 0, "empty result-template ID");
  check(isPlainObject(template), "result template is not an object: " + templateId);
  if (!isPlainObject(template)) {
    continue;
  }
  check(allowedStatuses.has(template.status), "invalid template status: " + templateId);
  check(
    typeof template.scope === "string" && template.scope.trim().length > 0,
    "template scope is empty: " + templateId,
  );
  check(
    typeof template.rationale === "string" && template.rationale.trim().length > 0,
    "template rationale is empty: " + templateId,
  );
  const evidenceRefs = Array.isArray(template.evidence_refs) ? template.evidence_refs : [];
  check(Array.isArray(template.evidence_refs), "template evidence_refs is not an array: " + templateId);
  check(
    new Set(evidenceRefs).size === evidenceRefs.length &&
      evidenceRefs.every((reference) => typeof reference === "string" && reference.length > 0),
    "template evidence_refs is not a unique nonempty string set: " + templateId,
  );
  for (const reference of evidenceRefs) {
    check(evidenceReferenceResolves(reference), "unresolved template evidence reference: " + reference);
  }
  if (template.status === "supported" || template.status === "refuted") {
    check(evidenceRefs.length > 0, "supported/refuted template lacks evidence: " + templateId);
    for (const reference of evidenceRefs) {
      check(
        sourceArtifactPaths.has(evidenceReferencePath(reference)),
        `supported/refuted template evidence is not SHA-pinned: ${templateId}: ${reference}`,
      );
    }
  }
  if (template.status === "supported") {
    check(
      template.obligation_id === undefined,
      "supported template unexpectedly carries an obligation: " + templateId,
    );
  } else if (template.status === "unknown") {
    check(
      typeof template.obligation_id === "string" &&
        obligationIds.has(template.obligation_id),
      "unknown template lacks a known obligation: " + templateId,
    );
    if (typeof template.obligation_id === "string") {
      referencedObligationIds.add(template.obligation_id);
    }
  } else if (template.status === "refuted") {
    check(
      template.obligation_id === undefined ||
        (typeof template.obligation_id === "string" &&
          obligationIds.has(template.obligation_id)),
      "refuted template has an unknown obligation: " + templateId,
    );
    if (typeof template.obligation_id === "string") {
      referencedObligationIds.add(template.obligation_id);
    }
  } else {
    check(
      template.obligation_id === undefined,
      "not-applicable template unexpectedly carries an obligation: " + templateId,
    );
  }
}

function requireTemplate(templateId, context) {
  check(
    typeof templateId === "string" && Object.hasOwn(resultTemplates, templateId),
    context + " references an unknown result template " + String(templateId),
  );
  if (typeof templateId === "string") {
    usedTemplateIds.add(templateId);
  }
}

const mapProfilesByPredicate = new Map();
const subgroupProfilesByPredicate = new Map();
check(
  canonicalJson(sortedObjectKeys(profiles)) ===
    canonicalJson(sortedObjectKeys(expectedProfiles)),
  "profile set changed from the authoritative verifier matrix",
);
for (const [profileId, expectedProfile] of Object.entries(expectedProfiles)) {
  const actualProfile = profiles[profileId];
  check(isPlainObject(actualProfile), `missing authoritative profile ${profileId}`);
  if (!isPlainObject(actualProfile)) {
    continue;
  }
  check(
    actualProfile.kind === expectedProfile.kind &&
      actualProfile.predicate_id === expectedProfile.predicate_id &&
      canonicalJson(actualProfile.claims) === canonicalJson(expectedProfile.claims) &&
      actualProfile.full_transfer_pipeline === expectedProfile.full_transfer_pipeline,
    `claim/template assignment changed for authoritative profile ${profileId}`,
  );
}
for (const [profileId, profile] of Object.entries(profiles)) {
  check(isPlainObject(profile), "profile is not an object: " + profileId);
  if (!isPlainObject(profile)) {
    continue;
  }
  const expectedClaimIds =
    profile.kind === "map_base"
      ? mapClaimIds
      : profile.kind === "subgroup_base"
        ? subgroupClaimIds
        : [];
  check(
    profile.kind === "map_base" || profile.kind === "subgroup_base",
    "invalid profile kind: " + profileId,
  );
  check(
    typeof profile.predicate_id === "string" && profile.predicate_id.length > 0,
    "invalid profile predicate: " + profileId,
  );
  check(
    Array.isArray(profile.hypotheses) &&
      profile.hypotheses.length > 0 &&
      new Set(profile.hypotheses).size === profile.hypotheses.length &&
      profile.hypotheses.every(
        (hypothesis) => typeof hypothesis === "string" && hypothesis.trim().length > 0,
      ),
    "profile hypotheses are not a unique nonempty string set: " + profileId,
  );
  check(
    typeof profile.boundary === "string" && profile.boundary.trim().length > 0,
    "profile boundary is empty: " + profileId,
  );
  check(isPlainObject(profile.claims), "profile claims is not an object: " + profileId);
  const actualClaimIds = isPlainObject(profile.claims) ? sortedObjectKeys(profile.claims) : [];
  check(
    canonicalJson(actualClaimIds) === canonicalJson(expectedClaimIds),
    "profile claim set is incomplete or has extras: " + profileId,
  );
  for (const [claimId, templateId] of Object.entries(profile.claims ?? {})) {
    check(expectedClaimIds.includes(claimId), "profile uses an unknown claim: " + profileId);
    requireTemplate(templateId, "profile " + profileId + " claim " + claimId);
  }
  if (profile.kind === "map_base") {
    check(
      !mapProfilesByPredicate.has(profile.predicate_id),
      "duplicate map predicate profile " + profile.predicate_id,
    );
    mapProfilesByPredicate.set(profile.predicate_id, profileId);
    requireTemplate(profile.full_transfer_pipeline, "map profile " + profileId + " full pipeline");
  } else if (profile.kind === "subgroup_base") {
    check(
      !subgroupProfilesByPredicate.has(profile.predicate_id),
      "duplicate subgroup predicate profile " + profile.predicate_id,
    );
    subgroupProfilesByPredicate.set(profile.predicate_id, profileId);
    check(
      profile.full_transfer_pipeline === undefined,
      "subgroup profile unexpectedly defines full_transfer_pipeline: " + profileId,
    );
  }
}
const fullPipelinePrerequisites = [
  "cover_existence",
  "curve_map_correctness",
  "jacobian_arithmetic_for_cover",
  "pullback_point_classes",
  "pushforward_point_divisors",
  "arbitrary_jacobian_pushforward",
  "composition_on_typed_domain",
];
for (const [profileId, profile] of Object.entries(profiles)) {
  if (profile.kind !== "map_base") {
    continue;
  }
  const pipelineStatus = resultTemplates[profile.full_transfer_pipeline]?.status;
  if (pipelineStatus === "supported") {
    for (const claimId of fullPipelinePrerequisites) {
      const prerequisiteStatus = resultTemplates[profile.claims?.[claimId]]?.status;
      check(
        prerequisiteStatus === "supported",
        `full transfer pipeline ${profileId} is supported while prerequisite ${claimId} is ${String(prerequisiteStatus)}`,
      );
    }
  }
}
check(
  canonicalJson([...mapProfilesByPredicate.keys()].sort()) === canonicalJson(closedMapPredicates),
  "map profiles do not implement the closed predicate set",
);
check(
  canonicalJson([...subgroupProfilesByPredicate.keys()].sort()) ===
    canonicalJson(closedSubgroupPredicates),
  "subgroup profiles do not implement the closed predicate set",
);

const mapContexts = new Map();
const subgroupContexts = new Map();
for (const graphCurve of graph.curves) {
  if (graphCurve.cover_links.length !== 1) {
    continue;
  }
  const link = graphCurve.cover_links[0];
  const graphMap = graph.maps[link.map_uid];
  const coverResult = reportByModelHash.get(graphCurve.model_sha256);
  if (graphMap === undefined || coverResult === undefined) {
    continue;
  }
  const targetModel = JSON.parse(graphCurve.model_json);
  let predicateId = null;
  if (
    coverResult.certificate.construction === "binary_cubic_pullback_v1" &&
    targetModel.field.startsWith("f2m-") &&
    targetModel.form === "y^2+xy=x^3+a*x^2+b" &&
    graphMap.record.target_model_map === null
  ) {
    predicateId = "binary_cubic_same_field";
  } else if (
    coverResult.certificate.construction === "prime_quadratic_pullback_v1" &&
    targetModel.field.startsWith("fp-") &&
    targetModel.form === "y^2=x^3+a*x+b" &&
    graphMap.record.target_model_map === null
  ) {
    predicateId = "odd_prime_quadratic_direct_same_field";
  } else if (
    coverResult.certificate.construction === "prime_quadratic_pullback_v1" &&
    targetModel.field.startsWith("fp-") &&
    targetModel.form !== "y^2=x^3+a*x+b" &&
    graphMap.record.target_model_map !== null
  ) {
    predicateId = "odd_prime_quadratic_normalized_same_field";
  } else if (
    coverResult.certificate.construction === "prime_quadratic_pullback_v1" &&
    targetModel.field.startsWith("fpk-") &&
    targetModel.form === "y^2=x^3+a*x+b" &&
    graphMap.record.target_model_map === null
  ) {
    predicateId = "odd_extension_quadratic_direct_same_field";
  }
  check(predicateId !== null, "map matches no closed predicate: " + link.map_uid);
  mapContexts.set(link.map_uid, {
    predicateId,
    graphCurve,
    graphMap,
    coverResult,
  });
  for (const representation of graphCurve.representations) {
    const triple = [representation.ec1, representation.curve_uid, link.map_uid];
    const tripleKey = JSON.stringify(triple);
    check(!subgroupContexts.has(tripleKey), "duplicate subgroup triple " + tripleKey);
    subgroupContexts.set(tripleKey, {
      triple,
      predicateId:
        graphMap.record.degree === 3
          ? "linked_subgroup_degree3"
          : graphMap.record.degree === 2
            ? "linked_subgroup_degree2"
            : null,
      representation,
      graphCurve,
      graphMap,
    });
    check(
      graphMap.record.degree === 2 || graphMap.record.degree === 3,
      "subgroup map has an unsupported degree: " + tripleKey,
    );
  }
}
check(mapContexts.size === Object.keys(graph.maps).length, "not every CV1 map has a predicate context");
check(subgroupContexts.size === representationCount, "not every EC1 representation has a subgroup context");

const typedBindings = isPlainObject(assessment.typed_bindings)
  ? assessment.typed_bindings
  : {};
const actualFieldBindings = isPlainObject(typedBindings.fields)
  ? typedBindings.fields
  : {};
const actualJacobianBindings = isPlainObject(typedBindings.jacobians)
  ? typedBindings.jacobians
  : {};
const actualTypedMapBindings = isPlainObject(typedBindings.maps)
  ? typedBindings.maps
  : {};
check(isPlainObject(assessment.typed_bindings), "typed_bindings is not an object");
check(isPlainObject(typedBindings.fields), "typed_bindings.fields is not an object");
check(isPlainObject(typedBindings.jacobians), "typed_bindings.jacobians is not an object");
check(isPlainObject(typedBindings.maps), "typed_bindings.maps is not an object");

const expectedFieldBindings = new Map();
const expectedJacobianBindings = new Map();
const expectedTypedMapBindings = new Map();
const fieldUidByLabel = new Map();
for (const [mapUid, context] of mapContexts) {
  const link = context.graphCurve.cover_links[0];
  const sourceCover = graph.covers[link.cover_uid];
  if (sourceCover === undefined) {
    continue;
  }
  const targetModel = JSON.parse(context.graphCurve.model_json);
  const targetFieldRecord = targetDeclaredFieldRecord(targetModel);
  check(
    canonicalJson(sourceCover.record.field) === canonicalJson(targetFieldRecord),
    "typed binding source/target field mismatch for " + mapUid,
  );
  const field = expectedDeclaredFieldBinding(sourceCover.record.field);
  const priorField = expectedFieldBindings.get(field.field_uid);
  check(
    priorField === undefined || canonicalJson(priorField) === canonicalJson(field),
    "declared-field UID collision for " + mapUid,
  );
  expectedFieldBindings.set(field.field_uid, field);
  const priorLabelUid = fieldUidByLabel.get(field.record.label);
  check(
    priorLabelUid === undefined || priorLabelUid === field.field_uid,
    "declared-field label names multiple records: " + field.record.label,
  );
  fieldUidByLabel.set(field.record.label, field.field_uid);

  if (field.representation === "prime") {
    check(
      field.record.label === `fp-${field.characteristic}` && field.extension_degree === 1,
      "invalid prime-field context for " + field.record.label,
    );
  } else if (field.representation === "odd_polynomial_extension") {
    check(
      Array.isArray(field.record.modulus) &&
        field.record.modulus.length === field.extension_degree &&
        field.record.label.startsWith(`fpk-${field.characteristic}-${field.extension_degree}-`),
      "invalid odd-extension context for " + field.record.label,
    );
  } else {
    check(
      field.representation === "binary_polynomial_extension" &&
        field.characteristic === "2" &&
        typeof field.record.modulus === "string" &&
        field.record.label.startsWith(`f2m-${field.extension_degree}-`),
      "invalid binary-extension context for " + field.record.label,
    );
  }

  const binary = context.predicateId === "binary_cubic_same_field";
  const infinityConfiguration = binary
    ? "one_rational_point"
    : "two_rational_points";
  check(
    binary
      ? sourceCover.record.genus === 3 && sourceCover.record.f.length === 8
      : sourceCover.record.genus === 2 && sourceCover.record.f.length === 7,
    "Jacobian context has the wrong source geometry for " + mapUid,
  );
  const jacobianRecord = {
    schema_version: "jacobian-context/v1",
    source_cover_uid: link.cover_uid,
    declared_field_uid: field.field_uid,
    dimension: sourceCover.record.genus,
    object: "Picard_degree_zero",
    infinity_configuration: infinityConfiguration,
  };
  const jacobianUid =
    `urn:jacobian-context:1:sha256:${sha256(canonicalJson(jacobianRecord))}`;
  const jacobian = { jacobian_uid: jacobianUid, record: jacobianRecord };
  const priorJacobian = expectedJacobianBindings.get(jacobianUid);
  check(
    priorJacobian === undefined || canonicalJson(priorJacobian) === canonicalJson(jacobian),
    "Jacobian-context UID collision for " + mapUid,
  );
  expectedJacobianBindings.set(jacobianUid, jacobian);

  check(!expectedTypedMapBindings.has(mapUid), "duplicate expected typed map " + mapUid);
  expectedTypedMapBindings.set(mapUid, {
    map_uid: mapUid,
    target_icv1_slug: context.graphCurve.icv1_identity.slug,
    target_model_sha256: context.graphCurve.model_sha256,
    source_cover_uid: link.cover_uid,
    source_jacobian_uid: jacobianUid,
    field_relationship: {
      kind: "same_declared_field_identity",
      source_field_uid: field.field_uid,
      target_field_uid: field.field_uid,
      degree: 1,
    },
  });
}

check(expectedFieldBindings.size === 318, "unexpected declared-field context count");
check(fieldUidByLabel.size === 318, "declared-field labels are not bijective with field contexts");
check(expectedJacobianBindings.size === 2098, "unexpected Jacobian-context count");
check(expectedTypedMapBindings.size === 2100, "unexpected typed-map binding count");
check(
  canonicalJson(actualFieldBindings) === canonicalJson(sortedRecord(expectedFieldBindings)),
  "typed declared-field bindings are not the exact independently derived set",
);
check(
  canonicalJson(actualJacobianBindings) === canonicalJson(sortedRecord(expectedJacobianBindings)),
  "typed Jacobian bindings are not the exact independently derived set",
);
check(
  canonicalJson(actualTypedMapBindings) === canonicalJson(sortedRecord(expectedTypedMapBindings)),
  "typed map bindings are not the exact independently derived set",
);

const memberships = isPlainObject(assessment.memberships) ? assessment.memberships : {};
const mapMemberships = Array.isArray(memberships.map_base) ? memberships.map_base : [];
const subgroupMemberships = Array.isArray(memberships.subgroup_base)
  ? memberships.subgroup_base
  : [];
check(isPlainObject(assessment.memberships), "memberships is not an object");
check(Array.isArray(memberships.map_base), "memberships.map_base is not an array");
check(Array.isArray(memberships.subgroup_base), "memberships.subgroup_base is not an array");

const mapProfileByMember = new Map();
const mapMembershipProfileIds = new Set();
for (const group of mapMemberships) {
  check(isPlainObject(group), "map membership group is not an object");
  if (!isPlainObject(group)) {
    continue;
  }
  const profile = profiles[group.profile_id];
  check(
    profile !== undefined &&
      profile.kind === "map_base" &&
      !mapMembershipProfileIds.has(group.profile_id),
    "map membership has an unknown, wrong-kind, or duplicate profile",
  );
  mapMembershipProfileIds.add(group.profile_id);
  const members = Array.isArray(group.members) ? group.members : [];
  check(Array.isArray(group.members), "map membership members is not an array");
  check(group.count === members.length, "map membership count mismatch for " + String(group.profile_id));
  check(
    members.every((member) => typeof member === "string" && member.length > 0),
    "map membership contains a non-string member for " + String(group.profile_id),
  );
  check(
    canonicalJson(members) === canonicalJson([...members].sort(strictStringOrder)),
    "map membership is not sorted for " + String(group.profile_id),
  );
  check(
    new Set(members).size === members.length,
    "map membership contains duplicates for " + String(group.profile_id),
  );
  check(
    group.members_sha256 ===
      memberDigest(`${mapMembershipDomain}/${group.profile_id}`, members),
    "map membership digest mismatch for " + String(group.profile_id),
  );
  for (const mapUid of members) {
    const context = mapContexts.get(mapUid);
    check(context !== undefined, "map membership references an unknown target " + String(mapUid));
    check(!mapProfileByMember.has(mapUid), "map appears in multiple base profiles " + String(mapUid));
    mapProfileByMember.set(mapUid, group.profile_id);
    if (context !== undefined && profile !== undefined) {
      check(
        context.predicateId === profile.predicate_id,
        "map does not satisfy its profile predicate " + String(mapUid),
      );
    }
  }
}
check(
  canonicalJson([...mapMembershipProfileIds].sort()) ===
    canonicalJson(
      Object.entries(profiles)
        .filter(([, profile]) => profile.kind === "map_base")
        .map(([profileId]) => profileId)
        .sort(),
    ),
  "map membership groups do not cover every map profile exactly once",
);
check(
  canonicalJson([...mapProfileByMember.keys()].sort(strictStringOrder)) ===
    canonicalJson([...mapContexts.keys()].sort(strictStringOrder)),
  "map memberships are not an exact partition of all CV1 UIDs",
);

const subgroupProfileByMember = new Map();
const subgroupMembershipProfileIds = new Set();
for (const group of subgroupMemberships) {
  check(isPlainObject(group), "subgroup membership group is not an object");
  if (!isPlainObject(group)) {
    continue;
  }
  const profile = profiles[group.profile_id];
  check(
    profile !== undefined &&
      profile.kind === "subgroup_base" &&
      !subgroupMembershipProfileIds.has(group.profile_id),
    "subgroup membership has an unknown, wrong-kind, or duplicate profile",
  );
  subgroupMembershipProfileIds.add(group.profile_id);
  const members = Array.isArray(group.members) ? group.members : [];
  check(Array.isArray(group.members), "subgroup membership members is not an array");
  check(
    group.count === members.length,
    "subgroup membership count mismatch for " + String(group.profile_id),
  );
  check(
    members.every(
      (member) =>
        Array.isArray(member) &&
        member.length === 3 &&
        member.every((part) => typeof part === "string" && part.length > 0),
    ),
    "subgroup membership contains a malformed triple for " + String(group.profile_id),
  );
  check(
    canonicalJson(members) === canonicalJson([...members].sort(compareTriples)),
    "subgroup membership is not sorted for " + String(group.profile_id),
  );
  const memberKeys = members.map((member) => JSON.stringify(member));
  check(
    new Set(memberKeys).size === memberKeys.length,
    "subgroup membership contains duplicate triples for " + String(group.profile_id),
  );
  check(
    group.members_sha256 ===
      memberDigest(`${subgroupMembershipDomain}/${group.profile_id}`, members),
    "subgroup membership digest mismatch for " + String(group.profile_id),
  );
  for (const triple of members) {
    const tripleKey = JSON.stringify(triple);
    const context = subgroupContexts.get(tripleKey);
    check(context !== undefined, "subgroup membership references an unknown triple " + tripleKey);
    check(
      !subgroupProfileByMember.has(tripleKey),
      "subgroup triple appears in multiple base profiles " + tripleKey,
    );
    subgroupProfileByMember.set(tripleKey, group.profile_id);
    if (context !== undefined && profile !== undefined) {
      check(
        context.predicateId === profile.predicate_id,
        "subgroup does not satisfy its profile predicate " + tripleKey,
      );
    }
  }
}
check(
  canonicalJson([...subgroupMembershipProfileIds].sort()) ===
    canonicalJson(
      Object.entries(profiles)
        .filter(([, profile]) => profile.kind === "subgroup_base")
        .map(([profileId]) => profileId)
        .sort(),
    ),
  "subgroup membership groups do not cover every subgroup profile exactly once",
);
check(
  canonicalJson([...subgroupProfileByMember.keys()].sort(strictStringOrder)) ===
    canonicalJson([...subgroupContexts.keys()].sort(strictStringOrder)),
  "subgroup memberships are not an exact partition of all EC1/map triples",
);

const expectedPredicateCounts = {
  binary_cubic_same_field: 160,
  odd_prime_quadratic_direct_same_field: 586,
  odd_prime_quadratic_normalized_same_field: 30,
  odd_extension_quadratic_direct_same_field: 1324,
  linked_subgroup_degree3: 134,
  linked_subgroup_degree2: 530,
};
const actualPredicateCounts = {};
for (const context of mapContexts.values()) {
  actualPredicateCounts[context.predicateId] =
    (actualPredicateCounts[context.predicateId] ?? 0) + 1;
}
for (const context of subgroupContexts.values()) {
  actualPredicateCounts[context.predicateId] =
    (actualPredicateCounts[context.predicateId] ?? 0) + 1;
}
check(
  canonicalJson(actualPredicateCounts) === canonicalJson(expectedPredicateCounts),
  "closed profile predicate counts changed",
);

const expectedExactOverlays = {
  "map.binary-f2m7-degree3": {
    kind: "map",
    target:
      "urn:curve-cover-map:1:sha256:0dc60016b03df7d55ee829880508251fe8261c92a241b76a78987eeaa9787a03",
    claim_overrides: {
      curve_map_correctness: "map.f2m7-geometry.supported",
    },
    evidence_additions: [
      "research/hyperelliptic_cover_infrastructure_20261006/evidence.json#icv1-f2m7-t13-616700dd",
      "research/hyperelliptic_cover_infrastructure_20261006/README.md",
    ],
  },
  "subgroup.binary-f2m7-order29": {
    kind: "subgroup",
    target: [
      "EC1N7Ce0hb6d297a2ca08",
      "urn:ec-record:1:sha256:b6d297a2ca085e9d7b8a6858df7421083e12e740c89b1647660b413e00ed02ff",
      "urn:curve-cover-map:1:sha256:0dc60016b03df7d55ee829880508251fe8261c92a241b76a78987eeaa9787a03",
    ],
    claim_overrides: {
      numeric_exact_source_order: "subgroup.exact-order29.supported",
      image_order_equals_catalog_order: "subgroup.exact-image-order29.supported",
    },
    evidence_additions: [
      "research/hyperelliptic_cover_infrastructure_20261006/evidence.json#subgroup_transfer",
      "src/binary_ecc/cover_transfer.rs#subgroup_status_and_arbitrary_pushforward_boundary",
    ],
  },
  "map.prime-fp827-degree2": {
    kind: "map",
    target:
      "urn:curve-cover-map:1:sha256:5768ce11909c0f56d84f05e70d71dfb6a46fb175b6f9cc0ed727eec14aef9c8e",
    claim_overrides: {
      curve_map_correctness: "map.fp827-geometry.supported",
    },
    evidence_additions: [
      "research/hyperelliptic_cover_infrastructure_20261006/evidence.json#icv1-fp10-t5-192cb216",
      "research/hyperelliptic_cover_infrastructure_20261006/README.md",
    ],
  },
  "subgroup.prime-fp827-order823": {
    kind: "subgroup",
    target: [
      "EC1P10Cfphcf15d87ab724",
      "urn:ec-record:1:sha256:cf15d87ab724c8d1fc74d88af320937e4c1cc7e706d4e0e8a0d11734831342a3",
      "urn:curve-cover-map:1:sha256:5768ce11909c0f56d84f05e70d71dfb6a46fb175b6f9cc0ed727eec14aef9c8e",
    ],
    claim_overrides: {},
    evidence_additions: [
      "docs/curves/cover-links.yaml#EC1P10Cfphcf15d87ab724",
      "research/hyperelliptic_cover_infrastructure_20261006/evidence.json#icv1-fp10-t5-192cb216",
    ],
  },
};
const exactOverlays = Array.isArray(assessment.exact_overlays)
  ? assessment.exact_overlays
  : [];
check(Array.isArray(assessment.exact_overlays), "exact_overlays is not an array");
check(
  canonicalJson(
    exactOverlays
      .map((overlay) => overlay?.overlay_id)
      .sort(strictStringOrder),
  ) === canonicalJson(sortedObjectKeys(expectedExactOverlays)),
  "exact-overlay ID set changed from the authoritative verifier matrix",
);
const overlayIds = new Set();
const overlayByMap = new Map();
const overlayBySubgroup = new Map();
for (const overlay of exactOverlays) {
  check(isPlainObject(overlay), "exact_overlays contains a non-object");
  if (!isPlainObject(overlay)) {
    continue;
  }
  check(
    typeof overlay.overlay_id === "string" &&
      overlay.overlay_id.length > 0 &&
      !overlayIds.has(overlay.overlay_id),
    "duplicate or invalid exact-overlay ID",
  );
  overlayIds.add(overlay.overlay_id);
  const expectedOverlay = expectedExactOverlays[overlay.overlay_id];
  check(expectedOverlay !== undefined, `unexpected exact overlay ${String(overlay.overlay_id)}`);
  if (expectedOverlay !== undefined) {
    check(
      overlay.kind === expectedOverlay.kind &&
        canonicalJson(overlay.target) === canonicalJson(expectedOverlay.target) &&
        canonicalJson(overlay.claim_overrides) ===
          canonicalJson(expectedOverlay.claim_overrides) &&
        canonicalJson(overlay.evidence_additions) ===
          canonicalJson(expectedOverlay.evidence_additions),
      `typed target/claim/evidence binding changed for exact overlay ${overlay.overlay_id}`,
    );
  }
  check(
    typeof overlay.boundary === "string" && overlay.boundary.trim().length > 0,
    "exact overlay boundary is empty: " + String(overlay.overlay_id),
  );
  check(
    isPlainObject(overlay.claim_overrides),
    "exact overlay claim_overrides is not an object: " + String(overlay.overlay_id),
  );
  const evidenceAdditions = Array.isArray(overlay.evidence_additions)
    ? overlay.evidence_additions
    : [];
  check(
    Array.isArray(overlay.evidence_additions) &&
      new Set(evidenceAdditions).size === evidenceAdditions.length &&
      evidenceAdditions.every(
        (reference) => typeof reference === "string" && reference.length > 0,
      ),
    "exact overlay evidence_additions is not a unique string set: " + String(overlay.overlay_id),
  );
  for (const reference of evidenceAdditions) {
    check(evidenceReferenceResolves(reference), "unresolved overlay evidence reference: " + reference);
    check(
      sourceArtifactPaths.has(evidenceReferencePath(reference)),
      `overlay evidence is not SHA-pinned: ${String(overlay.overlay_id)}: ${reference}`,
    );
  }

  let allowedClaims = [];
  let targetKey = null;
  if (overlay.kind === "map") {
    check(
      typeof overlay.target === "string" && mapContexts.has(overlay.target),
      "map overlay targets an unknown CV1 UID: " + String(overlay.overlay_id),
    );
    targetKey = overlay.target;
    allowedClaims = [...mapClaimIds, "full_transfer_pipeline"];
    check(
      typeof targetKey === "string" && !overlayByMap.has(targetKey),
      "multiple exact overlays target one map: " + String(targetKey),
    );
    if (typeof targetKey === "string") {
      overlayByMap.set(targetKey, overlay);
    }
  } else if (overlay.kind === "subgroup") {
    const validTriple =
      Array.isArray(overlay.target) &&
      overlay.target.length === 3 &&
      overlay.target.every((part) => typeof part === "string" && part.length > 0);
    targetKey = validTriple ? JSON.stringify(overlay.target) : null;
    check(
      validTriple && subgroupContexts.has(targetKey),
      "subgroup overlay targets an unknown exact triple: " + String(overlay.overlay_id),
    );
    allowedClaims = subgroupClaimIds;
    check(
      targetKey !== null && !overlayBySubgroup.has(targetKey),
      "multiple exact overlays target one subgroup triple: " + String(targetKey),
    );
    if (targetKey !== null) {
      overlayBySubgroup.set(targetKey, overlay);
    }
  } else {
    check(false, "invalid exact-overlay kind: " + String(overlay.overlay_id));
  }

  for (const [claimId, templateId] of Object.entries(overlay.claim_overrides ?? {})) {
    check(
      allowedClaims.includes(claimId),
      "exact overlay overrides an unknown claim: " + String(overlay.overlay_id),
    );
    requireTemplate(
      templateId,
      "exact overlay " + String(overlay.overlay_id) + " claim " + claimId,
    );
  }
}

function zeroStatusCounts() {
  return {
    unknown: 0,
    supported: 0,
    refuted: 0,
    not_applicable: 0,
  };
}

function addTemplateStatus(counts, templateId, context) {
  const template = resultTemplates[templateId];
  check(template !== undefined, context + " has no effective result template");
  if (template !== undefined && allowedStatuses.has(template.status)) {
    counts[template.status] += 1;
  }
}

const effectiveStatusCounts = {
  map_claims: zeroStatusCounts(),
  map_full_transfer_pipeline: zeroStatusCounts(),
  subgroup_claims: zeroStatusCounts(),
  all: zeroStatusCounts(),
};
for (const mapUid of [...mapContexts.keys()].sort(strictStringOrder)) {
  const profileId = mapProfileByMember.get(mapUid);
  const profile = profiles[profileId];
  const overlay = overlayByMap.get(mapUid);
  if (profile === undefined) {
    check(false, "map has no effective base profile: " + mapUid);
    continue;
  }
  for (const claimId of mapClaimIds) {
    const templateId = overlay?.claim_overrides?.[claimId] ?? profile.claims[claimId];
    addTemplateStatus(effectiveStatusCounts.map_claims, templateId, "map " + mapUid + " claim " + claimId);
    addTemplateStatus(effectiveStatusCounts.all, templateId, "map " + mapUid + " claim " + claimId);
  }
  const pipelineTemplate =
    overlay?.claim_overrides?.full_transfer_pipeline ?? profile.full_transfer_pipeline;
  addTemplateStatus(
    effectiveStatusCounts.map_full_transfer_pipeline,
    pipelineTemplate,
    "map " + mapUid + " full_transfer_pipeline",
  );
  addTemplateStatus(
    effectiveStatusCounts.all,
    pipelineTemplate,
    "map " + mapUid + " full_transfer_pipeline",
  );
}
for (const subgroupKey of [...subgroupContexts.keys()].sort(strictStringOrder)) {
  const profileId = subgroupProfileByMember.get(subgroupKey);
  const profile = profiles[profileId];
  const overlay = overlayBySubgroup.get(subgroupKey);
  if (profile === undefined) {
    check(false, "subgroup has no effective base profile: " + subgroupKey);
    continue;
  }
  for (const claimId of subgroupClaimIds) {
    const templateId = overlay?.claim_overrides?.[claimId] ?? profile.claims[claimId];
    addTemplateStatus(
      effectiveStatusCounts.subgroup_claims,
      templateId,
      "subgroup " + subgroupKey + " claim " + claimId,
    );
    addTemplateStatus(
      effectiveStatusCounts.all,
      templateId,
      "subgroup " + subgroupKey + " claim " + claimId,
    );
  }
}

const declaredEffectiveStatusCounts =
  assessment.effective_status_counts ?? assessment.coverage?.effective_status_counts;
if (declaredEffectiveStatusCounts !== undefined) {
  check(
    canonicalJson(declaredEffectiveStatusCounts) === canonicalJson(effectiveStatusCounts),
    "declared effective-status counts do not match all effective claim assignments",
  );
}
check(
  canonicalJson([...usedTemplateIds].sort()) === canonicalJson(sortedObjectKeys(resultTemplates)),
  "one or more result templates are unused",
);

function gcd(left, right) {
  let a = left;
  let b = right;
  while (b !== 0n) {
    const remainder = a % b;
    a = b;
    b = remainder;
  }
  return a;
}
for (const [subgroupKey, context] of subgroupContexts) {
  const subgroupOrder = BigInt(context.representation.curve.subgroup_order);
  const degree = BigInt(context.graphMap.record.degree);
  check(subgroupOrder > 1n, "subgroup order is not positive for " + subgroupKey);
  check(gcd(subgroupOrder, degree) === 1n, "subgroup order is not coprime to map degree for " + subgroupKey);
}

const coverage = isPlainObject(assessment.coverage) ? assessment.coverage : {};
check(isPlainObject(assessment.coverage), "coverage is not an object");
const expectedCoverageKeys = [
  "catalog_maps",
  "unique_cover_models",
  "geometry_verified_maps",
  "ec1_representations",
  "assessed_maps",
  "unassessed_maps",
  "assessed_subgroups",
  "unassessed_subgroups",
  "map_membership_entries",
  "subgroup_membership_entries",
  "map_profiles",
  "subgroup_profiles",
  "scope_note",
].sort();
check(
  canonicalJson(sortedObjectKeys(coverage)) === canonicalJson(expectedCoverageKeys),
  "coverage key set changed",
);
const geometryVerifiedMaps = report.curves.filter(
  (curve) => curve.status === "verified_over_declared_field",
).length;
const expectedCoverage = {
  catalog_maps: Object.keys(graph.maps).length,
  unique_cover_models: Object.keys(graph.covers).length,
  geometry_verified_maps: geometryVerifiedMaps,
  ec1_representations: representationCount,
  assessed_maps: mapProfileByMember.size,
  unassessed_maps: Object.keys(graph.maps).length - mapProfileByMember.size,
  assessed_subgroups: subgroupProfileByMember.size,
  unassessed_subgroups: representationCount - subgroupProfileByMember.size,
  map_membership_entries: mapMemberships.reduce(
    (sum, group) => sum + (Array.isArray(group.members) ? group.members.length : 0),
    0,
  ),
  subgroup_membership_entries: subgroupMemberships.reduce(
    (sum, group) => sum + (Array.isArray(group.members) ? group.members.length : 0),
    0,
  ),
  map_profiles: [...mapProfilesByPredicate.keys()].length,
  subgroup_profiles: [...subgroupProfilesByPredicate.keys()].length,
};
for (const [key, expectedValue] of Object.entries(expectedCoverage)) {
  check(coverage[key] === expectedValue, "coverage count mismatch for " + key);
}
check(
  typeof coverage.scope_note === "string" && coverage.scope_note.trim().length > 0,
  "coverage scope_note is empty",
);
check(coverage.assessed_maps === 2100 && coverage.unassessed_maps === 0, "map coverage is incomplete");
check(
  coverage.assessed_subgroups === 664 && coverage.unassessed_subgroups === 0,
  "subgroup coverage is incomplete",
);

const inlineAssessments = [];
for (const artifact of sourceArtifacts) {
  if (isPlainObject(artifact.freshness)) {
    inlineAssessments.push(["source artifact " + artifact.path, artifact.freshness]);
  }
}
const controls = Array.isArray(assessment.controls) ? assessment.controls : [];
check(Array.isArray(assessment.controls), "controls is not an array");
check(
  canonicalJson(
    controls.map((control) => control?.control).sort(strictStringOrder),
  ) === canonicalJson(sortedObjectKeys(expectedControlStatuses)),
  "control ID set changed from the authoritative verifier matrix",
);
for (const control of controls) {
  const expected = expectedControlStatuses[control?.control];
  check(expected !== undefined, `unexpected control ${String(control?.control)}`);
  if (expected !== undefined) {
    check(
      control.assessment?.status === expected[0] &&
        (control.assessment?.obligation_id ?? null) === expected[1],
      `authoritative status/obligation mismatch for control ${control.control}`,
    );
  }
  if (isPlainObject(control.assessment)) {
    inlineAssessments.push(["control " + String(control.control), control.assessment]);
  }
}
for (const [label, inlineAssessment] of inlineAssessments) {
  check(allowedStatuses.has(inlineAssessment.status), label + " has an invalid status");
  const refs = Array.isArray(inlineAssessment.evidence_refs)
    ? inlineAssessment.evidence_refs
    : [];
  check(Array.isArray(inlineAssessment.evidence_refs), label + " evidence_refs is not an array");
  check(
    typeof inlineAssessment.scope === "string" && inlineAssessment.scope.trim().length > 0,
    label + " scope is empty",
  );
  check(
    typeof inlineAssessment.rationale === "string" &&
      inlineAssessment.rationale.trim().length > 0,
    label + " rationale is empty",
  );
  check(
    new Set(refs).size === refs.length,
    label + " evidence_refs contains duplicates",
  );
  for (const reference of refs) {
    check(evidenceReferenceResolves(reference), "unresolved inline evidence reference: " + reference);
  }
  if (inlineAssessment.status === "supported" || inlineAssessment.status === "refuted") {
    check(refs.length > 0, label + " supported/refuted result lacks evidence");
    for (const reference of refs) {
      check(
        sourceArtifactPaths.has(evidenceReferencePath(reference)),
        `${label} supported/refuted evidence is not SHA-pinned: ${reference}`,
      );
    }
  }
  if (inlineAssessment.status === "supported") {
    check(
      inlineAssessment.obligation_id === undefined,
      label + " supported result unexpectedly carries an obligation",
    );
  } else if (inlineAssessment.status === "unknown") {
    check(
      typeof inlineAssessment.obligation_id === "string" &&
        obligationIds.has(inlineAssessment.obligation_id),
      label + " unknown result lacks a known obligation",
    );
    if (typeof inlineAssessment.obligation_id === "string") {
      referencedObligationIds.add(inlineAssessment.obligation_id);
    }
  } else if (inlineAssessment.status === "refuted") {
    check(
      inlineAssessment.obligation_id === undefined ||
        (typeof inlineAssessment.obligation_id === "string" &&
          obligationIds.has(inlineAssessment.obligation_id)),
      label + " refuted result has an unknown obligation",
    );
    if (typeof inlineAssessment.obligation_id === "string") {
      referencedObligationIds.add(inlineAssessment.obligation_id);
    }
  } else {
    check(
      inlineAssessment.obligation_id === undefined,
      label + " not-applicable result unexpectedly carries an obligation",
    );
  }
}
check(
  canonicalJson([...referencedObligationIds].sort(strictStringOrder)) ===
    canonicalJson([...obligationIds].sort(strictStringOrder)),
  "open obligations are orphaned or referenced without declarations",
);
const allowedCostStatuses = new Set(["measured", "modeled", "unknown", "not_applicable"]);
const costItems = Array.isArray(assessment.global_cost_accounting)
  ? assessment.global_cost_accounting
  : [];
check(
  canonicalJson(costItems.map((item) => item?.phase).sort(strictStringOrder)) ===
    canonicalJson(sortedObjectKeys(expectedCostStatuses)),
  "global cost phase set changed from the authoritative verifier matrix",
);
for (const item of costItems) {
  check(allowedCostStatuses.has(item.status), "global cost phase has an invalid status");
  check(
    expectedCostStatuses[item.phase] === item.status,
    `authoritative cost status mismatch for ${String(item.phase)}`,
  );
  if (item.status === "unknown" || item.status === "not_applicable") {
    check(item.value === null, "unknown/not-applicable cost phase has a numeric value");
  } else {
    check(
      typeof item.value === "number" && Number.isFinite(item.value),
      "measured/modeled cost phase lacks a finite numeric value",
    );
  }
  const sourceRefs = Array.isArray(item.source_refs) ? item.source_refs : [];
  check(Array.isArray(item.source_refs), "global cost source_refs is not an array");
  check(new Set(sourceRefs).size === sourceRefs.length, "global cost source_refs contains duplicates");
  for (const reference of sourceRefs) {
    check(evidenceReferenceResolves(reference), "unresolved global-cost evidence reference: " + reference);
    if (item.status === "measured" || item.status === "modeled") {
      check(
        sourceArtifactPaths.has(evidenceReferencePath(reference)),
        `measured/modeled cost evidence is not SHA-pinned: ${reference}`,
      );
    }
  }
  if (item.status === "measured" || item.status === "modeled") {
    check(sourceRefs.length > 0, "measured/modeled cost phase lacks evidence");
  }
}
const reportArtifacts = isPlainObject(assessment.report) ? assessment.report : {};
check(isPlainObject(assessment.report), "report is not an object");
check(
  canonicalJson(sortedObjectKeys(reportArtifacts)) ===
    canonicalJson(
      [
        "markdown",
        "pdf",
        "graph_source",
        "graph_svg",
        "graph_png",
        "context_graph_source",
        "context_graph_svg",
        "context_graph_png",
        "release_receipt",
      ].sort(),
    ),
  "report artifact key set changed",
);
const hashBoundReportKinds = new Set([
  "markdown",
  "pdf",
  "graph_source",
  "graph_svg",
  "graph_png",
  "context_graph_source",
  "context_graph_svg",
  "context_graph_png",
  "release_receipt",
]);
const reportCoverageMarker =
  "cover-transfer-assessments/v2; " +
  [
    "catalog_maps",
    "unique_cover_models",
    "geometry_verified_maps",
    "ec1_representations",
    "assessed_maps",
    "unassessed_maps",
    "assessed_subgroups",
    "unassessed_subgroups",
  ]
    .map((key) => `${key}=${coverage[key]}`)
    .join("; ");
const staleReportPatterns = [
  /the other 2,089 map rows and 662 EC1 representations remain unassessed/iu,
  /2,089 map rows unassessed here/iu,
  /662 other EC1 representations/iu,
  /assess the other 2,089 current map rows/iu,
  /assess the remaining 2,089 maps and 662 EC1 representations/iu,
];
const staleCurrentPopulationPatterns = [
  /current 2,091-row/iu,
  /2,091 current (?:map rows|models|inputs)/iu,
  /1,931-row odd catalog population/iu,
  /577 direct(?:-prime| short-Weierstrass)/iu,
  /737 (?:forms inside|models total)/iu,
  /737\/1,354 result/iu,
  /all 2,091 current inputs/iu,
  /all 2,091 model forms/iu,
];
for (const [kind, relativePath] of Object.entries(reportArtifacts)) {
  const absolute = resolveRepoPath(relativePath);
  check(absolute !== null, "unsafe report artifact path for " + kind);
  check(absolute !== null && fs.existsSync(absolute), "missing report artifact " + String(relativePath));
  if (hashBoundReportKinds.has(kind)) {
    check(
      sourceArtifactPaths.has(relativePath),
      `report artifact ${kind} is not bound by a source_artifacts SHA-256 pin`,
    );
  }
  if (absolute === null || !fs.existsSync(absolute)) {
    continue;
  }
  const bytes = fs.readFileSync(absolute);
  if (kind === "markdown" || kind === "graph_source") {
    const contents = bytes.toString("utf8");
    check(
      contents.includes(reportCoverageMarker),
      `${kind} report lacks the exact v2 coverage marker: ${reportCoverageMarker}`,
    );
    for (const stalePattern of staleReportPatterns) {
      check(!stalePattern.test(contents), `${kind} report retains stale partial-assessment language`);
    }
    for (const stalePattern of staleCurrentPopulationPatterns) {
      check(!stalePattern.test(contents), `${kind} report retains a stale current-population count`);
    }
  }
  if (kind === "context_graph_source") {
    const contents = bytes.toString("utf8");
    for (const marker of [
      "318 fields / 2,098 Jacobian contexts / 2,100 maps",
      "current-source test: 2,100 / 2,100",
      "746 inside / 1,354 outside old allowlist",
    ]) {
      check(contents.includes(marker), `context graph lacks current marker: ${marker}`);
    }
    for (const stalePattern of staleCurrentPopulationPatterns) {
      check(!stalePattern.test(contents), "context graph retains a stale current-population count");
    }
  }
  if (kind === "release_receipt") {
    const contents = bytes.toString("utf8");
    for (const marker of [
      "current 2,100-row",
      "586 direct short-Weierstrass prime-field models",
      "2,100 total",
      "746/1,354 result",
      "13 passed; 0 failed",
    ]) {
      check(contents.includes(marker), `release receipt lacks current marker: ${marker}`);
    }
    for (const stalePattern of staleCurrentPopulationPatterns) {
      check(!stalePattern.test(contents), "release receipt retains a stale current-population count");
    }
  }
  if (kind === "pdf") {
    check(bytes.subarray(0, 5).toString("ascii") === "%PDF-", "report PDF signature is invalid");
    check(
      pdfLiteralMetadata(bytes, "Subject") === reportCoverageMarker,
      "report PDF lacks the embedded exact v2 coverage marker",
    );
  } else if (kind === "graph_png" || kind === "context_graph_png") {
    check(
      bytes.subarray(0, 8).toString("hex") === "89504e470d0a1a0a",
      "report PNG signature is invalid",
    );
  } else if (kind === "graph_svg" || kind === "context_graph_svg") {
    check(/<svg\b/u.test(bytes.toString("utf8")), "report SVG is invalid");
  } else if (kind === "graph_source" || kind === "context_graph_source") {
    check(/\bdigraph\b/u.test(bytes.toString("utf8")), "report graph source is invalid");
  }
}

const summary = {
  root,
  problem_count: problems.length,
  problems: problems.slice(0, 100),
  counts: {
    registry_rows: registry.curves.length,
    cover_report_rows: report.curves.length,
    graph_curve_nodes: graph.curves.length,
    maps: Object.keys(graph.maps).length,
    unique_covers: Object.keys(graph.covers).length,
    representations: representationCount,
    unique_ec1: ec1Aliases.size,
    shared_covers: [...coverReferenceCounts.values()].filter((count) => count === 2).length,
  },
  construction_counts: constructionCounts,
  arithmetic_checks: arithmeticChecks,
  profile_counts: {
    map: mapProfilesByPredicate.size,
    subgroup: subgroupProfilesByPredicate.size,
    exact_overlays: exactOverlays.length,
  },
  typed_binding_counts: {
    fields: expectedFieldBindings.size,
    jacobians: expectedJacobianBindings.size,
    maps: expectedTypedMapBindings.size,
  },
  effective_status_counts: effectiveStatusCounts,
  scope_note:
    "Checks the complete v2 membership ledger, explicit field/Jacobian/map bindings, graph/content joins, and cover formulas " +
    "over every declared field. It does not independently prove odd p primality, " +
    "extension-modulus irreducibility, EC1 subgroup orders, or transfer advantage.",
};

console.log(JSON.stringify(summary, null, 2));
if (problems.length !== 0) {
  process.exitCode = 1;
}
