"""A dependency-free validator for the JSON Schema subset the ICMS schemas use.

Supported keywords: type (string or list), enum, const, required, properties,
additionalProperties (bool or schema), patternProperties, items, minItems,
maxItems, uniqueItems, minimum, maximum, exclusiveMinimum, minLength, pattern,
oneOf, anyOf, allOf, if/then/else, $ref to "#/$defs/<name>", and $defs.
Unsupported keywords are refused when the schema is loaded, so a schema can
never silently lose a constraint.
"""
from __future__ import annotations

import re

SUPPORTED = {
    "$schema", "$id", "$defs", "$ref", "$comment", "title", "description", "examples", "default",
    "type", "enum", "const", "required", "properties", "additionalProperties", "patternProperties",
    "items", "minItems", "maxItems", "uniqueItems", "minimum", "maximum", "exclusiveMinimum",
    "minLength", "pattern", "oneOf", "anyOf", "allOf", "if", "then", "else", "x-icms",
}

_TYPES = {
    "object": dict, "array": list, "string": str, "boolean": bool, "null": type(None),
}


def _canon(value):
    """A JSON-faithful key: booleans stay apart from numbers (True is not 1),
    1 and 1.0 are one number, and object key order does not matter."""
    if isinstance(value, bool):
        return ("bool", value)
    if isinstance(value, (int, float)):
        return ("num", float(value)) if isinstance(value, float) and not value.is_integer() else ("num", int(value))
    if isinstance(value, dict):
        return ("obj", tuple(sorted((k, _canon(v)) for k, v in value.items())))
    if isinstance(value, list):
        return ("arr", tuple(_canon(v) for v in value))
    return (type(value).__name__, value)


def json_equal(a, b) -> bool:
    return _canon(a) == _canon(b)


def _is_type(value, name: str) -> bool:
    if name == "integer":
        return isinstance(value, int) and not isinstance(value, bool)
    if name == "number":
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    return isinstance(value, _TYPES[name])


def audit_schema(schema, path: str = "#") -> list[str]:
    """Return every unsupported keyword in the schema tree."""
    bad = []
    if isinstance(schema, dict):
        for k, v in schema.items():
            if k not in SUPPORTED:
                bad.append(f"{path}: unsupported keyword {k!r}")
            if k in ("properties", "patternProperties", "$defs"):
                for name, sub in v.items():
                    bad += audit_schema(sub, f"{path}/{k}/{name}")
            elif k in ("items", "additionalProperties", "if", "then", "else") and isinstance(v, dict):
                bad += audit_schema(v, f"{path}/{k}")
            elif k == "items":
                bad.append(f"{path}: only the single-schema form of items is supported")
            elif k in ("oneOf", "anyOf", "allOf"):
                for i, sub in enumerate(v):
                    bad += audit_schema(sub, f"{path}/{k}/{i}")
    return bad


class Validator:
    def __init__(self, schema: dict):
        problems = audit_schema(schema)
        if problems:
            raise ValueError("schema uses unsupported keywords: " + "; ".join(problems))
        self.root = schema

    def _resolve(self, ref: str) -> dict:
        if not ref.startswith("#/$defs/"):
            raise ValueError(f"only local $defs refs are supported: {ref}")
        return self.root["$defs"][ref[len("#/$defs/"):]]

    def errors(self, value, schema: dict | None = None, path: str = "$") -> list[str]:
        s = self.root if schema is None else schema
        if "$ref" in s:
            return self.errors(value, self._resolve(s["$ref"]), path)
        out: list[str] = []
        if "type" in s:
            types = s["type"] if isinstance(s["type"], list) else [s["type"]]
            if not any(_is_type(value, t) for t in types):
                return [f"{path}: expected {'/'.join(types)}, got {type(value).__name__}"]
        if "const" in s and not json_equal(value, s["const"]):
            out.append(f"{path}: must equal {s['const']!r}")
        if "enum" in s and not any(json_equal(value, e) for e in s["enum"]):
            out.append(f"{path}: {value!r} not in {s['enum']}")
        if isinstance(value, str):
            if "minLength" in s and len(value) < s["minLength"]:
                out.append(f"{path}: shorter than {s['minLength']}")
            if "pattern" in s and not re.search(s["pattern"], value):
                out.append(f"{path}: {value!r} does not match {s['pattern']}")
        if _is_type(value, "number"):
            if "minimum" in s and value < s["minimum"]:
                out.append(f"{path}: {value} < minimum {s['minimum']}")
            if "maximum" in s and value > s["maximum"]:
                out.append(f"{path}: {value} > maximum {s['maximum']}")
            if "exclusiveMinimum" in s and value <= s["exclusiveMinimum"]:
                out.append(f"{path}: {value} <= exclusiveMinimum {s['exclusiveMinimum']}")
        if isinstance(value, list):
            if "minItems" in s and len(value) < s["minItems"]:
                out.append(f"{path}: fewer than {s['minItems']} items")
            if "maxItems" in s and len(value) > s["maxItems"]:
                out.append(f"{path}: more than {s['maxItems']} items")
            if s.get("uniqueItems"):
                seen = [_canon(v) for v in value]
                if len(seen) != len(set(seen)):
                    out.append(f"{path}: items are not unique")
            if isinstance(s.get("items"), dict):
                for i, item in enumerate(value):
                    out += self.errors(item, s["items"], f"{path}[{i}]")
        if isinstance(value, dict):
            for req in s.get("required", []):
                if req not in value:
                    out.append(f"{path}: missing required {req!r}")
            props = s.get("properties", {})
            pats = s.get("patternProperties", {})
            for k, v in value.items():
                # properties and every matching patternProperties both apply;
                # additionalProperties applies only when neither did.
                matched = k in props
                if matched:
                    out += self.errors(v, props[k], f"{path}.{k}")
                for pat, sub in pats.items():
                    if re.search(pat, k):
                        matched = True
                        out += self.errors(v, sub, f"{path}.{k}")
                if matched:
                    continue
                ap = s.get("additionalProperties", True)
                if ap is False:
                    out.append(f"{path}: unknown key {k!r}")
                elif isinstance(ap, dict):
                    out += self.errors(v, ap, f"{path}.{k}")
        for sub in s.get("allOf", []):
            out += self.errors(value, sub, path)
        if "anyOf" in s and not any(not self.errors(value, sub, path) for sub in s["anyOf"]):
            out.append(f"{path}: matches none of anyOf")
        if "oneOf" in s:
            n = sum(1 for sub in s["oneOf"] if not self.errors(value, sub, path))
            if n != 1:
                out.append(f"{path}: matches {n} of oneOf, expected exactly 1")
        if "if" in s and not self.errors(value, s["if"], path):
            if "then" in s:
                out += self.errors(value, s["then"], path)
        elif "if" in s and "else" in s:
            out += self.errors(value, s["else"], path)
        return out
