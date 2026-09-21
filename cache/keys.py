"""Key schema for cached index-calculus precomputation artifacts.

A cache key has one job: name the bytes.  If two computations can produce
different bytes under the same key, the cache is a correctness bug, not a
speedup.  So everything that can change the bytes goes in the key --

  * the curve, down to its *representation* and not just its isomorphism
    class (`curve_fingerprint`),
  * the artifact parameters (which summation polynomial, which monomial
    order, ...),
  * the serializer that produced the bytes (`codec`), because a caller
    that swaps its algebra backend swaps the encoding with it,
  * `SCHEMA_VERSION`, the kill switch for everything above that we got
    wrong or later change.

Representation, not isomorphism class.  `F_2[z]/(z^9+z+1)` and
`F_2[z]/(z^9+z^4+1)` are the same field up to isomorphism and wholly
different as byte strings, and this repository carries both kinds of
choice in the same tree: `scripts/ecc2k130_point_decomposition.py` works in
a polynomial basis with an explicit reduction polynomial, while
`ecc2k130/codegen/decomp.py` works in the permuted type-II optimal normal
basis.  `src/bin/ic/params.rs` distinguishes them as the `binary` and
`binary_abstract` field kinds, and this module consumes that same document
shape so the distinction survives into the key.

What is deliberately *not* in a curve fingerprint: the subgroup order, the
cofactor, the generator, the target point, the factor base, and the curve's
human-readable name.  A summation polynomial depends on the curve equation
and the field alone, so folding an instance into its key would silently
turn a fleet-wide artifact into a per-run one.  The name is excluded for
the same reason in reverse: `ECC2K-95` and `ecc2k-95` are one curve.

Key layout, with the Redis hash tag around the fingerprint so that one
curve's artifacts land on one shard and a warming worker does not fan out
across a cluster:

    redis   ic/v1/{f0e1...}/sumpoly/codec=gz-json/form=symmetric/n=5
    s3      indexcalc/v1/f0e1.../sumpoly/codec=gz-json/form=symmetric/n=5

Parameters are emitted in sorted order, so a caller that passes them in a
different order still hits the same entry.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import re
from typing import Any, Iterable, Mapping

# Bump on any change to a serialization, to a fingerprint input, or to an
# algorithm whose output is cached.  Old keys are then never read again and
# age out on their own.  Never mutate an artifact in place under a live key:
# a reader that already holds the old bytes has no way to find out.
SCHEMA_VERSION = 1

KEY_PREFIX = "ic"

# Namespaces.  One per artifact family; they do not share parameter sets.
NS_SUMMATION = "sumpoly"
NS_GROEBNER = "gbstruct"
NS_SPECIALIZATION = "spec"

FINGERPRINT_BITS = 128
_FINGERPRINT_CHARS = FINGERPRINT_BITS // 4

# Redis hash tags are delimited by the first `{` and the first following `}`,
# so neither may appear anywhere else in the key.  S3 tolerates both but asks
# for special handling, which is why `s3_key` drops the tag entirely.
_SEGMENT_OK = re.compile(r"\A[A-Za-z0-9._:+-]+\Z")


class CacheKeyError(ValueError):
    """A key could not be built from the inputs given."""


# -- integers -------------------------------------------------------------

def canonical_int(value: Any, *, what: str = "integer") -> int:
    """Read the repository's integer encoding: `0x`-prefixed hex, else decimal.

    This mirrors `number()` in `src/bin/ic/params.rs`.  It matters here
    because `docs/ic/binary-example.json` writes coordinates as `"0x160"`
    and coefficients as `"0"`; without a canonical form, `"0"` and `"0x0"`
    would fingerprint the same curve twice and simply never share a cache
    entry.
    """
    if isinstance(value, bool):
        raise CacheKeyError(f"{what}: bool is not an integer")
    if isinstance(value, int):
        parsed = value
    elif isinstance(value, str):
        text = value.strip()
        if not text:
            raise CacheKeyError(f"{what}: empty string")
        try:
            parsed = int(text[2:], 16) if text.startswith("0x") else int(text, 10)
        except ValueError as exc:
            raise CacheKeyError(f"{what}: {value!r} is not an integer") from exc
    else:
        raise CacheKeyError(f"{what}: unsupported type {type(value).__name__}")
    if parsed < 0:
        raise CacheKeyError(f"{what}: must be non-negative, got {parsed}")
    return parsed


# -- fields ---------------------------------------------------------------

def canonical_field(field: Mapping[str, Any]) -> tuple:
    """Canonicalise one `field` object from an ic parameter document.

    Accepts the three kinds `src/bin/ic/params.rs` defines:

        {"kind": "binary",          "degree": 9,   "polynomial_terms": [0, 1, 9]}
        {"kind": "binary_abstract", "degree": 131, "representation": "..."}
        {"kind": "prime",           "modulus": "97"}

    Returns a hashable tuple.  `polynomial_terms` is sorted, because
    `params.rs` builds it as `low_terms` with the degree pushed on the end
    (`[6, 0, 97]` for ECC2K-95) while `docs/ic/binary-example.json` writes
    it ascending (`[0, 1, 9]`).  Two spellings of one field must not read as
    two fields -- that costs a recomputed Groebner basis, not correctness,
    but it costs it every time.
    """
    if not isinstance(field, Mapping):
        raise CacheKeyError("field must be a mapping")
    kind = field.get("kind")
    if kind == "binary":
        degree = canonical_int(field["degree"], what="field degree")
        terms = sorted({canonical_int(t, what="polynomial term")
                        for t in field["polynomial_terms"]})
        if degree < 1:
            raise CacheKeyError("binary field degree must be positive")
        if not terms or terms[-1] != degree:
            raise CacheKeyError(
                f"polynomial_terms {terms} must include the degree {degree} "
                "as its leading term")
        if terms[0] != 0:
            raise CacheKeyError(
                f"polynomial_terms {terms} must include the constant term 0")
        return ("binary", degree, tuple(terms))
    if kind == "binary_abstract":
        degree = canonical_int(field["degree"], what="field degree")
        representation = field.get("representation")
        if not isinstance(representation, str) or not representation.strip():
            raise CacheKeyError("binary_abstract needs a non-empty representation")
        # Collapse whitespace only.  Case and wording are load-bearing here:
        # the string is the sole record of which basis the bytes are in, and
        # this module cannot tell two descriptions apart any other way.
        return ("binary_abstract", degree, " ".join(representation.split()))
    if kind == "prime":
        modulus = canonical_int(field["modulus"], what="field modulus")
        if modulus < 2:
            raise CacheKeyError("prime field modulus must be at least 2")
        return ("prime", modulus)
    raise CacheKeyError(f"unknown field kind {kind!r}")


def field_order_bound(field: tuple) -> int:
    """The exclusive upper bound on a canonical coefficient in `field`."""
    if field[0] in ("binary", "binary_abstract"):
        return 1 << field[1]
    return field[1]


def canonical_coefficient(value: Any, field: tuple, *, what: str) -> int:
    """Reduce a curve coefficient into its field's canonical range.

    Prime fields are reduced mod `p`, since `a = 2` and `a = 99` over
    `F_97` are one curve and should share one cache entry.  Binary fields
    are *not* reduced: an out-of-range value there means the caller is
    holding a different representation than it thinks it is, and a loud
    failure beats a key that quietly names the wrong bytes.
    """
    raw = canonical_int(value, what=what)
    bound = field_order_bound(field)
    if field[0] == "prime":
        return raw % bound
    if raw >= bound:
        raise CacheKeyError(
            f"{what}: {raw} is not reduced in a degree-{field[1]} binary field")
    return raw


# -- fingerprints ---------------------------------------------------------

def _digest(parts: Iterable[Any]) -> str:
    h = hashlib.sha256()
    for part in parts:
        h.update(repr(part).encode("utf-8"))
        h.update(b"\x1f")
    return h.hexdigest()[:_FINGERPRINT_CHARS]


def curve_fingerprint(curve: Mapping[str, Any]) -> str:
    """Fingerprint the curve equation and its representation.

    `curve` is an ic parameter document, or any mapping carrying its
    `field`, `a` and `b` entries -- `docs/ic/binary-example.json` and the
    profiles in `src/bin/ic/params.rs` both work as-is.  Everything
    instance-specific in such a document is ignored; see the module
    docstring for why.
    """
    if not isinstance(curve, Mapping):
        raise CacheKeyError("curve must be a mapping")
    for required in ("field", "a", "b"):
        if required not in curve:
            raise CacheKeyError(f"curve is missing {required!r}")
    field = canonical_field(curve["field"])
    a = canonical_coefficient(curve["a"], field, what="curve a")
    b = canonical_coefficient(curve["b"], field, what="curve b")
    return _digest((SCHEMA_VERSION, field, a, b))


# Scalars whose `repr` is a faithful, stable description of the value.
_ATOMIC = (bool, int, float, str, bytes, bytearray, type(None))


def inputs_fingerprint(*parts: Any) -> str:
    """Fingerprint arbitrary caller-supplied inputs, e.g. a specialization.

    Mappings are sorted by key so that two equal dicts fingerprint alike.
    Containers are tagged by kind, so that `{"a": 1}` and `[("a", 1)]` --
    which have the same canonical shape and are not the same input -- do not
    share an entry; that collision would hand back the wrong artifact, not
    merely a slow one.

    Anything that is not a scalar or a container must define its own
    `__repr__`.  A plain object's repr is its memory address, which is
    unstable across processes and, worse, *reused* within one: two distinct
    objects can land on the same address and so on the same key.
    """
    return _digest(_stable(p) for p in (SCHEMA_VERSION,) + parts)


def _stable(value: Any) -> Any:
    if isinstance(value, _ATOMIC):
        return value
    if isinstance(value, Mapping):
        return ("map", tuple(sorted((str(k), _stable(v)) for k, v in value.items())))
    if isinstance(value, (list, tuple)):
        return ("seq", tuple(_stable(v) for v in value))
    if isinstance(value, (set, frozenset)):
        return ("set", tuple(sorted(repr(_stable(v)) for v in value)))
    if type(value).__repr__ is object.__repr__:
        raise CacheKeyError(
            f"cannot fingerprint {type(value).__name__}: it has no __repr__ of "
            "its own, so it would be keyed by memory address")
    return value


# -- keys -----------------------------------------------------------------

def _segment(value: Any, *, what: str) -> str:
    text = str(value)
    if not _SEGMENT_OK.match(text):
        raise CacheKeyError(
            f"{what}: {text!r} must match {_SEGMENT_OK.pattern} "
            "(no braces, slashes or spaces -- they break Redis hash tags)")
    return text


@dataclass(frozen=True)
class CacheKey:
    """One artifact's address, rendered per tier.

    `params` is normalised to a sorted tuple of `(name, value)` pairs on
    construction, so callers may pass them in any order.
    """

    namespace: str
    curve_fp: str
    params: tuple[tuple[str, str], ...]

    @classmethod
    def build(cls, namespace: str, curve_fp: str, **params: Any) -> "CacheKey":
        ns = _segment(namespace, what="namespace")
        fp = _segment(curve_fp, what="curve fingerprint")
        items = tuple(sorted(
            (_segment(name, what="parameter name"),
             _segment(value, what=f"parameter {name}"))
            for name, value in params.items()))
        return cls(ns, fp, items)

    @property
    def suffix(self) -> str:
        """Everything below the fingerprint: namespace and parameters."""
        tail = "/".join(f"{name}={value}" for name, value in self.params)
        return f"{self.namespace}/{tail}" if tail else self.namespace

    def redis_key(self) -> str:
        """Hash-tagged so a curve's artifacts share a cluster-mode shard."""
        return f"{KEY_PREFIX}/v{SCHEMA_VERSION}/{{{self.curve_fp}}}/{self.suffix}"

    def s3_key(self, prefix: str) -> str:
        """The same address without the hash tag; S3 wants no braces."""
        head = prefix.strip("/")
        path = f"v{SCHEMA_VERSION}/{self.curve_fp}/{self.suffix}"
        return f"{head}/{path}" if head else path

    def l1_key(self) -> tuple:
        return (SCHEMA_VERSION, self.namespace, self.curve_fp, self.params)

    def __str__(self) -> str:
        return self.redis_key()
