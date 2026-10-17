"""Exact native Boolean-domain MQ search in characteristic two.

Coefficients may belong to any supported polynomial-basis GF(2**d), while
variables range over precisely {0, 1}. Positive powers therefore collapse
to exponent one. Equations must be quadratic after this collapse and exact
coefficient cancellation. Packed Gray-code evaluation enumerates every
assignment; this is a faster exhaustive kernel, not a new complexity bound.

The optional C++17 compiler runs only when the binary is absent or stale.
Call ``build_native_boolean`` before timing installed execution. A solve that
triggers a build includes it in end-to-end wall time and reports its duration.
``max_work`` counts assignments, not instructions, elapsed time, or bytes.
"""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from time import perf_counter

from .binaryfield import BinaryField
from .fields import ExtensionField, PrimeField


_NATIVE_DIR = Path(__file__).resolve().parent.parent / "native"
_SOURCE = _NATIVE_DIR / "boolean_mq.cpp"
_BINARY = _NATIVE_DIR / "boolean_mq"
_WORD_MASK = (1 << 64) - 1


def _width(field):
    if not isinstance(field, (BinaryField, PrimeField, ExtensionField)):
        raise ValueError("Native Boolean MQ requires a supported polynomial-basis finite field")
    if field.characteristic != 2:
        raise ValueError("Native Boolean MQ requires characteristic two")
    return field.degree


def _prepare(equations, field, n_variables=None):
    width = _width(field)
    if not isinstance(equations, (list, tuple)):
        raise ValueError("Equations must be a list or tuple of sparse equations")
    prepared = []
    for equation in equations:
        if not isinstance(equation, (list, tuple)):
            raise ValueError("Each equation must be a list or tuple of sparse terms")
        merged = {}
        for term in equation:
            if not isinstance(term, (list, tuple)) or len(term) != 2:
                raise ValueError("Each term must contain a coefficient and exponent tuple")
            coefficient, powers = term
            if type(coefficient) is not int:
                raise ValueError("Coefficients must be canonical integer field elements")
            field.validate(coefficient)
            if (not isinstance(powers, (list, tuple)) or not powers or
                    any(type(exponent) is not int or exponent < 0 for exponent in powers)):
                raise ValueError("Powers must contain one nonnegative integer per variable")
            if n_variables is None:
                n_variables = len(powers)
            if not 1 <= n_variables <= 62:
                raise ValueError("Native Boolean MQ supports 1..62 variables")
            if len(powers) != n_variables:
                raise ValueError("Exponent tuple length must match n_variables")
            support = tuple(i for i, exponent in enumerate(powers) if exponent)
            merged[support] = merged.get(support, 0) ^ coefficient
        retained = [(coefficient, support) for support, coefficient in merged.items() if coefficient]
        if any(len(support) > 2 for _, support in retained):
            raise ValueError("Native Boolean MQ supports degree at most two after Boolean exponent collapse")
        prepared.append(retained)
    return prepared, width


def supports_boolean_mq(equations, field, domain):
    """Return whether the exact canonical domain and equations fit this backend."""
    if (not isinstance(domain, (list, tuple)) or len(domain) != 2 or
            any(type(value) is not int for value in domain) or list(domain) != [0, 1]):
        return False
    try:
        _prepare(equations, field)
    except (TypeError, ValueError):
        return False
    return True


def _build():
    if _BINARY.is_file() and _BINARY.stat().st_mtime_ns >= _SOURCE.stat().st_mtime_ns:
        return _BINARY, 0.0
    compiler = shutil.which("g++") or shutil.which("c++")
    if compiler is None:
        raise RuntimeError("Native Boolean MQ requires a C++17 compiler (g++ or c++)")
    started = perf_counter()
    with tempfile.TemporaryDirectory(prefix=".boolean-build-", dir=_NATIVE_DIR) as temporary:
        target = Path(temporary) / "boolean_mq"
        command = [compiler, "-O3", "-DNDEBUG", "-std=c++17", str(_SOURCE), "-o", str(target)]
        process = subprocess.run(command, capture_output=True, text=True, check=False)
        if process.returncode:
            raise RuntimeError("Native Boolean MQ build failed: " + process.stderr.strip())
        os.replace(target, _BINARY)
    return _BINARY, perf_counter() - started


def build_native_boolean():
    """Build once and return the executable path; install work precedes timings."""
    return _build()[0]


def _serialize(prepared, n_variables, width, max_solutions):
    # A coordinate bit that is zero in every coefficient of an equation is
    # zero for every Boolean assignment. Omit precisely these identically
    # zero bit planes. Keep all active planes; no sampling or probabilistic
    # equation reduction is performed. In particular, GF(2**83) coefficients
    # in {0,1} need only one packed bit per equation, as they do over GF(2).
    masks = []
    for equation in prepared:
        active = 0
        for coefficient, _ in equation:
            active |= coefficient
        masks.append(active)
    packed_bits = sum(mask.bit_count() for mask in masks)
    words = max(1, (packed_bits + 63) // 64)
    constant = 0
    linear = [0] * n_variables
    quadratic = {}
    shift = 0
    full_mask = (1 << width) - 1
    for equation, active in zip(prepared, masks):
        indices = [i for i in range(width) if active & (1 << i)]
        for coefficient, support in equation:
            if active == full_mask:
                coefficient <<= shift
            else:
                coefficient = sum(((coefficient >> bit) & 1) << (shift + index)
                                  for index, bit in enumerate(indices))
            if not support:
                constant ^= coefficient
            elif len(support) == 1:
                linear[support[0]] ^= coefficient
            else:
                quadratic[support] = quadratic.get(support, 0) ^ coefficient
        shift += len(indices)

    def encoded(vector):
        return " ".join(str((vector >> (64 * i)) & _WORD_MASK) for i in range(words))

    lines = [f"GRBOOLMQ1 {n_variables} {words} {max_solutions} {1 << n_variables}",
             encoded(constant)]
    lines.extend(encoded(vector) for vector in linear)
    lines.extend(encoded(quadratic.get((j, k), 0))
                 for k in range(n_variables) for j in range(k))
    return "\n".join(lines) + "\n", words, packed_bits


def solve_boolean_mq(equations, n_variables, field, *, max_work=1_000_000,
                     max_solutions=10_000):
    """Exactly solve Boolean-domain quadratic equations without a Groebner basis.

    Work exceeding the cap is refused before compilation or process startup.
    Output truncation still returns an exact solution count. Retained roots
    follow Gray enumeration and are sorted lexicographically for reporting.
    This backend supports 1..62 variables; larger inputs are explicitly refused.
    """
    if type(n_variables) is not int or not 1 <= n_variables <= 62:
        raise ValueError("Native Boolean MQ n_variables must be between 1 and 62")
    if type(max_work) is not int or max_work < 0:
        raise ValueError("max_work must be a nonnegative integer")
    if type(max_solutions) is not int or max_solutions < 0:
        raise ValueError("max_solutions must be a nonnegative integer")
    prepared, width = _prepare(equations, field, n_variables)
    total = 1 << n_variables
    report = {
        "method": "native_boolean_gray", "status": "budget_exceeded",
        "solutions": [], "solution_count_exact": None,
        "work_estimate": total, "table_assignment_estimate": 0,
        "assignments_visited": 0, "native_compilation_seconds": 0.0,
        "reason": f"The estimate of {total} assignments exceeds max_work={max_work}; search was not started.",
    }
    if total > max_work:
        return report
    # Never ask the kernel to retain more roots than exist or than uint64 can encode.
    payload, words, packed_bits = _serialize(prepared, n_variables, width, min(total, max_solutions))
    binary, compilation_seconds = _build()
    started = perf_counter()
    process = subprocess.run([str(binary)], input=payload, capture_output=True,
                             text=True, check=False)
    process_seconds = perf_counter() - started
    if process.returncode:
        raise RuntimeError("Native Boolean MQ failed: " + process.stderr.strip())
    parts = process.stdout.split()
    if len(parts) < 4 or parts[0] != "GRBOOLMQ1":
        raise RuntimeError("Native Boolean MQ returned an invalid header")
    try:
        count, visited, retained = map(int, parts[1:4])
        masks = [int(mask) for mask in parts[4:]]
    except ValueError as error:
        raise RuntimeError("Native Boolean MQ returned malformed integers") from error
    if (not 0 <= count <= total or visited != total or retained != len(masks) or
            retained != min(max_solutions, count) or len(set(masks)) != len(masks) or
            any(not 0 <= mask < total for mask in masks)):
        raise RuntimeError("Native Boolean MQ returned an inconsistent certificate")
    solutions = sorted([[int(bool(mask & (1 << i))) for i in range(n_variables)]
                        for mask in masks])
    # Recheck every retained tuple against the collapsed, exact input equations.
    # Exhaustive completeness follows from the independently tested Gray invariant.
    for mask in masks:
        for equation in prepared:
            value = 0
            for coefficient, support in equation:
                if all(mask & (1 << i) for i in support):
                    value ^= coefficient
            if value:
                raise RuntimeError("Native Boolean MQ returned a tuple that fails exact verification")
    report.update(status="complete" if count == retained else "output_truncated",
                  solutions=solutions, solution_count_exact=count,
                  assignments_visited=total, packed_words=words,
                  packed_bits=packed_bits,
                  omitted_zero_bitplanes=len(prepared) * width - packed_bits,
                  native_compilation_seconds=compilation_seconds,
                  native_process_seconds=process_seconds,
                  reason="Every Boolean assignment was enumerated by exact packed Gray-code evaluation; retained tuples were verified.")
    return report
