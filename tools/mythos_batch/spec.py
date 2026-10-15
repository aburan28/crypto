"""Request files: parsing, defaults and local validation.

Everything here runs without network or credentials.  The checks encode the
constraints of ``claude-mythos-5-1`` so that a bad row fails before a billable
batch is created: a batch is accepted or rejected as a whole only after a long
asynchronous round trip, and individual invalid rows surface hours later.
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field
from typing import Any, Iterable

# Claude Mythos 5.1 (Project Glasswing). Override with MYTHOS_MODEL or --model.
DEFAULT_MODEL = "claude-mythos-5-1"

# Thinking is always on and counts toward max_tokens, so leave room for it.
# Batches are asynchronous: no HTTP timeout pushes this down.
DEFAULT_MAX_TOKENS = 32_000
MAX_OUTPUT_TOKENS = 128_000
EFFORTS = ("low", "medium", "high", "xhigh", "max")

MAX_REQUESTS_PER_BATCH = 100_000
# API limit is 256 MB per batch; stay under it with headroom for envelope bytes.
MAX_BATCH_BYTES = 200 * 1024 * 1024

CUSTOM_ID_RE = re.compile(r"^[A-Za-z0-9_-]{1,64}$")
ROW_KEYS = {"custom_id", "prompt", "messages", "system", "max_tokens", "effort", "params"}
# Rejected on Mythos 5.1 (400 from the API) or silently ignored; fail fast.
FORBIDDEN_PARAMS = {
    "temperature": "sampling parameters are not accepted by Claude Mythos 5.1",
    "top_p": "sampling parameters are not accepted by Claude Mythos 5.1",
    "top_k": "sampling parameters are not accepted by Claude Mythos 5.1",
    "stream": "batch requests are non-streaming",
    "model": "set the model with --model / MYTHOS_MODEL, not per row",
    "max_tokens": "use the row-level max_tokens key",
}


class RequestError(ValueError):
    """One or more rows failed validation. ``problems`` is a list of strings."""

    def __init__(self, problems: list[str]):
        self.problems = problems
        super().__init__("; ".join(problems[:5]) + (" ..." if len(problems) > 5 else ""))


@dataclass
class Defaults:
    model: str = field(default_factory=lambda: os.environ.get("MYTHOS_MODEL") or DEFAULT_MODEL)
    max_tokens: int = field(
        default_factory=lambda: int(os.environ.get("MYTHOS_MAX_TOKENS") or DEFAULT_MAX_TOKENS)
    )
    effort: str | None = None
    system: str | None = None


def _check_params(extra: dict[str, Any], where: str, problems: list[str]) -> None:
    for key, why in FORBIDDEN_PARAMS.items():
        if key in extra:
            problems.append(f"{where}: params.{key}: {why}")
    thinking = extra.get("thinking")
    if thinking is not None:
        kind = thinking.get("type") if isinstance(thinking, dict) else None
        if kind != "adaptive":
            problems.append(
                f"{where}: params.thinking must be omitted or {{'type': 'adaptive'}} "
                "(always-on thinking; disabled/enabled+budget_tokens are a 400)"
            )
    choice = extra.get("tool_choice")
    if isinstance(choice, dict) and choice.get("type") in ("any", "tool"):
        problems.append(
            f"{where}: params.tool_choice type {choice.get('type')!r} is a 400 on Claude "
            "Mythos 5.1; use 'auto' (or 'none') and name the tool in the prompt"
        )


def normalize_row(row: dict[str, Any], defaults: Defaults, where: str) -> dict[str, Any]:
    """Validate one input row and return ``{"custom_id", "params"}`` for the API."""
    problems: list[str] = []
    if not isinstance(row, dict):
        raise RequestError([f"{where}: row must be a JSON object"])
    unknown = sorted(set(row) - ROW_KEYS)
    if unknown:
        problems.append(f"{where}: unknown keys {unknown} (allowed: {sorted(ROW_KEYS)})")

    custom_id = row.get("custom_id")
    if not isinstance(custom_id, str) or not CUSTOM_ID_RE.match(custom_id):
        problems.append(f"{where}: custom_id must match [A-Za-z0-9_-]{{1,64}}, got {custom_id!r}")

    has_prompt, has_messages = "prompt" in row, "messages" in row
    if has_prompt == has_messages:
        problems.append(f"{where}: give exactly one of 'prompt' or 'messages'")
        messages: list[Any] = []
    elif has_prompt:
        if not isinstance(row["prompt"], str) or not row["prompt"].strip():
            problems.append(f"{where}: prompt must be a non-empty string")
        messages = [{"role": "user", "content": row["prompt"]}]
    else:
        messages = row["messages"]
        if not isinstance(messages, list) or not messages:
            problems.append(f"{where}: messages must be a non-empty list")
            messages = []
        elif not isinstance(messages[-1], dict) or messages[-1].get("role") != "user":
            problems.append(
                f"{where}: last message must have role 'user' (assistant prefill is a 400 "
                "on Claude Mythos 5.1)"
            )

    max_tokens = row.get("max_tokens", defaults.max_tokens)
    if not isinstance(max_tokens, int) or isinstance(max_tokens, bool) or not (
        1 <= max_tokens <= MAX_OUTPUT_TOKENS
    ):
        problems.append(f"{where}: max_tokens must be an integer in 1..{MAX_OUTPUT_TOKENS}")

    effort = row.get("effort", defaults.effort)
    if effort is not None and effort not in EFFORTS:
        problems.append(f"{where}: effort must be one of {EFFORTS}, got {effort!r}")

    extra = row.get("params", {})
    if not isinstance(extra, dict):
        problems.append(f"{where}: params must be an object")
        extra = {}
    _check_params(extra, where, problems)

    system = row.get("system", defaults.system)
    if system is not None and not isinstance(system, (str, list)):
        problems.append(f"{where}: system must be a string or a list of text blocks")

    if problems:
        raise RequestError(problems)

    params: dict[str, Any] = {
        "model": defaults.model,
        "max_tokens": max_tokens,
        "messages": messages,
    }
    if system:
        params["system"] = system
    if effort:
        params["output_config"] = {"effort": effort}
    for key, value in extra.items():
        if key == "output_config" and effort:
            params["output_config"] = {**value, **params["output_config"]}
        else:
            params[key] = value
    return {"custom_id": custom_id, "params": params}


def load_requests(lines: Iterable[str], defaults: Defaults) -> list[dict[str, Any]]:
    """Parse JSONL ``lines``; raise :class:`RequestError` listing every problem."""
    out: list[dict[str, Any]] = []
    problems: list[str] = []
    seen: dict[str, int] = {}
    for number, line in enumerate(lines, 1):
        if not line.strip():
            continue
        where = f"line {number}"
        try:
            row = json.loads(line)
        except json.JSONDecodeError as exc:
            problems.append(f"{where}: invalid JSON ({exc.msg})")
            continue
        try:
            item = normalize_row(row, defaults, where)
        except RequestError as exc:
            problems.extend(exc.problems)
            continue
        cid = item["custom_id"]
        if cid in seen:
            problems.append(f"{where}: duplicate custom_id {cid!r} (first on line {seen[cid]})")
            continue
        seen[cid] = number
        out.append(item)
    if not out and not problems:
        problems.append("no requests found")
    if problems:
        raise RequestError(problems)
    return out


def single_prompt(prompt: str, defaults: Defaults) -> list[dict[str, Any]]:
    return [normalize_row({"custom_id": "prompt-1", "prompt": prompt}, defaults, "--prompt")]


def chunk_requests(requests: list[dict[str, Any]]) -> list[list[dict[str, Any]]]:
    """Split into groups that respect the per-batch request-count and byte limits."""
    chunks: list[list[dict[str, Any]]] = []
    current: list[dict[str, Any]] = []
    size = 0
    for item in requests:
        item_size = len(json.dumps(item, separators=(",", ":")).encode()) + 1
        if item_size > MAX_BATCH_BYTES:
            raise RequestError([f"request {item['custom_id']!r} alone exceeds the batch size limit"])
        if current and (len(current) >= MAX_REQUESTS_PER_BATCH or size + item_size > MAX_BATCH_BYTES):
            chunks.append(current)
            current, size = [], 0
        current.append(item)
        size += item_size
    if current:
        chunks.append(current)
    return chunks


def build_params(prompt: str, **overrides: Any) -> dict[str, Any]:
    """Convenience for library users: the API params for a single prompt."""
    return single_prompt(prompt, Defaults(**overrides))[0]["params"]
