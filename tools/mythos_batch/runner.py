"""Submit, wait for and collect Message Batches; state lives in an output directory.

The output directory is the unit of work::

    manifest.json        model, input hash, batch ids (no credentials)
    requests.jsonl       the normalized requests that were submitted
    results.jsonl        one row per request, in input order (after ``collect``)
    retryable.jsonl      input rows worth resubmitting (errored/expired/canceled/missing)
    summary.json         counts, token usage and estimated cost

``submit`` refuses to run twice into the same directory, so a re-run cannot bill a
second batch; ``run`` resumes from the manifest instead.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

from .spec import chunk_requests

# Batch list price per MTok for Claude Mythos 5.1 (50% of $10 / $50). Estimate only:
# cache reads/writes are not priced here.
BATCH_PRICE_PER_MTOK = {"input": 5.00, "output": 25.00}
MANIFEST = "manifest.json"
SCHEMA = 1


class BatchStateError(RuntimeError):
    pass


def make_client(**kwargs: Any) -> Any:
    """Anthropic client from the environment (``ANTHROPIC_API_KEY``, ``ANTHROPIC_BASE_URL``)."""
    import anthropic  # imported lazily: validation and tests need no SDK

    return anthropic.Anthropic(max_retries=kwargs.pop("max_retries", 5), **kwargs)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _atomic_write(path: Path, text: str) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


def _jsonl(rows: list[dict[str, Any]]) -> str:
    return "".join(json.dumps(r, ensure_ascii=False, sort_keys=True) + "\n" for r in rows)


def load_manifest(out: Path) -> dict[str, Any]:
    path = out / MANIFEST
    if not path.exists():
        raise BatchStateError(f"{path} not found; run `submit` first")
    return json.loads(path.read_text())


def request_digest(requests: list[dict[str, Any]]) -> str:
    return hashlib.sha256(_jsonl(requests).encode()).hexdigest()


def _counts(batch: Any) -> dict[str, int]:
    c = batch.request_counts
    return {k: int(getattr(c, k)) for k in ("processing", "succeeded", "errored", "canceled", "expired")}


def submit(client: Any, requests: list[dict[str, Any]], out: Path, *, model: str, log=print) -> dict[str, Any]:
    """Create the batch(es) and record them. Raises if ``out`` already holds a submission."""
    out.mkdir(parents=True, exist_ok=True)
    if (out / MANIFEST).exists():
        raise BatchStateError(
            f"{out / MANIFEST} exists: this directory already has a submission "
            "(use `wait`/`collect`, or pick a new --out)"
        )
    chunks = chunk_requests(requests)
    manifest: dict[str, Any] = {
        "schema": SCHEMA,
        "model": model,
        "created_at": _now(),
        "input_sha256": request_digest(requests),
        "n_requests": len(requests),
        "batches": [],
    }
    _atomic_write(out / "requests.jsonl", _jsonl(requests))
    for i, chunk in enumerate(chunks, 1):
        batch = client.messages.batches.create(requests=chunk)
        manifest["batches"].append(
            {"id": batch.id, "n_requests": len(chunk), "first_custom_id": chunk[0]["custom_id"]}
        )
        # Persist after every batch: an interrupted multi-batch submit stays recoverable.
        _atomic_write(out / MANIFEST, json.dumps(manifest, indent=2) + "\n")
        log(f"submitted batch {i}/{len(chunks)}: {batch.id} ({len(chunk)} requests)")
    return manifest


def status(client: Any, out: Path) -> list[dict[str, Any]]:
    rows = []
    for entry in load_manifest(out)["batches"]:
        batch = client.messages.batches.retrieve(entry["id"])
        rows.append(
            {"id": batch.id, "processing_status": batch.processing_status, **_counts(batch)}
        )
    return rows


def wait(
    client: Any,
    out: Path,
    *,
    poll_seconds: float = 60.0,
    timeout_seconds: float = 25 * 3600,
    sleep: Callable[[float], None] = time.sleep,
    clock: Callable[[], float] = time.monotonic,
    log=print,
) -> list[dict[str, Any]]:
    """Block until every batch has ended. Batches finish in 24 h at most."""
    start = clock()
    while True:
        rows = status(client, out)
        done = sum(r["processing_status"] == "ended" for r in rows)
        log(
            f"[{_now()}] {done}/{len(rows)} batches ended; "
            + ", ".join(f"{r['id']}: {r['processing_status']} "
                        f"(ok {r['succeeded']}, err {r['errored']}, pending {r['processing']})"
                        for r in rows)
        )
        if done == len(rows):
            return rows
        if clock() - start > timeout_seconds:
            raise TimeoutError(
                f"batches still running after {timeout_seconds / 3600:.1f} h; they keep running "
                f"server-side - resume with `wait --out {out}`"
            )
        sleep(poll_seconds)


def _error_fields(err: Any) -> tuple[str, str]:
    inner = getattr(err, "error", err)  # ErrorResponse wraps the actual error
    return str(getattr(inner, "type", "unknown")), str(getattr(inner, "message", ""))


def _usage(msg: Any) -> dict[str, int]:
    u = msg.usage
    return {
        k: int(getattr(u, k, 0) or 0)
        for k in ("input_tokens", "output_tokens", "cache_creation_input_tokens", "cache_read_input_tokens")
    }


def result_row(item: Any) -> dict[str, Any]:
    """One output row from a ``MessageBatchIndividualResponse``."""
    res = item.result
    row: dict[str, Any] = {"custom_id": item.custom_id, "status": res.type}
    if res.type == "succeeded":
        msg = res.message
        row.update(
            stop_reason=msg.stop_reason,
            model=msg.model,
            text="".join(b.text for b in msg.content if b.type == "text"),
            usage=_usage(msg),
        )
        if msg.stop_reason == "refusal":
            # Safeguards declined the request; HTTP 200, so it counts as "succeeded" upstream.
            row["status"] = "refused"
            details = getattr(msg, "stop_details", None)
            row["stop_details"] = (
                {"category": getattr(details, "category", None),
                 "explanation": getattr(details, "explanation", None)}
                if details else None
            )
        elif msg.stop_reason == "max_tokens":
            row["truncated"] = True
    elif res.type == "errored":
        etype, emsg = _error_fields(res.error)
        row["error"] = {"type": etype, "message": emsg,
                        "request_id": getattr(res.error, "request_id", None)}
    return row


def collect(client: Any, out: Path, *, log=print) -> dict[str, Any]:
    """Download results of ended batches and write ``results.jsonl`` / ``summary.json``."""
    manifest = load_manifest(out)
    inputs = [json.loads(l) for l in (out / "requests.jsonl").read_text().splitlines() if l.strip()]
    order = [r["custom_id"] for r in inputs]
    by_id: dict[str, dict[str, Any]] = {}
    pending = []
    for entry in manifest["batches"]:
        batch = client.messages.batches.retrieve(entry["id"])
        if batch.processing_status != "ended":
            pending.append(entry["id"])
            continue
        # Results arrive in any order: key by custom_id, never by position.
        for item in client.messages.batches.results(entry["id"]):
            by_id[item.custom_id] = result_row(item)
    if pending:
        raise BatchStateError(f"batches not ended yet: {', '.join(pending)}; run `wait` first")

    rows = [by_id.get(cid) or {"custom_id": cid, "status": "missing"} for cid in order]
    _atomic_write(out / "results.jsonl", _jsonl(rows))

    counts: dict[str, int] = {}
    tokens = {"input_tokens": 0, "output_tokens": 0, "cache_creation_input_tokens": 0, "cache_read_input_tokens": 0}
    truncated = 0
    for row in rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
        truncated += bool(row.get("truncated"))
        for k, v in row.get("usage", {}).items():
            tokens[k] += v
    # invalid_request_error rows will fail again unchanged; everything else may be retried.
    retry_ids = {
        r["custom_id"] for r in rows
        if r["status"] in ("expired", "canceled", "missing")
        or (r["status"] == "errored" and r["error"]["type"] != "invalid_request_error")
    }
    retry_rows = [r for r in inputs if r["custom_id"] in retry_ids]
    _atomic_write(out / "retryable.jsonl", _jsonl(
        [{"custom_id": r["custom_id"], **_retry_row(r["params"])} for r in retry_rows]))
    cost = (tokens["input_tokens"] * BATCH_PRICE_PER_MTOK["input"]
            + tokens["output_tokens"] * BATCH_PRICE_PER_MTOK["output"]) / 1e6
    summary = {
        "model": manifest["model"],
        "collected_at": _now(),
        "batches": [e["id"] for e in manifest["batches"]],
        "n_requests": len(rows),
        "status_counts": counts,
        "truncated_max_tokens": truncated,
        "retryable": len(retry_ids),
        "tokens": tokens,
        "estimated_cost_usd_input_output_only": round(cost, 4),
        "cost_note": f"batch list price ${BATCH_PRICE_PER_MTOK['input']}/${BATCH_PRICE_PER_MTOK['output']} "
                     "per MTok; cache tokens not priced; confirm against the Console usage page",
    }
    _atomic_write(out / "summary.json", json.dumps(summary, indent=2) + "\n")
    log(f"collected {len(rows)} results: {counts}; estimated ${cost:.2f}")
    return summary


def _retry_row(params: dict[str, Any]) -> dict[str, Any]:
    """Inverse of ``normalize_row`` for resubmission (the model comes from the next run's flags)."""
    row: dict[str, Any] = {"messages": params["messages"], "max_tokens": params["max_tokens"]}
    if "system" in params:
        row["system"] = params["system"]
    if "output_config" in params and "effort" in params["output_config"]:
        row["effort"] = params["output_config"]["effort"]
    rest = {k: v for k, v in params.items()
            if k not in ("model", "messages", "max_tokens", "system", "output_config")}
    if rest:
        row["params"] = rest
    return row


def cancel(client: Any, out: Path, *, log=print) -> None:
    for entry in load_manifest(out)["batches"]:
        batch = client.messages.batches.cancel(entry["id"])
        log(f"{entry['id']}: {batch.processing_status}")


def describe_api_error(exc: BaseException) -> str:
    """Actionable text for the failures this integration is likely to meet."""
    name = type(exc).__name__
    code = getattr(exc, "status_code", None)
    text = str(exc)
    hints = {
        "AuthenticationError": "ANTHROPIC_API_KEY is missing, malformed or revoked.",
        "PermissionDeniedError": "this key may lack access to claude-mythos-5-1 (Project Glasswing "
                                 "participants only) or to the Batches API.",
        "NotFoundError": "model or batch not found; claude-mythos-5-1 is hidden from keys outside "
                         "Project Glasswing, and batch results expire after 29 days.",
        "RateLimitError": "rate limited; the Mythos models share their own rate-limit pool. "
                          "Retry later or submit fewer requests.",
        "RequestTooLargeError": "batch exceeds 256 MB; split the input file.",
    }
    hint = hints.get(name, "")
    if name == "BadRequestError" and "data retention" in text.lower():
        hint = ("Mythos requires 30-day data retention on the organization/workspace "
                "(not available under zero data retention).")
    elif name == "BadRequestError":
        hint = "the API rejected the request body; fix the input rows."
    return f"{name}" + (f" (HTTP {code})" if code else "") + f": {text}" + (f"\n  hint: {hint}" if hint else "")
