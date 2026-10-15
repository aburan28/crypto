"""Command line: ``python3 -m tools.mythos_batch <command>``."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import runner
from .spec import (
    DEFAULT_MODEL,
    EFFORTS,
    Defaults,
    RequestError,
    chunk_requests,
    load_requests,
    single_prompt,
)

EXIT_OK, EXIT_USAGE, EXIT_API, EXIT_FAILURES, EXIT_INTERRUPTED = 0, 2, 3, 4, 130


def _add_input_args(p: argparse.ArgumentParser) -> None:
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--input", type=Path, help="JSONL: one request per line (see README)")
    src.add_argument("--prompt", help="single prompt (custom_id prompt-1), for smoke tests")
    p.add_argument("--model", help=f"default $MYTHOS_MODEL or {DEFAULT_MODEL}")
    p.add_argument("--max-tokens", type=int, help="default for rows without max_tokens (thinking counts)")
    p.add_argument("--effort", choices=EFFORTS, help="default output_config.effort for rows without one")
    p.add_argument("--system", help="default system prompt for rows without one")


def _add_wait_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--poll-seconds", type=float, default=60.0)
    p.add_argument("--timeout-hours", type=float, default=25.0)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="mythos_batch", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("validate", help="check the input locally; no API call, no key needed")
    _add_input_args(p)
    p = sub.add_parser("submit", help="create the batch and return")
    _add_input_args(p)
    p.add_argument("--out", type=Path, required=True)
    p = sub.add_parser("run", help="submit (or resume) + wait + collect")
    _add_input_args(p)
    p.add_argument("--out", type=Path, required=True)
    _add_wait_args(p)
    for name, text in (("status", "show batch progress"), ("wait", "block until batches end"),
                       ("collect", "download results"), ("cancel", "cancel the batches")):
        p = sub.add_parser(name, help=text)
        p.add_argument("--out", type=Path, required=True)
        if name == "wait":
            _add_wait_args(p)
    p = sub.add_parser("list", help="recent batches on this workspace")
    p.add_argument("--limit", type=int, default=20)
    return ap


def _load(args: argparse.Namespace) -> tuple[Defaults, list[dict]]:
    d = Defaults()
    for attr in ("model", "max_tokens", "effort", "system"):
        if getattr(args, attr, None) is not None:
            setattr(d, attr, getattr(args, attr))
    if args.prompt is not None:
        return d, single_prompt(args.prompt, d)
    with args.input.open() as fh:
        return d, load_requests(fh, d)


def _results_exit(summary: dict) -> int:
    bad = {k: v for k, v in summary["status_counts"].items() if k != "succeeded"}
    if bad:
        print(f"not all requests succeeded: {bad} (see results.jsonl; retryable.jsonl lists "
              "resubmittable rows)", file=sys.stderr)
        return EXIT_FAILURES
    if summary["truncated_max_tokens"]:
        print(f"warning: {summary['truncated_max_tokens']} responses hit max_tokens "
              "(thinking counts toward it); raise max_tokens and resubmit those", file=sys.stderr)
    return EXIT_OK


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.cmd == "validate":
            defaults, reqs = _load(args)
            chunks = chunk_requests(reqs)
            print(f"ok: {len(reqs)} requests, model {defaults.model}, {len(chunks)} batch(es)")
            return EXIT_OK

        client = runner.make_client()
        if args.cmd in ("submit", "run"):
            out: Path = args.out
            defaults, reqs = _load(args)
            existing = out / runner.MANIFEST
            if args.cmd == "run" and existing.exists():
                if runner.load_manifest(out)["input_sha256"] != runner.request_digest(reqs):
                    raise runner.BatchStateError(
                        f"{out} holds a submission for different requests; pick a new --out")
                print(f"resuming {out} (already submitted)")
            else:
                runner.submit(client, reqs, out, model=defaults.model)
            if args.cmd == "submit":
                print(f"submitted; continue with: wait --out {out} && collect --out {out}")
                return EXIT_OK
        if args.cmd in ("run", "wait"):
            runner.wait(client, args.out, poll_seconds=args.poll_seconds,
                        timeout_seconds=args.timeout_hours * 3600)
        if args.cmd in ("run", "collect"):
            return _results_exit(runner.collect(client, args.out))
        if args.cmd == "status":
            for row in runner.status(client, args.out):
                print(json.dumps(row))
        elif args.cmd == "cancel":
            runner.cancel(client, args.out)
        elif args.cmd == "list":
            for i, b in enumerate(client.messages.batches.list(limit=args.limit)):
                if i >= args.limit:
                    break
                print(b.id, b.processing_status, b.created_at)
        return EXIT_OK
    except RequestError as exc:
        print("invalid input:\n  " + "\n  ".join(exc.problems), file=sys.stderr)
        return EXIT_USAGE
    except (runner.BatchStateError, TimeoutError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE
    except TypeError as exc:
        # The SDK raises TypeError when no credential source resolves.
        if "authentication method" not in str(exc):
            raise
        print("error: no Anthropic credentials found. Export ANTHROPIC_API_KEY "
              "(or run `ant auth login`).", file=sys.stderr)
        return EXIT_USAGE
    except KeyboardInterrupt:
        print("\ninterrupted; submitted batches keep running server-side - resume with "
              "`wait`/`collect` (or re-run `run`) on the same --out", file=sys.stderr)
        return EXIT_INTERRUPTED
    except Exception as exc:  # noqa: BLE001 - anthropic.* errors are mapped, others re-raised
        if type(exc).__module__.startswith("anthropic"):
            print(runner.describe_api_error(exc), file=sys.stderr)
            return EXIT_API
        raise
