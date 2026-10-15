"""Tests for tools/mythos_batch. Run: python3 -m pytest tools/test_mythos_batch.py

The end-to-end test drives the real ``anthropic`` SDK against a local fake of the
Message Batches endpoints (via ANTHROPIC_BASE_URL) and is skipped without the SDK.
"""

import json
import sys
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
from types import SimpleNamespace as NS

import pytest

sys.path.insert(0, str(Path(__file__).parent))
import mythos_batch as mb  # noqa: E402
from mythos_batch import cli, runner, spec  # noqa: E402

D = spec.Defaults


def rows(*items):
    return [json.dumps(i) for i in items]


# ---- validation -------------------------------------------------------------


def test_prompt_row_becomes_mythos_request():
    (item,) = spec.load_requests(rows({"custom_id": "a", "prompt": "hi", "effort": "high"}), D())
    p = item["params"]
    assert p["model"] == "claude-mythos-5-1"
    assert p["messages"] == [{"role": "user", "content": "hi"}]
    assert p["output_config"] == {"effort": "high"}
    # thinking is always on for Mythos: never send a thinking block unless asked
    assert "thinking" not in p and "temperature" not in p


def test_env_model_override(monkeypatch):
    monkeypatch.setenv("MYTHOS_MODEL", "claude-mythos-5-2")
    assert D().model == "claude-mythos-5-2"


@pytest.mark.parametrize(
    "row, needle",
    [
        ({"custom_id": "bad id!", "prompt": "x"}, "custom_id"),
        ({"custom_id": "a"}, "exactly one"),
        ({"custom_id": "a", "prompt": "x", "messages": []}, "exactly one"),
        ({"custom_id": "a", "prompt": "x", "bogus": 1}, "unknown keys"),
        ({"custom_id": "a", "prompt": "x", "max_tokens": 10**6}, "max_tokens"),
        ({"custom_id": "a", "prompt": "x", "effort": "turbo"}, "effort"),
        ({"custom_id": "a", "prompt": "x", "params": {"temperature": 0}}, "temperature"),
        ({"custom_id": "a", "prompt": "x", "params": {"thinking": {"type": "disabled"}}}, "thinking"),
        ({"custom_id": "a", "prompt": "x", "params": {"thinking": {"type": "enabled", "budget_tokens": 5000}}}, "thinking"),
        ({"custom_id": "a", "prompt": "x", "params": {"tool_choice": {"type": "any"}}}, "tool_choice"),
        ({"custom_id": "a", "prompt": "x", "params": {"model": "other"}}, "params.model"),
        ({"custom_id": "a", "messages": [{"role": "user", "content": "q"}, {"role": "assistant", "content": "{"}]}, "prefill"),
    ],
)
def test_rejects_invalid_rows(row, needle):
    with pytest.raises(spec.RequestError) as e:
        spec.load_requests(rows(row), D())
    assert needle in str(e.value.problems)


def test_reports_every_problem_and_duplicates():
    lines = rows({"custom_id": "a", "prompt": "x"}, {"custom_id": "a", "prompt": "y"}) + ["{nope", ""]
    with pytest.raises(spec.RequestError) as e:
        spec.load_requests(lines, D())
    assert len(e.value.problems) == 2
    assert any("duplicate" in p for p in e.value.problems)
    assert any("invalid JSON" in p for p in e.value.problems)


def test_adaptive_thinking_and_auto_tool_choice_allowed():
    ok = {"custom_id": "a", "prompt": "x",
          "params": {"thinking": {"type": "adaptive", "display": "summarized"}, "tool_choice": {"type": "auto"}}}
    (item,) = spec.load_requests(rows(ok), D())
    assert item["params"]["thinking"]["display"] == "summarized"


def test_chunking_by_count_and_bytes(monkeypatch):
    reqs = [{"custom_id": f"r{i}", "params": {"messages": []}} for i in range(5)]
    monkeypatch.setattr(spec, "MAX_REQUESTS_PER_BATCH", 2)
    assert [len(c) for c in spec.chunk_requests(reqs)] == [2, 2, 1]
    monkeypatch.setattr(spec, "MAX_REQUESTS_PER_BATCH", 100)
    monkeypatch.setattr(spec, "MAX_BATCH_BYTES", 130)
    assert sum(len(c) for c in spec.chunk_requests(reqs)) == 5
    assert len(spec.chunk_requests(reqs)) > 1


# ---- runner with fakes ------------------------------------------------------


def message(text, stop="end_turn", **kw):
    return NS(
        type="succeeded",
        message=NS(
            stop_reason=stop, model="claude-mythos-5-1",
            content=[NS(type="thinking"), NS(type="text", text=text)],
            usage=NS(input_tokens=10, output_tokens=20, cache_creation_input_tokens=0, cache_read_input_tokens=0),
            stop_details=kw.get("stop_details"),
        ),
    )


class FakeBatches:
    def __init__(self, results, ended=True):
        self.created, self.results_by_id, self.ended = [], results, ended

    def create(self, requests):
        self.created.append(requests)
        return NS(id=f"msgbatch_{len(self.created)}")

    def retrieve(self, bid):
        c = NS(processing=0, succeeded=0, errored=0, canceled=0, expired=0)
        return NS(id=bid, processing_status="ended" if self.ended else "in_progress", request_counts=c)

    def results(self, bid):
        return iter(self.results_by_id)

    def cancel(self, bid):
        return NS(processing_status="canceling")


def fake_client(results, ended=True):
    return NS(messages=NS(batches=FakeBatches(results, ended)))


def reqs(*ids):
    return spec.load_requests(rows(*({"custom_id": i, "prompt": f"p {i}"} for i in ids)), D())


def test_submit_collect_orders_by_input_and_classifies(tmp_path):
    r = reqs("a", "b", "c", "d", "e")
    results = [  # deliberately out of order
        NS(custom_id="e", result=NS(type="expired")),
        NS(custom_id="c", result=NS(type="errored", error=NS(request_id="req_1", error=NS(type="overloaded_error", message="busy")))),
        NS(custom_id="b", result=message("declined", stop="refusal", stop_details=NS(category="cyber", explanation="x"))),
        NS(custom_id="a", result=message("answer")),
        # "d" never returns -> missing
    ]
    client = fake_client(results)
    m = runner.submit(client, r, tmp_path, model="claude-mythos-5-1", log=lambda *_: None)
    assert m["n_requests"] == 5 and m["batches"][0]["id"] == "msgbatch_1"
    s = runner.collect(client, tmp_path, log=lambda *_: None)
    out = [json.loads(l) for l in (tmp_path / "results.jsonl").read_text().splitlines()]
    assert [o["custom_id"] for o in out] == ["a", "b", "c", "d", "e"]
    assert [o["status"] for o in out] == ["succeeded", "refused", "errored", "missing", "expired"]
    assert out[0]["text"] == "answer" and out[1]["stop_details"]["category"] == "cyber"
    assert out[2]["error"] == {"type": "overloaded_error", "message": "busy", "request_id": "req_1"}
    assert s["tokens"]["input_tokens"] == 20 and s["tokens"]["output_tokens"] == 40
    assert s["estimated_cost_usd_input_output_only"] == round((20 * 5 + 40 * 25) / 1e6, 4)
    retry = [json.loads(l)["custom_id"] for l in (tmp_path / "retryable.jsonl").read_text().splitlines()]
    assert retry == ["c", "d", "e"]  # the refusal is not retryable
    # the retry file is itself valid input for another submission
    assert len(spec.load_requests((tmp_path / "retryable.jsonl").open(), D())) == 3


def test_invalid_request_errors_are_not_retryable(tmp_path):
    client = fake_client([NS(custom_id="a", result=NS(type="errored", error=NS(request_id=None, error=NS(type="invalid_request_error", message="bad"))))])
    runner.submit(client, reqs("a"), tmp_path, model="m", log=lambda *_: None)
    s = runner.collect(client, tmp_path, log=lambda *_: None)
    assert s["retryable"] == 0 and s["status_counts"] == {"errored": 1}


def test_truncation_flagged(tmp_path):
    client = fake_client([NS(custom_id="a", result=message("cut", stop="max_tokens"))])
    runner.submit(client, reqs("a"), tmp_path, model="m", log=lambda *_: None)
    assert runner.collect(client, tmp_path, log=lambda *_: None)["truncated_max_tokens"] == 1


def test_second_submit_into_same_dir_refused(tmp_path):
    client = fake_client([])
    runner.submit(client, reqs("a"), tmp_path, model="m", log=lambda *_: None)
    with pytest.raises(runner.BatchStateError):
        runner.submit(client, reqs("a"), tmp_path, model="m", log=lambda *_: None)
    assert len(client.messages.batches.created) == 1  # no second billable batch


def test_collect_before_end_raises(tmp_path):
    client = fake_client([], ended=False)
    runner.submit(client, reqs("a"), tmp_path, model="m", log=lambda *_: None)
    with pytest.raises(runner.BatchStateError):
        runner.collect(client, tmp_path, log=lambda *_: None)


def test_wait_polls_then_times_out(tmp_path):
    client = fake_client([], ended=False)
    runner.submit(client, reqs("a"), tmp_path, model="m", log=lambda *_: None)
    t = iter(range(0, 1000, 100))
    with pytest.raises(TimeoutError):
        runner.wait(client, tmp_path, poll_seconds=1, timeout_seconds=250, sleep=lambda s: None,
                    clock=lambda: next(t), log=lambda *_: None)


def test_multi_batch_manifest(tmp_path, monkeypatch):
    monkeypatch.setattr(spec, "MAX_REQUESTS_PER_BATCH", 2)
    client = fake_client([])
    m = runner.submit(client, reqs("a", "b", "c"), tmp_path, model="m", log=lambda *_: None)
    assert [b["n_requests"] for b in m["batches"]] == [2, 1]


def test_manifest_has_no_credentials(tmp_path, monkeypatch):
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-ant-secret-value")
    runner.submit(fake_client([]), reqs("a"), tmp_path, model="m", log=lambda *_: None)
    for f in tmp_path.iterdir():
        assert "sk-ant-secret-value" not in f.read_text()


def test_describe_api_error_hints():
    class BadRequestError(Exception):
        status_code = 400

    msg = runner.describe_api_error(BadRequestError("organization must have data retention enabled"))
    assert "30-day" in msg


def test_missing_key_is_a_clear_error(tmp_path, monkeypatch, capsys):
    pytest.importorskip("anthropic")
    for var in ("ANTHROPIC_API_KEY", "ANTHROPIC_AUTH_TOKEN", "ANTHROPIC_PROFILE"):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("HOME", str(tmp_path))  # no ~/.config/anthropic profile
    monkeypatch.setenv("ANTHROPIC_BASE_URL", "http://127.0.0.1:9")
    code = cli.main(["submit", "--prompt", "hi", "--out", str(tmp_path / "o")])
    assert code == cli.EXIT_USAGE
    assert "ANTHROPIC_API_KEY" in capsys.readouterr().err
    assert not (tmp_path / "o" / runner.MANIFEST).exists()


def test_cli_validate_needs_no_key(tmp_path, monkeypatch, capsys):
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    f = tmp_path / "in.jsonl"
    f.write_text(json.dumps({"custom_id": "a", "prompt": "x"}) + "\n")
    assert cli.main(["validate", "--input", str(f)]) == 0
    assert "claude-mythos-5-1" in capsys.readouterr().out
    f.write_text(json.dumps({"custom_id": "a", "prompt": "x", "params": {"temperature": 1}}) + "\n")
    assert cli.main(["validate", "--input", str(f)]) == cli.EXIT_USAGE


# ---- end to end: real SDK against a local fake Batches API ------------------


class FakeAPI(BaseHTTPRequestHandler):
    batches: dict = {}
    seen: list = []

    def log_message(self, *a):
        pass

    def _send(self, obj, ctype="application/json"):
        body = obj if isinstance(obj, bytes) else json.dumps(obj).encode()
        self.send_response(200)
        self.send_header("content-type", ctype)
        self.send_header("content-length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _batch(self, bid):
        reqs = self.batches[bid]
        return {"id": bid, "type": "message_batch", "processing_status": "ended",
                "request_counts": {"processing": 0, "succeeded": len(reqs), "errored": 0, "canceled": 0, "expired": 0},
                "created_at": "2026-10-15T00:00:00Z", "expires_at": "2026-10-16T00:00:00Z",
                "ended_at": "2026-10-15T00:05:00Z", "cancel_initiated_at": None, "archived_at": None,
                "results_url": f"http://{self.headers['host']}/v1/messages/batches/{bid}/results"}

    def do_POST(self):
        n = int(self.headers["content-length"])
        body = json.loads(self.rfile.read(n))
        FakeAPI.seen.append((self.path, dict(self.headers), body))
        bid = f"msgbatch_{len(self.batches) + 1}"
        self.batches[bid] = body["requests"]
        self._send(self._batch(bid))

    def do_GET(self):
        parts = self.path.split("?")[0].strip("/").split("/")
        bid = parts[3]
        if parts[-1] == "results":
            lines = []
            for r in reversed(self.batches[bid]):  # out of order on purpose
                lines.append(json.dumps({"custom_id": r["custom_id"], "result": {"type": "succeeded", "message": {
                    "id": "msg_1", "type": "message", "role": "assistant", "model": r["params"]["model"],
                    "content": [{"type": "thinking", "thinking": "", "signature": "s"},
                                {"type": "text", "text": "echo: " + r["params"]["messages"][-1]["content"]}],
                    "stop_reason": "end_turn", "stop_sequence": None,
                    "usage": {"input_tokens": 7, "output_tokens": 9}}}}))
            self._send(("\n".join(lines) + "\n").encode(), "application/binary")
        else:
            self._send(self._batch(bid))


def test_end_to_end_with_real_sdk(tmp_path, monkeypatch, capsys):
    pytest.importorskip("anthropic")
    FakeAPI.batches, FakeAPI.seen = {}, []
    srv = HTTPServer(("127.0.0.1", 0), FakeAPI)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-ant-test")
    monkeypatch.setenv("ANTHROPIC_BASE_URL", f"http://127.0.0.1:{srv.server_port}")
    inp, out = tmp_path / "in.jsonl", tmp_path / "run"
    inp.write_text("".join(json.dumps(r) + "\n" for r in (
        {"custom_id": "one", "prompt": "alpha", "effort": "high"},
        {"custom_id": "two", "prompt": "beta", "system": "be brief"})))
    try:
        code = cli.main(["run", "--input", str(inp), "--out", str(out), "--poll-seconds", "0"])
        # re-running resumes instead of creating a second batch
        code2 = cli.main(["run", "--input", str(inp), "--out", str(out), "--poll-seconds", "0"])
    finally:
        srv.shutdown()
    assert code == 0 and code2 == 0
    creates = [s for s in FakeAPI.seen if s[0].endswith("/v1/messages/batches")]
    assert len(creates) == 1
    path, headers, body = creates[0]
    assert {k.lower(): v for k, v in headers.items()}["x-api-key"] == "sk-ant-test"
    first = body["requests"][0]["params"]
    assert first["model"] == "claude-mythos-5-1" and first["output_config"] == {"effort": "high"}
    got = [json.loads(l) for l in (out / "results.jsonl").read_text().splitlines()]
    assert [(g["custom_id"], g["text"]) for g in got] == [("one", "echo: alpha"), ("two", "echo: beta")]
    summary = json.loads((out / "summary.json").read_text())
    assert summary["status_counts"] == {"succeeded": 2} and summary["tokens"]["output_tokens"] == 18
