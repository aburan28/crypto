import pytest

from taskq import protocol
from conftest import make_spec

SHA = "a" * 40


def test_defaults_are_explicit_and_hash_is_stable():
    a = protocol.normalize_spec(make_spec(SHA, ["ok.py"]))
    b = protocol.normalize_spec({**make_spec(SHA, ["ok.py"]),
                                 "limits": {"timeout_seconds": 3600}})
    assert a["retry"]["max_attempts"] == 3 and a["command"]["cwd"] == "."
    assert protocol.spec_sha256(a) == protocol.spec_sha256(b)
    assert "benchmark" not in a


@pytest.mark.parametrize("bad, where", [
    ({"source": {"repo": "crypto", "commit": "main"}}, "source/commit"),
    ({"source": {"repo": "https://evil", "commit": SHA}}, "source/repo"),
    ({"kind": "shell"}, "kind"),
    ({"surprise": 1}, "<root>"),
])
def test_rejects(bad, where):
    with pytest.raises(protocol.SpecError, match=where):
        protocol.normalize_spec({**make_spec(SHA, ["ok.py"]), **bad})


def test_cwd_may_not_escape():
    spec = make_spec(SHA, ["x"])
    spec["command"]["cwd"] = "../.."
    with pytest.raises(protocol.SpecError, match="escapes"):
        protocol.normalize_spec(spec)


def test_task_ids_sort_by_time():
    import time
    a = protocol.new_task_id()
    time.sleep(0.002)
    b = protocol.new_task_id()
    assert a < b and a.startswith("T-")


def test_sparse_paths_may_not_escape():
    spec = make_spec(SHA, ["x"])
    spec["source"]["sparse_paths"] = ["research/../../etc"]
    with pytest.raises(protocol.SpecError, match="escapes"):
        protocol.normalize_spec(spec)


def test_shipped_example_is_valid():
    import json
    import pathlib
    ex = pathlib.Path(__file__).parent.parent / "examples" / "benchmark.json"
    protocol.normalize_spec(json.loads(ex.read_text()))
