import os
import shutil
import socket
import subprocess
import time


import pytest
import redis

from taskq.execute import RepoCache
from taskq.store import Store


def _free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


@pytest.fixture(scope="session")
def redis_url(tmp_path_factory):
    if os.environ.get("TASKQ_TEST_REDIS_URL"):
        yield os.environ["TASKQ_TEST_REDIS_URL"]
        return
    exe = shutil.which("redis-server")
    if not exe:
        pytest.skip("redis-server not installed and TASKQ_TEST_REDIS_URL unset")
    port = _free_port()
    d = tmp_path_factory.mktemp("redis")
    proc = subprocess.Popen([exe, "--port", str(port), "--save", "", "--appendonly", "no",
                             "--dir", str(d)], stdout=subprocess.DEVNULL)
    url = f"redis://127.0.0.1:{port}/0"
    for _ in range(100):
        try:
            redis.Redis.from_url(url).ping()
            break
        except redis.ConnectionError:
            time.sleep(0.05)
    yield url
    proc.terminate()
    proc.wait()


@pytest.fixture
def store(redis_url, request):
    ns = f"t{abs(hash(request.node.nodeid)) % 10**8}"
    s = Store.from_url(redis_url, ns)
    for k in s.r.scan_iter(f"{ns}:*"):
        s.r.delete(k)
    return s


@pytest.fixture
def repo(tmp_path):
    """A local git repo standing in for `crypto`, with a few test programs."""
    root = tmp_path / "src-repo"
    root.mkdir()
    (root / "ok.py").write_text(
        "import json, os, sys\n"
        "out = os.environ['TASKQ_OUTPUT_DIR']\n"
        "json.dump({'answer': 42, 'rep': os.environ['TASKQ_REPETITION']},\n"
        "          open(os.path.join(out, 'metrics.json'), 'w'))\n"
        "open(os.path.join(out, 'table.csv'), 'w').write('a,b\\n1,2\\n')\n"
        "print('hello'); sum(range(200000))\n")
    (root / "fail.py").write_text("import sys; print('boom', file=sys.stderr); sys.exit(3)\n")
    (root / "sleep.py").write_text("import time, sys; time.sleep(float(sys.argv[1]))\n")
    # toy "solver": brute-forces k = 17 for P = (0, 10) of order 50 on
    # y^2 = x^3 + 2x + 3 over F_97, or lies
    (root / "solve.py").write_text(
        "import json, os, sys\n"
        "p, a, b = 97, 2, 3\n"
        "def add(P, Q):\n"
        "    if P is None: return Q\n"
        "    if Q is None: return P\n"
        "    if P[0] == Q[0] and (P[1] + Q[1]) % p == 0: return None\n"
        "    l = ((3*P[0]*P[0]+a) * pow(2*P[1], -1, p) if P == Q else (Q[1]-P[1]) * pow(Q[0]-P[0], -1, p)) % p\n"
        "    x = (l*l - P[0] - Q[0]) % p\n"
        "    return (x, (l*(P[0]-x) - P[1]) % p)\n"
        "P, Q = (0, 10), None\n"
        "for _ in range(17): Q = add(Q, P)\n"
        "mode = sys.argv[1]\n"
        "R, k = None, 0\n"
        "while R != Q: R, k = add(R, P), k + 1\n"
        "if mode == 'lie': k += 1\n"
        "cert = {'kind': 'discrete_log', 'curve': {'field': 'prime', 'p': p, 'a': a, 'b': b},\n"
        "        'statement': {'P': list(P), 'Q': list(Q), 'k': k}}\n"
        "if mode != 'silent':\n"
        "    json.dump(cert, open(os.path.join(os.environ['TASKQ_OUTPUT_DIR'], 'certificate.json'), 'w'))\n")
    (root / "check.py").write_text(
        "import json, os, sys\n"
        "c = json.load(open(os.environ['TASKQ_CERTIFICATE']))\n"
        "print(json.dumps({'k': c['statement']['k']}))\n"
        "sys.exit(0 if c['statement']['k'] == 17 else 1)\n")
    (root / "sub").mkdir()
    (root / "sub" / "where.py").write_text("import os; print(os.getcwd())\n")
    env = {**os.environ, "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t",
           "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@t"}
    for cmd in (["init", "-q"], ["add", "."], ["commit", "-qm", "init"]):
        subprocess.run(["git", *cmd], cwd=root, check=True, env=env)
    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, capture_output=True,
                         text=True, check=True).stdout.strip()
    return root, sha


@pytest.fixture
def repos(tmp_path, repo):
    return RepoCache(tmp_path / "cache", {"crypto": str(repo[0])})


def make_spec(sha, argv, **over):
    import sys
    spec = {"schema": "taskq.task-spec/v1", "queue": "cpu", "kind": "command",
            "source": {"repo": "crypto", "commit": sha},
            "command": {"argv": [sys.executable, *argv]}}
    spec.update(over)
    return spec
