#!/usr/bin/env bash
# Opt-in Linux acceptance using approved, preinstalled native binaries only.
# Needs Python 3, /proc, and a real UI-enabled Cairn on PATH. No installs/builds.
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
python3 -I - "$repo" <<'PY'
import json
import hashlib
from html.parser import HTMLParser
import os
from pathlib import Path
import signal
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from urllib.parse import urljoin, urlsplit

repo = Path(sys.argv[1])
assert Path('/proc/self/environ').exists(), 'Linux /proc required for credential assertion'
wrapper = ['bash', str(repo / 'scripts/private-workspace.sh')]
# Only this loopback acceptance client bypasses system proxies.
http = urllib.request.build_opener(urllib.request.ProxyHandler({}))

class Scripts(HTMLParser):
    def __init__(self):
        super().__init__()
        self.sources = []

    def handle_starttag(self, tag, attrs):
        if tag == 'script' and dict(attrs).get('src'):
            self.sources.append(dict(attrs)['src'])

def port():
    with socket.socket() as sock:
        sock.bind(('127.0.0.1', 0))
        return sock.getsockname()[1]

with tempfile.TemporaryDirectory(prefix='private-workspace-native-') as temp:
    root = Path(temp)
    state = root / 'state'
    # Isolate user configuration; public settings deliberately tempt the wrapper.
    env = {k: v for k, v in os.environ.items() if not k.startswith('CAIRN_')}
    env.update(HOME=str(root), CAIRN_PRIVATE_STATE=str(state),
               CAIRN_KEY=str(root / 'absent-key'), WORKSPACE_OFFLINE='1',
               CONDUCTOR_ENDPOINT='http://127.0.0.1:8080',
               CONDUCTOR_TOKEN='acceptance-sentinel-not-a-real-token',
               CAIRN_SEEDS='builtin', CAIRN_PEERS='open', CAIRN_PORTMAP='on',
               CAIRN_PORTMAP_HTTP='on', CAIRN_BEACON_PORT='47396')
    subprocess.run(wrapper + ['check'], env=env, cwd=root, check=True, timeout=15)
    binary = Path(shutil.which('cairn', path=env.get('PATH')))
    print('tested Cairn binary SHA-256: ' + hashlib.sha256(binary.read_bytes()).hexdigest(), flush=True)
    assert not state.exists(), 'check created node state'
    original_identity = None
    for iteration in range(2):
        p2p, web = port(), port()
        while web == p2p:
            web = port()
        env.update(CAIRN_NODE_P2P=f'127.0.0.1:{p2p}', CAIRN_NODE_HTTP=f'127.0.0.1:{web}')
        log = root / f'launch-{iteration}.log'
        with log.open('wb') as output:
            process = subprocess.Popen(wrapper + ['run'], env=env, cwd=root,
                                       stdout=output, stderr=output, start_new_session=True)
            try:
                def get(path):
                    with http.open(f'http://127.0.0.1:{web}{path}', timeout=2) as response:
                        assert response.status == 200
                        return response.read()

                deadline = time.monotonic() + 30
                while True:
                    assert process.poll() is None, 'native process exited before readiness'
                    try:
                        assert get('/health') == b'ok\n'
                        break
                    except (OSError, urllib.error.URLError):
                        assert time.monotonic() < deadline, 'native readiness deadline exceeded'
                        time.sleep(.1)
                network = json.loads(get('/network'))['node']
                assert network['peers_policy']['policy'] == 'allowlist'
                assert network['peers_policy']['allowed'] == 0
                assert network['reach']['external']['mode'] == 'off'
                assert 'tailnet' not in network['reach']
                assert json.loads(get('/ui/seeds.json'))['seeds'] == []
                # Feature flags alone do not prove the release contains UI assets.
                html = get('/ui/')
                assert b'<html' in html.lower(), 'packaged reader HTML missing'
                scripts = Scripts()
                scripts.feed(html.decode())
                assert scripts.sources, 'packaged reader has no JavaScript bundles'
                for source in scripts.sources:
                    resolved = urlsplit(urljoin(f'http://127.0.0.1:{web}/ui/', source))
                    assert resolved.scheme == 'http' and resolved.netloc == f'127.0.0.1:{web}', \
                        'reader references a script outside its loopback HTTP origin'
                    bundle = get(resolved.path + ('?' + resolved.query if resolved.query else ''))
                    assert bundle and not bundle.lstrip().lower().startswith((b'<!doctype', b'<html'))
                child_env = Path(f'/proc/{process.pid}/environ').read_bytes().split(b'\0')
                assert not any(v.startswith(b'CONDUCTOR_TOKEN=') for v in child_env)
                for setting in (b'CAIRN_NETWORK=private', b'CAIRN_ROLES=relay',
                                b'CAIRN_MCP_MAX_SPEND=0'):
                    assert setting in child_env, f'missing private child setting: {setting!r}'
                assert state.stat().st_mode & 0o777 == 0o700
                # Compare only the public identity; never print private key material.
                identity = json.loads((state / 'node.identity.json').read_text())['public']
                assert identity
                if original_identity is None:
                    original_identity = identity
                else:
                    assert identity == original_identity, 'restart changed transport identity'
                print(f'native launch {iteration + 1}: health, packaged reader, private discovery, '
                      'credential isolation and identity persistence passed', flush=True)
            except BaseException:
                # Native diagnostics contain no supplied bearer token; redact defensively.
                print(log.read_text(errors='replace').replace(env['CONDUCTOR_TOKEN'], '[redacted]'),
                      file=sys.stderr)
                raise
            finally:
                if process.poll() is None:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait(timeout=5)
    # A real bind failure must reach the supervisor, without replacing identity
    # or switching to another port/public launch path behind the operator's back.
    with socket.socket() as occupied:
        occupied.bind(('127.0.0.1', 0))
        occupied.listen()
        env['CAIRN_NODE_HTTP'] = f'127.0.0.1:{occupied.getsockname()[1]}'
        env['CAIRN_NODE_P2P'] = f'127.0.0.1:{port()}'
        failed = subprocess.run(wrapper + ['run'], env=env, cwd=root,
                                capture_output=True, timeout=15)
        assert failed.returncode != 0, 'occupied HTTP port was silently bypassed'
        diagnostics = (failed.stdout + failed.stderr).lower()
        assert b'bind' in diagnostics or b'address already in use' in diagnostics
        assert json.loads((state / 'node.identity.json').read_text())['public'] == original_identity
    print('native occupied-port failure propagated; transport identity preserved', flush=True)
    assert env['CONDUCTOR_TOKEN'] == 'acceptance-sentinel-not-a-real-token'
print('native workspace acceptance passed (local Linux process; no control-plane authentication claim)')
PY
