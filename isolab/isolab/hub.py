"""`isolab hub`: run a NATS JetStream server with a configuration isolab writes.

A hub is one nats-server. Several hubs cluster over routes; a hub can also
hang off another lab as a leaf node. ``nats_download`` fetches the official
release binary for this platform and checks it against the published
SHA256SUMS, so a lab host needs nothing but this package to start.
"""
from __future__ import annotations

import hashlib
import io
import os
import platform
import secrets
import shutil
import string
import sys
import tarfile
import urllib.request
import zipfile
from importlib import resources
from pathlib import Path

NATS_VERSION = "2.15.0"
RELEASES = "https://github.com/nats-io/nats-server/releases/download"


def platform_asset(version: str = NATS_VERSION) -> tuple[str, str]:
    system = platform.system().lower()
    arch = {"x86_64": "amd64", "amd64": "amd64", "aarch64": "arm64", "arm64": "arm64"}.get(platform.machine(), platform.machine())
    ext = "zip" if system == "windows" else "tar.gz"
    name = f"nats-server-v{version}-{system}-{arch}"
    return name, f"{name}.{ext}"


def nats_download(dest_dir: Path, version: str = NATS_VERSION) -> Path:
    name, asset = platform_asset(version)
    dest_dir.mkdir(parents=True, exist_ok=True)
    target = dest_dir / "nats-server"
    url = f"{RELEASES}/v{version}/{asset}"
    with urllib.request.urlopen(url, timeout=120) as resp:
        data = resp.read()
    sums = urllib.request.urlopen(f"{RELEASES}/v{version}/SHA256SUMS", timeout=60).read().decode()
    want = next((line.split()[0] for line in sums.splitlines() if line.strip().endswith(asset)), None)
    got = hashlib.sha256(data).hexdigest()
    if want and got != want:
        raise RuntimeError(f"{asset}: sha256 {got} does not match published {want}")
    if asset.endswith(".zip"):
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            member = next(m for m in z.namelist() if m.endswith("nats-server") or m.endswith("nats-server.exe"))
            target.write_bytes(z.read(member))
    else:
        with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as t:
            member = next(m for m in t.getmembers() if m.name.endswith("/nats-server"))
            fh = t.extractfile(member)
            assert fh is not None
            target.write_bytes(fh.read())
    target.chmod(0o755)
    return target


def find_nats_server(extra: list[str] | None = None) -> str | None:
    for cand in [*(extra or []), "/opt/isolab/bin/nats-server", str(Path.home() / ".isolab/bin/nats-server")]:
        if cand and Path(cand).is_file():
            return cand
    return shutil.which("nats-server")


def generate_token() -> str:
    return secrets.token_urlsafe(32)


def render_config(*, server_name: str, listen: str, store_dir: str, token: str, http: str | None = None,
                  max_memory: str = "1GB", max_file: str = "50GB", cluster_name: str | None = None,
                  cluster_listen: str | None = None, routes: list[str] | None = None,
                  leaf_remote: str | None = None, leaf_listen: str | None = None) -> str:
    tmpl = string.Template(resources.files("isolab").joinpath("nats.conf.tmpl").read_text())
    cluster_block = ""
    if cluster_name:
        route_lines = "\n".join(f'    "{r}"' for r in routes or [])
        cluster_block = (f"\ncluster {{\n  name: \"{cluster_name}\"\n  listen: \"{cluster_listen or '0.0.0.0:6222'}\"\n"
                         f"  authorization {{ token: \"{token}\" }}\n  routes = [\n{route_lines}\n  ]\n}}\n")
    leaf_block = ""
    if leaf_remote or leaf_listen:
        leaf_block = "\nleafnodes {\n"
        if leaf_listen:
            leaf_block += f"  listen: \"{leaf_listen}\"\n"
        if leaf_remote:
            leaf_block += f"  remotes = [ {{ url: \"{leaf_remote}\" }} ]\n"
        leaf_block += "}\n"
    return tmpl.substitute(SERVER_NAME=server_name, LISTEN=listen, STORE_DIR=store_dir, TOKEN=token,
                           HTTP_LINE=f'http: "{http}"' if http else "", MAX_MEMORY=max_memory, MAX_FILE=max_file,
                           CLUSTER_BLOCK=cluster_block, LEAF_BLOCK=leaf_block)


def run_hub(args) -> int:
    store = Path(args.store).expanduser()
    store.mkdir(parents=True, exist_ok=True)
    token = args.token or os.environ.get("ISOLAB_TOKEN")
    if not token:
        token_file = store / "token"
        if token_file.exists():
            token = token_file.read_text().strip()
        else:
            token = generate_token()
            token_file.write_text(token + "\n")
            token_file.chmod(0o600)
            print(f"isolab hub: generated token, saved to {token_file}", file=sys.stderr)
    conf = render_config(server_name=args.name or platform.node(), listen=args.listen, store_dir=str(store / "jetstream"),
                         token=token, http=args.http, max_memory=args.max_memory, max_file=args.max_file,
                         cluster_name=args.cluster_name, cluster_listen=args.cluster_listen,
                         routes=[r for r in (args.routes or "").split(",") if r], leaf_remote=args.leaf_remote,
                         leaf_listen=args.leaf_listen)
    conf_path = store / "nats.conf"
    conf_path.write_text(conf)
    conf_path.chmod(0o600)
    exe = find_nats_server([args.nats_server] if args.nats_server else None)
    if not exe and args.download:
        exe = str(nats_download(store / "bin"))
    if not exe:
        print("isolab hub: nats-server not found; install it, pass --nats-server, or use --download", file=sys.stderr)
        return 2
    host = args.listen.split(":")[0] or "HOST"
    port = args.listen.rsplit(":", 1)[-1]
    print(f"isolab hub: {exe} -c {conf_path}", file=sys.stderr)
    print(f"isolab hub: clients connect with ISOLAB_URL=nats://{token}@{host if host not in ('0.0.0.0', '') else '<this host>'}:{port}", file=sys.stderr)
    if args.print_config:
        print(conf)
        return 0
    os.execv(exe, [exe, "-c", str(conf_path)])
    return 0
