#!/usr/bin/env python3
"""Start a broker-authenticated worker without AWS tools or credentials."""

import hashlib
import json
import os
import shutil
import sys
import urllib.error
import urllib.request


def request(path, body):
    endpoint = os.environ["ECC_UPLOAD_API"].rstrip("/")
    token = os.environ["ECC_UPLOAD_TOKEN"]
    workerId = os.environ["ECC_WORKER_ID"]
    data = json.dumps(dict(body, workerId=workerId), separators=(",", ":")).encode()
    call = urllib.request.Request(
        endpoint + path, data=data, method="POST",
        headers={"authorization": "Bearer " + token,
                 "content-type": "application/json", "accept": "application/json",
                 "user-agent": "ecc2k130-worker/1"})
    try:
        with urllib.request.urlopen(call, timeout=30) as response:
            return json.loads(response.read())
    except urllib.error.HTTPError as error:
        try:
            detail = json.loads(error.read()).get("error", "request refused")
        except Exception:
            detail = "request refused"
        raise RuntimeError("upload API refused bootstrap object: " + detail)


def download(key, path, required=True):
    result = request("/v1/objects/access", {"key": key, "download": True})
    if not result.get("exists"):
        if required:
            raise RuntimeError("campaign object is missing: " + key)
        return False
    url = result.get("downloadUrl")
    if not url:
        raise RuntimeError("upload API omitted a download URL")
    temporary = path + ".part"
    with urllib.request.urlopen(url, timeout=300) as source, open(temporary, "wb") as output:
        shutil.copyfileobj(source, output, 1024 * 1024)
    os.replace(temporary, path)
    return True


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    root = os.environ.get("ECC_ROOT", "/work/ecc2k130")
    os.makedirs(root, exist_ok=True)
    campaignPath = os.path.join(root, "campaign.bootstrap.json")
    download("campaign.json", campaignPath)
    with open(campaignPath) as source:
        campaign = json.load(source)
    key = campaign.get("binaryKey")
    expected = campaign.get("binarySha256")
    if not key or not expected:
        raise RuntimeError("campaign does not name a pinned GPU binary")
    binary = os.path.join(root, "ecc2k130")
    download(key, binary)
    if sha256(binary) != expected:
        raise RuntimeError("downloaded GPU binary does not match campaign.json")
    os.chmod(binary, 0o755)
    prefix = key.rsplit("/", 1)[0]
    library = os.path.join(root, "libgomp.so.1")
    if download(prefix + "/libgomp.so.1", library, required=False):
        os.environ["LD_LIBRARY_PATH"] = root + ":" + os.environ.get("LD_LIBRARY_PATH", "")
    os.environ["ECC_ROOT"] = root
    os.environ["ECC_CLIENT"] = binary
    os.environ.setdefault("ECC_GPU", "0")
    os.execv(sys.executable, [sys.executable, "-u", "/opt/ecc2k130/worker.py"])


if __name__ == "__main__":
    main()
