"""S3: the corpus, and the only place the campaign's bytes exist.

Two rules, both inherited from the strict storage protocol in `protocol.py`
rather than invented here:

  * **Objects are immutable.**  A key is written once.  The record objects are
    content-addressed (`<stream>-<offset>-<sha256>.bin`), so a retry after a
    connection reset rewrites identical bytes to an identical key, and an
    object that differs from one already at that key is a bug loud enough to
    stop for.  `putBytes` refuses the overwrite by default rather than letting
    it through.
  * **The bytes land before the ledger row.**  An object with no row is
    invisible to the dashboard and will be picked up by the next ingest pass.
    A row with no object is a record of points that do not exist, which the
    merge would then go looking for.  One of those is a delay and the other is
    a corruption, so the order is never the other way round.

boto3 is imported lazily.  `LocalObjectStore` is the same interface over a
directory, which is what the rehearsal scripts and the suite run against.
"""

from __future__ import annotations

import hashlib
import io
import os
import threading

from .config import Config


class ObjectExists(RuntimeError):
    """A key already holds bytes that differ from the ones offered."""


def sha256Bytes(data):
    return hashlib.sha256(data).hexdigest()


class ObjectStore:
    """S3 under one bucket and prefix."""

    def __init__(self, bucket, prefix="", client=None, endpoint="", region=""):
        self.bucket = bucket
        self.prefix = prefix.strip("/")
        self.endpoint = endpoint or ""
        self.region = region or ""
        self._client = client
        self._lock = threading.Lock()

    @classmethod
    def fromEnv(cls, config=None):
        config = config or Config.fromEnv()
        if not config.bucket:
            raise RuntimeError("ECC_BUCKET is not set")
        return cls(config.bucket, config.prefix, endpoint=config.s3Endpoint,
                   region=config.s3Region)

    @property
    def client(self):
        with self._lock:
            if self._client is None:
                import boto3  # lazy

                kw = {}
                if self.region:
                    kw["region_name"] = self.region
                if self.endpoint:
                    # A non-AWS endpoint is addressed by path: a bucket name
                    # in the hostname needs DNS that a MinIO on a private
                    # network does not have.  Real S3 keeps boto3's default.
                    from botocore.config import Config as BotoConfig

                    kw["endpoint_url"] = self.endpoint
                    kw["config"] = BotoConfig(s3={"addressing_style": "path"})
                self._client = boto3.client("s3", **kw)
            return self._client

    def ensureBucket(self):
        """Create the bucket if it is missing.  For an S3-compatible endpoint.

        Never called against real S3 by anything in this package: a deployed
        bucket is made by Terraform with versioning and a lifecycle policy,
        and a control plane that can conjure its own corpus bucket is a
        control plane that will silently write a campaign into the wrong one.
        """
        if not self.endpoint:
            raise RuntimeError("refusing to create a bucket outside a configured endpoint")
        try:
            self.client.head_bucket(Bucket=self.bucket)
            return False
        except Exception as exc:
            if not _isNotFound(exc):
                raise
        self.client.create_bucket(Bucket=self.bucket)
        return True

    def key(self, name):
        name = str(name).lstrip("/")
        return "%s/%s" % (self.prefix, name) if self.prefix else name

    # -- reads --------------------------------------------------------------
    def head(self, name):
        """`{bytes, sha256, uploadedAt, etag}` or None."""
        try:
            r = self.client.head_object(Bucket=self.bucket, Key=self.key(name))
        except Exception as exc:
            if _isNotFound(exc):
                return None
            raise
        meta = r.get("Metadata") or {}
        last = r.get("LastModified")
        return {
            "bytes": int(r.get("ContentLength") or 0),
            "sha256": meta.get("sha256"),
            "etag": (r.get("ETag") or "").strip('"'),
            "uploadedAt": int(last.timestamp()) if last is not None else 0,
        }

    def getBytes(self, name):
        r = self.client.get_object(Bucket=self.bucket, Key=self.key(name))
        return r["Body"].read()

    def list(self, prefix="", limit=None):
        """Every object under a prefix, oldest first, paginated.

        Paginated by hand rather than with a paginator so this works against
        an injected fake, and so the 1,000-key page limit cannot silently
        truncate the corpus the way an unpaginated `list_objects_v2` would.
        """
        out, token = [], None
        full = self.key(prefix) if prefix else self.prefix
        while True:
            kw = {"Bucket": self.bucket}
            if full:
                kw["Prefix"] = full
            if token:
                kw["ContinuationToken"] = token
            page = self.client.list_objects_v2(**kw)
            for item in page.get("Contents") or []:
                last = item.get("LastModified")
                out.append({
                    "key": item["Key"],
                    "bytes": int(item.get("Size") or 0),
                    "uploadedAt": int(last.timestamp()) if last is not None else 0,
                })
                if limit and len(out) >= limit:
                    return out
            if not page.get("IsTruncated"):
                break
            token = page.get("NextContinuationToken")
            if not token:
                break
        out.sort(key=lambda o: (o["uploadedAt"], o["key"]))
        return out

    # -- writes -------------------------------------------------------------
    def putBytes(self, name, data, overwrite=False, metadata=None):
        """Write bytes once.  Returns `{key, bytes, sha256, existed}`.

        An identical object already at the key is success, not a conflict:
        that is what a retried upload looks like, and it is the reason the
        supervisor can retry an upload without first asking whether the last
        attempt got through.
        """
        digest = sha256Bytes(data)
        key = self.key(name)
        if not overwrite:
            existing = self.head(name)
            if existing is not None:
                if existing.get("sha256") and existing["sha256"] != digest:
                    raise ObjectExists("%s holds different bytes (%s != %s)"
                                       % (key, existing["sha256"], digest))
                if not existing.get("sha256") and existing.get("bytes") != len(data):
                    raise ObjectExists("%s holds %d bytes, offered %d"
                                       % (key, existing.get("bytes"), len(data)))
                return {"key": key, "bytes": len(data), "sha256": digest, "existed": True}
        meta = dict(metadata or {})
        meta["sha256"] = digest
        self.client.put_object(Bucket=self.bucket, Key=key, Body=io.BytesIO(data),
                               ContentLength=len(data), Metadata=meta)
        return {"key": key, "bytes": len(data), "sha256": digest, "existed": False}

    def putFile(self, name, path, overwrite=False, metadata=None):
        with open(path, "rb") as fh:
            return self.putBytes(name, fh.read(), overwrite=overwrite, metadata=metadata)


class LocalObjectStore:
    """The same interface over a directory, for rehearsals and tests."""

    def __init__(self, root, prefix=""):
        self.root = root
        self.prefix = prefix.strip("/")
        os.makedirs(root, exist_ok=True)

    def key(self, name):
        name = str(name).lstrip("/")
        return "%s/%s" % (self.prefix, name) if self.prefix else name

    def _path(self, name):
        return os.path.join(self.root, self.key(name))

    def head(self, name):
        path = self._path(name)
        if not os.path.exists(path):
            return None
        with open(path, "rb") as fh:
            data = fh.read()
        return {"bytes": len(data), "sha256": sha256Bytes(data),
                "etag": sha256Bytes(data)[:32], "uploadedAt": int(os.path.getmtime(path))}

    def getBytes(self, name):
        with open(self._path(name), "rb") as fh:
            return fh.read()

    def list(self, prefix="", limit=None):
        base = os.path.join(self.root, self.key(prefix) if prefix else self.prefix)
        out = []
        for dirpath, _, names in os.walk(base if os.path.isdir(base) else self.root):
            for name in names:
                path = os.path.join(dirpath, name)
                key = os.path.relpath(path, self.root)
                if prefix and not key.startswith(self.key(prefix)):
                    continue
                out.append({"key": key, "bytes": os.path.getsize(path),
                            "uploadedAt": int(os.path.getmtime(path))})
        out.sort(key=lambda o: (o["uploadedAt"], o["key"]))
        return out[:limit] if limit else out

    def putBytes(self, name, data, overwrite=False, metadata=None):
        digest = sha256Bytes(data)
        path = self._path(name)
        if not overwrite and os.path.exists(path):
            existing = self.head(name)
            if existing["sha256"] != digest:
                raise ObjectExists("%s holds different bytes" % self.key(name))
            return {"key": self.key(name), "bytes": len(data), "sha256": digest, "existed": True}
        os.makedirs(os.path.dirname(path), exist_ok=True)
        tmp = path + ".tmp"
        with open(tmp, "wb") as fh:
            fh.write(data)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
        return {"key": self.key(name), "bytes": len(data), "sha256": digest, "existed": False}

    def putFile(self, name, path, overwrite=False, metadata=None):
        with open(path, "rb") as fh:
            return self.putBytes(name, fh.read(), overwrite=overwrite, metadata=metadata)


def _isNotFound(exc):
    response = getattr(exc, "response", None) or {}
    code = str((response.get("Error") or {}).get("Code") or "")
    if code in ("404", "NoSuchKey", "NotFound"):
        return True
    return type(exc).__name__ in ("NoSuchKey", "ClientError404") or "404" in str(exc)


def storeFromEnv(config=None):
    """The object store a CLI should use: the local directory when one is set."""
    config = config or Config.fromEnv()
    if config.localStore:
        return LocalObjectStore(os.path.join(config.localStore, "s3"), config.prefix)
    return ObjectStore.fromEnv(config)


def copyTree(src, dst, prefix=""):
    """Copy every object under a prefix between two stores (rehearsal aid)."""
    moved = 0
    for item in src.list(prefix):
        dst.putBytes(item["key"], src.getBytes(item["key"]))
        moved += 1
    return moved


__all__ = ["ObjectStore", "LocalObjectStore", "ObjectExists", "sha256Bytes",
           "storeFromEnv", "copyTree"]
