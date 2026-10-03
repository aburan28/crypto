"""Content-addressed inputs: hashing, a local cache, ELF sniffing, local-file resolution."""
from __future__ import annotations

import base64
import hashlib
import os
import shutil
import struct
from pathlib import Path
from typing import Any, Callable


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class LocalBlobCache:
    """``<root>/<sha256>`` files; verified on store and on use."""

    def __init__(self, root: str | Path):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    def path_for(self, digest: str) -> Path:
        return self.root / digest

    def has(self, digest: str) -> bool:
        p = self.path_for(digest)
        return p.exists() and sha256_file(p) == digest

    def store_bytes(self, data: bytes) -> str:
        digest = sha256_bytes(data)
        tmp = self.path_for(digest + ".tmp")
        tmp.write_bytes(data)
        os.replace(tmp, self.path_for(digest))
        return digest

    def store_file(self, path: str | Path) -> str:
        digest = sha256_file(path)
        target = self.path_for(digest)
        if not target.exists():
            tmp = self.path_for(digest + ".tmp")
            shutil.copyfile(path, tmp)
            os.replace(tmp, target)
        return digest


_ELF_MACHINES = {0x3E: "x86_64", 0xB7: "aarch64", 0x03: "i386", 0x28: "arm",
                 0xF3: "riscv64", 0x15: "ppc64le", 0x16: "s390x"}


def elf_info(path: str | Path) -> dict[str, Any] | None:
    """Architecture and linkage of an ELF file, or None for anything else."""
    try:
        with open(path, "rb") as fh:
            head = fh.read(64)
            if len(head) < 52 or head[:4] != b"\x7fELF":
                return None
            is64 = head[4] == 2
            little = head[5] == 1
            endian = "<" if little else ">"
            e_machine = struct.unpack(endian + "H", head[18:20])[0]
            if is64:
                e_phoff, = struct.unpack(endian + "Q", head[32:40])
                e_phentsize, e_phnum = struct.unpack(endian + "HH", head[54:58])
            else:
                e_phoff, = struct.unpack(endian + "I", head[28:32])
                e_phentsize, e_phnum = struct.unpack(endian + "HH", head[42:46])
            interp = None
            fh.seek(e_phoff)
            table = fh.read(e_phentsize * min(e_phnum, 64))
            for i in range(min(e_phnum, 64)):
                ph = table[i * e_phentsize:(i + 1) * e_phentsize]
                p_type = struct.unpack(endian + "I", ph[:4])[0]
                if p_type == 3:  # PT_INTERP
                    if is64:
                        off, = struct.unpack(endian + "Q", ph[8:16])
                        size, = struct.unpack(endian + "Q", ph[32:40])
                    else:
                        off, = struct.unpack(endian + "I", ph[4:8])
                        size, = struct.unpack(endian + "I", ph[16:20])
                    fh.seek(off)
                    interp = fh.read(min(size, 256)).split(b"\0", 1)[0].decode(errors="replace")
                    break
    except OSError:
        return None
    return {"arch": _ELF_MACHINES.get(e_machine, f"em{e_machine}"),
            "bits": 64 if is64 else 32, "static": interp is None, "interpreter": interp}


def decode_content(inp: dict[str, Any]) -> bytes:
    if inp.get("encoding") == "base64":
        return base64.b64decode(inp["content"])
    return inp["content"].encode()


def resolve_local_inputs(spec: dict[str, Any],
                         upload: Callable[[str, Path], None],
                         has_blob: Callable[[str], bool] | None = None,
                         base_dir: str | Path | None = None) -> list[dict[str, Any]]:
    """Replace every ``local_file`` input by a ``sha256`` blob reference.

    ``upload(digest, path)`` is called for each file the hub does not have
    yet. Returns the list of uploads performed, for the submitter's record.
    """
    done = []
    for i, inp in enumerate(spec.get("inputs") or []):
        if "local_file" not in inp:
            continue
        src = Path(inp["local_file"]).expanduser()
        if base_dir and not src.is_absolute():
            src = Path(base_dir) / src
        if not src.is_file():
            raise FileNotFoundError(f"inputs/{i}: {src} is not a file")
        digest = sha256_file(src)
        if has_blob is None or not has_blob(digest):
            upload(digest, src)
            done.append({"path": inp["path"], "sha256": digest, "bytes": src.stat().st_size})
        mode = inp.get("mode") or ("0755" if os.access(src, os.X_OK) else "0644")
        del inp["local_file"]
        inp.update({"sha256": digest, "bytes": src.stat().st_size, "mode": mode})
    return done
