"""Read access to docs/ic/measurement/registry.json."""
from __future__ import annotations

import json
import os
from typing import Any

from .canonical import sha256_file

HERE = os.path.dirname(os.path.abspath(__file__))
REGISTRY_PATH = os.path.abspath(os.path.join(HERE, "..", "..", "docs", "ic", "measurement", "registry.json"))


class Registry:
    def __init__(self, doc: dict[str, Any], path: str = REGISTRY_PATH):
        self.doc = doc
        self.path = path
        self.sha256 = sha256_file(path) if os.path.exists(path) else None
        self._families = {f["family"]: f for f in doc["factor_base_families"]["list"]}
        self._units = {u["id"]: u for u in doc["units"]["list"]}
        self._windows = {w["name"]: w for w in doc["windows"]["list"]}
        self._rho = {r["id"]: r for r in doc["references"]["rho"]}
        self._floor = {f["id"]: f for f in doc["references"]["floor"]}
        self._phases = [p["name"] for p in doc["phases"]["list"]]

    @classmethod
    def load(cls, path: str = REGISTRY_PATH) -> "Registry":
        with open(path, encoding="utf-8") as fh:
            return cls(json.load(fh), path)

    def family(self, name: str) -> dict[str, Any] | None:
        return self._families.get(name)

    def has_unit(self, uid: str) -> bool:
        return uid in self._units

    def unit(self, uid: str) -> dict[str, Any] | None:
        return self._units.get(uid)

    def window(self, name: str) -> dict[str, Any] | None:
        return self._windows.get(name)

    def has_rho(self, rid: str) -> bool:
        return rid in self._rho

    def has_floor(self, fid: str) -> bool:
        return fid in self._floor

    def solver_options(self, solver: str) -> list[str]:
        return list(self.doc["solver_options_that_change_search"]["list"].get(solver, []))

    @property
    def phases(self) -> list[str]:
        return list(self._phases)

    def phase(self, name: str) -> dict[str, Any] | None:
        for p in self.doc["phases"]["list"]:
            if p["name"] == name:
                return p
        return None
