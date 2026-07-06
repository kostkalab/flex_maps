"""Resolve release-hosted upstream artifacts into repo-local caches."""

from __future__ import annotations

import hashlib
import json
import re
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class ResolvedArtifact:
    """A cached upstream release artifact."""

    name: str
    path: Path
    asset_name: str
    url: str
    size: int | None
    sha256: str
    role: str
    representation: str


class UpstreamSourceResolver:
    """Resolve configured GitHub release assets and cache them under this repo."""

    def __init__(self, config_path: str | Path, repo_root: str | Path = ".") -> None:
        self.config_path = Path(config_path)
        self.repo_root = Path(repo_root).resolve()
        with self.config_path.open() as fh:
            self.config = yaml.safe_load(fh)
        self._release_json: dict[str, Any] | None = None

    def resolve_many(
        self,
        names: list[str],
        lock_path: str | Path | None = None,
        offline: bool = False,
        force_download: bool = False,
    ) -> dict[str, ResolvedArtifact]:
        resolved = {
            name: self.resolve(
                name,
                offline=offline,
                force_download=force_download,
            )
            for name in names
        }
        if lock_path is not None:
            self.write_lock(resolved, lock_path)
        return resolved

    def resolve(
        self,
        name: str,
        offline: bool = False,
        force_download: bool = False,
    ) -> ResolvedArtifact:
        artifacts = self.config.get("artifacts", {})
        if name not in artifacts:
            known = ", ".join(sorted(artifacts))
            raise KeyError(f"Unknown artifact {name!r}. Known artifacts: {known}")

        artifact_cfg = artifacts[name]
        cache_root_key = artifact_cfg["cache_root"]
        cache_root = self.repo_root / self.config["cache_roots"][cache_root_key]
        cache_root.mkdir(parents=True, exist_ok=True)

        if offline:
            asset = self._resolve_cached_asset(name, artifact_cfg, cache_root)
        else:
            asset = self._resolve_release_asset(name, artifact_cfg)

        path = cache_root / asset["name"]
        if not offline and (force_download or not path.exists()):
            self._download(asset["browser_download_url"], path)
        if not path.exists():
            raise FileNotFoundError(
                f"Resolved {name} to {path}, but the file is not cached"
            )

        return ResolvedArtifact(
            name=name,
            path=path,
            asset_name=asset["name"],
            url=asset.get("browser_download_url", ""),
            size=asset.get("size"),
            sha256=_sha256(path),
            role=str(artifact_cfg.get("role", "")),
            representation=str(artifact_cfg.get("representation", "")),
        )

    def write_lock(
        self,
        resolved: dict[str, ResolvedArtifact],
        lock_path: str | Path,
    ) -> None:
        lock = {
            "upstream": self.config["upstream"],
            "artifacts": {
                name: {
                    "asset_name": artifact.asset_name,
                    "path": str(artifact.path.relative_to(self.repo_root)),
                    "url": artifact.url,
                    "size": artifact.size,
                    "sha256": artifact.sha256,
                    "role": artifact.role,
                    "representation": artifact.representation,
                }
                for name, artifact in sorted(resolved.items())
            },
        }
        path = self.repo_root / lock_path
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w") as fh:
            yaml.safe_dump(lock, fh, sort_keys=False)

    def _resolve_release_asset(self, name: str, artifact_cfg: dict) -> dict:
        pattern = re.compile(artifact_cfg["asset_pattern"])
        matches = [
            asset
            for asset in self._release()["assets"]
            if pattern.fullmatch(asset["name"])
        ]
        if len(matches) != 1:
            matched_names = ", ".join(asset["name"] for asset in matches) or "none"
            raise RuntimeError(
                f"Expected exactly one release asset for {name}; matched {matched_names}"
            )
        return matches[0]

    def _resolve_cached_asset(
        self,
        name: str,
        artifact_cfg: dict,
        cache_root: Path,
    ) -> dict:
        pattern = re.compile(artifact_cfg["asset_pattern"])
        matches = [path for path in cache_root.iterdir() if pattern.fullmatch(path.name)]
        if len(matches) != 1:
            matched_names = ", ".join(path.name for path in matches) or "none"
            raise RuntimeError(
                f"Expected exactly one cached file for {name}; matched {matched_names}"
            )
        path = matches[0]
        return {
            "name": path.name,
            "browser_download_url": self._download_url(path.name),
            "size": path.stat().st_size,
        }

    def _release(self) -> dict:
        if self._release_json is None:
            upstream = self.config["upstream"]
            api_base = upstream.get("api_base", "https://api.github.com").rstrip("/")
            repo = upstream["repo"]
            release = upstream["release"]
            url = f"{api_base}/repos/{repo}/releases/tags/{release}"
            with urllib.request.urlopen(url) as response:
                self._release_json = json.load(response)
        return self._release_json

    def _download(self, url: str, path: Path) -> None:
        tmp_path = path.with_suffix(path.suffix + ".tmp")
        with urllib.request.urlopen(url) as response, tmp_path.open("wb") as out:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)
        tmp_path.replace(path)

    def _download_url(self, asset_name: str) -> str:
        upstream = self.config["upstream"]
        return (
            f"https://github.com/{upstream['repo']}/releases/download/"
            f"{upstream['release']}/{asset_name}"
        )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()
