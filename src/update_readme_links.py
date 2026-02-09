#!/usr/bin/env python3
"""Update README download links based on current maps/ artifacts."""

from __future__ import annotations

import re
import sys
from pathlib import Path


def parse_species_config(path: Path) -> dict[str, str]:
    name = ""
    code = ""
    prefix = ""
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line.startswith("name:"):
            name = line.split(":", 1)[1].strip()
        elif line.startswith("species_code:"):
            code = line.split(":", 1)[1].strip()
        elif line.startswith("output_prefix:"):
            prefix = line.split(":", 1)[1].strip()
    if not (name and code and prefix):
        raise ValueError(f"Missing required fields in {path}")
    return {"name": name, "code": code, "prefix": prefix}


def latest_match(maps_dir: Path, pattern: str) -> Path:
    matches = list(maps_dir.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No files matching {pattern} in {maps_dir}")
    return max(matches, key=lambda p: p.stat().st_mtime)


def build_links(root: Path) -> list[str]:
    species_dir = root / "species"
    maps_dir = root / "maps"
    lines: list[str] = []
    for cfg in sorted(species_dir.glob("*.yaml")):
        spec = parse_species_config(cfg)
        prefix = spec["prefix"]
        graph = latest_match(maps_dir, f"{prefix}*.graphml.gz")
        pdf = latest_match(maps_dir, f"{prefix}*.pdf")
        lines.append(
            f"- {spec['name']} ([map](./maps/{graph.name}?raw=1), [stats](./maps/{pdf.name}))"
        )
    return lines


def replace_section(text: str, lines: list[str]) -> str:
    new_block = "\n".join(lines) + "\n"
    pattern = re.compile(r"(## Download links for maps\n)(?:\n)?(.*?)(?=\n## |\Z)", re.S)
    match = pattern.search(text)
    if not match:
        raise ValueError("README missing 'Download links for maps' section")
    start, end = match.span(2)
    return text[:start] + new_block + text[end:]


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    readme = root / "README.md"
    text = readme.read_text()
    lines = build_links(root)
    updated = replace_section(text, lines)
    if updated != text:
        readme.write_text(updated)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
