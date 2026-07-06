#!/usr/bin/env python3
"""Package flex_maps release artifacts from completed pipeline outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import tarfile
from dataclasses import dataclass
from pathlib import Path


REACTION_EMBEDDINGS = [
    ("reaction_embeddings", "reaction_embeddings_side_mean_fingerprint_raw.parquet"),
    ("reaction_embeddings", "reaction_embeddings_side_mean_fingerprint_pca256.parquet"),
    ("reaction_embeddings", "reaction_embeddings_drfp.parquet"),
    ("reaction_embeddings_dl", "reaction_embeddings_reactiont5_4vec.parquet"),
    ("reaction_embeddings_dl", "reaction_embeddings_rxngraphormer_2vec.parquet"),
]


@dataclass(frozen=True)
class Species:
    code: str
    name: str
    output_prefix: str


@dataclass(frozen=True)
class ReleaseFile:
    archive: str
    source: Path
    path_in_archive: str


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--timestamp",
        help="Pipeline timestamp to package; defaults to latest complete timestamp",
    )
    parser.add_argument(
        "--output-dir",
        default="release",
        help="Directory for release archives",
    )
    args = parser.parse_args()

    root = Path.cwd()
    species = load_species(root / "species")
    timestamp = args.timestamp or detect_latest_complete_timestamp(root, species)
    output_dir = root / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    archives = {
        "maps": f"flex_maps_maps.{timestamp}.tar.gz",
        "reaction_tables": f"flex_maps_reaction_tables.{timestamp}.tar.gz",
        "reaction_embeddings": f"flex_maps_reaction_embeddings.{timestamp}.tar.gz",
        "validation": f"flex_maps_validation.{timestamp}.tar.gz",
    }

    files = collect_release_files(root, species, timestamp, archives)
    manifest_path = output_dir / f"MANIFEST.{timestamp}.tsv"
    write_manifest(files, manifest_path)
    files.append(
        ReleaseFile(
            archive=archives["validation"],
            source=manifest_path,
            path_in_archive=f"validation/{manifest_path.name}",
        )
    )

    for archive in archives.values():
        path = output_dir / archive
        if path.exists():
            path.unlink()
        write_archive(path, [entry for entry in files if entry.archive == archive])

    latest_path = output_dir / "LATEST_RELEASE"
    latest_path.write_text(timestamp + "\n")

    print(f"Packaged release timestamp: {timestamp}")
    for archive in archives.values():
        path = output_dir / archive
        print(f"{path}\t{path.stat().st_size} bytes")
    print(f"{manifest_path}\t{manifest_path.stat().st_size} bytes")
    print(f"{latest_path}\t{latest_path.stat().st_size} bytes")
    return 0


def load_species(species_dir: Path) -> list[Species]:
    species = []
    for path in sorted(species_dir.glob("*.yaml")):
        fields = parse_simple_yaml_fields(path, {"species_code", "name", "output_prefix"})
        species.append(
            Species(
                code=fields["species_code"],
                name=fields["name"],
                output_prefix=fields["output_prefix"],
            )
        )
    if not species:
        raise RuntimeError(f"No species configs found in {species_dir}")
    return species


def parse_simple_yaml_fields(path: Path, names: set[str]) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw in path.read_text().splitlines():
        if ":" not in raw:
            continue
        key, value = raw.split(":", 1)
        key = key.strip()
        if key in names:
            values[key] = value.strip().strip("'\"")
    missing = sorted(names - set(values))
    if missing:
        raise ValueError(f"Missing fields in {path}: {', '.join(missing)}")
    return values


def detect_latest_complete_timestamp(root: Path, species: list[Species]) -> str:
    candidates = set()
    pattern = re.compile(r"\.reactions\.(\d{8}_\d{6})\.tsv$")
    for path in (root / "results").glob("*/*.reactions.*.tsv"):
        match = pattern.search(path.name)
        if match:
            candidates.add(match.group(1))

    for timestamp in sorted(candidates, reverse=True):
        try:
            collect_release_files(
                root,
                species,
                timestamp,
                {
                    "maps": "maps.tar.gz",
                    "reaction_tables": "reaction_tables.tar.gz",
                    "reaction_embeddings": "reaction_embeddings.tar.gz",
                    "validation": "validation.tar.gz",
                },
            )
        except FileNotFoundError:
            continue
        return timestamp

    raise RuntimeError("Could not find a complete release timestamp")


def collect_release_files(
    root: Path,
    species: list[Species],
    timestamp: str,
    archives: dict[str, str],
) -> list[ReleaseFile]:
    files: list[ReleaseFile] = []

    for spec in species:
        result_dir = root / "results" / spec.code
        map_graph = root / "maps" / f"{spec.output_prefix}.{timestamp}.graphml.gz"
        map_pdf = root / "maps" / f"{spec.output_prefix}.{timestamp}.pdf"
        add(files, archives["maps"], map_graph, f"maps/{map_graph.name}")
        add(files, archives["maps"], map_pdf, f"maps/{map_pdf.name}")

        table_files = [
            result_dir / f"{spec.output_prefix}.reactions.{timestamp}.tsv",
            result_dir / f"{spec.output_prefix}.smiles_dropped_reactions.{timestamp}.tsv",
            result_dir / f"{spec.output_prefix}.smiles_dropped_compounds.{timestamp}.tsv",
            result_dir / "reaction_dg0.tsv",
        ]
        for path in table_files:
            add(
                files,
                archives["reaction_tables"],
                path,
                f"reaction_tables/{spec.code}/{path.name}",
            )

        for subdir, filename in REACTION_EMBEDDINGS:
            path = result_dir / subdir / filename
            add(
                files,
                archives["reaction_embeddings"],
                path,
                f"reaction_embeddings/{spec.code}/{subdir}/{filename}",
            )

    validation_dir = root / "results" / "validation"
    validation_files = [
        validation_dir / f"update_validation_summary.{timestamp}.md",
        validation_dir / f"species_graph_validation.{timestamp}.csv",
        validation_dir / f"reaction_embedding_validation.{timestamp}.csv",
    ]
    for path in validation_files:
        add(files, archives["validation"], path, f"validation/{path.name}")

    return files


def add(files: list[ReleaseFile], archive: str, source: Path, path_in_archive: str) -> None:
    if not source.exists():
        raise FileNotFoundError(source)
    files.append(ReleaseFile(archive=archive, source=source, path_in_archive=path_in_archive))


def write_manifest(files: list[ReleaseFile], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["archive", "path_in_archive", "size_bytes", "sha256"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for entry in files:
            writer.writerow(
                {
                    "archive": entry.archive,
                    "path_in_archive": entry.path_in_archive,
                    "size_bytes": entry.source.stat().st_size,
                    "sha256": sha256(entry.source),
                }
            )


def write_archive(path: Path, files: list[ReleaseFile]) -> None:
    with tarfile.open(path, "w:gz") as tar:
        for entry in files:
            tar.add(entry.source, arcname=entry.path_in_archive)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


if __name__ == "__main__":
    raise SystemExit(main())
