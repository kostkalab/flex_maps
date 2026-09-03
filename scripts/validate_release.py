#!/usr/bin/env python3
"""Validate graph, sidecar, DG0, and embedding release artifacts."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import networkx as nx
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.validation import (
    validate_embedding_table,
    validate_full_graph_sidecar,
    validate_graph_sidecar,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timestamp", required=True, help="Artifact timestamp to validate")
    parser.add_argument("--results-dir", default="results")
    parser.add_argument("--species-dir", default="species")
    return parser.parse_args()


def load_species(species_dir: Path) -> list[dict]:
    species = []
    for path in sorted(species_dir.glob("*.yaml")):
        data = yaml.safe_load(path.read_text())
        species.append(
            {
                "species": data["species_code"],
                "output_prefix": data.get(
                    "output_prefix", f"{data['species_code']}_metabolic_graph"
                ),
            }
        )
    if not species:
        raise RuntimeError(f"No species configs found in {species_dir}")
    return species


def require(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def graph_row(root: Path, spec: dict, timestamp: str) -> tuple[dict, pd.DataFrame]:
    code = spec["species"]
    prefix = spec["output_prefix"]
    result_dir = root / code
    graph_path = require(result_dir / f"{prefix}.{timestamp}.graphml")
    full_graph_path = require(result_dir / f"{prefix}.full.{timestamp}.graphml")
    sidecar_path = require(result_dir / f"{prefix}.reactions.{timestamp}.tsv")
    dropped_reactions = require(
        result_dir / f"{prefix}.smiles_dropped_reactions.{timestamp}.tsv"
    )
    dropped_compounds = require(
        result_dir / f"{prefix}.smiles_dropped_compounds.{timestamp}.tsv"
    )
    timestamped_dg0 = result_dir / f"reaction_dg0.{timestamp}.tsv"
    dg0_path = require(timestamped_dg0 if timestamped_dg0.exists() else result_dir / "reaction_dg0.tsv")

    graph = nx.read_graphml(graph_path)
    full_graph = nx.read_graphml(full_graph_path)
    sidecar = pd.read_csv(sidecar_path, sep="\t", dtype=str, keep_default_na=False)
    dg0 = pd.read_csv(dg0_path, sep="\t", dtype=str, keep_default_na=False)
    semantic_errors = validate_graph_sidecar(graph, sidecar)
    semantic_errors.extend(
        f"full graph: {error}"
        for error in validate_full_graph_sidecar(full_graph, sidecar)
    )
    graph_nodes = {
        str(node)
        for node, data in graph.nodes(data=True)
        if data.get("group") == "Reaction"
    }
    sidecar_nodes = set(sidecar["reaction_node"].astype(str))
    full_count = sum(
        data.get("group") == "Reaction" for _, data in full_graph.nodes(data=True)
    )
    statuses = dg0.get("delta_g_status", pd.Series(dtype=str)).value_counts().to_dict()
    return {
        "species": code,
        "timestamp": timestamp,
        "graph_reaction_nodes": len(graph_nodes),
        "full_graph_reaction_nodes": int(full_count),
        "reaction_sidecar_rows": len(sidecar),
        "reaction_sidecar_unique_nodes": sidecar["reaction_node"].nunique(),
        "sidecar_matches_graph_reactions": graph_nodes == sidecar_nodes,
        "smiles_dropped_reactions": len(pd.read_csv(dropped_reactions, sep="\t")),
        "smiles_dropped_compounds": len(pd.read_csv(dropped_compounds, sep="\t")),
        "dg0_rows": len(dg0),
        "dg0_calculated_rows": int((dg0.get("delta_g_status") == "equilibrator_calculated").sum()),
        "dg0_status_counts": json.dumps(statuses),
        "semantic_validation_errors": len(semantic_errors),
        "semantic_validation_details": "|".join(semantic_errors),
        "passes": graph_nodes == sidecar_nodes and not semantic_errors,
        "sidecar_path": str(sidecar_path),
        "graphml_path": str(graph_path),
        "full_graphml_path": str(full_graph_path),
    }, sidecar


def embedding_paths(result_dir: Path) -> dict[str, Path]:
    return {
        "side_mean_fingerprint_raw": result_dir / "reaction_embeddings/reaction_embeddings_side_mean_fingerprint_raw.parquet",
        "side_mean_fingerprint_pca256": result_dir / "reaction_embeddings/reaction_embeddings_side_mean_fingerprint_pca256.parquet",
        "drfp": result_dir / "reaction_embeddings/reaction_embeddings_drfp.parquet",
        "reactiont5_4vec": result_dir / "reaction_embeddings_dl/reaction_embeddings_reactiont5_4vec.parquet",
        "rxngraphormer_2vec": result_dir / "reaction_embeddings_dl/reaction_embeddings_rxngraphormer_2vec.parquet",
    }


def markdown_table(table: pd.DataFrame) -> str:
    """Render a small DataFrame without adding a tabulate dependency."""
    headers = [str(column) for column in table.columns]
    rows = []
    for values in table.itertuples(index=False, name=None):
        rows.append([str(value).replace("|", "\\|") for value in values])
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def write_summary(path: Path, timestamp: str, graphs: pd.DataFrame, embeddings: pd.DataFrame) -> None:
    graph_columns = [
        "species", "graph_reaction_nodes", "full_graph_reaction_nodes",
        "reaction_sidecar_rows", "sidecar_matches_graph_reactions",
        "semantic_validation_errors", "passes", "smiles_dropped_reactions",
        "smiles_dropped_compounds", "dg0_rows", "dg0_calculated_rows",
    ]
    embedding_columns = [
        "species", "representation", "rows", "unique_reaction_nodes",
        "embedding_dims", "finite_rows", "status_counts",
        "missing_reaction_nodes", "extra_reaction_nodes", "passes",
    ]
    text = [
        f"# Cross-Species Update Validation ({timestamp})",
        "", "## Graphs and Sidecars", "",
        markdown_table(graphs[graph_columns]),
        "", "## Reaction Embeddings", "",
        markdown_table(embeddings[embedding_columns]),
        "", "## Files", "",
        f"- species graph validation CSV: `results/validation/species_graph_validation.{timestamp}.csv`",
        f"- reaction embedding validation CSV: `results/validation/reaction_embedding_validation.{timestamp}.csv`",
        "",
    ]
    path.write_text("\n".join(text))


def main() -> int:
    args = parse_args()
    results_dir = Path(args.results_dir)
    graph_rows = []
    embedding_rows = []
    for spec in load_species(Path(args.species_dir)):
        row, sidecar = graph_row(results_dir, spec, args.timestamp)
        graph_rows.append(row)
        for representation, path in embedding_paths(results_dir / spec["species"]).items():
            require(path)
            embedding_row = validate_embedding_table(sidecar, path, representation)
            embedding_rows.append({"species": spec["species"], **embedding_row})

    graph_table = pd.DataFrame(graph_rows)
    embedding_table = pd.DataFrame(embedding_rows)
    output_dir = results_dir / "validation"
    output_dir.mkdir(parents=True, exist_ok=True)
    graph_path = output_dir / f"species_graph_validation.{args.timestamp}.csv"
    embedding_path = output_dir / f"reaction_embedding_validation.{args.timestamp}.csv"
    summary_path = output_dir / f"update_validation_summary.{args.timestamp}.md"
    graph_table.to_csv(graph_path, index=False)
    embedding_table.to_csv(embedding_path, index=False)
    write_summary(summary_path, args.timestamp, graph_table, embedding_table)

    failed_graphs = int((~graph_table["passes"]).sum())
    failed_embeddings = int((~embedding_table["passes"]).sum())
    print(f"Wrote {graph_path}")
    print(f"Wrote {embedding_path}")
    print(f"Wrote {summary_path}")
    if failed_graphs or failed_embeddings:
        print(f"Validation failed: {failed_graphs} graph rows, {failed_embeddings} embedding rows")
        return 1
    print("Validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
