"""Validation helpers for oriented graphs, reaction sidecars, and embeddings."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd


EXPECTED_EMBEDDING_DIMS = {
    "side_mean_fingerprint_raw": 2048,
    "side_mean_fingerprint_pca256": 512,
    "drfp": 2048,
    "reactiont5_4vec": 3072,
    "rxngraphormer_2vec": 1536,
}


class ValidationError(ValueError):
    """Raised when graph/sidecar semantic validation fails."""


def _compound_id(node: str, data: dict) -> str:
    return str(data.get("name") or str(node).replace("kegg:", "", 1))


def _split_compounds(value: object) -> set[str]:
    if value is None or pd.isna(value):
        return set()
    return {part for part in str(value).split("|") if part}


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes"}


def _dg0_value(value: object) -> float | None:
    if value is None or pd.isna(value) or str(value).strip() == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def validate_graph_sidecar(
    graph: nx.DiGraph,
    sidecar: pd.DataFrame,
    *,
    raise_on_error: bool = False,
) -> list[str]:
    """Validate reaction coverage, edge direction, and orientation semantics."""
    errors: list[str] = []
    required = {
        "reaction_node",
        "mnx_left_compounds",
        "mnx_right_compounds",
        "input_compounds",
        "output_compounds",
        "input_smiles",
        "output_smiles",
        "I_to_O_string",
        "O_to_I_string",
        "dg0",
        "orientation",
        "orientation_source",
        "thermo_agrees_with_mnx",
        "orientation_flip_from_mnx",
    }
    missing_columns = sorted(required - set(sidecar.columns))
    if missing_columns:
        errors.append("missing sidecar columns: " + ", ".join(missing_columns))
        if raise_on_error:
            raise ValidationError("; ".join(errors))
        return errors

    reaction_nodes = {
        str(node)
        for node, data in graph.nodes(data=True)
        if data.get("group") == "Reaction"
    }
    sidecar_nodes = sidecar["reaction_node"].astype(str)
    sidecar_node_set = set(sidecar_nodes)
    duplicates = sorted(sidecar_nodes[sidecar_nodes.duplicated()].unique())
    if duplicates:
        errors.append(f"duplicate reaction_node rows: {duplicates[:10]}")
    missing_reactions = sorted(reaction_nodes - sidecar_node_set)
    extra_reactions = sorted(sidecar_node_set - reaction_nodes)
    if missing_reactions:
        errors.append(f"graph reactions missing from sidecar: {missing_reactions[:10]}")
    if extra_reactions:
        errors.append(f"sidecar reactions missing from graph: {extra_reactions[:10]}")

    compound_nodes = {
        _compound_id(node, data): str(node)
        for node, data in graph.nodes(data=True)
        if data.get("group") == "Compound"
    }

    for _, row in sidecar.iterrows():
        reaction_node = str(row["reaction_node"])
        if reaction_node not in reaction_nodes:
            continue
        prefix = f"{reaction_node}: "
        left = _split_compounds(row["mnx_left_compounds"])
        right = _split_compounds(row["mnx_right_compounds"])
        inputs = _split_compounds(row["input_compounds"])
        outputs = _split_compounds(row["output_compounds"])
        dg0 = _dg0_value(row["dg0"])
        if dg0 is None and str(row["dg0"]).strip():
            errors.append(prefix + f"invalid dg0 value {row['dg0']!r}")

        graph_inputs = {
            _compound_id(node, graph.nodes[node])
            for node in graph.predecessors(reaction_node)
            if graph.nodes[node].get("group") == "Compound"
        }
        graph_outputs = {
            _compound_id(node, graph.nodes[node])
            for node in graph.successors(reaction_node)
            if graph.nodes[node].get("group") == "Compound"
        }
        surviving_compounds = set(compound_nodes)
        if graph_inputs != inputs & surviving_compounds:
            errors.append(prefix + "graph input compounds do not match sidecar")
        if graph_outputs != outputs & surviving_compounds:
            errors.append(prefix + "graph output compounds do not match sidecar")

        if dg0 is None:
            expected_inputs, expected_outputs = left, right
            expected_orientation = "unknown"
            expected_source = "mnx_fallback"
            expected_flip = False
            expected_agreement = ""
        elif dg0 > 0:
            expected_inputs, expected_outputs = right, left
            expected_orientation = "reverse"
            expected_source = "dg0"
            expected_flip = True
            expected_agreement = False
        else:
            expected_inputs, expected_outputs = left, right
            expected_orientation = "forward"
            expected_source = "dg0"
            expected_flip = False
            expected_agreement = True

        if inputs != expected_inputs or outputs != expected_outputs:
            errors.append(prefix + "oriented compounds do not match DG0/MetaNetX sides")
        if str(row["orientation"]) != expected_orientation:
            errors.append(prefix + f"orientation should be {expected_orientation}")
        if str(row["orientation_source"]) != expected_source:
            errors.append(prefix + f"orientation_source should be {expected_source}")
        if _as_bool(row["orientation_flip_from_mnx"]) != expected_flip:
            errors.append(prefix + "orientation_flip_from_mnx is inconsistent")
        agreement = row["thermo_agrees_with_mnx"]
        if expected_agreement == "":
            if not (agreement is None or pd.isna(agreement) or str(agreement) == ""):
                errors.append(prefix + "thermo agreement must be blank without DG0")
        elif _as_bool(agreement) != expected_agreement:
            errors.append(prefix + "thermo_agrees_with_mnx is inconsistent")

        for compound in inputs:
            compound_node = compound_nodes.get(compound)
            if compound_node is None:
                continue
            if not graph.has_edge(compound_node, reaction_node):
                errors.append(prefix + f"missing input edge {compound_node} -> {reaction_node}")
            if compound not in outputs and graph.has_edge(
                reaction_node, compound_node
            ):
                errors.append(prefix + f"unexpected reverse input edge to {compound_node}")
        for compound in outputs:
            compound_node = compound_nodes.get(compound)
            if compound_node is None:
                continue
            if not graph.has_edge(reaction_node, compound_node):
                errors.append(prefix + f"missing output edge {reaction_node} -> {compound_node}")
            if compound not in inputs and graph.has_edge(
                compound_node, reaction_node
            ):
                errors.append(prefix + f"unexpected reverse output edge from {compound_node}")

        input_smiles = "" if pd.isna(row["input_smiles"]) else str(row["input_smiles"])
        output_smiles = "" if pd.isna(row["output_smiles"]) else str(row["output_smiles"])
        if str(row["I_to_O_string"]) != f"{input_smiles}>>{output_smiles}":
            errors.append(prefix + "I_to_O_string does not match side SMILES")
        if str(row["O_to_I_string"]) != f"{output_smiles}>>{input_smiles}":
            errors.append(prefix + "O_to_I_string does not match side SMILES")

    if raise_on_error and errors:
        preview = "; ".join(errors[:10])
        if len(errors) > 10:
            preview += f"; ... ({len(errors)} errors total)"
        raise ValidationError(preview)
    return errors


def validate_full_graph_sidecar(
    full_graph: nx.DiGraph,
    sidecar: pd.DataFrame,
) -> list[str]:
    """Validate surviving sidecar reactions against their full-graph orientation."""
    surviving = set(sidecar["reaction_node"].astype(str))
    scoped = full_graph.copy()
    scoped.remove_nodes_from(
        node
        for node, data in full_graph.nodes(data=True)
        if data.get("group") == "Reaction" and str(node) not in surviving
    )
    return validate_graph_sidecar(scoped, sidecar)


def validate_embedding_table(
    sidecar: pd.DataFrame,
    path: str | Path,
    representation: str,
) -> dict:
    """Validate one embedding parquet against an oriented reaction sidecar."""
    path = Path(path)
    expected_dim = EXPECTED_EMBEDDING_DIMS[representation]
    table = pd.read_parquet(path)
    required = {"reaction_node", "embedding_status", "embedding"}
    missing_columns = sorted(required - set(table.columns))
    sidecar_nodes = set(sidecar["reaction_node"].astype(str))
    if missing_columns:
        return {
            "representation": representation,
            "expected_dim": expected_dim,
            "rows": len(table),
            "unique_reaction_nodes": 0,
            "embedding_dims": "{}",
            "finite_rows": 0,
            "status_counts": "{}",
            "missing_reaction_nodes": len(sidecar_nodes),
            "extra_reaction_nodes": 0,
            "passes": False,
            "path": str(path),
            "validation_errors": "missing columns: " + ", ".join(missing_columns),
        }

    nodes = table["reaction_node"].astype(str)
    node_set = set(nodes)
    dims: Counter[int] = Counter()
    finite_rows = 0
    malformed = 0
    for value in table["embedding"]:
        try:
            vector = np.asarray(value, dtype=float).reshape(-1)
        except (TypeError, ValueError):
            malformed += 1
            continue
        dims[len(vector)] += 1
        finite_rows += int(np.isfinite(vector).all())
    status_counts = table["embedding_status"].fillna("").astype(str).value_counts()
    missing_nodes = sidecar_nodes - node_set
    extra_nodes = node_set - sidecar_nodes
    errors = []
    if nodes.duplicated().any():
        errors.append("duplicate reaction_node values")
    if missing_nodes:
        errors.append("missing reaction nodes")
    if extra_nodes:
        errors.append("extra reaction nodes")
    if dims != Counter({expected_dim: len(table)}):
        errors.append("unexpected embedding dimensions")
    if finite_rows != len(table):
        errors.append("non-finite or malformed embeddings")
    if len(table) != len(sidecar):
        errors.append("row count differs from sidecar")
    invalid_statuses = set(status_counts.index) - {"ok", "partial"}
    if invalid_statuses:
        errors.append("invalid embedding statuses: " + ", ".join(sorted(invalid_statuses)))
    return {
        "representation": representation,
        "expected_dim": expected_dim,
        "rows": len(table),
        "unique_reaction_nodes": nodes.nunique(),
        "embedding_dims": json.dumps({str(k): v for k, v in sorted(dims.items())}),
        "finite_rows": finite_rows,
        "status_counts": json.dumps(status_counts.to_dict()),
        "missing_reaction_nodes": len(missing_nodes),
        "extra_reaction_nodes": len(extra_nodes),
        "passes": not errors and malformed == 0,
        "path": str(path),
        "validation_errors": "|".join(errors),
    }
