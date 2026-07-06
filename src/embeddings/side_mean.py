"""Side-specific mean compound-vector reaction embeddings."""

from __future__ import annotations

import gzip
import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd


COMMON_METADATA_COLUMNS = [
    "reaction_node",
    "embedding_status",
    "error_message",
]


def load_compound_vectors(path: str | Path, representation: str) -> tuple[dict[str, np.ndarray], int]:
    """Load KEGG compound vectors from a release fingerprint artifact."""
    with gzip.open(path, "rb") as fh:
        obj = pickle.load(fh)
    if not isinstance(obj, pd.DataFrame):
        raise TypeError(f"Expected a pandas DataFrame in {path}, found {type(obj)}")
    if "kegg_id" not in obj.columns:
        raise ValueError(f"Fingerprint artifact {path} has no kegg_id column")

    if representation == "pca256":
        vector_columns = sorted(
            [col for col in obj.columns if re.fullmatch(r"PC\d+", str(col))],
            key=lambda col: int(str(col)[2:]),
        )
        if not vector_columns:
            raise ValueError(f"No PC columns found in {path}")
        matrix = obj[vector_columns].to_numpy(dtype=np.float32)
    elif representation == "raw":
        if "fingerprint" not in obj.columns:
            raise ValueError(f"Raw fingerprint artifact {path} has no fingerprint column")
        matrix = np.vstack(
            obj["fingerprint"].map(lambda value: np.asarray(value, dtype=np.float32))
        )
    else:
        raise ValueError(f"Unknown side-mean representation {representation!r}")

    ids = obj["kegg_id"].astype(str).tolist()
    vectors = {cid: matrix[idx] for idx, cid in enumerate(ids)}
    return vectors, int(matrix.shape[1])


def build_side_mean_embeddings(
    reactions: pd.DataFrame,
    compound_vectors: dict[str, np.ndarray],
    vector_dim: int,
) -> tuple[np.ndarray, pd.DataFrame, dict]:
    """Build concat(mean(inputs), mean(outputs)) for every reaction row."""
    embedding_rows = []
    metadata_rows = []
    status_counts: dict[str, int] = {}

    for _, row in reactions.reset_index(drop=True).iterrows():
        inputs = _split_compounds(row.get("input_compounds", ""))
        outputs = _split_compounds(row.get("output_compounds", ""))
        status, error, vector = _reaction_vector(
            inputs, outputs, compound_vectors, vector_dim
        )
        status_counts[status] = status_counts.get(status, 0) + 1

        embedding_rows.append(vector)
        metadata_rows.append(_metadata_row(row, inputs, outputs, status, error))

    embeddings = np.vstack(embedding_rows).astype(np.float32)
    metadata = pd.DataFrame(metadata_rows, columns=COMMON_METADATA_COLUMNS)
    return embeddings, metadata, status_counts


def _reaction_vector(
    inputs: list[str],
    outputs: list[str],
    compound_vectors: dict[str, np.ndarray],
    vector_dim: int,
) -> tuple[str, str, np.ndarray]:
    if not inputs or not outputs:
        error = "empty_input_side" if not inputs else "empty_output_side"
        return "failed", error, np.full(vector_dim * 2, np.nan, dtype=np.float32)

    available_inputs = [cid for cid in inputs if cid in compound_vectors]
    available_outputs = [cid for cid in outputs if cid in compound_vectors]
    missing_inputs = [cid for cid in inputs if cid not in compound_vectors]
    missing_outputs = [cid for cid in outputs if cid not in compound_vectors]

    if not available_inputs or not available_outputs:
        missing = sorted(set(missing_inputs + missing_outputs))
        error = "missing_side_compound_embeddings:" + "|".join(missing)
        return "failed", error, np.full(vector_dim * 2, np.nan, dtype=np.float32)

    input_mean = np.mean([compound_vectors[cid] for cid in available_inputs], axis=0)
    output_mean = np.mean([compound_vectors[cid] for cid in available_outputs], axis=0)
    missing = sorted(set(missing_inputs + missing_outputs))
    if missing:
        status = "partial"
        error = "missing_compound_embeddings:" + "|".join(missing)
    else:
        status = "ok"
        error = ""
    return status, error, np.concatenate([input_mean, output_mean]).astype(np.float32)


def _metadata_row(
    row: pd.Series,
    inputs: list[str],
    outputs: list[str],
    status: str,
    error: str,
) -> dict:
    return {
        "reaction_node": row.get("reaction_node", ""),
        "embedding_status": status,
        "error_message": error,
    }


def _split_compounds(value: object) -> list[str]:
    if value is None or pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    return [item for item in text.split("|") if item]
