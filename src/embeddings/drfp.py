"""DRFP reaction fingerprints."""

from __future__ import annotations

import numpy as np
import pandas as pd


DRFP_DIM = 2048


def build_drfp_embeddings(
    reactions: pd.DataFrame,
    n_folded_length: int = DRFP_DIM,
) -> tuple[np.ndarray, pd.DataFrame, dict]:
    """Encode oriented reaction SMILES as DRFP vectors."""
    embedding_rows = []
    metadata_rows = []
    status_counts: dict[str, int] = {}

    for _, row in reactions.reset_index(drop=True).iterrows():
        reaction_smiles = _reaction_smiles(row)
        status, error, vector = _encode_one(reaction_smiles, n_folded_length)
        status_counts[status] = status_counts.get(status, 0) + 1
        embedding_rows.append(vector)
        metadata_rows.append(
            {
                "reaction_node": row.get("reaction_node", ""),
                "embedding_status": status,
                "error_message": error,
            }
        )

    embeddings = np.vstack(embedding_rows).astype(np.uint8)
    metadata = pd.DataFrame(
        metadata_rows,
        columns=["reaction_node", "embedding_status", "error_message"],
    )
    return embeddings, metadata, status_counts


def _encode_one(
    reaction_smiles: str,
    n_folded_length: int,
) -> tuple[str, str, np.ndarray]:
    if not reaction_smiles or reaction_smiles == ">>":
        return (
            "failed",
            "missing_reaction_smiles",
            np.zeros(n_folded_length, dtype=np.uint8),
        )
    try:
        from drfp import DrfpEncoder

        vector = DrfpEncoder.encode(
            reaction_smiles,
            n_folded_length=n_folded_length,
            show_progress_bar=False,
        )[0]
    except Exception as exc:  # DRFP can raise parser-specific exceptions.
        return (
            "failed",
            f"drfp_encode_error:{type(exc).__name__}:{exc}",
            np.zeros(n_folded_length, dtype=np.uint8),
        )
    return "ok", "", np.asarray(vector, dtype=np.uint8)


def _reaction_smiles(row: pd.Series) -> str:
    value = row.get("I_to_O_string", "")
    if value is None or pd.isna(value):
        return ""
    return str(value)
