"""Shared output helpers for reaction embeddings."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def write_embedding_parquet(
    output_dir: str | Path,
    name: str,
    embeddings: np.ndarray,
    metadata: pd.DataFrame,
) -> Path:
    """Write a self-contained reaction embedding table."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"reaction_embeddings_{name}.parquet"
    table = metadata.copy()
    table["embedding"] = [row.astype(float).tolist() for row in embeddings]
    table.to_parquet(path, index=False)
    return path
