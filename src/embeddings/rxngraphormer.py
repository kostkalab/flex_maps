"""RXNGraphormer pretrained reaction embeddings."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_MODEL_PATH = (
    "data/rxngraphormer/RXNGraphormer/model_path/pretrained_classification_model"
)


def build_rxngraphormer_embeddings(
    reactions: pd.DataFrame,
    model_path: str | Path = DEFAULT_MODEL_PATH,
    batch_size: int = 64,
    include_reverse: bool = False,
) -> tuple[np.ndarray, pd.DataFrame, dict, dict]:
    """Build E_graphormer(I_to_O), optionally concat E_graphormer(O_to_I)."""
    from rxngraphormer.rxn_emb import RXNEMB

    model_path = Path(model_path)
    if not (model_path / "parameters.json").exists():
        raise FileNotFoundError(f"Missing RXNGraphormer parameters.json: {model_path}")
    if not (model_path / "model" / "valid_checkpoint.pt").exists():
        raise FileNotFoundError(
            f"Missing RXNGraphormer checkpoint: {model_path / 'model' / 'valid_checkpoint.pt'}"
        )

    calculator = RXNEMB(pretrained_model_path=str(model_path), model_type="classifier")
    i_to_o = reactions["I_to_O_string"].fillna("").astype(str).tolist()
    forward = _as_numpy(calculator.gen_rxn_emb(i_to_o, batch_size=batch_size))

    if include_reverse:
        o_to_i = reactions["O_to_I_string"].fillna("").astype(str).tolist()
        reverse = _as_numpy(calculator.gen_rxn_emb(o_to_i, batch_size=batch_size))
        embeddings = np.concatenate([forward, reverse], axis=1).astype(np.float32)
        order = "E_graphormer(I_to_O)|E_graphormer(O_to_I)"
    else:
        embeddings = forward.astype(np.float32)
        order = "E_graphormer(I_to_O)"

    metadata = _ok_metadata(reactions)
    status_counts = {"ok": int(len(reactions))}
    run_metadata = {
        "rxngraphormer_model_path": str(model_path),
        "rxngraphormer_order": order,
        "batch_size": batch_size,
        "forward_dim": int(forward.shape[1]),
    }
    return embeddings, metadata, status_counts, run_metadata


def _as_numpy(value) -> np.ndarray:
    if hasattr(value, "detach"):
        value = value.detach().cpu().numpy()
    return np.asarray(value, dtype=np.float32)


def _ok_metadata(reactions: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "reaction_node": reactions["reaction_node"].astype(str),
            "embedding_status": "ok",
            "error_message": "",
        }
    )
