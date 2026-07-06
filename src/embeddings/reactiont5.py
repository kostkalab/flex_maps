"""ReactionT5 encoder embeddings."""

from __future__ import annotations

import numpy as np
import pandas as pd


FORWARD_MODEL = "sagawa/ReactionT5v2-forward"
RETRO_MODEL = "sagawa/ReactionT5v2-retrosynthesis"


def build_reactiont5_4vec_embeddings(
    reactions: pd.DataFrame,
    forward_model_name: str = FORWARD_MODEL,
    retro_model_name: str = RETRO_MODEL,
    batch_size: int = 32,
    device: str | None = None,
    local_files_only: bool = False,
) -> tuple[np.ndarray, pd.DataFrame, dict, dict]:
    """Build [Ef(I->O), Er(O->I), Er(I->O), Ef(O->I)] embeddings."""
    import torch
    from transformers import AutoTokenizer, T5EncoderModel

    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"

    i_to_o = reactions["I_to_O_string"].fillna("").astype(str).tolist()
    o_to_i = reactions["O_to_I_string"].fillna("").astype(str).tolist()

    forward_tokenizer = AutoTokenizer.from_pretrained(
        forward_model_name,
        local_files_only=local_files_only,
    )
    forward_model = (
        T5EncoderModel.from_pretrained(
            forward_model_name,
            local_files_only=local_files_only,
        )
        .to(device)
        .eval()
    )
    retro_tokenizer = AutoTokenizer.from_pretrained(
        retro_model_name,
        local_files_only=local_files_only,
    )
    retro_model = (
        T5EncoderModel.from_pretrained(
            retro_model_name,
            local_files_only=local_files_only,
        )
        .to(device)
        .eval()
    )

    ef_i_to_o = _encode(forward_tokenizer, forward_model, i_to_o, batch_size, device)
    er_o_to_i = _encode(retro_tokenizer, retro_model, o_to_i, batch_size, device)
    er_i_to_o = _encode(retro_tokenizer, retro_model, i_to_o, batch_size, device)
    ef_o_to_i = _encode(forward_tokenizer, forward_model, o_to_i, batch_size, device)
    embeddings = np.concatenate(
        [ef_i_to_o, er_o_to_i, er_i_to_o, ef_o_to_i],
        axis=1,
    ).astype(np.float32)

    metadata = _ok_metadata(reactions)
    status_counts = {"ok": int(len(reactions))}
    run_metadata = {
        "reactiont5_forward_model": forward_model_name,
        "reactiont5_retro_model": retro_model_name,
        "reactiont5_order": "E_forward(I_to_O)|E_retro(O_to_I)|E_retro(I_to_O)|E_forward(O_to_I)",
        "device": device,
        "batch_size": batch_size,
        "hidden_dim_each": int(ef_i_to_o.shape[1]),
        "local_files_only": local_files_only,
    }
    return embeddings, metadata, status_counts, run_metadata


def _encode(tokenizer, model, strings: list[str], batch_size: int, device: str) -> np.ndarray:
    import torch

    chunks = []
    for start in range(0, len(strings), batch_size):
        batch_strings = strings[start : start + batch_size]
        batch = tokenizer(batch_strings, padding=True, return_tensors="pt").to(device)
        with torch.no_grad():
            hidden = model(**batch).last_hidden_state
            mask = batch["attention_mask"].unsqueeze(-1)
            pooled = (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1)
            pooled = torch.nn.functional.normalize(pooled, p=2, dim=1)
        chunks.append(pooled.detach().cpu().numpy())
    return np.vstack(chunks).astype(np.float32)


def _ok_metadata(reactions: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "reaction_node": reactions["reaction_node"].astype(str),
            "embedding_status": "ok",
            "error_message": "",
        }
    )
