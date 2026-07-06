#!/usr/bin/env python
"""Generate reaction embeddings from an oriented reaction sidecar table."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

def _early_configure_cuda_from_argv() -> None:
    gpu = "0"
    for idx, arg in enumerate(sys.argv):
        if arg == "--gpu" and idx + 1 < len(sys.argv):
            gpu = sys.argv[idx + 1]
            break
        if arg.startswith("--gpu="):
            gpu = arg.split("=", 1)[1]
            break
    if gpu != "cpu":
        os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu)


_early_configure_cuda_from_argv()

import numpy as np
import pandas as pd

from src.embeddings.io import write_embedding_parquet
from src.embeddings.drfp import build_drfp_embeddings
from src.embeddings.reactiont5 import build_reactiont5_4vec_embeddings
from src.embeddings.rxngraphormer import (
    DEFAULT_MODEL_PATH as DEFAULT_RXNGRAPHORMER_MODEL_PATH,
    build_rxngraphormer_embeddings,
)
from src.embeddings.side_mean import build_side_mean_embeddings, load_compound_vectors
from src.upstream_sources import UpstreamSourceResolver


SIDE_MEAN_REPRESENTATIONS = {
    "side_mean_fingerprint_raw": ("fingerprints_raw", "raw"),
    "side_mean_fingerprint_pca256": ("fingerprints_pca256", "pca256"),
}

REPRESENTATIONS = [
    *SIDE_MEAN_REPRESENTATIONS.keys(),
    "drfp",
    "reactiont5_4vec",
    "rxngraphormer_2vec",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reaction_table", help="Oriented reaction sidecar TSV")
    parser.add_argument(
        "--representation",
        default="all",
        choices=[
            "all",
            "dl",
            *REPRESENTATIONS,
            "raw",
            "pca256",
            "rxngraphormer",
            "rxngraphormer_2",
        ],
        help="Embedding representation to generate",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory; defaults to <reaction table dir>/reaction_embeddings",
    )
    parser.add_argument(
        "--sources-config",
        default="config/upstream_sources.yaml",
        help="Upstream artifact configuration",
    )
    parser.add_argument(
        "--lock-path",
        default="data/upstream_sources.lock.yaml",
        help="Resolved artifact lock/provenance file",
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Use already cached artifacts and do not query GitHub",
    )
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Download matched release assets even if cached",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Generate embeddings for the first N reaction rows only",
    )
    parser.add_argument(
        "--gpu",
        default="0",
        help="GPU index for CUDA DL representations; use 'cpu' to force CPU",
    )
    parser.add_argument(
        "--dl-batch-size",
        type=int,
        default=32,
        help="Batch size for DL embedding models",
    )
    parser.add_argument(
        "--rxngraphormer-model-path",
        default=DEFAULT_RXNGRAPHORMER_MODEL_PATH,
        help="Path to RXNGraphormer pretrained_classification_model directory",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    selected = _selected_representations(args.representation)

    table_path = Path(args.reaction_table)
    output_dir = Path(args.output_dir) if args.output_dir else table_path.parent / "reaction_embeddings"

    reactions = pd.read_csv(table_path, sep="\t")
    if args.limit is not None:
        reactions = reactions.head(args.limit).copy()

    needed_artifacts = ["joint_smiles"]
    needed_artifacts.extend(
        SIDE_MEAN_REPRESENTATIONS[name][0]
        for name in selected
        if name in SIDE_MEAN_REPRESENTATIONS
    )

    resolver = UpstreamSourceResolver(args.sources_config)
    resolved = resolver.resolve_many(
        sorted(set(needed_artifacts)),
        lock_path=args.lock_path,
        offline=args.offline,
        force_download=args.force_download,
    )

    report_blocks = []
    for name in selected:
        artifact = None
        extra_metadata = {}
        if name in SIDE_MEAN_REPRESENTATIONS:
            artifact_name, vector_representation = SIDE_MEAN_REPRESENTATIONS[name]
            artifact = resolved[artifact_name]
            vectors, vector_dim = load_compound_vectors(
                artifact.path,
                vector_representation,
            )
            embeddings, metadata, status_counts = build_side_mean_embeddings(
                reactions,
                vectors,
                vector_dim,
            )
        elif name == "drfp":
            embeddings, metadata, status_counts = build_drfp_embeddings(reactions)
        elif name == "reactiont5_4vec":
            _configure_gpu(args.gpu)
            embeddings, metadata, status_counts, extra_metadata = (
                build_reactiont5_4vec_embeddings(
                    reactions,
                    batch_size=args.dl_batch_size,
                    device=_torch_device(args.gpu),
                    local_files_only=args.offline,
                )
            )
        elif name == "rxngraphormer_2vec":
            _configure_gpu(args.gpu)
            embeddings, metadata, status_counts, extra_metadata = (
                build_rxngraphormer_embeddings(
                    reactions,
                    model_path=args.rxngraphormer_model_path,
                    batch_size=args.dl_batch_size,
                    include_reverse=True,
                )
            )
        else:
            raise ValueError(f"Unsupported representation {name!r}")

        output_path = write_embedding_parquet(output_dir, name, embeddings, metadata)
        report_blocks.append(
            _report_block(
                name=name,
                reaction_table=table_path,
                artifact=artifact,
                joint_smiles=resolved["joint_smiles"],
                embeddings=embeddings,
                status_counts=status_counts,
                metadata=metadata,
                extra_metadata=extra_metadata,
                output_path=output_path,
            )
        )

    report_path = output_dir / "reaction_embedding_generation_report.md"
    report_path.write_text("\n\n".join(report_blocks) + "\n")
    print(f"Wrote report: {report_path}")


def _selected_representations(value: str) -> list[str]:
    if value == "all":
        return list(REPRESENTATIONS)
    if value == "dl":
        return ["reactiont5_4vec", "rxngraphormer_2vec"]
    if value == "raw":
        return ["side_mean_fingerprint_raw"]
    if value == "pca256":
        return ["side_mean_fingerprint_pca256"]
    if value in {"rxngraphormer", "rxngraphormer_2"}:
        return ["rxngraphormer_2vec"]
    return [value]


def _configure_gpu(gpu: str) -> None:
    if gpu != "cpu":
        os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu)
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")
    os.environ.setdefault("DGLBACKEND", "pytorch")
    os.environ.setdefault("DGL_DOWNLOAD_DIR", "/tmp/dgl")
    os.environ.setdefault("XDG_CACHE_HOME", "/tmp/cache")


def _torch_device(gpu: str) -> str | None:
    if gpu == "cpu":
        return "cpu"
    return None


def _report_block(
    name: str,
    reaction_table: Path,
    artifact,
    joint_smiles,
    embeddings,
    status_counts: dict[str, int],
    metadata: pd.DataFrame,
    extra_metadata: dict,
    output_path: Path,
) -> str:
    status_lines = "\n".join(
        f"- {status}: {count}" for status, count in sorted(status_counts.items())
    )
    non_ok = metadata[metadata["embedding_status"] != "ok"].head(10)
    if non_ok.empty:
        non_ok_lines = "- none"
    else:
        non_ok_lines = "\n".join(
            f"- {row.reaction_node}: {row.embedding_status}; {row.error_message}"
            for row in non_ok.itertuples(index=False)
        )
    source_lines = [
        f"- source_joint_smiles: {joint_smiles.asset_name}",
        f"- source_joint_smiles_sha256: {joint_smiles.sha256}",
    ]
    if artifact is not None:
        source_lines.extend(
            [
                f"- source_compound_vectors: {artifact.asset_name}",
                f"- source_compound_vectors_sha256: {artifact.sha256}",
            ]
        )
    for key, value in sorted(extra_metadata.items()):
        source_lines.append(f"- {key}: {value}")
    source_text = "\n".join(source_lines)
    return f"""# {name}

- reaction_table: {reaction_table}
- n_reactions_total: {embeddings.shape[0]}
- embedding_dim: {embeddings.shape[1]}
- finite_rows: {int(np.isfinite(embeddings).all(axis=1).sum())}
{source_text}
- output_parquet: {output_path}

embedding_status:
{status_lines}

non_ok_examples:
{non_ok_lines}
"""


if __name__ == "__main__":
    main()
