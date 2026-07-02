"""Chemical identifier and SMILES helpers."""

from __future__ import annotations

import gzip
import pickle
import re
from pathlib import Path
from urllib.parse import urlparse

import requests


DEFAULT_KEGG_SMILES_URL = (
    "https://github.com/kostkalab/flex_embeddings/releases/download/v1.0.3/"
    "kegg_smiles_joint.20260701.tsv.gz"
)

RE_KEGG_COMPOUND = re.compile(r"(C\d{5})")


def normalize_kegg_compound_id(value: object) -> str | None:
    """Normalize a KEGG compound identifier to ``C00031`` form."""
    if value is None:
        return None
    match = RE_KEGG_COMPOUND.search(str(value))
    if not match:
        return None
    return match.group(1)


def _cache_path_for_url(url: str, cache_dir: Path) -> Path:
    parsed = urlparse(url)
    filename = Path(parsed.path).name
    return cache_dir / filename


def resolve_smiles_source(source: str | Path | None, cache_dir: Path) -> Path | None:
    """Resolve a local or URL SMILES source, downloading URLs into ``cache_dir``."""
    if source is None:
        source = DEFAULT_KEGG_SMILES_URL

    source_str = str(source)
    if source_str.startswith(("http://", "https://")):
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path = _cache_path_for_url(source_str, cache_dir)
        if cache_path.exists():
            return cache_path

        print(f"  Downloading KEGG SMILES map to {cache_path}...")
        try:
            response = requests.get(source_str, timeout=120)
            response.raise_for_status()
        except requests.RequestException as exc:
            print(f"  Warning: failed to download KEGG SMILES map: {exc}")
            return None

        cache_path.write_bytes(response.content)
        return cache_path

    path = Path(source_str)
    if not path.exists():
        print(f"  Warning: KEGG SMILES map not found: {path}")
        return None
    return path


def _truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def _bool_text(value: bool) -> str:
    return "true" if value else "false"


def load_kegg_smiles_records(
    source: str | Path | None,
    cache_dir: Path,
) -> dict[str, dict[str, str | bool]]:
    """Load KEGG compound SMILES records from the joint TSV or legacy pickle files."""
    path = resolve_smiles_source(source, cache_dir)
    if path is None:
        return {}

    if "".join(path.suffixes[-2:]) == ".tsv.gz" or path.suffix == ".tsv":
        import pandas as pd

        df = pd.read_csv(path, sep="\t")
        if not {"kegg_id", "selected_smiles"}.issubset(set(df.columns)):
            raise ValueError(f"Unsupported SMILES TSV schema in {path}")

        records = {}
        for _, row in df.iterrows():
            cid = normalize_kegg_compound_id(row.get("kegg_id"))
            smiles = row.get("selected_smiles")
            if cid is None or smiles is None:
                continue
            smiles_value = str(smiles).strip()
            if not smiles_value or smiles_value.lower() == "nan":
                continue
            records[cid] = {
                "smiles": smiles_value,
                "smiles_source": str(row.get("selected_source", "joint")).strip(),
                "smiles_has_wildcard": _bool_text("*" in smiles_value),
                "smiles_rdkit_parse_ok": _bool_text(
                    _truthy(row.get("pubchem_parse_ok"))
                    or _truthy(row.get("metanetx_parse_ok"))
                ),
                "smiles_selection_reason": str(
                    row.get("selection_reason", "")
                ).strip(),
            }
        print(f"  Loaded {len(records):,} KEGG compound SMILES from {path}")
        return records

    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rb") as handle:
            raw = pickle.load(handle)
    except ModuleNotFoundError as exc:
        if exc.name in {"pyarrow", "rdkit"}:
            raise ModuleNotFoundError(
                "Loading the KEGG/PubChem SMILES pickle requires "
                f"{exc.name}. Install pyarrow and rdkit in the pipeline "
                "environment or provide a plain dict pickle with KEGG "
                "compound IDs mapped to SMILES."
            ) from exc
        raise

    records: dict[str, dict[str, str | bool]] = {}
    if hasattr(raw, "columns") and {"kegg", "smiles"}.issubset(set(raw.columns)):
        items = zip(raw["kegg"], raw["smiles"], strict=False)
    elif isinstance(raw, dict):
        items = raw.items()
    else:
        try:
            items = dict(raw).items()
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Unsupported SMILES map format in {path}") from exc

    for key, value in items:
        cid = normalize_kegg_compound_id(key)
        if cid is None or value is None:
            continue
        smiles_value = str(value).strip()
        if smiles_value:
            records[cid] = {
                "smiles": smiles_value,
                "smiles_source": "kegg_pubchem",
                "smiles_has_wildcard": _bool_text("*" in smiles_value),
                "smiles_rdkit_parse_ok": "true",
                "smiles_selection_reason": "legacy_smiles_map",
            }

    print(f"  Loaded {len(records):,} KEGG compound SMILES from {path}")
    return records


def load_kegg_smiles_map(source: str | Path | None, cache_dir: Path) -> dict[str, str]:
    """Load the KEGG compound to SMILES mapping."""
    records = load_kegg_smiles_records(source, cache_dir)
    return {cid: str(record["smiles"]) for cid, record in records.items()}


def canonicalize_smiles(smiles: str) -> str:
    """Canonicalize a molecule SMILES when RDKit is installed."""
    try:
        from rdkit import Chem, RDLogger
    except ImportError:
        return smiles

    RDLogger.DisableLog("rdApp.warning")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles
    return Chem.MolToSmiles(mol, canonical=True)


def reaction_side_smiles(compounds: list[str], smiles_map: dict[str, str]) -> tuple[str, list[str]]:
    """Return deterministic dot-separated SMILES and missing compound IDs."""
    side_smiles = []
    missing = []
    for compound in compounds:
        cid = normalize_kegg_compound_id(compound)
        if cid is None:
            continue
        value = smiles_map.get(cid)
        if value:
            side_smiles.append(canonicalize_smiles(value))
        else:
            missing.append(cid)
    return ".".join(sorted(side_smiles)), sorted(set(missing))
