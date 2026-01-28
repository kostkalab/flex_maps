"""KEGG data download and parsing utilities."""

import csv
import math
import re
from datetime import datetime
from pathlib import Path

import pandas as pd
import requests


def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with snake_case column names."""
    df = df.copy()
    df.columns = (
        df.columns.astype(str)
        .str.strip()
        .str.lower()
        .str.replace(r"[^0-9a-z]+", "_", regex=True)
        .str.strip("_")
    )
    return df


def make_kegg_header() -> str:
    """Fetch KEGG version info and format as a comment header."""
    url = "https://rest.kegg.jp/info/kegg"
    response = requests.get(url, timeout=30)
    kegg_info = response.text.strip()
    now = datetime.now()
    date_str = now.strftime("%Y-%m-%d %H:%M:%S")
    header = f"KEGG data version:\n{kegg_info}\nRetrieved on: {date_str}"
    header = re.sub(r"^", "# ", header, flags=re.MULTILINE)
    return header


def _download_or_load(
    path: Path,
    url: str,
    columns: list[str],
    prefix_strip: dict[str, str] | None = None,
    header: str | None = None,
    quoting: int | None = None,
) -> pd.DataFrame:
    """Load cached CSV or download from KEGG REST API."""
    if path.exists():
        return pd.read_csv(path, comment="#")

    df = pd.read_csv(url, sep="\t", header=None).dropna(how="all")
    df.columns = columns

    if prefix_strip:
        for col, prefix in prefix_strip.items():
            df[col] = df[col].str.replace(prefix, "", regex=False)

    with open(path, "w") as f:
        f.write((header or "") + "\n")

    kwargs = {"index": False, "mode": "a"}
    if quoting is not None:
        kwargs["quoting"] = quoting
    df.to_csv(path, **kwargs)

    return df


def load_ko_to_reaction(kegg_dir: Path, header: str) -> pd.DataFrame:
    """Load KO to reaction mapping."""
    return _download_or_load(
        path=kegg_dir / "dn_ko_to_rn.csv",
        url="https://rest.kegg.jp/link/reaction/ko",
        columns=["ko", "rn"],
        prefix_strip={"ko": "ko:", "rn": "rn:"},
        header=header,
    )


def load_all_reactions(kegg_dir: Path, header: str) -> pd.DataFrame:
    """Load all KEGG reactions."""
    return _download_or_load(
        path=kegg_dir / "dn_reactions.csv",
        url="https://rest.kegg.jp/list/reaction",
        columns=["rn", "description"],
        header=header,
    )


def load_gene_to_ko(kegg_dir: Path, species_code: str, header: str) -> pd.DataFrame:
    """Load species gene to KO mapping."""
    return _download_or_load(
        path=kegg_dir / f"dn_{species_code}gene_to_ko.csv",
        url=f"https://rest.kegg.jp/link/{species_code}/ko",
        columns=["ko", "gene"],
        prefix_strip={"ko": "ko:", "gene": f"{species_code}:"},
        header=header,
    )


def load_reaction_to_module(kegg_dir: Path, header: str) -> pd.DataFrame:
    """Load reaction to module mapping."""
    return _download_or_load(
        path=kegg_dir / "tb_rn_to_md.csv",
        url="https://rest.kegg.jp/link/reaction/module",
        columns=["md", "rn"],
        prefix_strip={"md": "md:", "rn": "rn:"},
        header=header,
    )


def load_module_names(kegg_dir: Path, header: str) -> pd.DataFrame:
    """Load module names."""
    return _download_or_load(
        path=kegg_dir / "dn_modulenames.csv",
        url="https://rest.kegg.jp/list/module",
        columns=["md", "module_name"],
        prefix_strip={"md": "md:"},
        header=header,
        quoting=csv.QUOTE_ALL,
    )


def load_species_modules(kegg_dir: Path, species_code: str, header: str) -> pd.DataFrame:
    """Load species-specific module mappings."""
    df = _download_or_load(
        path=kegg_dir / f"dn_{species_code}_to_modules.csv",
        url=f"https://rest.kegg.jp/link/{species_code}/module",
        columns=["md", "gene"],
        prefix_strip={"md": "md:", "gene": f"{species_code}:"},
        header=header,
    )
    # Entries look like "mmu_M00001" after removing "md:"
    df["md"] = df["md"].str.replace(f"{species_code}_", "", regex=False)
    return df


def load_species_pathways(kegg_dir: Path, species_code: str, header: str) -> pd.DataFrame:
    """Load species pathway mappings."""
    return _download_or_load(
        path=kegg_dir / f"dn_{species_code}_pathways.csv",
        url=f"https://rest.kegg.jp/link/{species_code}/pathway",
        columns=["pathway", "gene"],
        prefix_strip={"pathway": "path:", "gene": f"{species_code}:"},
        header=header,
    )


def load_pathway_to_reaction(kegg_dir: Path, header: str) -> pd.DataFrame:
    """Load pathway to reaction mapping."""
    return _download_or_load(
        path=kegg_dir / "dn_pathway_to_rn.csv",
        url="https://rest.kegg.jp/link/reaction/pathway",
        columns=["path", "rn"],
        prefix_strip={"path": "path:", "rn": "rn:"},
        header=header,
    )


def load_pathway_names(kegg_dir: Path, header: str) -> pd.DataFrame:
    """Load pathway names."""
    return _download_or_load(
        path=kegg_dir / "dn_pathway_names.csv",
        url="https://rest.kegg.jp/list/pathway",
        columns=["pathway", "pathway_name"],
        prefix_strip={"pathway": "path:"},
        header=header,
    )


def load_gene_annotations(kegg_dir: Path, species_code: str, header: str) -> pd.DataFrame:
    """Load or fetch gene annotations (Ensembl, NCBI, EC numbers)."""
    path = kegg_dir / f"{species_code}_ensembl_ec.tsv"

    if path.exists():
        df = clean_columns(pd.read_csv(path, sep="\t", comment="#"))
        df["kegg_gene_id"] = df["kegg_gene_id"].astype(str).str.replace(
            f"{species_code}:", "", regex=False
        )
        ncbi_series = df["ncbi_gene_id"]
        mask = ncbi_series.notna()
        if mask.any() and not ncbi_series[mask].astype(str).str.fullmatch(r"\d+").all():
            print("  Cached gene annotations appear invalid; refreshing cache...")
            path.unlink(missing_ok=True)
            return load_gene_annotations(kegg_dir, species_code, header)
        df["ncbi_gene_id"] = ncbi_series.astype(str)
        return df

    print(f"Fetching gene annotations for {species_code}...")

    # Step 1: KEGG gene → enzyme mappings
    print("  Fetching KEGG enzyme links...")
    ec_url = f"https://rest.kegg.jp/link/{species_code}/enzyme"
    ec_df = pd.read_csv(ec_url, sep="\t", header=None, names=["ec_number", "kegg_gene_id"])
    ec_df["ec_number"] = ec_df["ec_number"].str.replace("ec:", "", regex=False)
    ec_df["kegg_gene_id"] = ec_df["kegg_gene_id"].str.replace(
        f"{species_code}:", "", regex=False
    )

    # Step 2: KEGG gene → NCBI gene ID conversion
    print("  Fetching KEGG → NCBI gene mappings...")
    conv_df = _fetch_kegg_to_ncbi(species_code)

    gene_df = pd.merge(ec_df, conv_df, on="kegg_gene_id", how="left")
    gene_df = gene_df.dropna(subset=["ncbi_gene_id"])
    gene_df["ncbi_gene_id"] = gene_df["ncbi_gene_id"].astype(str)

    # Step 3: Query MyGene for Ensembl IDs and gene symbols
    print("  Fetching gene symbols and Ensembl IDs from MyGene...")
    meta_map = _fetch_gene_metadata(
        gene_df["ncbi_gene_id"].unique().tolist(),
        species_code,
    )
    meta_df = pd.DataFrame.from_dict(meta_map, orient="index").reset_index()
    if meta_df.empty:
        meta_df = pd.DataFrame(columns=["ncbi_gene_id", "gene_symbol", "ensembl_id"])
    else:
        meta_df = meta_df.rename(columns={"index": "ncbi_gene_id"})

    result = pd.merge(
        gene_df,
        meta_df,
        on="ncbi_gene_id",
        how="left",
    )

    result["kegg_enzyme_id"] = result["ec_number"].apply(
        lambda x: f"ec:{x}" if pd.notna(x) and x != "" else None
    )

    result = result[
        [
            "kegg_gene_id",
            "ncbi_gene_id",
            "gene_symbol",
            "ensembl_id",
            "ec_number",
            "kegg_enzyme_id",
        ]
    ]

    # Save human-readable TSV with KEGG prefixes
    file_df = result.copy()
    file_df["kegg_gene_id"] = (
        file_df["kegg_gene_id"].astype(str).apply(lambda g: f"{species_code}:{g}")
    )
    file_df = file_df.fillna("")
    file_df = file_df.rename(
        columns={
            "kegg_gene_id": "KEGG_Gene_ID",
            "ncbi_gene_id": "NCBI_Gene_ID",
            "gene_symbol": "Gene_Symbol",
            "ensembl_id": "Ensembl_ID",
            "ec_number": "EC_Number",
            "kegg_enzyme_id": "KEGG_Enzyme_ID",
        }
    )

    with open(path, "w") as f:
        f.write(header + "\n")
    file_df.to_csv(path, sep="\t", index=False, mode="a")
    print(f"  Saved to {path}")

    return result


def _fetch_kegg_to_ncbi(species_code: str) -> pd.DataFrame:
    """Fetch KEGG gene to NCBI Gene ID conversion table."""
    conv_url = f"https://rest.kegg.jp/conv/ncbi-geneid/{species_code}"
    conv_df = pd.read_csv(conv_url, sep="\t", header=None, names=["kegg_gene", "ncbi_gene"])
    conv_df["kegg_gene_id"] = conv_df["kegg_gene"].str.replace(
        f"{species_code}:", "", regex=False
    )
    conv_df["ncbi_gene_id"] = conv_df["ncbi_gene"].str.replace(
        "ncbi-geneid:", "", regex=False
    )
    return conv_df[["kegg_gene_id", "ncbi_gene_id"]]


def _fetch_gene_metadata(ncbi_ids: list[str], species_code: str) -> dict[str, dict[str, str]]:
    """Fetch gene symbols and Ensembl IDs from MyGene.info."""
    import time

    cleaned_ids: list[str] = []
    for value in ncbi_ids:
        if value is None:
            continue
        if isinstance(value, float):
            if math.isnan(value):
                continue
            value = str(int(value)) if value.is_integer() else str(value)
        else:
            value = str(value).strip()
        if value:
            cleaned_ids.append(value)

    if not cleaned_ids:
        return {}

    species_map = {
        "mmu": "mouse",
        "hsa": "human",
        "dre": "zebrafish",
        "xla": "8355", #- https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=info&id=8355
    }
    species = species_map.get(species_code, species_code)

    meta_map: dict[str, dict[str, str]] = {}
    batch_size = 1000

    for i in range(0, len(cleaned_ids), batch_size):
        batch = cleaned_ids[i : i + batch_size]

        response = requests.post(
            "https://mygene.info/v3/gene",
            json={
                "ids": batch,
                "fields": "ensembl.gene,symbol",
                "species": species,
            },
            timeout=60,
        )
        response.raise_for_status()

        for item in response.json():
            query_id = str(item.get("query", ""))
            symbol = item.get("symbol")
            ensembl = item.get("ensembl")

            if isinstance(ensembl, list):
                ensembl_gene = next(
                    (entry.get("gene") for entry in ensembl if isinstance(entry, dict) and entry.get("gene")),
                    None,
                )
            elif isinstance(ensembl, dict):
                ensembl_gene = ensembl.get("gene")
            else:
                ensembl_gene = None

            entry = meta_map.setdefault(query_id, {})
            if symbol:
                entry["gene_symbol"] = symbol
            if ensembl_gene:
                entry["ensembl_id"] = ensembl_gene

        if i + batch_size < len(cleaned_ids):
            time.sleep(0.5)

    for entry in meta_map.values():
        entry.setdefault("gene_symbol", None)
        entry.setdefault("ensembl_id", None)

    return meta_map
