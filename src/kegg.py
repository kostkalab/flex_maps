"""KEGG data download and parsing utilities."""

import csv
import gzip
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


def _get_timestamped_cache_path(cache_dir: Path, base_name: str, max_age_days: int = 7) -> Path:
    """Get path for timestamped cache file, reusing recent files if available.
    
    Args:
        cache_dir: Directory for cache files
        base_name: Base filename (e.g., 'gene_info.gz')
        max_age_days: Maximum age in days to reuse existing cache (default: 7)
    
    Returns:
        Path to cache file (existing or new with today's date)
    """
    from datetime import datetime, timedelta
    
    cache_dir.mkdir(parents=True, exist_ok=True)
    today = datetime.now().strftime("%Y%m%d")
    
    # Check for today's file first
    today_file = cache_dir / f"{base_name.replace('.gz', '')}.{today}.gz"
    if today_file.exists():
        return today_file
    
    # Look for recent files within max_age_days
    name_base = base_name.replace('.gz', '')
    existing = sorted(cache_dir.glob(f"{name_base}.*.gz"), reverse=True)
    
    for cached_file in existing:
        try:
            # Extract date from filename: gene_info.20260206.gz -> 20260206
            parts = cached_file.stem.split('.')
            if len(parts) >= 2:
                date_str = parts[-1]
                file_date = datetime.strptime(date_str, "%Y%m%d")
                age_days = (datetime.now() - file_date).days
                
                if age_days <= max_age_days:
                    print(f"  Using cached {cached_file.name} ({age_days} days old)")
                    return cached_file
        except (ValueError, IndexError):
            continue
    
    # No recent file found, return path for new download
    return today_file


def _download_ncbi_file(url: str, cache_path: Path) -> None:
    """Download NCBI flat file with progress reporting."""
    if cache_path.exists():
        return
    
    print(f"  Downloading {cache_path.name}...")
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    
    response = requests.get(url, stream=True, timeout=120)
    response.raise_for_status()
    
    total_size = int(response.headers.get('content-length', 0))
    downloaded = 0
    
    with open(cache_path, 'wb') as f:
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
                downloaded += len(chunk)
                if total_size > 0:
                    pct = (downloaded / total_size) * 100
                    print(f"    {pct:.1f}% ({downloaded // 1024 // 1024}MB / {total_size // 1024 // 1024}MB)", end="\r")
    
    print(f"\n    Downloaded {cache_path.name}")


def _load_ncbi_gene_info(cache_dir: Path, tax_id: str) -> dict[str, str]:
    """Load gene symbols from NCBI gene_info.gz file."""
    cache_path = _get_timestamped_cache_path(cache_dir, "gene_info.gz")
    url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz"
    
    _download_ncbi_file(url, cache_path)
    
    print(f"  Parsing gene_info.gz for tax_id {tax_id}...")
    gene_symbols = {}
    
    with gzip.open(cache_path, 'rt', encoding='utf-8') as f:
        # Skip header line
        header = f.readline()
        
        for line_num, line in enumerate(f, start=2):
            if line_num % 100000 == 0:
                print(f"    Processed {line_num:,} lines, found {len(gene_symbols):,} genes for species", end="\r")
            
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            
            if fields[0] == tax_id:  # tax_id is first column
                ncbi_id = fields[1]  # GeneID is second column
                symbol = fields[2]   # Symbol is third column
                if symbol and symbol != '-':
                    gene_symbols[ncbi_id] = symbol
    
    print(f"\n    Loaded {len(gene_symbols):,} gene symbols for tax_id {tax_id}")
    return gene_symbols


def _load_ncbi_gene2ensembl(cache_dir: Path, tax_id: str) -> dict[str, str]:
    """Load Ensembl IDs from NCBI gene2ensembl.gz file."""
    cache_path = _get_timestamped_cache_path(cache_dir, "gene2ensembl.gz")
    url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"
    
    _download_ncbi_file(url, cache_path)
    
    print(f"  Parsing gene2ensembl.gz for tax_id {tax_id}...")
    ensembl_ids = {}
    
    with gzip.open(cache_path, 'rt', encoding='utf-8') as f:
        # Skip header line
        header = f.readline()
        
        for line_num, line in enumerate(f, start=2):
            if line_num % 100000 == 0:
                print(f"    Processed {line_num:,} lines, found {len(ensembl_ids):,} mappings for species", end="\r")
            
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            
            if fields[0] == tax_id:  # tax_id is first column
                ncbi_id = fields[1]       # GeneID is second column
                ensembl_id = fields[2]    # Ensembl_gene_identifier is third column
                if ensembl_id and ensembl_id != '-':
                    ensembl_ids[ncbi_id] = ensembl_id
    
    print(f"\n    Loaded {len(ensembl_ids):,} Ensembl mappings for tax_id {tax_id}")
    return ensembl_ids


def _get_metadata_from_ncbi(data_dir: Path, tax_id: str, ncbi_ids: list[str]) -> dict[str, dict[str, str]]:
    """Fetch gene symbols and Ensembl IDs from NCBI flat files.
    
    Args:
        data_dir: Root data directory (not kegg-specific)
        tax_id: NCBI taxonomy ID (e.g., '10090' for mouse)
        ncbi_ids: List of NCBI gene IDs to lookup
    """
    cache_dir = data_dir / "ncbi"
    cache_dir.mkdir(parents=True, exist_ok=True)
    
    # Load both mappings
    gene_symbols = _load_ncbi_gene_info(cache_dir, tax_id)
    ensembl_ids = _load_ncbi_gene2ensembl(cache_dir, tax_id)
    
    # Build result dictionary
    print(f"  Building metadata dictionary for {len(ncbi_ids):,} genes...")
    result = {}
    for ncbi_id in ncbi_ids:
        meta = {}
        if ncbi_id in gene_symbols:
            meta['gene_symbol'] = gene_symbols[ncbi_id]
        if ncbi_id in ensembl_ids:
            meta['ensembl_id'] = ensembl_ids[ncbi_id]
        if meta:  # Only add if we have at least one field
            result[ncbi_id] = meta
    
    print(f"  Found metadata for {len(result):,} / {len(ncbi_ids):,} genes ({100*len(result)/len(ncbi_ids):.1f}%)")
    return result


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
    df = _download_or_load(
        path=kegg_dir / f"dn_{species_code}gene_to_ko.csv",
        url=f"https://rest.kegg.jp/link/{species_code}/ko",
        columns=["ko", "gene"],
        prefix_strip={"ko": "ko:", "gene": f"{species_code}:"},
        header=header,
    )
    df["gene"] = df["gene"].astype(str)
    return df


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
    df["gene"] = df["gene"].astype(str)
    return df


def load_species_pathways(kegg_dir: Path, species_code: str, header: str) -> pd.DataFrame:
    """Load species pathway mappings."""
    df = _download_or_load(
        path=kegg_dir / f"dn_{species_code}_pathways.csv",
        url=f"https://rest.kegg.jp/link/{species_code}/pathway",
        columns=["pathway", "gene"],
        prefix_strip={"pathway": "path:", "gene": f"{species_code}:"},
        header=header,
    )
    df["gene"] = df["gene"].astype(str)
    return df


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


def load_gene_annotations(kegg_dir: Path, species_code: str, header: str, ncbi_taxonomy_id: str = "") -> pd.DataFrame:
    """Load or fetch gene annotations (Ensembl, NCBI, EC numbers).
    
    Args:
        kegg_dir: Path to KEGG data directory
        species_code: KEGG species code (e.g., 'mmu')
        header: Header comment for TSV file
        ncbi_taxonomy_id: NCBI taxonomy ID (e.g., '10090' for mouse). Required for new fetches.
    """
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

    # Step 1: KEGG gene → NCBI gene ID conversion
    print("  Fetching KEGG → NCBI gene mappings...")
    conv_df = _fetch_kegg_to_ncbi(species_code)

    # Step 2: KEGG gene → enzyme mappings
    print("  Fetching KEGG enzyme links...")
    ec_url = f"https://rest.kegg.jp/link/{species_code}/enzyme"
    ec_df = pd.read_csv(ec_url, sep="\t", header=None, names=["ec_number", "kegg_gene_id"])
    ec_df = ec_df.astype(str)
    ec_df["ec_number"] = ec_df["ec_number"].str.replace("ec:", "", regex=False)
    ec_df["kegg_gene_id"] = ec_df["kegg_gene_id"].str.replace(
        f"{species_code}:", "", regex=False
    )

    # Base the annotations on all genes that have NCBI mappings, then join enzymes if they exist
    gene_df = pd.merge(conv_df, ec_df, on="kegg_gene_id", how="left")
    gene_df["ncbi_gene_id"] = gene_df["ncbi_gene_id"].astype(str)

    # Step 3: Get gene symbols and Ensembl IDs from NCBI flat files
    print("  Fetching gene symbols and Ensembl IDs from NCBI flat files...")
    if not ncbi_taxonomy_id:
        raise ValueError(f"ncbi_taxonomy_id is required for species {species_code}. Add it to the species YAML config.")
    
    data_dir = kegg_dir.parent  # Go up from data/kegg to data/
    meta_map = _get_metadata_from_ncbi(
        data_dir,
        ncbi_taxonomy_id,
        gene_df["ncbi_gene_id"].unique().tolist()
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

    result["kegg_gene_id"] = result["kegg_gene_id"].astype(str)
    return result


def _fetch_kegg_to_ncbi(species_code: str) -> pd.DataFrame:
    """Fetch KEGG gene to NCBI Gene ID conversion table."""
    conv_url = f"https://rest.kegg.jp/conv/ncbi-geneid/{species_code}"
    conv_df = pd.read_csv(conv_url, sep="\t", header=None, names=["kegg_gene", "ncbi_gene"])
    conv_df = conv_df.astype(str)
    conv_df["kegg_gene_id"] = conv_df["kegg_gene"].str.replace(
        f"{species_code}:", "", regex=False
    )
    conv_df = conv_df.astype(str)
    conv_df["ncbi_gene_id"] = conv_df["ncbi_gene"].str.replace(
        "ncbi-geneid:", "", regex=False
    )
    return conv_df[["kegg_gene_id", "ncbi_gene_id"]]
