# Upstream MetaNetX SMILES Supplement

## Goal

Build a lightweight KEGG compound to MetaNetX SMILES artifact for use as a
supplemental SMILES source. This should avoid requiring downstream projects to
parse the full `MNXref.ttl` RDF file.

The intended downstream precedence is:

```text
KEGG/PubChem v1.0.2 SMILES first
MetaNetX SMILES fallback
```

MetaNetX SMILES should be treated as a coverage supplement, not as a silent
replacement for PubChem-derived concrete structures.

## MetaNetX Inputs

Use the tab-delimited MNXref release files, not the Turtle RDF distribution:

```text
chem_xref.tsv
chem_prop.tsv
```

Required fields:

`chem_xref.tsv`

```text
XREF      external compound identifier, e.g. kegg.compound:C00031
MNX_ID    MetaNetX compound identifier
STRING    description from the external resource
```

`chem_prop.tsv`

```text
MNX_ID
name
reference
formula
charge
mass
inchi
inchikey
smiles
```

MetaNetX notes that formula/protonation/R-group handling may differ from the
external reference resource, so preserve provenance and validation flags.

## Extraction Logic

1. Load `chem_xref.tsv`, ignoring comment lines beginning with `#`.
2. Keep rows where `XREF` is a KEGG compound:

```text
kegg.compound:C00031
keggC:C00031
```

3. Normalize KEGG IDs to `Cxxxxx`.
4. Join to `chem_prop.tsv` on `MNX_ID`.
5. Keep rows with non-empty `smiles`.
6. If multiple MetaNetX rows map to the same KEGG compound, choose deterministically:

Recommended precedence:

```text
1. reference starts with kegg.compound: or keggC:
2. RDKit-parseable SMILES
3. no wildcard atom "*"
4. formula count matches RDKit formula after normalizing element order and charge suffix
5. lexicographically smallest MNX_ID as final tie-breaker
```

Do not drop wildcard SMILES globally. Keep them and flag them. They provide
substantial coverage for generic biochemical compounds, and RDKit can parse and
fingerprint them.

## Output Artifact

Write a compressed TSV:

```text
kegg_metanetx_smiles.<YYYYMMDD>.tsv.gz
```

Required columns:

```text
kegg_compound
smiles
smiles_source
mnx_id
mnx_reference
mnx_formula
mnx_charge
mnx_mass
mnx_inchi
mnx_inchikey
rdkit_parse_ok
has_wildcard
rdkit_formula
formula_count_matches_mnx
rdkit_inchikey
inchikey_matches_mnx
```

Use:

```text
smiles_source = metanetx
```

## Validation

Run these checks before publishing:

```text
n_rows
n_unique_kegg_compounds
n_rdkit_parse_ok
n_has_wildcard
n_formula_count_matches_mnx
n_inchikey_matches_mnx
```

Also report current-map coverage if the flex maps compound list is available:

```text
n_current_map_compounds
n_current_map_covered_by_pubchem
n_current_map_covered_by_metanetx
n_current_map_covered_by_pubchem_or_metanetx
n_current_map_uncovered
```

Expected behavior from the proof-of-concept on the current maps:

```text
MetaNetX fallback recovers many reactions lost by PubChem-only coverage.
Most recovered MetaNetX-only structures contain wildcard atoms.
RDKit can parse and Morgan-fingerprint those wildcard structures.
Neural reaction encoders may still treat wildcard SMILES as out-of-distribution.
```

## Minimal Python Skeleton

```python
from pathlib import Path
import re

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


RE_KEGG = re.compile(r"(?:kegg\\.compound:|keggC:)(C\\d{5})$")


def formula_counts(formula: str) -> dict[str, int]:
    formula = str(formula or "").replace("+", "").replace("-", "")
    counts = {}
    for elem, num in re.findall(r"(\\*|[A-Z][a-z]?)(\\d*)", formula):
        counts[elem] = counts.get(elem, 0) + (int(num) if num else 1)
    return counts


def normalize_kegg(xref: str) -> str | None:
    match = RE_KEGG.search(str(xref))
    if match:
        return match.group(1)
    return None


def load_mnxref_tsv(path: Path, columns: list[str]) -> pd.DataFrame:
    return pd.read_csv(path, sep="\\t", comment="#", header=None, names=columns)


def build_kegg_metanetx_smiles(
    chem_xref_path: Path,
    chem_prop_path: Path,
    output_path: Path,
) -> pd.DataFrame:
    xref = load_mnxref_tsv(
        chem_xref_path,
        ["xref", "mnx_id", "xref_description"],
    )
    prop = load_mnxref_tsv(
        chem_prop_path,
        [
            "mnx_id",
            "name",
            "reference",
            "mnx_formula",
            "mnx_charge",
            "mnx_mass",
            "mnx_inchi",
            "mnx_inchikey",
            "smiles",
        ],
    )

    xref["kegg_compound"] = xref["xref"].map(normalize_kegg)
    xref = xref[xref["kegg_compound"].notna()].copy()

    merged = xref.merge(prop, on="mnx_id", how="inner")
    merged = merged[merged["smiles"].notna() & (merged["smiles"].astype(str) != "")]
    merged["smiles_source"] = "metanetx"

    records = []
    for _, row in merged.iterrows():
        smiles = str(row["smiles"])
        mol = Chem.MolFromSmiles(smiles)
        rdkit_parse_ok = mol is not None
        rdkit_formula = ""
        rdkit_inchikey = ""
        if mol is not None:
            rdkit_formula = rdMolDescriptors.CalcMolFormula(mol)
            try:
                rdkit_inchikey = Chem.MolToInchiKey(mol)
            except Exception:
                rdkit_inchikey = ""

        formula_match = (
            bool(row["mnx_formula"])
            and formula_counts(row["mnx_formula"]) == formula_counts(rdkit_formula)
        )
        inchikey_match = (
            bool(row["mnx_inchikey"])
            and bool(rdkit_inchikey)
            and row["mnx_inchikey"] == rdkit_inchikey
        )

        records.append(
            {
                "kegg_compound": row["kegg_compound"],
                "smiles": smiles,
                "smiles_source": "metanetx",
                "mnx_id": row["mnx_id"],
                "mnx_reference": row["reference"],
                "mnx_formula": row["mnx_formula"],
                "mnx_charge": row["mnx_charge"],
                "mnx_mass": row["mnx_mass"],
                "mnx_inchi": row["mnx_inchi"],
                "mnx_inchikey": row["mnx_inchikey"],
                "rdkit_parse_ok": rdkit_parse_ok,
                "has_wildcard": "*" in smiles,
                "rdkit_formula": rdkit_formula,
                "formula_count_matches_mnx": formula_match,
                "rdkit_inchikey": rdkit_inchikey,
                "inchikey_matches_mnx": inchikey_match,
            }
        )

    out = pd.DataFrame.from_records(records)
    out["kegg_reference"] = out["mnx_reference"].astype(str).str.startswith(
        ("kegg.compound:", "keggC:")
    )
    out = out.sort_values(
        [
            "kegg_compound",
            "kegg_reference",
            "rdkit_parse_ok",
            "has_wildcard",
            "formula_count_matches_mnx",
            "mnx_id",
        ],
        ascending=[True, False, False, True, False, True],
    )
    out = out.drop_duplicates("kegg_compound", keep="first")
    out = out.drop(columns=["kegg_reference"])

    out.to_csv(output_path, sep="\\t", index=False)
    return out
```

## Downstream Merge Contract

Downstream code should merge this with the PubChem-derived source as:

```text
if kegg_pubchem_smiles exists:
    use kegg_pubchem_smiles, smiles_source=kegg_pubchem
elif kegg_metanetx_smiles exists:
    use metanetx_smiles, smiles_source=metanetx
else:
    mark missing
```

Do not silently neutralize or charge-standardize the stored SMILES. If needed,
neutralized canonical SMILES can be used for equivalence checks, but the emitted
reaction strings should preserve the selected source SMILES.
