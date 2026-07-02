"""Main pipeline orchestration."""

from dataclasses import dataclass, field
from pathlib import Path

import networkx as nx
import pandas as pd
from weasyprint import HTML, CSS
from markdown import markdown

from . import chem, graph, kegg, metanetx, reactions, thermo


@dataclass
class SpeciesConfig:
    """Configuration for a species-specific metabolic map."""

    species_code: str
    name: str
    ncbi_taxonomy_id: str = ""
    min_enzyme_coverage: float = 0.66
    rescue_modules: list[str] = field(default_factory=list)
    excluded_modules: list[str] = field(default_factory=list)
    output_prefix: str = "metabolic_graph"


@dataclass
class PipelineResult:
    """Results from pipeline execution."""

    graph: nx.DiGraph
    metrics: dict
    dropped_modules: pd.DataFrame


def build_reaction_sets(
    ko_to_rn: pd.DataFrame,
    gene_to_ko: pd.DataFrame,
    all_reactions: pd.DataFrame,
) -> dict[str, set]:
    """
    Build sets of reaction IDs for filtering.
    
    Returns dict with keys:
        - all: all KEGG reactions
        - with_enzyme: reactions with KO assignments
        - without_enzyme: reactions without KO assignments
        - species: reactions with species-specific genes
    """
    all_rn = set(all_reactions["rn"])
    with_enzyme = set(ko_to_rn["rn"])
    without_enzyme = all_rn - with_enzyme

    # Merge gene -> KO -> reaction
    gene_to_rn = pd.merge(gene_to_ko, ko_to_rn, on="ko", how="inner")
    species_rn = set(gene_to_rn["rn"])

    return {
        "all": all_rn,
        "with_enzyme": with_enzyme,
        "without_enzyme": without_enzyme,
        "species": species_rn,
        "gene_to_rn": gene_to_rn,
    }


def filter_modules(
    rn_to_module: pd.DataFrame,
    module_names: pd.DataFrame,
    species_modules: pd.DataFrame,
    species_reactions: set,
    enzyme_reactions: set,
    min_coverage: float,
    rescue_modules: list[str],
    excluded_modules: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Filter modules to those well-represented in the species.
    
    Returns:
        Tuple of (selected_modules_df, dropped_modules_df)
    """
    # Calculate coverage stats per module
    stats = (
        rn_to_module.assign(
            has_ko=lambda df: df["rn"].isin(enzyme_reactions),
            has_species=lambda df: df["rn"].isin(species_reactions),
        )
        .groupby("md", as_index=False)
        .agg(
            reactions=("rn", "nunique"),
            reactions_with_ko=("has_ko", "sum"),
            reactions_with_species=("has_species", "sum"),
        )
    )

    stats["coverage"] = stats["reactions_with_species"] / stats["reactions_with_ko"]
    stats = stats.fillna(0)

    # Add module names
    stats = pd.merge(stats, module_names, on="md", how="left")

    # Track inclusion/exclusion flags
    excluded_modules = excluded_modules or []
    rescue_modules = rescue_modules or []

    stats["has_species_version"] = stats["md"].isin(set(species_modules["md"]))
    stats["meets_coverage"] = stats["coverage"] >= min_coverage
    stats["is_rescued"] = stats["md"].isin(rescue_modules)
    stats["is_excluded"] = stats["md"].isin(excluded_modules)

    selected_mask = (
        stats["meets_coverage"]
        & (stats["has_species_version"] | stats["is_rescued"])
        & ~stats["is_excluded"]
    )
    selected = stats[selected_mask].copy()

    dropped_mask = stats["meets_coverage"] & ~selected_mask
    dropped = stats[dropped_mask].copy()

    def _drop_reason(row: pd.Series) -> str:
        reasons = []
        if row["is_excluded"]:
            reasons.append("explicitly excluded")
        if not row["has_species_version"] and not row["is_rescued"]:
            reasons.append("no species-specific KEGG submap")
        if not reasons:
            reasons.append("other")
        return "; ".join(reasons)

    if not dropped.empty:
        dropped["drop_reason"] = dropped.apply(_drop_reason, axis=1)
    else:
        dropped["drop_reason"] = pd.Series(dtype=str)

    return selected, dropped


def collect_final_reactions(
    species_reactions: set,
    module_stats: pd.DataFrame,
    rn_to_module: pd.DataFrame,
    species_pathways: pd.DataFrame,
    pathway_to_rn: pd.DataFrame,
    enzyme_reactions: set,
) -> set:
    """Collect the final set of reaction IDs to include."""
    # Reactions from selected modules
    module_ids = set(module_stats["md"])
    module_rns = set(
        rn_to_module[rn_to_module["md"].isin(module_ids)]["rn"]
    )

    # Reactions from species pathways (convert species code to map)
    species_pathways = species_pathways.copy()
    species_pathways["pathway"] = species_pathways["pathway"].str.replace(
        r"^[a-z]+", "map", regex=True
    )
    pathway_ids = set(species_pathways["pathway"])

    pathway_rns = set(
        pathway_to_rn[pathway_to_rn["path"].isin(pathway_ids)]["rn"]
    )

    # Include pathway reactions without enzymes
    pathway_no_enzyme = pathway_rns - enzyme_reactions

    # Combine: species enzymes + module reactions + pathway non-enzyme reactions
    return species_reactions | module_rns | pathway_no_enzyme


def run_pipeline(
    config: SpeciesConfig,
    data_dir: Path,
    output_dir: Path,
    metanetx_ttl: Path | None,
    timestamp: str | None = None,
    kegg_only: bool = False,
    smiles_source: str | Path | None = None,
    dg0_table: str | Path | None = None,
) -> PipelineResult:
    """
    Run the full pipeline for a species.
    
    Args:
        config: Species configuration
        data_dir: Directory for cached KEGG data
        output_dir: Directory for output files
        metanetx_ttl: Path to MetaNetX TTL file
        timestamp: Optional timestamp for output filename
    
    Returns:
        PipelineResult with graph (empty if kegg_only) and metrics
    """
    if not kegg_only and metanetx_ttl is None:
        raise ValueError("metanetx_ttl is required unless kegg_only=True")

    kegg_dir = data_dir / "kegg"
    kegg_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    species_output_dir = output_dir / config.species_code
    species_output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Building metabolic map for {config.name} ({config.species_code})")

    # Download/load KEGG data
    print("Loading KEGG data...")
    kegg_header = kegg.make_kegg_header()

    ko_to_rn = kegg.load_ko_to_reaction(kegg_dir, kegg_header)
    all_reactions = kegg.load_all_reactions(kegg_dir, kegg_header)
    gene_to_ko = kegg.load_gene_to_ko(kegg_dir, config.species_code, kegg_header)
    rn_to_module = kegg.load_reaction_to_module(kegg_dir, kegg_header)
    module_names = kegg.load_module_names(kegg_dir, kegg_header)
    species_modules = kegg.load_species_modules(kegg_dir, config.species_code, kegg_header)
    species_pathways = kegg.load_species_pathways(kegg_dir, config.species_code, kegg_header)
    pathway_to_rn = kegg.load_pathway_to_reaction(kegg_dir, kegg_header)
    pathway_names = kegg.load_pathway_names(kegg_dir, kegg_header)
    compound_names = kegg.load_compound_names(kegg_dir, kegg_header)
    gene_annot = kegg.load_gene_annotations(
        kegg_dir, config.species_code, kegg_header, config.ncbi_taxonomy_id
    )

    # Build reaction sets
    print("Building reaction sets...")
    rxn_sets = build_reaction_sets(ko_to_rn, gene_to_ko, all_reactions)

    # Filter modules
    print("Filtering modules...")
    module_stats, dropped_modules = filter_modules(
        rn_to_module=rn_to_module,
        module_names=module_names,
        species_modules=species_modules,
        species_reactions=rxn_sets["species"],
        enzyme_reactions=rxn_sets["with_enzyme"],
        min_coverage=config.min_enzyme_coverage,
        rescue_modules=config.rescue_modules,
        excluded_modules=config.excluded_modules,
    )
    print(f"  Selected: {len(module_stats)} modules")
    print(
        "  Dropped: "
        f"{len(dropped_modules)} modules (met coverage but excluded or lacked species submap)"
    )

    # Collect final reactions
    print("Collecting final reaction set...")
    final_reactions = collect_final_reactions(
        species_reactions=rxn_sets["species"],
        module_stats=module_stats,
        rn_to_module=rn_to_module,
        species_pathways=species_pathways,
        pathway_to_rn=pathway_to_rn,
        enzyme_reactions=rxn_sets["with_enzyme"],
    )

    if kegg_only:
        suffix = f".{timestamp}" if timestamp else ""
        dropped_path = (
            species_output_dir / f"{config.output_prefix}.dropped_modules{suffix}.csv"
        )
        retained_path = (
            species_output_dir / f"{config.output_prefix}.retained_modules{suffix}.csv"
        )
        dropped_modules.to_csv(dropped_path, index=False)
        module_stats.to_csv(retained_path, index=False)

        metrics = {
            "total_kegg_reactions": len(rxn_sets["all"]),
            "species_reactions": len(rxn_sets["species"]),
            "modules_selected": len(module_stats),
            "modules_dropped": len(dropped_modules),
            "final_reactions_input": len(final_reactions),
            "kegg_only": True,
        }
        print(f"Dropped modules written to {dropped_path}")
        print(f"Retained modules written to {retained_path}")
        return PipelineResult(
            graph=nx.DiGraph(), metrics=metrics, dropped_modules=dropped_modules
        )

    # Load MetaNetX and build graph
    print("Loading MetaNetX...")
    mnx_indices = metanetx.load_metanetx_indices(metanetx_ttl)

    print("Building reaction graph...")
    G = metanetx.build_kegg_digraph(final_reactions, mnx_indices)

    # Extract largest component and prune
    print("Extracting largest component...")
    G, num_components = graph.get_largest_component(G)

    print("Pruning dead-end metabolites...")
    G, removed_mets, removed_rxns, removed_gens, removed_mdls, removed_pwys = (
        graph.prune_dead_ends(G)
    )

    # Annotate compound names (KEGG common name = first entry)
    print("Adding compound names...")
    compound_names = compound_names.copy()
    compound_names["compound_name"] = compound_names["compound_name"].astype(str)
    compound_names["compound_name"] = (
        compound_names["compound_name"].str.split(";").str[0].str.strip()
    )
    compound_name_map = dict(
        zip(compound_names["compound"], compound_names["compound_name"], strict=False)
    )
    for node, data in G.nodes(data=True):
        if data.get("group") != "Compound":
            continue
        cid = data.get("name") or node.replace("kegg:", "", 1)
        common_name = compound_name_map.get(cid, "NA")
        if pd.isna(common_name) or common_name == "":
            common_name = "NA"
        G.nodes[node]["compound_name"] = common_name

    # Annotate with genes
    print("Adding gene annotations...")
    gene_to_rn = rxn_sets["gene_to_rn"].copy()

    # Merge gene annotations
    cols = ["ncbi_gene_id", "gene_symbol", "ensembl_id"]
    gene_agg = (
        gene_annot.groupby("kegg_gene_id", as_index=False)
        .agg({col: lambda s: "|".join(s.dropna().astype(str).unique()) for col in cols})
    )
    gene_to_rn = pd.merge(
        gene_to_rn,
        gene_agg,
        left_on="gene",
        right_on="kegg_gene_id",
        how="left",
    )
    graph.add_gene_annotations(G, gene_to_rn)

    # Annotate with modules
    print("Adding module annotations...")
    module_rns = pd.merge(
        rn_to_module[rn_to_module["md"].isin(module_stats["md"])],
        module_stats[["md", "module_name"]],
        on="md",
        how="left",
    )
    graph.add_module_annotations(G, module_rns)

    # Annotate with pathways
    print("Adding pathway annotations...")
    pathway_rns = pd.merge(
        pathway_to_rn,
        pathway_names,
        left_on="path",
        right_on="pathway",
        how="left",
    )
    pathway_rns = pathway_rns[~pathway_rns["path"].str.startswith("rn")]
    graph.add_pathway_annotations(G, pathway_rns)

    # Add orientation/SMILES annotations, save the full graph, then filter to the
    # SMILES-supported graph used by downstream embeddings.
    print("Loading joint KEGG SMILES...")
    smiles_records = chem.load_kegg_smiles_records(smiles_source, kegg_dir)
    print("Loading reaction DG0 annotations...")
    dg0_records = thermo.load_dg0_table(dg0_table)
    thermo.annotate_graph_dg0(G, dg0_records)
    print(f"  Loaded DG0 annotations for {len(dg0_records):,} reactions")
    suffix = f".{timestamp}" if timestamp else ""

    print("Annotating full graph with reaction orientation and SMILES...")
    reactions.build_oriented_reaction_table(G, smiles_map=smiles_records)

    full_graph_path = (
        species_output_dir / f"{config.output_prefix}.full{suffix}.graphml"
    )
    nx.write_graphml(G, full_graph_path)
    print(f"Full graph written to {full_graph_path}")

    print("Filtering graph to SMILES-supported reactions...")
    G, smiles_dropped_reactions = reactions.filter_smiles_supported_graph(G)
    smiles_initial_drop_count = len(smiles_dropped_reactions)

    print("Pruning dead ends after SMILES filtering...")
    smiles_prune_input_graph = G.copy()
    (
        G,
        smiles_removed_mets,
        smiles_removed_rxns,
        smiles_removed_gens,
        smiles_removed_mdls,
        smiles_removed_pwys,
    ) = graph.prune_dead_ends(G)

    if smiles_removed_rxns:
        extra_rxn_rows = []
        for node in sorted(smiles_removed_rxns):
            node_data = (
                smiles_prune_input_graph.nodes[node]
                if node in smiles_prune_input_graph
                else {}
            )
            status = (
                reactions.reaction_side_status(smiles_prune_input_graph, node)
                if node in smiles_prune_input_graph
                else {
                    "input_compounds": [],
                    "output_compounds": [],
                    "input_smiles_present": [],
                    "output_smiles_present": [],
                    "missing_input_compounds": [],
                    "missing_output_compounds": [],
                }
            )
            extra_rxn_rows.append(
                {
                    "reaction_node": node,
                    "reaction_id": node_data.get("name") or node.replace("kegg:", "", 1),
                    "drop_reason": "dead_end_after_smiles_filter",
                    "input_compounds": reactions.join_compounds(
                        status["input_compounds"]
                    ),
                    "output_compounds": reactions.join_compounds(
                        status["output_compounds"]
                    ),
                    "input_smiles_present": reactions.join_compounds(
                        status["input_smiles_present"]
                    ),
                    "output_smiles_present": reactions.join_compounds(
                        status["output_smiles_present"]
                    ),
                    "missing_input_compounds": reactions.join_compounds(
                        status["missing_input_compounds"]
                    ),
                    "missing_output_compounds": reactions.join_compounds(
                        status["missing_output_compounds"]
                    ),
                }
            )
        smiles_dropped_reactions = pd.concat(
            [smiles_dropped_reactions, pd.DataFrame(extra_rxn_rows)],
            ignore_index=True,
        )

    smiles_dropped_compounds = pd.DataFrame(
        [
            {
                "compound_node": node,
                "compound_id": (
                    smiles_prune_input_graph.nodes[node].get("name")
                    if node in smiles_prune_input_graph
                    else node.replace("kegg:", "", 1)
                ),
                "drop_reason": "dead_end_after_smiles_filter",
            }
            for node in sorted(smiles_removed_mets)
        ],
        columns=["compound_node", "compound_id", "drop_reason"],
    )

    dropped_reactions_path = (
        species_output_dir
        / f"{config.output_prefix}.smiles_dropped_reactions{suffix}.tsv"
    )
    dropped_compounds_path = (
        species_output_dir
        / f"{config.output_prefix}.smiles_dropped_compounds{suffix}.tsv"
    )
    smiles_dropped_reactions.to_csv(dropped_reactions_path, sep="\t", index=False)
    smiles_dropped_compounds.to_csv(dropped_compounds_path, sep="\t", index=False)
    print(f"SMILES dropped reactions written to {dropped_reactions_path}")
    print(f"SMILES dropped compounds written to {dropped_compounds_path}")

    reaction_table_path = (
        species_output_dir / f"{config.output_prefix}.reactions{suffix}.tsv"
    )
    reaction_table = reactions.build_oriented_reaction_table(
        G,
        smiles_map=smiles_records,
        output_path=reaction_table_path,
    )
    print(f"Reaction table written to {reaction_table_path}")

    # Save output
    if timestamp:
        filename = f"{config.output_prefix}.{timestamp}.graphml"
    else:
        filename = f"{config.output_prefix}.graphml"

    output_path = species_output_dir / filename
    nx.write_graphml(G, output_path)
    print(f"Graph written to {output_path}")

    # Collect metrics
    counts = graph.count_by_group(G)
    metrics = {
        "total_kegg_reactions": len(rxn_sets["all"]),
        "species_reactions": len(rxn_sets["species"]),
        "modules_selected": len(module_stats),
        "modules_dropped": len(dropped_modules),
        "final_reactions_input": len(final_reactions),
        "num_components": num_components,
        "pruned_metabolites": len(removed_mets),
        "pruned_reactions": len(removed_rxns),
        "pruned_genes": len(removed_gens),
        "pruned_modules": len(removed_mdls),
        "pruned_pathways": len(removed_pwys),
        "smiles_filter_reactions_dropped_missing_side": smiles_initial_drop_count,
        "smiles_filter_reactions_dropped_dead_end": len(smiles_removed_rxns),
        "smiles_filter_compounds_dropped_dead_end": len(smiles_removed_mets),
        "smiles_filter_genes_dropped_orphan": len(smiles_removed_gens),
        "smiles_filter_modules_dropped_orphan": len(smiles_removed_mdls),
        "smiles_filter_pathways_dropped_orphan": len(smiles_removed_pwys),
        "reaction_table_rows": len(reaction_table),
        "reaction_table_complete_smiles": int(
            (reaction_table["missing_smiles_compounds"] == "").sum()
        ),
        **{f"nodes_{k.lower()}": v for k, v in counts.items()},
        "edges": G.number_of_edges(),
    }

    print("\n" + graph.graph_summary(G))

    return PipelineResult(graph=G, metrics=metrics, dropped_modules=dropped_modules)


def generate_report(
    config: SpeciesConfig,
    result: PipelineResult,
    output_path: Path,
) -> None:
    """Generate a PDF summary report."""
    counts = graph.count_by_group(result.graph)
    m = result.metrics

    md_content = f"""
# Metabolic Map: {config.name}

**Species code:** `{config.species_code}`  
**Enzyme coverage threshold:** {config.min_enzyme_coverage:.0%}  
**Rescued modules:** {len(config.rescue_modules)}

---

## Summary Statistics

### Input Data

| Metric | Count |
|--------|-------|
| Total KEGG reactions | {m.get('total_kegg_reactions', 'N/A'):,} |
| Species enzyme-backed reactions | {m.get('species_reactions', 'N/A'):,} |
| Modules selected | {m.get('modules_selected', 'N/A'):,} |
| Reactions before graph construction | {m.get('final_reactions_input', 'N/A'):,} |

### Final Graph

| Node Type | Count |
|-----------|-------|
| Reactions | {counts.get('Reaction', 0):,} |
| Compounds | {counts.get('Compound', 0):,} |
| Genes | {counts.get('Gene', 0):,} |
| Modules | {counts.get('Module', 0):,} |
| Pathways | {counts.get('Pathway', 0):,} |

**Total edges:** {m.get('edges', 0):,}

---

## Graph Schema

### Node Types

Each node has a `group` attribute: `Compound`, `Reaction`, `Gene`, `Module`, or `Pathway`.

### Node Identifiers

| Group | ID Format | Example |
|-------|-----------|---------|
| Compound | `kegg:CXXXXX` | `kegg:C00001` |
| Reaction | `kegg:RXXXXX` | `kegg:R00286` |
| Module | `kegg:MXXXXX` | `kegg:M00001` |
| Pathway | `kegg:mapXXXXX` | `kegg:map00010` |
| Gene | `kegg:{config.species_code}:XXXXXXXX` | `kegg:{config.species_code}:100037283` |

### Additional Node Attributes

| Group | Attributes |
|-------|------------|
| Compound | `compound_name`, `smiles`, `smiles_source` |
| Reaction | `raw_reaction`, `dg0`, `dg0_reaction`, `delta_g_status`, `balance_status`, `balance_compound`, `orientation`, `orientation_source`, `thermo_agrees_with_mnx`, `orientation_flip_from_mnx` |
| Gene | `gene_symbol`, `ensembl_id`, `ncbi_gene_id` |
| Module | `module_name` |
| Pathway | `pathway_name` |

### Edge Attributes

- **Compound ↔ Reaction:** `coef` (stoichiometric coefficient), `compartment`
- **Gene/Module/Pathway → Reaction:** no additional attributes

---

## Rescued Modules

The following modules were manually included despite lacking species-specific KEGG submaps:

"""
    for mod in config.rescue_modules:
        md_content += f"- `{mod}`\n"

    if not result.dropped_modules.empty:
        md_content += f"""

---

## Dropped Modules

The following {len(result.dropped_modules)} modules met the coverage threshold but were not included:

| Module | Coverage | Reactions | Reason | Name |
|--------|----------|-----------|--------|------|
"""
        for _, row in result.dropped_modules.iterrows():
            coverage_pct = (
                f"{row['coverage']:.0%}" if pd.notna(row.get("coverage")) else "N/A"
            )
            name = row.get("module_name", "N/A") or "N/A"
            if len(name) > 40:
                name = name[:37] + "..."
            reason = row.get("drop_reason", "N/A")
            md_content += (
                f"| `{row['md']}` | {coverage_pct} | "
                f"{int(row['reactions'])} | {reason} | {name} |\n"
            )

    html_body = markdown(md_content, extensions=["tables"])

    css = CSS(string="""
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            font-size: 11pt;
            line-height: 1.5;
            max-width: 800px;
            margin: 40px auto;
            padding: 0 20px;
            color: #333;
        }
        h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
        h2 { color: #34495e; margin-top: 30px; }
        h3 { color: #7f8c8d; }
        table { border-collapse: collapse; width: 100%; margin: 15px 0; }
        th, td { border: 1px solid #ddd; padding: 8px 12px; text-align: left; }
        th { background-color: #f8f9fa; font-weight: 600; }
        tr:nth-child(even) { background-color: #f8f9fa; }
        code { background-color: #f4f4f4; padding: 2px 6px; border-radius: 3px; font-size: 10pt; }
        hr { border: none; border-top: 1px solid #eee; margin: 30px 0; }
    """)

    html_full = f"<html><body>{html_body}</body></html>"
    HTML(string=html_full).write_pdf(output_path, stylesheets=[css])
    print(f"Report written to {output_path}")
