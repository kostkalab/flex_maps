"""Graph manipulation and annotation utilities."""

import networkx as nx
import pandas as pd


def get_largest_component(G: nx.DiGraph) -> tuple[nx.DiGraph, int]:
    """Return the largest weakly connected component and total component count."""
    ccs = list(nx.weakly_connected_components(G))
    largest = max(ccs, key=len)
    return G.subgraph(largest).copy(), len(ccs)


def prune_dead_ends(
    G: nx.DiGraph, reversible: bool = False
) -> tuple[nx.DiGraph, list[str], list[str]]:
    """
    Recursively remove dead-end metabolites and their associated reactions.
    
    Dead-end metabolites are compounds with no producers or no consumers.
    If reversible=True, treats edges as undirected (degree <= 1).
    
    Returns:
        Tuple of (pruned_graph, removed_metabolites, removed_reactions)
    """
    g = G.copy()
    removed_mets = []
    removed_rxns = []

    while True:
        if reversible:
            dead_mets = [
                m for m, data in g.nodes(data=True)
                if data.get("group") == "Compound" and g.degree(m) <= 1
            ]
        else:
            dead_mets = [
                m for m, data in g.nodes(data=True)
                if data.get("group") == "Compound"
                and (g.in_degree(m) == 0 or g.out_degree(m) == 0)
            ]

        if not dead_mets:
            break

        dead_rxns = set()
        for m in dead_mets:
            for nbr in list(g.predecessors(m)) + list(g.successors(m)):
                if g.nodes[nbr].get("group") == "Reaction":
                    dead_rxns.add(nbr)

        g.remove_nodes_from(dead_mets)
        g.remove_nodes_from(dead_rxns)
        removed_mets.extend(dead_mets)
        removed_rxns.extend(dead_rxns)

    return g, removed_mets, removed_rxns


def add_gene_annotations(
    G: nx.DiGraph,
    gene_to_reaction: pd.DataFrame,
) -> None:
    """
    Add gene nodes and edges to reaction nodes in-place.
    
    Expects gene_to_reaction to have columns:
        gene, rn, gene_symbol, ensembl_id, ncbi_gene_id
    """
    for _, row in gene_to_reaction.iterrows():
        rn = row["rn"]
        gene = row["gene"]
        symbol = row.get("gene_symbol") or "NA"
        ensid = row.get("ensembl_id") or "NA"
        ncbiid = row.get("ncbi_gene_id") or "NA"

        # Handle NaN values
        if pd.isna(symbol):
            symbol = "NA"
        if pd.isna(ensid):
            ensid = "NA"
        if pd.isna(ncbiid):
            ncbiid = "NA"

        rxn_node = f"kegg:{rn}"
        gene_node = f"kegg:{gene}"

        if rxn_node in G:
            if gene_node not in G:
                G.add_node(
                    gene_node,
                    group="Gene",
                    name=gene,
                    gene_symbol=symbol,
                    ensembl_id=ensid,
                    ncbi_gene_id=ncbiid,
                )
            G.add_edge(gene_node, rxn_node)


def add_module_annotations(
    G: nx.DiGraph,
    module_reactions: pd.DataFrame,
) -> None:
    """
    Add module nodes and edges to reaction nodes in-place.
    
    Expects module_reactions to have columns: md, rn, module_name
    """
    for _, row in module_reactions.iterrows():
        rn = row["rn"]
        md = row["md"]
        module_name = row.get("module_name") or "NA"

        if pd.isna(module_name):
            module_name = "NA"

        rxn_node = f"kegg:{rn}"
        module_node = f"kegg:{md}"

        if rxn_node in G:
            if module_node not in G:
                G.add_node(module_node, group="Module", name=md, module_name=module_name)
            G.add_edge(module_node, rxn_node)


def add_pathway_annotations(
    G: nx.DiGraph,
    pathway_reactions: pd.DataFrame,
) -> None:
    """
    Add pathway nodes and edges to reaction nodes in-place.
    
    Expects pathway_reactions to have columns: path, rn, pathway_name
    """
    for _, row in pathway_reactions.iterrows():
        rn = row["rn"]
        path = row["path"]
        pathway_name = row.get("pathway_name") or "NA"

        if pd.isna(pathway_name):
            pathway_name = "NA"

        rxn_node = f"kegg:{rn}"
        pathway_node = f"kegg:{path}"

        if rxn_node in G:
            if pathway_node not in G:
                G.add_node(
                    pathway_node, group="Pathway", name=path, pathway_name=pathway_name
                )
            G.add_edge(pathway_node, rxn_node)


def count_by_group(G: nx.DiGraph) -> dict[str, int]:
    """Count nodes by group attribute."""
    counts = {}
    for _, data in G.nodes(data=True):
        group = data.get("group", "Unknown")
        counts[group] = counts.get(group, 0) + 1
    return counts


def graph_summary(G: nx.DiGraph) -> str:
    """Return a formatted summary of graph contents."""
    counts = count_by_group(G)
    lines = [
        f"Nodes: {G.number_of_nodes()}",
        f"Edges: {G.number_of_edges()}",
        "By type:",
    ]
    for group in ["Reaction", "Compound", "Gene", "Module", "Pathway"]:
        if group in counts:
            lines.append(f"  {group}: {counts[group]}")
    return "\n".join(lines)
