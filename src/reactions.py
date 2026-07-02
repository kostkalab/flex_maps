"""Oriented reaction table construction and graph annotations."""

from __future__ import annotations

from pathlib import Path

import networkx as nx
import pandas as pd

from . import chem


REACTION_TABLE_COLUMNS = [
    "reaction_node",
    "reaction_id",
    "raw_reaction",
    "mnx_left_compounds",
    "mnx_right_compounds",
    "input_compounds",
    "output_compounds",
    "input_smiles",
    "output_smiles",
    "I_to_O_string",
    "O_to_I_string",
    "dg0",
    "dg0_reaction",
    "delta_g_status",
    "balance_status",
    "balance_compound",
    "orientation",
    "orientation_source",
    "thermo_agrees_with_mnx",
    "orientation_flip_from_mnx",
    "missing_smiles_compounds",
]


REACTION_NODE_DEFAULTS = {
    "dg0": "",
    "dg0_reaction": "",
    "delta_g_status": "missing",
    "balance_status": "missing",
    "balance_compound": "",
    "orientation": "unknown",
    "orientation_source": "mnx_fallback",
    "thermo_agrees_with_mnx": "",
    "orientation_flip_from_mnx": False,
}


def _compound_id(node: str, data: dict) -> str:
    return str(data.get("name") or node.replace("kegg:", "", 1))


def _reaction_id(node: str, data: dict) -> str:
    return str(data.get("name") or node.replace("kegg:", "", 1))


def _side_compounds(G: nx.DiGraph, reaction_node: str, predecessors: bool) -> list[str]:
    nodes = G.predecessors(reaction_node) if predecessors else G.successors(reaction_node)
    compounds = []
    for node in nodes:
        data = G.nodes[node]
        if data.get("group") == "Compound":
            compounds.append(_compound_id(node, data))
    return sorted(set(compounds))


def join_compounds(compounds: list[str]) -> str:
    """Serialize a compound side for TSV/GraphML attributes."""
    return "|".join(sorted(compounds))


def _smiles_value(record: object) -> str:
    if isinstance(record, dict):
        return str(record.get("smiles", "") or "")
    return str(record or "")


def _smiles_records(smiles_data: dict) -> dict[str, dict]:
    records = {}
    for cid, record in smiles_data.items():
        if isinstance(record, dict):
            records[cid] = record
        else:
            records[cid] = {
                "smiles": str(record),
                "smiles_source": "unknown",
                "smiles_has_wildcard": "true" if "*" in str(record) else "false",
                "smiles_rdkit_parse_ok": "",
                "smiles_selection_reason": "",
            }
    return records


def annotate_compound_smiles(G: nx.DiGraph, smiles_data: dict) -> None:
    """Add ``smiles`` and ``smiles_source`` attributes to compound nodes."""
    smiles_records = _smiles_records(smiles_data)
    for node, data in G.nodes(data=True):
        if data.get("group") != "Compound":
            continue
        cid = _compound_id(node, data)
        record = smiles_records.get(cid)
        if record:
            data["smiles"] = _smiles_value(record)
            data["smiles_source"] = str(record.get("smiles_source", ""))
            data["smiles_has_wildcard"] = str(
                record.get("smiles_has_wildcard", "")
            )
            data["smiles_rdkit_parse_ok"] = str(
                record.get("smiles_rdkit_parse_ok", "")
            )
            data["smiles_selection_reason"] = str(
                record.get("smiles_selection_reason", "")
            )
        else:
            data["smiles"] = ""
            data["smiles_source"] = ""
            data["smiles_has_wildcard"] = ""
            data["smiles_rdkit_parse_ok"] = ""
            data["smiles_selection_reason"] = ""


def graph_smiles_map(G: nx.DiGraph) -> dict[str, str]:
    """Return compound SMILES already attached to the graph by upstream loaders."""
    values = {}
    for node, data in G.nodes(data=True):
        if data.get("group") != "Compound":
            continue
        smiles = data.get("smiles")
        if smiles:
            values[_compound_id(node, data)] = str(smiles)
    return values


def reaction_side_status(
    G: nx.DiGraph,
    reaction_node: str,
) -> dict:
    """Return input/output SMILES availability for an already oriented reaction."""
    inputs = _side_compounds(G, reaction_node, predecessors=True)
    outputs = _side_compounds(G, reaction_node, predecessors=False)
    smiles_map = graph_smiles_map(G)
    input_present = [cid for cid in inputs if smiles_map.get(cid)]
    output_present = [cid for cid in outputs if smiles_map.get(cid)]
    missing_input = [cid for cid in inputs if not smiles_map.get(cid)]
    missing_output = [cid for cid in outputs if not smiles_map.get(cid)]
    return {
        "input_compounds": inputs,
        "output_compounds": outputs,
        "input_smiles_present": input_present,
        "output_smiles_present": output_present,
        "missing_input_compounds": missing_input,
        "missing_output_compounds": missing_output,
        "all_inputs_missing": bool(inputs) and not input_present,
        "all_outputs_missing": bool(outputs) and not output_present,
    }


def find_smiles_unsupported_reactions(G: nx.DiGraph) -> pd.DataFrame:
    """Find reactions where a complete input or output side lacks SMILES."""
    columns = [
        "reaction_node",
        "reaction_id",
        "drop_reason",
        "input_compounds",
        "output_compounds",
        "input_smiles_present",
        "output_smiles_present",
        "missing_input_compounds",
        "missing_output_compounds",
    ]
    rows = []
    for node, data in sorted(G.nodes(data=True)):
        if data.get("group") != "Reaction":
            continue
        status = reaction_side_status(G, node)
        if not status["all_inputs_missing"] and not status["all_outputs_missing"]:
            continue
        if status["all_inputs_missing"] and status["all_outputs_missing"]:
            reason = "all_inputs_and_outputs_missing_smiles"
        elif status["all_inputs_missing"]:
            reason = "all_inputs_missing_smiles"
        else:
            reason = "all_outputs_missing_smiles"
        rows.append(
            {
                "reaction_node": node,
                "reaction_id": _reaction_id(node, data),
                "drop_reason": reason,
                "input_compounds": join_compounds(status["input_compounds"]),
                "output_compounds": join_compounds(status["output_compounds"]),
                "input_smiles_present": join_compounds(status["input_smiles_present"]),
                "output_smiles_present": join_compounds(
                    status["output_smiles_present"]
                ),
                "missing_input_compounds": join_compounds(
                    status["missing_input_compounds"]
                ),
                "missing_output_compounds": join_compounds(
                    status["missing_output_compounds"]
                ),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def filter_smiles_supported_graph(
    G: nx.DiGraph,
) -> tuple[nx.DiGraph, pd.DataFrame]:
    """Drop reactions where one full side has no SMILES support."""
    dropped = find_smiles_unsupported_reactions(G)
    filtered = G.copy()
    if not dropped.empty:
        filtered.remove_nodes_from(dropped["reaction_node"].tolist())
    return filtered, dropped


def _orient_sides(
    left_compounds: list[str],
    right_compounds: list[str],
    dg0: object,
) -> tuple[list[str], list[str], dict]:
    if dg0 in {"", None} or pd.isna(dg0):
        return left_compounds, right_compounds, {
            "orientation": "unknown",
            "orientation_source": "mnx_fallback",
            "thermo_agrees_with_mnx": "",
            "orientation_flip_from_mnx": False,
        }

    dg0_value = float(dg0)
    if dg0_value > 0:
        return right_compounds, left_compounds, {
            "orientation": "reverse",
            "orientation_source": "dg0",
            "thermo_agrees_with_mnx": False,
            "orientation_flip_from_mnx": True,
        }

    return left_compounds, right_compounds, {
        "orientation": "forward",
        "orientation_source": "dg0",
        "thermo_agrees_with_mnx": True,
        "orientation_flip_from_mnx": False,
    }


def orient_reaction_edges(G: nx.DiGraph, reaction_node: str, inputs: list[str], outputs: list[str]) -> None:
    """Physically orient compound-reaction edges as input -> reaction -> output."""
    input_nodes = {f"kegg:{cid}" for cid in inputs}
    output_nodes = {f"kegg:{cid}" for cid in outputs}
    compound_nodes = input_nodes | output_nodes

    existing_attrs = {}
    for node in compound_nodes:
        if G.has_edge(node, reaction_node):
            existing_attrs[node] = dict(G.edges[node, reaction_node])
        if G.has_edge(reaction_node, node):
            existing_attrs.setdefault(node, dict(G.edges[reaction_node, node]))

    for node in compound_nodes:
        if G.has_edge(node, reaction_node):
            G.remove_edge(node, reaction_node)
        if G.has_edge(reaction_node, node):
            G.remove_edge(reaction_node, node)

    for node in input_nodes:
        if node in G:
            G.add_edge(node, reaction_node, **existing_attrs.get(node, {}))
    for node in output_nodes:
        if node in G:
            G.add_edge(reaction_node, node, **existing_attrs.get(node, {}))


def build_oriented_reaction_table(
    G: nx.DiGraph,
    smiles_map: dict,
    output_path: Path | None = None,
) -> pd.DataFrame:
    """
    Annotate graph nodes/edges and return the oriented reaction sidecar table.

    Phase 1/2 proof of concept uses missing DG0 fallback orientation. Later DG0
    loaders can populate reaction node DG0 fields before this function is called.
    """
    annotate_compound_smiles(G, smiles_map)
    normalized_smiles_map = {
        cid: _smiles_value(record) for cid, record in smiles_map.items()
    }
    reaction_smiles_map = {**normalized_smiles_map, **graph_smiles_map(G)}
    rows = []

    reaction_nodes = sorted(
        node for node, data in G.nodes(data=True) if data.get("group") == "Reaction"
    )
    for node in reaction_nodes:
        data = G.nodes[node]
        for key, value in REACTION_NODE_DEFAULTS.items():
            data.setdefault(key, value)

        reaction_id = _reaction_id(node, data)
        left_compounds = _side_compounds(G, node, predecessors=True)
        right_compounds = _side_compounds(G, node, predecessors=False)
        inputs, outputs, orientation_data = _orient_sides(
            left_compounds, right_compounds, data.get("dg0")
        )

        data.update(orientation_data)
        orient_reaction_edges(G, node, inputs, outputs)

        input_smiles, missing_input = chem.reaction_side_smiles(
            inputs, reaction_smiles_map
        )
        output_smiles, missing_output = chem.reaction_side_smiles(
            outputs, reaction_smiles_map
        )
        missing_smiles = sorted(set(missing_input + missing_output))

        row = {
            "reaction_node": node,
            "reaction_id": reaction_id,
            "raw_reaction": data.get("raw_reaction", ""),
            "mnx_left_compounds": join_compounds(left_compounds),
            "mnx_right_compounds": join_compounds(right_compounds),
            "input_compounds": join_compounds(inputs),
            "output_compounds": join_compounds(outputs),
            "input_smiles": input_smiles,
            "output_smiles": output_smiles,
            "I_to_O_string": f"{input_smiles}>>{output_smiles}",
            "O_to_I_string": f"{output_smiles}>>{input_smiles}",
            "dg0": data.get("dg0", ""),
            "dg0_reaction": data.get("dg0_reaction", ""),
            "delta_g_status": data.get("delta_g_status", "missing"),
            "balance_status": data.get("balance_status", "missing"),
            "balance_compound": data.get("balance_compound", ""),
            "orientation": data["orientation"],
            "orientation_source": data["orientation_source"],
            "thermo_agrees_with_mnx": data["thermo_agrees_with_mnx"],
            "orientation_flip_from_mnx": data["orientation_flip_from_mnx"],
            "missing_smiles_compounds": join_compounds(missing_smiles),
        }
        rows.append(row)

    table = pd.DataFrame(rows, columns=REACTION_TABLE_COLUMNS)
    if output_path is not None:
        table.to_csv(output_path, sep="\t", index=False)
    return table
