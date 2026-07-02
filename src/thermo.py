"""DG0 calculation helpers for graph-derived reactions."""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
from typing import Iterable

import networkx as nx
import pandas as pd


DG0_COLUMNS = [
    "reaction_node",
    "reaction_id",
    "source_graphs",
    "dg0",
    "dg0_std",
    "dg0_reaction",
    "delta_g_status",
    "balance_status",
    "balance_compound",
    "error_message",
]

H_H2O_COMPOUNDS = {"C00001", "C00080"}
INORGANIC_BALANCING_COMPOUNDS = ["C00009", "C00011"]
CURRENCY_BALANCING_COMPOUNDS = [
    "C00010",  # CoA
    "C00002",  # ATP
    "C00008",  # ADP
    "C00020",  # AMP
    "C00003",  # NAD+
    "C00004",  # NADH
    "C00006",  # NADP+
    "C00005",  # NADPH
    "C00016",  # FAD
    "C01352",  # FADH2
    "C00035",  # GDP
    "C00044",  # GTP
]
BALANCING_COMPOUNDS = [
    *sorted(H_H2O_COMPOUNDS),
    *INORGANIC_BALANCING_COMPOUNDS,
    *CURRENCY_BALANCING_COMPOUNDS,
]


@dataclass
class ReactionFormula:
    """Graph-derived reaction formula in current left-to-right graph direction."""

    reaction_node: str
    reaction_id: str
    left: list[tuple[str, float]]
    right: list[tuple[str, float]]

    def to_equilibrator_formula(self) -> str:
        return f"{format_side(self.left)} <=> {format_side(self.right)}"


def compound_id(node: str, data: dict) -> str:
    return str(data.get("name") or node.replace("kegg:", "", 1))


def reaction_id(node: str, data: dict) -> str:
    return str(data.get("name") or node.replace("kegg:", "", 1))


def _coef(value: object) -> float:
    if value in {"", None}:
        return 1.0
    return float(value)


def _compound_side(
    G: nx.DiGraph,
    reaction_node: str,
    incoming: bool,
) -> list[tuple[str, float]]:
    edges = G.in_edges(reaction_node, data=True) if incoming else G.out_edges(
        reaction_node, data=True
    )
    side = []
    for edge in edges:
        compound_node = edge[0] if incoming else edge[1]
        edge_data = edge[2]
        node_data = G.nodes[compound_node]
        if node_data.get("group") != "Compound":
            continue
        side.append((compound_id(compound_node, node_data), _coef(edge_data.get("coef"))))
    return sorted(side)


def graph_reaction_formula(G: nx.DiGraph, reaction_node: str) -> ReactionFormula:
    data = G.nodes[reaction_node]
    return ReactionFormula(
        reaction_node=reaction_node,
        reaction_id=reaction_id(reaction_node, data),
        left=_compound_side(G, reaction_node, incoming=True),
        right=_compound_side(G, reaction_node, incoming=False),
    )


def format_side(side: Iterable[tuple[str, float]]) -> str:
    parts = []
    for cid, coef in sorted(side):
        if coef == 1:
            parts.append(cid)
        else:
            parts.append(f"{coef:g} {cid}")
    return " + ".join(parts)


def balance_status_for_compound(compound_id_value: str) -> str:
    if compound_id_value in H_H2O_COMPOUNDS:
        return "balanced_H_H2O"
    if compound_id_value in INORGANIC_BALANCING_COMPOUNDS:
        return "balanced_inorganic"
    return "balanced_currency"


def _component_contribution():
    cache_home = Path("data/equilibrator").resolve()
    cache_home.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_home))
    try:
        from equilibrator_api import ComponentContribution
    except ImportError as exc:
        raise RuntimeError(
            "DG0 calculation requires the optional equilibrator-api package. "
            "Install it in the DG0 calculation environment before running "
            "build_reaction_dg0.py."
        ) from exc
    return ComponentContribution()


def calculate_dg0_record(
    formula: ReactionFormula,
    component_contribution=None,
) -> dict:
    """Calculate DG0 for one graph-derived reaction formula."""
    cc = component_contribution or _component_contribution()
    dg0_reaction = formula.to_equilibrator_formula()
    base = {
        "reaction_node": formula.reaction_node,
        "reaction_id": formula.reaction_id,
        "source_graphs": "",
        "dg0": "",
        "dg0_std": "",
        "dg0_reaction": dg0_reaction,
        "delta_g_status": "failed",
        "balance_status": "failed",
        "balance_compound": "",
        "error_message": "",
    }

    try:
        reaction = cc.parse_reaction_formula(dg0_reaction)
        balance_status = "already_balanced"
        balance_compound = ""

        if not reaction.is_balanced():
            for cid in BALANCING_COMPOUNDS:
                try:
                    candidate = reaction.copy()
                    candidate.balance_with_compound(cc.get_compound(cid))
                except Exception:
                    continue
                if candidate.is_balanced():
                    reaction = candidate
                    balance_status = balance_status_for_compound(cid)
                    balance_compound = cid
                    break

        if not reaction.is_balanced():
            try:
                candidate = cc.balance_by_oxidation(reaction)
            except Exception:
                candidate = None
            if candidate is not None and candidate.is_balanced():
                reaction = candidate
                balance_status = "oxidation_balance"

        if not reaction.is_balanced():
            base["error_message"] = "reaction_not_balanced"
            return base

        result = cc.dg_prime(reaction)
        magnitude = result.magnitude
        base.update(
            {
                "dg0": magnitude.nominal_value,
                "dg0_std": magnitude.std_dev,
                "delta_g_status": "equilibrator_calculated",
                "balance_status": balance_status,
                "balance_compound": balance_compound,
            }
        )
        return base
    except Exception as exc:
        base["error_message"] = str(exc)
        return base


def calculate_graph_dg0_table(
    G: nx.DiGraph,
    limit: int | None = None,
    reaction_ids: set[str] | None = None,
) -> pd.DataFrame:
    """Calculate DG0 rows for reactions in a graph."""
    cc = _component_contribution()
    reaction_nodes = sorted(
        node for node, data in G.nodes(data=True) if data.get("group") == "Reaction"
    )
    if reaction_ids is not None:
        reaction_nodes = [
            node
            for node in reaction_nodes
            if reaction_id(node, G.nodes[node]) in reaction_ids
        ]
    if limit is not None:
        reaction_nodes = reaction_nodes[:limit]
    rows = [
        calculate_dg0_record(graph_reaction_formula(G, node), cc)
        for node in reaction_nodes
    ]
    return pd.DataFrame(rows, columns=DG0_COLUMNS)


def graph_reaction_formulas(
    G: nx.DiGraph,
    source_graph: str = "",
) -> dict[str, tuple[ReactionFormula, set[str]]]:
    """Collect graph-derived reaction formulas keyed by KEGG reaction ID."""
    records = {}
    for node, data in G.nodes(data=True):
        if data.get("group") != "Reaction":
            continue
        formula = graph_reaction_formula(G, node)
        records[formula.reaction_id] = (formula, {source_graph} if source_graph else set())
    return records


def merge_graph_reaction_formulas(
    formula_maps: Iterable[dict[str, tuple[ReactionFormula, set[str]]]],
) -> tuple[dict[str, tuple[ReactionFormula, set[str]]], dict[str, set[str]]]:
    """Merge formulas and report conflicting formulas for the same reaction ID."""
    merged: dict[str, tuple[ReactionFormula, set[str]]] = {}
    conflicts: dict[str, set[str]] = {}
    for formula_map in formula_maps:
        for rid, (formula, sources) in formula_map.items():
            formula_text = formula.to_equilibrator_formula()
            if rid not in merged:
                merged[rid] = (formula, set(sources))
                continue
            existing, existing_sources = merged[rid]
            existing_text = existing.to_equilibrator_formula()
            existing_sources.update(sources)
            if existing_text != formula_text:
                conflicts.setdefault(rid, set()).update({existing_text, formula_text})
    return merged, conflicts


def calculate_formula_records_dg0_table(
    formulas: dict[str, tuple[ReactionFormula, set[str]]],
    limit: int | None = None,
) -> pd.DataFrame:
    """Calculate DG0 rows for pre-collected unique reaction formulas."""
    cc = _component_contribution()
    items = sorted(formulas.items())
    if limit is not None:
        items = items[:limit]
    rows = []
    try:
        from tqdm import tqdm
    except ImportError:
        iterator = items
    else:
        iterator = tqdm(items, desc="Calculating DG0")
    for _, (formula, sources) in iterator:
        row = calculate_dg0_record(formula, cc)
        row["source_graphs"] = "|".join(sorted(s for s in sources if s))
        rows.append(row)
    return pd.DataFrame(rows, columns=DG0_COLUMNS)


def calculate_formula_records_dg0_rows(
    formulas: dict[str, tuple[ReactionFormula, set[str]]],
    checkpoint_path: str | Path | None = None,
    checkpoint_every: int = 25,
    initial_rows: list[dict] | None = None,
) -> list[dict]:
    """Calculate DG0 rows with optional periodic TSV checkpoints."""
    cc = _component_contribution()
    items = sorted(formulas.items())
    rows = list(initial_rows or [])

    try:
        from tqdm import tqdm
    except ImportError:
        iterator = items
    else:
        iterator = tqdm(items, desc="Calculating DG0")

    checkpoint = Path(checkpoint_path) if checkpoint_path is not None else None
    if checkpoint is not None:
        checkpoint.parent.mkdir(parents=True, exist_ok=True)

    for index, (_, (formula, sources)) in enumerate(iterator, start=1):
        row = calculate_dg0_record(formula, cc)
        row["source_graphs"] = "|".join(sorted(s for s in sources if s))
        rows.append(row)
        if checkpoint is not None and checkpoint_every > 0 and index % checkpoint_every == 0:
            pd.DataFrame(rows, columns=DG0_COLUMNS).to_csv(
                checkpoint,
                sep="\t",
                index=False,
            )

    if checkpoint is not None:
        pd.DataFrame(rows, columns=DG0_COLUMNS).to_csv(
            checkpoint,
            sep="\t",
            index=False,
        )
    return rows


def load_dg0_table(path: str | Path | None) -> dict[str, dict]:
    """Load DG0 rows keyed by KEGG reaction ID."""
    if path is None:
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    records = {}
    for _, row in df.iterrows():
        rid = str(row.get("reaction_id", ""))
        if not rid:
            continue
        records[rid] = {col: row.get(col, "") for col in df.columns}
    return records


def annotate_graph_dg0(G: nx.DiGraph, dg0_records: dict[str, dict]) -> None:
    """Attach DG0 fields to reaction nodes from a ``reaction_id`` keyed table."""
    for node, data in G.nodes(data=True):
        if data.get("group") != "Reaction":
            continue
        rid = reaction_id(node, data)
        record = dg0_records.get(rid)
        if not record:
            continue
        data["dg0"] = str(record.get("dg0", ""))
        data["dg0_reaction"] = str(record.get("dg0_reaction", ""))
        data["delta_g_status"] = str(record.get("delta_g_status", "missing"))
        data["balance_status"] = str(record.get("balance_status", "missing"))
        data["balance_compound"] = str(record.get("balance_compound", ""))
