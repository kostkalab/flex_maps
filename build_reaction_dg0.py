#!/usr/bin/env python3
"""Calculate DG0 values from graph-derived reaction formulas."""

from __future__ import annotations

import argparse
from pathlib import Path

import networkx as nx
import pandas as pd

from src import thermo


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate DG0 for reactions using the exact graph side order."
    )
    parser.add_argument(
        "graphml",
        type=Path,
        nargs="+",
        help="Input full GraphML file(s)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output TSV path for DG0 records",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Calculate only the first N reactions, for smoke tests",
    )
    parser.add_argument(
        "--reaction-id",
        action="append",
        default=None,
        help="Restrict calculation to a KEGG reaction ID; may be provided multiple times",
    )
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=25,
        help="Write partial output every N calculated reactions",
    )
    parser.add_argument(
        "--seed-table",
        type=Path,
        action="append",
        default=[],
        help=(
            "Reuse rows from an existing DG0 TSV when reaction_id and "
            "graph-derived formula match exactly; may be provided multiple times"
        ),
    )
    return parser.parse_args()


def load_seed_rows(seed_tables: list[Path]) -> dict[tuple[str, str], dict]:
    """Load reusable rows keyed by (reaction_id, dg0_reaction)."""
    rows = {}
    for seed_table in seed_tables:
        if not seed_table.exists():
            raise FileNotFoundError(seed_table)
        df = pd.read_csv(seed_table, sep="\t", dtype=str).fillna("")
        for _, row in df.iterrows():
            rid = str(row.get("reaction_id", ""))
            reaction = str(row.get("dg0_reaction", ""))
            if not rid or not reaction:
                continue
            rows.setdefault((rid, reaction), row.to_dict())
    return rows


def apply_seed_rows(
    formulas: dict[str, tuple[thermo.ReactionFormula, set[str]]],
    seed_rows: dict[tuple[str, str], dict],
) -> tuple[list[dict], dict[str, tuple[thermo.ReactionFormula, set[str]]]]:
    """Copy seed rows for formulas already calculated elsewhere."""
    reused = []
    remaining = {}
    for rid, (formula, sources) in sorted(formulas.items()):
        formula_text = formula.to_equilibrator_formula()
        seed = seed_rows.get((rid, formula_text))
        if seed is None:
            remaining[rid] = (formula, sources)
            continue
        row = {column: seed.get(column, "") for column in thermo.DG0_COLUMNS}
        row["reaction_node"] = formula.reaction_node
        row["reaction_id"] = formula.reaction_id
        row["source_graphs"] = "|".join(sorted(s for s in sources if s))
        row["dg0_reaction"] = formula_text
        reused.append(row)
    return reused, remaining


def main() -> int:
    args = parse_args()
    formula_maps = []
    for graphml in args.graphml:
        G = nx.read_graphml(graphml)
        formula_maps.append(
            thermo.graph_reaction_formulas(G, source_graph=str(graphml))
        )

    formulas, conflicts = thermo.merge_graph_reaction_formulas(formula_maps)
    if args.reaction_id:
        keep = set(args.reaction_id)
        formulas = {rid: value for rid, value in formulas.items() if rid in keep}
        conflicts = {rid: value for rid, value in conflicts.items() if rid in keep}
    if conflicts:
        print("ERROR: conflicting graph-derived formulas for reaction IDs:")
        for rid, values in sorted(conflicts.items()):
            print(f"  {rid}")
            for value in sorted(values):
                print(f"    {value}")
        return 1

    if args.output.exists():
        existing = pd.read_csv(args.output, sep="\t", dtype=str).fillna("")
        done = set(existing["reaction_id"]) if "reaction_id" in existing.columns else set()
        formulas = {rid: value for rid, value in formulas.items() if rid not in done}
        rows = existing.to_dict("records")
        print(f"Resuming from {args.output}: {len(done):,} rows already present")
    else:
        rows = []

    if args.limit is not None:
        formulas = dict(list(sorted(formulas.items()))[: args.limit])

    if args.seed_table:
        seed_rows = load_seed_rows(args.seed_table)
        reused_rows, formulas = apply_seed_rows(formulas, seed_rows)
        rows.extend(reused_rows)
        print(
            f"Reused {len(reused_rows):,} rows from {len(args.seed_table):,} seed table(s); "
            f"{len(formulas):,} rows require calculation"
        )

    try:
        new_rows = thermo.calculate_formula_records_dg0_rows(
            formulas,
            checkpoint_path=args.output,
            checkpoint_every=args.checkpoint_every,
            initial_rows=rows,
        )
        table = pd.DataFrame(new_rows, columns=thermo.DG0_COLUMNS)
    except RuntimeError as exc:
        print(f"ERROR: {exc}")
        return 1
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, sep="\t", index=False)

    status_counts = table["delta_g_status"].value_counts(dropna=False).to_dict()
    balance_counts = table["balance_status"].value_counts(dropna=False).to_dict()
    print(f"Wrote {len(table):,} DG0 rows to {args.output}")
    print(f"delta_g_status: {status_counts}")
    print(f"balance_status: {balance_counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
