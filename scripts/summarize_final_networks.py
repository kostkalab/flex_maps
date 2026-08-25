#!/usr/bin/env python3
"""Summarize the final flex metabolic networks as TSV and LaTeX."""

from __future__ import annotations

import argparse
import csv
import re
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


GRAPHML_NS = "{http://graphml.graphdrawing.org/xmlns}"
TIMESTAMP_RE = re.compile(r"\.(\d{8}_\d{6})\.graphml$")


@dataclass(frozen=True)
class SpeciesConfig:
    code: str
    name: str
    output_prefix: str
    min_enzyme_coverage: float
    rescue_modules: frozenset[str]
    excluded_modules: frozenset[str]


@dataclass(frozen=True)
class NetworkSummary:
    species: str
    genes: int
    reactions: int
    compounds: int
    direct: int
    module_added: int
    pathway_only: int
    balanced_as_represented: int
    balanced_after_repair: int
    no_thermodynamic_estimate: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--timestamp",
        help="Pipeline timestamp; defaults to the latest timestamp present for every species",
    )
    parser.add_argument("--data-dir", type=Path, default=Path("data/kegg"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    parser.add_argument("--species-dir", type=Path, default=Path("species"))
    parser.add_argument(
        "--tsv",
        type=Path,
        help="Write the compact six-column table as TSV",
    )
    parser.add_argument(
        "--latex",
        type=Path,
        help="Write LaTeX to this file instead of standard output",
    )
    return parser.parse_args()


def read_commented_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        lines = (line for line in handle if not line.startswith("#"))
        return list(csv.DictReader(lines))


def parse_species_config(path: Path) -> SpeciesConfig:
    scalar_keys = {"species_code", "name", "output_prefix", "min_enzyme_coverage"}
    scalars: dict[str, str] = {}
    lists: dict[str, list[str]] = {
        "rescue_modules": [],
        "excluded_modules": [],
    }
    active_list: str | None = None

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        content = raw_line.split("#", 1)[0].rstrip()
        stripped = content.strip()
        if not stripped:
            continue
        if stripped.startswith("-") and active_list is not None:
            value = stripped[1:].strip().strip("'\"")
            if value:
                lists[active_list].append(value)
            continue
        active_list = None
        if ":" not in stripped:
            continue
        key, value = stripped.split(":", 1)
        key, value = key.strip(), value.strip()
        if key in lists:
            active_list = key
        elif key in scalar_keys:
            scalars[key] = value.strip("'\"")

    missing = sorted(scalar_keys - set(scalars))
    if missing:
        raise ValueError(f"Missing fields in {path}: {', '.join(missing)}")
    return SpeciesConfig(
        code=scalars["species_code"],
        name=scalars["name"],
        output_prefix=scalars["output_prefix"],
        min_enzyme_coverage=float(scalars["min_enzyme_coverage"]),
        rescue_modules=frozenset(lists["rescue_modules"]),
        excluded_modules=frozenset(lists["excluded_modules"]),
    )


def load_species(species_dir: Path) -> list[SpeciesConfig]:
    configs = [parse_species_config(path) for path in sorted(species_dir.glob("*.yaml"))]
    if not configs:
        raise RuntimeError(f"No species configurations found in {species_dir}")
    return configs


def available_timestamps(results_dir: Path, config: SpeciesConfig) -> set[str]:
    timestamps = set()
    result_dir = results_dir / config.code
    for path in result_dir.glob(f"{config.output_prefix}.*.graphml"):
        if ".full." in path.name:
            continue
        match = TIMESTAMP_RE.search(path.name)
        if match:
            timestamps.add(match.group(1))
    return timestamps


def detect_latest_common_timestamp(
    results_dir: Path, configs: list[SpeciesConfig]
) -> str:
    common: set[str] | None = None
    for config in configs:
        timestamps = available_timestamps(results_dir, config)
        common = timestamps if common is None else common & timestamps
    if not common:
        raise RuntimeError("No common final-network timestamp found for all species")
    return max(common)


def graph_nodes(path: Path) -> list[dict[str, str]]:
    root = ET.parse(path).getroot()
    keys = {
        key.attrib["id"]: key.attrib.get("attr.name", key.attrib["id"])
        for key in root.findall(GRAPHML_NS + "key")
        if key.attrib.get("for") == "node"
    }
    nodes = []
    for node in root.findall(".//" + GRAPHML_NS + "node"):
        nodes.append(
            {
                keys.get(data.attrib["key"], data.attrib["key"]): data.text or ""
                for data in node.findall(GRAPHML_NS + "data")
            }
        )
    return nodes


def summarize_networks(
    configs: list[SpeciesConfig], data_dir: Path, results_dir: Path, timestamp: str
) -> list[NetworkSummary]:
    ko_to_reactions: dict[str, set[str]] = defaultdict(set)
    for row in read_commented_csv(data_dir / "dn_ko_to_rn.csv"):
        ko_to_reactions[row["ko"]].add(row["rn"])
    enzyme_reactions = set().union(*ko_to_reactions.values())

    module_reactions: dict[str, set[str]] = defaultdict(set)
    for row in read_commented_csv(data_dir / "tb_rn_to_md.csv"):
        module_reactions[row["md"]].add(row["rn"])

    pathway_reactions: dict[str, set[str]] = defaultdict(set)
    for row in read_commented_csv(data_dir / "dn_pathway_to_rn.csv"):
        pathway_reactions[row["path"]].add(row["rn"])

    summaries = []
    for config in configs:
        species_reactions: set[str] = set()
        for row in read_commented_csv(data_dir / f"dn_{config.code}gene_to_ko.csv"):
            species_reactions.update(ko_to_reactions.get(row["ko"], ()))

        species_modules = {
            row["md"].removeprefix(f"{config.code}_")
            for row in read_commented_csv(
                data_dir / f"dn_{config.code}_to_modules.csv"
            )
        }
        selected_modules = set()
        for module, reactions in module_reactions.items():
            denominator = len(reactions & enzyme_reactions)
            coverage = (
                len(reactions & species_reactions) / denominator if denominator else 0.0
            )
            if (
                coverage >= config.min_enzyme_coverage
                and (module in species_modules or module in config.rescue_modules)
                and module not in config.excluded_modules
            ):
                selected_modules.add(module)
        selected_module_reactions = set().union(
            *(module_reactions[module] for module in selected_modules)
        )

        species_pathways = {
            re.sub(r"^[a-z]+", "map", row["pathway"])
            for row in read_commented_csv(
                data_dir / f"dn_{config.code}_pathways.csv"
            )
        }
        pathway_context = set().union(
            *(pathway_reactions[pathway] for pathway in species_pathways)
        ) - enzyme_reactions

        graph_path = (
            results_dir
            / config.code
            / f"{config.output_prefix}.{timestamp}.graphml"
        )
        if not graph_path.exists():
            raise FileNotFoundError(graph_path)
        nodes = graph_nodes(graph_path)
        group_counts = Counter(node.get("group") for node in nodes)
        reaction_nodes = [node for node in nodes if node.get("group") == "Reaction"]
        final_reactions = {node["name"] for node in reaction_nodes}

        direct = final_reactions & species_reactions
        module_added = (final_reactions - direct) & selected_module_reactions
        pathway_only = (final_reactions - direct - module_added) & pathway_context
        unclassified = final_reactions - direct - module_added - pathway_only
        if unclassified:
            raise RuntimeError(
                f"{config.code}: {len(unclassified)} final reactions have no inclusion "
                f"route: {', '.join(sorted(unclassified))}"
            )

        thermo = Counter()
        for node in reaction_nodes:
            if node.get("delta_g_status") != "equilibrator_calculated":
                thermo["none"] += 1
            elif node.get("balance_status") == "already_balanced":
                thermo["represented"] += 1
            else:
                thermo["repaired"] += 1

        summaries.append(
            NetworkSummary(
                species=config.name,
                genes=group_counts["Gene"],
                reactions=group_counts["Reaction"],
                compounds=group_counts["Compound"],
                direct=len(direct),
                module_added=len(module_added),
                pathway_only=len(pathway_only),
                balanced_as_represented=thermo["represented"],
                balanced_after_repair=thermo["repaired"],
                no_thermodynamic_estimate=thermo["none"],
            )
        )
    return summaries


def percent(count: int, total: int) -> str:
    return f"{100 * count / total:.1f}"


def inclusion_route(summary: NetworkSummary) -> str:
    return " / ".join(
        percent(count, summary.reactions)
        for count in (summary.direct, summary.module_added, summary.pathway_only)
    )


def thermodynamic_annotation(summary: NetworkSummary) -> str:
    return " / ".join(
        percent(count, summary.reactions)
        for count in (
            summary.balanced_as_represented,
            summary.balanced_after_repair,
            summary.no_thermodynamic_estimate,
        )
    )


def abbreviated_species(name: str) -> str:
    genus, *rest = name.split()
    return f"{genus[0]}. {' '.join(rest)}"


def latex_table(summaries: list[NetworkSummary], timestamp: str) -> str:
    rows = []
    for summary in summaries:
        rows.append(
            " & ".join(
                [
                    rf"\textit{{{abbreviated_species(summary.species)}}}",
                    f"{summary.genes:,}",
                    f"{summary.reactions:,}",
                    f"{summary.compounds:,}",
                    inclusion_route(summary),
                    thermodynamic_annotation(summary),
                ]
            )
            + r" \\"
        )
    body = "\n".join(rows)
    return rf"""% Generated by scripts/summarize_final_networks.py --timestamp {timestamp}
\begin{{table}}[ht]
\centering
\caption{{Final structure-supported metabolic networks used by flex. Inclusion
routes are reported as direct gene support / module-added / pathway-only, with
routes assigned hierarchically in that order. Thermodynamic annotation is
reported as balanced as represented / balanced after repair / no thermodynamic
estimate. Values in the final two columns are percentages of final reactions.}}
\label{{tab:metabolic_networks}}
\small
\begin{{tabular}}{{lrrrrl}}
\toprule
Species & Genes & Reactions & Compounds & Inclusion (\%) & Thermo. (\%) \\
\midrule
{body}
\bottomrule
\end{{tabular}}
\end{{table}}
"""


def write_tsv(summaries: list[NetworkSummary], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "Species",
                "Genes",
                "Reactions",
                "Compounds",
                "Inclusion (%) (direct / module / pathway)",
                "Thermodynamic annotation (%) (represented / repaired / none)",
            ]
        )
        for summary in summaries:
            writer.writerow(
                [
                    summary.species,
                    summary.genes,
                    summary.reactions,
                    summary.compounds,
                    inclusion_route(summary),
                    thermodynamic_annotation(summary),
                ]
            )


def main() -> int:
    args = parse_args()
    configs = load_species(args.species_dir)
    timestamp = args.timestamp or detect_latest_common_timestamp(
        args.results_dir, configs
    )
    summaries = summarize_networks(
        configs, args.data_dir, args.results_dir, timestamp
    )
    rendered = latex_table(summaries, timestamp)
    if args.tsv:
        write_tsv(summaries, args.tsv)
    if args.latex:
        args.latex.parent.mkdir(parents=True, exist_ok=True)
        args.latex.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
