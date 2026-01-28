#!/usr/bin/env python3
"""CLI entry point for metabolic map generation."""

import argparse
from datetime import datetime
from pathlib import Path

import yaml

from src.pipeline import SpeciesConfig, run_pipeline, generate_report


def load_config(config_path: Path) -> SpeciesConfig:
    """Load species config from YAML file."""
    with open(config_path) as f:
        data = yaml.safe_load(f)

    return SpeciesConfig(
        species_code=data["species_code"],
        name=data["name"],
        min_enzyme_coverage=data.get("min_enzyme_coverage", 0.66),
        rescue_modules=data.get("rescue_modules") or [],
        excluded_modules=data.get("excluded_modules") or [],
        output_prefix=data.get(
            "output_prefix", f"{data['species_code']}_metabolic_graph"
        ),
    )


def main():
    parser = argparse.ArgumentParser(
        description="Generate species-specific metabolic maps from KEGG and MetaNetX"
    )
    parser.add_argument(
        "config",
        type=Path,
        help="Path to species config YAML file",
    )
    parser.add_argument(
        "--metanetx",
        type=Path,
        default=None,
        help="Path to MetaNetX TTL file (required unless --kegg-only)",
    )
    parser.add_argument(
        "--kegg-only",
        action="store_true",
        help="Skip MetaNetX intersection; output KEGG reaction set only",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("data"),
        help="Directory for cached data (default: data/)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results"),
        help="Directory for output files (default: results/)",
    )
    parser.add_argument(
        "--timestamp",
        type=str,
        default=None,
        help="Timestamp for output filename (default: current time)",
    )

    parser.add_argument(
        "--no-report",
        action="store_true",
        help="Skip PDF report generation",
    )

    args = parser.parse_args()

    if args.timestamp is None:
        args.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    if not args.kegg_only and args.metanetx is None:
        parser.error("--metanetx is required unless --kegg-only is specified")

    config = load_config(args.config)
    species_output_dir = args.output_dir / config.species_code

    result = run_pipeline(
        config=config,
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        metanetx_ttl=args.metanetx,
        timestamp=args.timestamp,
        kegg_only=args.kegg_only,
    )

    print(f"\nPipeline complete. Metrics:")
    for k, v in result.metrics.items():
        print(f"  {k}: {v}")

    # Generate PDF report (skip for kegg-only mode)
    if not args.no_report and not args.kegg_only:
        report_filename = f"{config.output_prefix}.{args.timestamp}.pdf"
        report_path = species_output_dir / report_filename
        generate_report(config, result, report_path)


if __name__ == "__main__":
    main()
