# metabolic-maps

This repository constructs species-specific metabolic reaction graphs by integrating KEGG and MetaNetX data. The pipeline is configurable via YAML files, making it straightforward to generate maps for different organisms.

## Download links for maps

- Human ([map](./maps/human_kegg_graph_pruned_annotated.20260127_203125.graphml.gz?raw=1), [stats](./maps/human_kegg_graph_pruned_annotated.20260127_203125.pdf))
- Mouse ([map](./maps/mouse_kegg_graph_pruned_annotated.20260127_203125.graphml.gz?raw=1), [stats](./maps/mouse_kegg_graph_pruned_annotated.20260127_203125.pdf))
- Xenopus laevis ([map](./maps/xenopus-laevis_kegg_graph_pruned_annotated.20260127_203125.graphml.gz?raw=1), [stats](./maps/xenopus-laevis_kegg_graph_pruned_annotated.20260127_203125.pdf))
- Zebrafish ([map](./maps/zebrafish_kegg_graph_pruned_annotated.20260127_203125.graphml.gz?raw=1), [stats](./maps/zebrafish_kegg_graph_pruned_annotated.20260127_203125.pdf))

## Strategy used to create maps

1. Start with all KEGG reactions
2. Keep reactions with a species-specific KEGG:KO for an annotated enzyme
3. Find KEGG modules likely active in the species:
   - All modules where more than 66% of KO-covered reactions have a species KO
   - Remove modules that don't have a species-specific KEGG incarnation
   - Rescue manually specified modules (see `species/*.yaml` for per-species lists with rationale)
4. Find species KEGG pathways; add reactions in those pathways that have no enzymes (reactions with species enzymes are already included)
5. Intersect this reaction set with MetaNetX to build a KEGG × MetaNetX graph
6. Prune dead-end compounds
7. Export GraphML and PDF report

## Structure of the map/graph

- **Node types**: Each node has a `group` attribute with one of: `Compound`, `Reaction`, `Gene`, `Module`, `Pathway`
- **Node identifiers** (`name` attribute):
  - `C00001` — Compound (KEGG)
  - `R00286` — Reaction (KEGG)
  - `M00001` — Module (KEGG)
  - `map00010` — Pathway (KEGG)
  - `mmu:100037283` — Gene (KEGG/NCBI)
- **Additional node attributes**:
  - Gene: `gene_symbol`, `ensembl_id`, `ncbi_gene_id`
  - Module: `module_name`
  - Pathway: `pathway_name`
- **Edge attributes**: Edges between compounds and reactions have `coef` (stoichiometric coefficient) and `compartment`

## Prerequisites

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate flex_map_env
```

Required external data:
- `data/metanetx/MNXref.ttl` — MetaNetX RDF dump (download from [MetaNetX](https://www.metanetx.org/))

First run requires outbound HTTPS access for KEGG and MyGene downloads; subsequent runs use cached data in `data/kegg/`.

## Project structure

```
metabolic-maps/
├── src/
│   ├── kegg.py          # KEGG download/parsing, gene annotations
│   ├── metanetx.py      # MetaNetX graph building
│   ├── graph.py         # Pruning, annotation logic
│   └── pipeline.py      # Orchestration, PDF report
├── species/
│   ├── mmu.yaml         # Mouse config (includes rescued modules)
│   ├── ...
│   └── hsa.yaml         # Human config
├── data/                # Cached downloads (gitignored)
│   ├── kegg/            # KEGG tables and gene annotations
│   └── metanetx/
├── results/             # Output per-species graphs/reports (e.g., results/mmu/)
├── run.py               # CLI entry point
├── Makefile
└── environment.yml
```

## Usage

### Using Make (recommended)

```bash
# List available species
make list

# Build mouse metabolic map (full pipeline)
make mmu

# Build KEGG reaction set only (no MetaNetX required)
make mmu-kegg

# Build human metabolic map
make hsa

# Build all species (full maps)
make all

# Package latest PDFs and GraphML.gz into maps/
make dist

# Clean all outputs
make clean

# Clean only mouse outputs
make clean-mmu
```

The Makefile tracks dependencies so that:
- Editing `species/mmu.yaml` rebuilds only the mouse graph
- Editing `src/*.py` or `MNXref.ttl` rebuilds all species

### Using the CLI directly

```bash
# Generate mouse metabolic map (full pipeline)
python run.py species/mmu.yaml --metanetx data/metanetx/MNXref.ttl

# KEGG-only mode (no MetaNetX required)
python run.py species/mmu.yaml --kegg-only

# Skip PDF report
python run.py species/mmu.yaml --metanetx data/metanetx/MNXref.ttl --no-report

# Custom output directory
python run.py species/mmu.yaml --metanetx data/metanetx/MNXref.ttl --output-dir my_results/
```

### Outputs

Full pipeline outputs are stored under `results/<species_code>/`:
- `results/<species>/<prefix>.<timestamp>.graphml` — Metabolic graph
- `results/<species>/<prefix>.<timestamp>.pdf` — Summary report with schema documentation

Run `make dist` to gzip the latest GraphML per species and copy it together with its PDF into `maps/` for publication.

KEGG-only mode outputs (also under `results/<species_code>/`):
- `results/<species>/<prefix>.dropped_modules.<timestamp>.csv` — Modules that met coverage before being excluded
- `results/<species>/<prefix>.retained_modules.<timestamp>.csv` — Modules that satisfied all filters

## Species configuration

Species-specific parameters are defined in YAML files under `species/`:

```yaml
# species/mmu.yaml
species_code: mmu
name: Mus musculus
min_enzyme_coverage: 0.66

# Modules to include even without species-specific KEGG submaps
# Each entry should include a comment explaining the rationale
rescue_modules:
  - M00088  # Ketone body biosynthesis - liver expressed
  - M00074  # N-glycan biosynthesis, high-mannose type
  # ...

# Modules explicitly excluded despite meeting coverage threshold
excluded_modules:
  - M00165  # Calvin cycle - plant only
  - M00144  # NADH:quinone oxidoreductase, prokaryotes
  # ...

output_prefix: mouse_kegg_graph_pruned_annotated
```

See `species/mmu.yaml` for the complete mouse configuration including all rescued modules with rationale.

### Adding a new species

1. Create `species/<code>.yaml` with the appropriate config
2. Run `make <code>`

Gene annotations (NCBI IDs, Ensembl IDs, EC numbers) are fetched automatically from KEGG and MyGene.info on first run and cached in `data/kegg/<code>_ensembl_ec.tsv`.

### Supported species

The pipeline supports any species in KEGG. The MyGene.info lookup for Ensembl IDs has built-in mappings for:
- `mmu` — Mouse
- `hsa` — Human
- `dre` — Zebrafish
- `xla` — Xenopus laevis
