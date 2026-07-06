# metabolic-maps

This repository constructs species-specific metabolic reaction graphs by integrating KEGG and MetaNetX data. The pipeline is configurable via YAML files, making it straightforward to generate maps for different organisms.

## Download links for maps

- Danio rerio ([map](./maps/zebrafish_kegg_graph_pruned_annotated.20260702_163441.graphml.gz?raw=1), [stats](./maps/zebrafish_kegg_graph_pruned_annotated.20260702_163441.pdf))
- Homo sapiens ([map](./maps/human_kegg_graph_pruned_annotated.20260702_163441.graphml.gz?raw=1), [stats](./maps/human_kegg_graph_pruned_annotated.20260702_163441.pdf))
- Mus musculus ([map](./maps/mouse_kegg_graph_pruned_annotated.20260702_163441.graphml.gz?raw=1), [stats](./maps/mouse_kegg_graph_pruned_annotated.20260702_163441.pdf))
- Xenopus laevis ([map](./maps/xenopus-laevis_kegg_graph_pruned_annotated.20260702_163441.graphml.gz?raw=1), [stats](./maps/xenopus-laevis_kegg_graph_pruned_annotated.20260702_163441.pdf))

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
7. Annotate compounds with joint KEGG SMILES support and orient reaction sides
8. Save an internal full GraphML, then drop reactions where all input compounds or all output compounds lack SMILES support
9. Prune any dead-end compounds and orphan annotations created by the SMILES filter
10. Export the SMILES-supported GraphML, reaction side table, drop reports, and PDF report

## Structure of the map/graph

- **Node types**: Each node has a `group` attribute with one of: `Compound`, `Reaction`, `Gene`, `Module`, `Pathway`
- **Node identifiers** (`name` attribute):
  - `C00001` — Compound (KEGG)
  - `R00286` — Reaction (KEGG)
  - `M00001` — Module (KEGG)
  - `map00010` — Pathway (KEGG)
  - `mmu:100037283` — Gene (KEGG/NCBI)
- **Additional node attributes**:
  - Compound: `compound_name`, `smiles`, `smiles_source`
  - Reaction: `raw_reaction`, `dg0`, `dg0_reaction`, `delta_g_status`, `balance_status`, `orientation`, `orientation_source`
  - Gene: `gene_symbol`, `ensembl_id`, `ncbi_gene_id`
  - Module: `module_name`
  - Pathway: `pathway_name`
- **Edge attributes**: Edges between compounds and reactions have `coef` (stoichiometric coefficient) and `compartment`

## Prerequisites

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate flex_maps_env
```

Reaction embeddings that use deep learning models should use a separate
environment. The main map environment intentionally avoids PyTorch, Hugging
Face, PyG, and RXNGraphormer dependencies so normal graph construction remains
lightweight and compatible with the Python 3.11 / NumPy 2 stack.

```bash
conda env create -f extras/dl_embeddings_environment.yml
conda activate flex_maps_dl_env
```

This optional environment is anchored on the RXNGraphormer PyTorch 2 support
stack: Python 3.10, PyTorch 2.2.1, CUDA 12.1, PyG 2.6.1, NumPy 1.26.4, and
RDKit 2024.03. It uses a newer Hugging Face stack than RXNGraphormer documents
because ReactionT5v2's tokenizer requires a newer `tokenizers` parser than
`transformers==4.30.0` provides. Use the CUDA version only on machines with a
compatible NVIDIA driver; CPU-only installs can replace `pytorch-cuda=12.1`
with `cpuonly` and use the PyG CPU wheel index.

On systems where `$HOME` is read-only or shared, set writable cache locations
before importing RXNGraphormer/localmapper/DGL components:

```bash
export MPLCONFIGDIR=/tmp/mplconfig
export DGLBACKEND=pytorch
export DGL_DOWNLOAD_DIR=/tmp/dgl
export XDG_CACHE_HOME=/tmp/cache
```

RXNGraphormer embedding extraction also requires pretrained model weights from
the RXNGraphormer Figshare archive. For reaction embeddings, download and
extract `pretrained_classification_model.zip` from
<https://doi.org/10.6084/m9.figshare.28356077> under:

```text
data/rxngraphormer/RXNGraphormer/model_path/pretrained_classification_model/
```

The expected files are:

```text
parameters.json
model/valid_checkpoint.pt
```

The larger task-specific model and dataset archives from the same Figshare
record are only needed for reproducing RXNGraphormer training/evaluation tasks,
not for pretrained embedding extraction.

DL reaction embeddings can then be generated from the oriented reaction sidecar
table. GPU `0` is used by default; pass `--gpu 1` for another GPU or
`--gpu cpu` to force CPU execution.

```bash
python generate_reaction_embeddings.py \
  results/mmu/mouse_kegg_graph_pruned_annotated.reactions.20260702_163441.tsv \
  --representation dl \
  --output-dir results/mmu/reaction_embeddings_dl \
  --offline \
  --gpu 0
```

The DL parquet files use one row per `reaction_node`; the `embedding` column is
the complete numeric vector for that representation. Reaction orientation comes
from the sidecar table:

- `I_to_O_string`: `input_smiles>>output_smiles`
- `O_to_I_string`: `output_smiles>>input_smiles`

`reaction_embeddings_reactiont5_4vec.parquet` is a concatenation of four
L2-normalized ReactionT5 encoder vectors:

```text
[
  E_forward(I_to_O),
  E_retro(O_to_I),
  E_retro(I_to_O),
  E_forward(O_to_I)
]
```

The forward model is `sagawa/ReactionT5v2-forward`; the retrosynthesis model is
`sagawa/ReactionT5v2-retrosynthesis`. Embeddings are encoder hidden states, not
generated sequences, pooled by masked mean over the encoder tokens.

`reaction_embeddings_rxngraphormer_2vec.parquet` is a concatenation of two
RXNGraphormer vectors:

```text
[
  E_graphormer(I_to_O),
  E_graphormer(O_to_I)
]
```

The first half is the one-direction RXNGraphormer representation, so a separate
one-direction RXNGraphormer file is not published by default.

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
- `results/<species>/<prefix>.<timestamp>.graphml` — Published SMILES-supported metabolic graph
- `results/<species>/<prefix>.full.<timestamp>.graphml` — Internal unfiltered graph before SMILES support filtering
- `results/<species>/<prefix>.reactions.<timestamp>.tsv` — Oriented reaction side table with input/output compounds and SMILES strings
- `results/<species>/<prefix>.smiles_dropped_reactions.<timestamp>.tsv` — Reactions removed because an entire side lacked SMILES, plus dead-end cascade removals
- `results/<species>/<prefix>.smiles_dropped_compounds.<timestamp>.tsv` — Compounds removed during post-SMILES dead-end pruning
- `results/<species>/<prefix>.<timestamp>.pdf` — Summary report with schema documentation

Run `make dist` to gzip the latest SMILES-supported GraphML per species and copy it together with its PDF into `maps/` for publication. The internal `*.full.*.graphml` files are retained in `results/` but are not published.

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

## Rationale for excluding non-SMILES network nodes

The published maps are filtered to the small-molecule reaction network that can support structural embeddings and cheminformatics workflows. Reactions are removed only when an entire input side or output side lacks supported SMILES, because that side cannot be represented structurally.

- **Mass-balance validation**: SMILES strings provide explicit molecular structure from which formulas can be derived. Generic placeholders, redox carriers without selected structures, and macromolecules lack defined small-molecule stoichiometry, making automated element and mass-balance checks unreliable.
- **Property prediction**: Structural-activity and embedding pipelines depend on fixed molecular descriptors such as Morgan fingerprints. Variable fatty acid pools, generic classes, and protein-linked entities cannot be converted into uniform numerical vectors without adding unsupported assumptions.
- **Scope definition**: Molecular graph tools such as RDKit are designed for small-molecule chemical spaces. Large protein carriers and enzyme-linked entities belong in the annotation or enzyme layer rather than the soluble metabolite layer used for graph embeddings.
