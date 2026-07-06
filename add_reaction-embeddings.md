# Add Reaction Embeddings

## Goal

Add reaction orientation, reaction SMILES, and reaction embedding generation while keeping the current graph pipeline as the authoritative source of map structure.

The core pipeline should produce:

- an oriented GraphML file
- a reaction sidecar table with deterministic reaction strings
- existing summary/report files

Embedding generation should be a separate step that consumes the reaction sidecar table. It should not add Hugging Face, Torch, RXNGraphormer, or other heavy model dependencies to normal map generation.

## Scope

Implement four reaction representation families:

1. `side_mean`: current side-specific biochemical compound-pooling baseline.
2. `reactiont5`: ReactionT5 forward/retro encoder representation.
3. `rxngraphormer_2vec`: graph-transformer reaction encoder baseline.
4. `drfp`: DRFP transformation fingerprint baseline.

Do not modify gene/enzyme embedding code or train/test split logic unless strictly necessary.

## Phase 1: Oriented Reaction Table From The Map Pipeline

Add a per-run reaction table:

```text
results/<species>/<prefix>.reactions.<timestamp>.tsv
```

Required columns:

```text
reaction_node
reaction_id
raw_reaction
mnx_left_compounds
mnx_right_compounds
input_compounds
output_compounds
input_smiles
output_smiles
I_to_O_string
O_to_I_string
dg0
dg0_reaction
delta_g_status
balance_status
orientation
orientation_source
thermo_agrees_with_mnx
orientation_flip_from_mnx
missing_smiles_compounds
```

Definitions:

- `reaction_node`: exact graph node key, e.g. `kegg:R00286`.
- `reaction_id`: KEGG reaction ID, e.g. `R00286`.
- `raw_reaction`: original un-oriented KEGG/MetaNetX reaction formula where available.
- `mnx_left_compounds`: original MetaNetX left side before thermodynamic orientation.
- `mnx_right_compounds`: original MetaNetX right side before thermodynamic orientation.
- `input_compounds`: orientor-chosen input side `I`.
- `output_compounds`: orientor-chosen output side `O`.
- `input_smiles`: deterministic dot-separated canonical SMILES for `I`.
- `output_smiles`: deterministic dot-separated canonical SMILES for `O`.
- `I_to_O_string`: `input_smiles>>output_smiles`.
- `O_to_I_string`: `output_smiles>>input_smiles`.
- `dg0`: reaction DG0 used for orientation.
- `dg0_reaction`: reaction formula associated with the DG0 lookup/calculation.
- `delta_g_status`: source/status for the DG0 value.
- `balance_status`: balancing status for DG0 values calculated or recalculated with equilibrator.
- `orientation`: `forward`, `reverse`, or `unknown`.
- `orientation_source`: `dg0` or `mnx_fallback`.
- `thermo_agrees_with_mnx`: boolean indicating whether thermodynamic orientation keeps MetaNetX left-to-right direction.
- `orientation_flip_from_mnx`: boolean indicating whether the final input/output sides are flipped relative to MetaNetX.
- `missing_smiles_compounds`: compounds lacking SMILES, if any.

Orientation rule:

```python
if dg0 is missing:
    I = MetaNetX left
    O = MetaNetX right
    orientation = "unknown"
    orientation_source = "mnx_fallback"
elif dg0 > 0:
    I = MetaNetX right
    O = MetaNetX left
    orientation = "reverse"
    orientation_source = "dg0"
    thermo_agrees_with_mnx = False
    orientation_flip_from_mnx = True
else:
    I = MetaNetX left
    O = MetaNetX right
    orientation = "forward"
    orientation_source = "dg0"
    thermo_agrees_with_mnx = True
    orientation_flip_from_mnx = False
```

Use no uncertainty threshold for orientation.

For missing DG0, set:

```python
thermo_agrees_with_mnx = ""
orientation_flip_from_mnx = False
```

The blank agreement value avoids implying that the thermodynamic orientor agrees with MetaNetX when no thermodynamic evidence exists.

## DG0 Status And Trust Fields

Track DG0 provenance and balancing explicitly. These fields should be present in both the reaction table and reaction node attributes where available.

`delta_g_status` values:

```text
existing_modelseed
equilibrator_calculated
failed
missing
```

`balance_status` values:

```text
already_balanced
balanced_H_H2O
balanced_inorganic
balanced_currency
oxidation_balance
failed
not_applicable
missing
```

Recommended trust interpretation:

- `existing_modelseed` with known ModelSEED DG0 should be treated separately from recalculated equilibrator values.
- `already_balanced` and `balanced_H_H2O` are higher-confidence than broad currency-compound balancing.
- `balanced_inorganic` is intermediate and should list the compound used where possible.
- `balanced_currency` includes ATP, ADP, AMP, NAD, NADH, NADP, NADPH, FAD, FADH2, CoA, GDP, and GTP balancing; downstream code should be able to filter these out.
- `oxidation_balance` is useful but should be considered lower trust than directly balanced reactions.
- `failed` or `missing` DG0 values must not drive thermodynamic orientation.

If a specific balancing compound is used, add an optional `balance_compound` column/attribute.

## Phase 2: GraphML Updates

The GraphML remains the authoritative graph output.

Add reaction node attributes:

```text
dg0
dg0_reaction
delta_g_status
balance_status
balance_compound
orientation
orientation_source
thermo_agrees_with_mnx
orientation_flip_from_mnx
```

Add compound node attributes:

```text
smiles
smiles_source
```

Orient compound-reaction edges physically according to the orientor:

```text
input compound -> reaction -> output compound
```

If DG0 is missing, keep the current MetaNetX orientation and mark the reaction as `unknown`.

## Phase 3: SMILES Source

Use the KEGG compound to SMILES mapping:

```text
https://github.com/kostkalab/flex_embeddings/releases/download/v1.0.1/kegg_pubchem_smiles.20260129_222600.pkl.gz
```

Cache locally under `data/kegg/`.

Implementation notes:

- Load with `gzip` and `pickle`.
- Normalize keys to KEGG compound IDs like `C00031`.
- Canonicalize molecule SMILES before constructing reaction strings if RDKit or another project-approved canonicalizer is available.
- Sort molecules within each side for deterministic reaction strings.
- If any compound lacks SMILES, still write the reaction row, but mark representation generation for that row as unavailable where needed.

## Phase 4: Embedding Generation CLI

Add a separate command:

```bash
python generate_reaction_embeddings.py \
  results/mmu/<prefix>.reactions.<timestamp>.tsv \
  --representation all \
  --output-dir results/mmu/reaction_embeddings
```

Supported representations:

```text
side_mean
reactiont5
rxngraphormer_2vec
drfp
all
```

Do not expose separate CLI flags for ReactionT5 sum/concat variants or RXNGraphormer forward/reverse variants. Each representation family should generate its useful internal variants automatically.

## Phase 5: Common Embedding Output Contract

Every representation should save one self-contained Parquet table:

```text
reaction_embeddings_<name>.parquet
```

Common fields:

```text
reaction_node
embedding_status
error_message
embedding
```

The `embedding` column stores the fixed-length numeric vector for that
representation. Reaction metadata such as compounds, DG0, orientation, and
reaction strings remain in the reaction sidecar table and should be joined by
`reaction_node` when needed.

Also write a run-level report:

```text
reaction_embedding_generation_report.md
```

## Phase 6: Representation 1, `side_mean`

Preserve the current side-specific biochemical baseline.

For oriented reaction `r = (I, O)`:

```python
z_input = mean(compound_embedding[c] for c in I)
z_output = mean(compound_embedding[c] for c in O)
z_side = concat(z_input, z_output)
```

Output:

```text
reaction_embeddings_side_mean_fingerprint_raw.parquet
reaction_embeddings_side_mean_fingerprint_pca256.parquet
```

Optional non-default variant:

```python
z_diff = z_output - z_input
z_side_plus_diff = concat(z_input, z_output, z_diff)
```

Do not make the difference variant the default.

## Phase 7: Representation 2, `reactiont5`

Use Hugging Face ReactionT5 models from the `sagawa/reactiont5` family.

Preferred models:

```text
sagawa/ReactionT5v2-forward
sagawa/ReactionT5v2-retrosynthesis
```

If exact model names differ, inspect the Hugging Face collection and use the closest non-dataset-specific v2 forward and retrosynthesis models. Record the actual model names in metadata.

Input strings:

```text
I_to_O_string
O_to_I_string
```

Do not use generation output. Use encoder hidden states as embeddings.

Pooling:

```python
hidden = encoder_last_hidden_state
mask = attention_mask
embedding = (hidden * mask[..., None]).sum(axis=1) / mask.sum(axis=1, keepdims=True)
```

Normalize extracted encoder vectors consistently, for example with L2 normalization before concatenating or summing.

Generate the canonical four-vector output:

```text
[
  E_forward(I_to_O),
  E_retro(O_to_I),
  E_retro(I_to_O),
  E_forward(O_to_I)
]
```

Output:

```text
reaction_embeddings_reactiont5_4vec.parquet
```

Document that the symmetric sum representation can be derived from the four-vector:

```python
z_top = E_forward(I_to_O) + E_retro(O_to_I)
z_bottom = E_retro(I_to_O) + E_forward(O_to_I)
z_sum = concat(z_top, z_bottom)
```

This construction satisfies:

```text
r_T5(O, I) = swap(r_T5(I, O))
```

## Phase 8: Representation 3, `rxngraphormer_2vec`

Add RXNGraphormer as a modern graph-based reaction encoder baseline:

```text
https://github.com/licheng-xu-echo/RXNGraphormer
```

Use the same deterministic oriented reaction strings as all other representations.

Extraction strategy, in order:

1. Use an official embedding or feature extraction script if provided.
2. If no explicit embedding script exists, load the pretrained model and extract the pooled reaction representation immediately before the prediction head.
3. If the model only exposes task predictions and no clean embedding hook, add a small wrapper that returns the final hidden reaction representation.
4. If pretrained weights or required preprocessing are unavailable, document this clearly and mark RXNGraphormer as `not_available`.

Do not run full model training.

Required output if extraction succeeds:

```text
reaction_embeddings_rxngraphormer_2vec.parquet
```

where:

```text
[
  E_graphormer(I_to_O),
  E_graphormer(O_to_I)
]
```

This includes the one-direction representation as its first block, so do not
publish a separate one-direction RXNGraphormer file unless needed for a specific
debugging comparison. If `rxngraphormer_2vec` cannot be generated, record the
blocker in the report.

## Phase 9: Representation 4, `drfp`

Use DRFP as the transformation-fingerprint baseline.

For each oriented reaction:

```text
z_drfp = DRFP(I_to_O_string)
```

Output:

```text
reaction_embeddings_drfp.parquet
```

DRFP is not side-specific. It is included as a cheap transformation control, not as the main biochemical representation.

## Phase 10: RXNFP Decision

Do not require RXNFP in the main implementation.

Reason:

- `rxnfp` is not available from conda-forge.
- The PyPI package is old and pins dependencies incompatible with the current Python 3.11 environment.

If RXNFP is needed later, support a precomputed table input rather than making it part of the normal environment:

```bash
--rxnfp-table path/to/precomputed_rxnfp.tsv
```

## Dependency Strategy

Keep normal map generation lightweight.

Core environment:

- existing project dependencies
- optionally `drfp` if it installs cleanly in the current environment

Model-heavy dependencies should be optional and documented separately:

- ReactionT5: Hugging Face `transformers`, `torch`, model weights
- RXNGraphormer: repository-specific environment/checkpoint requirements

Prefer separate optional environment files if needed:

```text
extras/reactiont5_environment.yml
extras/rxngraphormer_environment.yml
```

## Smoke Tests

Core reaction table:

- Verify every row has `reaction_id`, `reaction_node`, orientation fields, and deterministic reaction strings where SMILES are complete.
- Verify raw MetaNetX sides are preserved as `mnx_left_compounds` and `mnx_right_compounds`.
- Verify sampled graph edge directions match `input_compounds -> reaction -> output_compounds`.
- Verify `thermo_agrees_with_mnx` and `orientation_flip_from_mnx` match the DG0 orientation rule.
- Verify `delta_g_status` and `balance_status` are populated for every row.
- Verify missing DG0 reactions retain MetaNetX orientation and are marked `mnx_fallback`.

DRFP:

- Verify finite numeric vectors for a small sample.
- Verify deterministic output for repeated runs.

ReactionT5:

- Verify finite numeric vectors for a small sample.
- Verify deterministic output in eval mode.
- Verify the documented swap property for the derived symmetric sum representation.

RXNGraphormer:

- Verify finite numeric vectors if extraction succeeds.
- Verify deterministic output in eval mode.
- If `_2` is generated, verify forward/reverse block swapping behaves as expected.
- If extraction fails, verify the report names the exact blocker.

## Report Contents

`reaction_embedding_generation_report.md` should include:

```text
n_reactions_total
n_reactions_with_complete_smiles
n_reactions_missing_smiles
n_reactions_dg0_oriented
n_reactions_mnx_fallback
n_reactions_thermo_agrees_with_mnx
n_reactions_orientation_flip_from_mnx
n_reactions_delta_g_existing_modelseed
n_reactions_delta_g_equilibrator_calculated
n_reactions_delta_g_failed
n_reactions_balance_already_balanced
n_reactions_balance_H_H2O
n_reactions_balance_inorganic
n_reactions_balance_currency
n_reactions_balance_oxidation_balance
n_reactions_balance_failed
n_reactions_embedded_side_mean
n_reactions_embedded_reactiont5
n_reactions_embedded_rxngraphormer_2vec
n_reactions_embedded_drfp
reactiont5_forward_model
reactiont5_retro_model
rxngraphormer_checkpoint
rxngraphormer_repo_commit
rxngraphormer_environment
failures_by_representation
```

Also include a short interpretation of each representation:

- `side_mean`: side-identity biochemical compound pooling baseline.
- `reactiont5`: modern reaction-string encoder baseline with forward/retro half-space structure.
- `rxngraphormer_2vec`: modern graph-transformer reaction encoder baseline.
- `drfp`: reaction-change fingerprint baseline.

## Implementation Notes

Recommended module split:

```text
src/reactions.py       # oriented reaction table construction
src/chem.py            # SMILES loading, canonicalization, reaction string construction
src/thermo.py          # DG0 loading/orientation helpers
src/embeddings/
  side_mean.py
  reactiont5.py
  rxngraphormer.py
  drfp.py
  io.py
generate_reaction_embeddings.py
```

The existing `process_kegg_reactions.py` should be treated as prototype/reference code. Reusable pieces should be moved into `src/` functions before being called by the pipeline.
