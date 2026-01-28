# Metabolic Maps - Build System
#
# Usage:
#   make mmu           # Build mouse map
#   make hsa           # Build human map
#   make all           # Build all species
#   make clean         # Remove all outputs
#   make clean-mmu     # Remove mouse outputs only

SHELL := /bin/bash
PYTHON := python
TIMESTAMP := $(shell date +%Y%m%d_%H%M%S)

# Directories
DATA_DIR := data
RESULTS_DIR := results
SPECIES_DIR := species
SRC_DIR := src
MAPS_DIR := maps

# External dependencies
METANETX_TTL := $(DATA_DIR)/metanetx/MNXref.ttl

# Source files (changes trigger full rebuild)
SRC_FILES := $(wildcard $(SRC_DIR)/*.py) run.py

# Species configs
SPECIES_CONFIGS := $(wildcard $(SPECIES_DIR)/*.yaml)
SPECIES_CODES := $(basename $(notdir $(SPECIES_CONFIGS)))

# Default target
.PHONY: all
all: $(SPECIES_CODES)

# KEGG-only targets (no MetaNetX required)
KEGG_ONLY_TARGETS := $(addsuffix -kegg,$(SPECIES_CODES))
.PHONY: $(KEGG_ONLY_TARGETS)
$(KEGG_ONLY_TARGETS): %-kegg: $(RESULTS_DIR)/.built.%-kegg

$(RESULTS_DIR)/.built.%-kegg: $(SPECIES_DIR)/%.yaml $(SRC_FILES)
	@echo "=== Building $* KEGG reaction set (no MetaNetX) ==="
	@mkdir -p $(RESULTS_DIR)
	$(PYTHON) run.py $< \
		--kegg-only \
		--data-dir $(DATA_DIR) \
		--output-dir $(RESULTS_DIR) \
		--timestamp $(TIMESTAMP)
	@touch $@
	@echo "=== $* KEGG-only complete ==="

# Pattern rule for species builds
# Depends on: config, source files, metanetx
.PHONY: $(SPECIES_CODES)
$(SPECIES_CODES): %: $(RESULTS_DIR)/.built.%

$(RESULTS_DIR)/.built.%: $(SPECIES_DIR)/%.yaml $(SRC_FILES) $(METANETX_TTL)
	@echo "=== Building $* metabolic map ==="
	@mkdir -p $(RESULTS_DIR)
	$(PYTHON) run.py $< \
		--metanetx $(METANETX_TTL) \
		--data-dir $(DATA_DIR) \
		--output-dir $(RESULTS_DIR) \
		--timestamp $(TIMESTAMP)
	@touch $@
	@echo "=== $* complete ==="

# MetaNetX dependency check
$(METANETX_TTL):
	@echo "ERROR: Missing MetaNetX file: $@"
	@echo "Download from https://www.metanetx.org/mnxdoc/mnxref.html"
	@exit 1

# Cleaning
.PHONY: clean
clean:
	rm -f $(RESULTS_DIR)/.built.*
	rm -f $(RESULTS_DIR)/*.graphml
	rm -f $(RESULTS_DIR)/*.pdf
	rm -f $(RESULTS_DIR)/*.csv
	rm -rf $(RESULTS_DIR)/*/
	rm -f $(DATA_DIR)/kegg/*.csv
	rm -f $(DATA_DIR)/kegg/*.tsv

.PHONY: clean-%
clean-%:
	rm -f $(RESULTS_DIR)/.built.$*
	rm -f $(RESULTS_DIR)/.built.$*-kegg
	rm -rf $(RESULTS_DIR)/$*/
	rm -f $(DATA_DIR)/kegg/$*_*.csv
	rm -f $(DATA_DIR)/kegg/$*_*.tsv

# List available species
.PHONY: list
list:
	@echo "Available species:"
	@for cfg in $(SPECIES_CONFIGS); do \
		code=$$(basename $$cfg .yaml); \
		name=$$(grep '^name:' $$cfg | cut -d: -f2 | xargs); \
		echo "  $$code - $$name"; \
	done

# Help
.PHONY: help
help:
	@echo "Metabolic Maps Build System"
	@echo ""
	@echo "Targets:"
	@echo "  make <species>       Build full metabolic map (e.g., make mmu)"
	@echo "  make <species>-kegg  Build KEGG reaction set only (e.g., make mmu-kegg)"
	@echo "  make all             Build all species (full maps)"
	@echo "  make list            List available species"
	@echo "  make clean           Remove all outputs"
	@echo "  make clean-<species> Remove outputs for one species"
	@echo "  make dist            Package latest PDFs + GraphML.gz into $(MAPS_DIR)/"
	@echo ""
	@echo "Available species: $(SPECIES_CODES)"

.PHONY: dist
dist: $(SPECIES_CODES)
	@mkdir -p $(MAPS_DIR)
	@for code in $(SPECIES_CODES); do \
		sp_dir=$(RESULTS_DIR)/$$code; \
		graphml=$$(ls -t $$sp_dir/*.graphml 2>/dev/null | head -n 1); \
		pdf=$$(ls -t $$sp_dir/*.pdf 2>/dev/null | head -n 1); \
		if [ -z "$$graphml" ]; then \
			echo "ERROR: No GraphML found for $$code. Run make $$code first."; exit 1; \
		fi; \
		if [ -z "$$pdf" ]; then \
			echo "ERROR: No PDF found for $$code. Run make $$code first."; exit 1; \
		fi; \
		out_graph=$(MAPS_DIR)/$$(basename $$graphml).gz; \
		out_pdf=$(MAPS_DIR)/$$(basename $$pdf); \
		$(PYTHON) -c "import gzip, pathlib, shutil, sys; src=pathlib.Path(sys.argv[1]); dst=pathlib.Path(sys.argv[2]); fin=src.open(\"rb\"); fout=gzip.open(dst, \"wb\"); shutil.copyfileobj(fin, fout); fin.close(); fout.close()" "$$graphml" "$$out_graph"; \
		cp "$$pdf" "$$out_pdf"; \
		echo "Packaged $$code -> $(MAPS_DIR)/"; \
	done
