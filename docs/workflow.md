# scMicroglia project workflow

## Overview

This project analyzes microglial states in Alzheimer’s disease using single-cell data, differential expression, gene-set enrichment, and regulatory network analysis. The workflow is organized around reproducible Python and R figure-generation scripts and a clean project folder structure.

## Recommended folder layout

```text
scMicroglia/
├── data/
│   ├── raw/
│   └── processed/
├── analysis/
├── results/
├── scripts/
├── R_scripts/
├── docs/
├── out/
└── README.md
```

## Workflow

1. Prepare raw files under `data/raw/`.
2. Create processed tables and QC outputs in `data/processed/`.
3. Run analysis scripts from `scripts/` and `R_scripts/`.
4. Save generated figures to `results/figures/` and tables to `results/tables/`.
5. Keep runtime exports in `out/` when needed for local inspection.
6. Document any major changes in this workflow file.

## Tools and languages

- Python: plotting, summary tables, annotation logic, figure generation
- R: publication-ready panels and shared style configuration

## Notes

- Raw data should remain separate from tracked results whenever it is large or sensitive.
- Reproducible scripts should be preferred over ad hoc notebook-only analysis.
- If a new analysis pipeline is introduced, create a short documentation note in this folder.
