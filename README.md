# scMicroglia

Single-cell analysis of microglial states in Alzheimer’s disease, with a focus on sex-specific and cross-region molecular differences between hippocampus (HPC) and prefrontal cortex (PFC).

This repository contains the reproducible analysis and figure-generation workflow used to study microglial state annotations, differential expression, pathway enrichment, and transcription-factor regulatory networks in Alzheimer’s disease. The project combines Python-based figure generation with R-based publication plotting and a shared manuscript style system.

## Study summary

The codebase implements a single-cell analysis pipeline that:

- profiles microglial states across brain regions
- compares disease-associated inflammatory states against homeostatic reference states
- tests sex-shared, female-specific, male-specific, and sex-dimorphic gene programs
- examines overlap between HPC and PFC signatures
- performs functional enrichment and network-level interpretation
- generates publication-quality figures and summary tables

## Repository layout

The repository currently preserves the original analysis scripts while adding a cleaner project organization for ongoing work.

```text
scMicroglia/
├── README.md
├── requirements.txt
├── scripts/
│   ├── fig1.py
│   ├── fig2.py
│   ├── fig3.py
│   ├── fig4_submission.py
│   ├── build_supp_figs_S1-S5.py
│   ├── figdata.py
│   ├── figstyle.py
│   ├── fig3_enrichGO_recompute.R
│   ├── tables_T1_3.py
│   ├── tables_T4.py
│   └── README.md
├── R_scripts/
│   ├── Figure1.R
│   ├── Figure2.R
│   ├── Figure3.R
│   ├── Figure4.R
│   ├── Supplementary_Figure_S1.R
│   ├── Supplementary_Figure_S2.R
│   ├── Supplementary_Figure_S3.R
│   ├── Supplementary_Figure_S4.R
│   ├── Supplementary_Figure_S5.R
│   ├── Supplementary_Figure_S6.R
│   ├── theme_scm.R
│   └── README.md
├── analysis/
│   ├── python/
│   │   └── README.md
│   ├── r/
│   │   └── README.md
│   └── README.md
├── data/
│   ├── raw/
│   │   └── README.md
│   ├── processed/
│   │   └── README.md
│   └── README.md
├── results/
│   ├── figures/
│   │   └── README.md
│   ├── tables/
│   │   └── README.md
│   └── README.md
├── docs/
│   ├── workflow.md
│   └── README.md
└── out/
    └── generated figures and exports from local runs
```

## Primary analysis modules

### Python figure generation

The `scripts/` directory contains the primary Python-based plotting pipeline:

- `fig1.py`: study design, DEG categories, overlap analysis, and core gene heatmap
- `fig2.py`: region- and sex-specific biological interpretation
- `fig3.py`: enrichment and mechanistic analyses
- `fig4_submission.py`: final submission-ready figures
- `build_supp_figs_S1-S5.py`: supplementary analysis figures
- `figdata.py` and `figstyle.py`: shared data and plotting utilities
- `tables_T1_3.py` and `tables_T4.py`: summary tables

### R graphics workflow

The `R_scripts/` directory contains the R plotting workflow used for consistent manuscript styling and figure rendering:

- `theme_scm.R`: shared color palette, typography, and export helpers
- `Figure1.R` to `Figure4.R`: main paper figures
- `Supplementary_Figure_S*.R`: supplementary panels

## Recommended data organization

For reproducibility and long-term maintenance, raw inputs should live under `data/raw/`, while processed and analysis-ready tables should be placed in `data/processed/`.

Suggested structure:

- `data/raw/`: original single-cell matrices, metadata, and external annotations
- `data/processed/`: cleaned cell-level objects, DE tables, marker lists, and intermediate datasets
- `analysis/`: scripts and notebooks for the computational workflow
- `results/`: final generated figures, summary tables, and manuscript outputs
- `out/`: runtime-generated figures from local script execution

## Setup and dependencies

The project relies on Python and R packages for data processing, plotting, and export. The Python dependencies are listed in `requirements.txt`.

### Python environment

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

### R dependencies

The R scripts use `ggplot2`, `patchwork`, `ggrepel`, `scales`, and `ragg`.

```r
install.packages(c("ggplot2", "patchwork", "ggrepel", "scales", "ragg"))
```

## Execution workflow

1. Place raw input data in `data/raw/`.
2. Run the preprocessing and analysis steps required for your datasets.
3. Store cleaned/intermediate outputs in `data/processed/`.
4. Use scripts in `scripts/` and `R_scripts/` to reproduce figures.
5. Save final figures and tables in `results/figures/` and `results/tables/`.
6. Preserve raw data separately from generated outputs to maintain a clean, auditable pipeline.

## Project conventions

- Keep raw data outside of tracked outputs when file sizes are large.
- Store processed tables and reusable analytical objects in `data/processed/`.
- Prefer script-driven reproducibility over notebook-only workflows.
- Update `docs/workflow.md` whenever the analysis pipeline changes.
- Keep figure-generation scripts modular and version-aware.

## Scientific context

The repository focuses on microglial states in Alzheimer’s disease and examines how inflammatory disease-associated signatures vary by sex and brain region. The analysis integrates differential expression, gene-set enrichment, and regulatory network interpretation to identify core microglial programs relevant to disease biology.

## Notes

- The repository currently contains the analytical code and figure-generation pipeline but not necessarily the full raw single-cell dataset.
- The project structure is designed to support both reproducible manuscript analysis and transparent collaboration.
- The existing `scripts/` and `R_scripts/` folders remain in place for compatibility with the current workflow.

## License

This repository does not currently include a license file. If the project will be shared publicly or submitted for publication, add an explicit license before release.

## Contact and maintenance

This project is best maintained by keeping the analytical workflow organized, documenting generated outputs, and updating the README whenever major analysis steps are added or revised.

---

For detailed execution steps and project conventions, see `docs/workflow.md`.
