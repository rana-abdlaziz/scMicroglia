# scMicroglia

Single-cell analysis of microglial states in Alzheimer’s disease, with a focus on sex-specific and cross-region differences between hippocampus (HPC) and prefrontal cortex (PFC).

This repository contains the analysis code, plotting scripts, and manuscript figure generation pipeline used to study microglial state markers, differential expression, enriched pathways, and transcription-factor activity.

## Project overview

The analysis workflow implemented in this codebase is built around:

- single-cell / microglia state annotation across brain regions
- within-sex differential expression contrasts comparing inflammatory AD-associated states vs homeostatic references
- sex-shared, sex-specific, and sex-dimorphic gene sets
- cross-region overlap between HPC and PFC signatures
- functional enrichment and network inference for biological interpretation
- publication-ready figure generation in Python and R

## Repository structure

```text
scMicroglia/
├── README.md
├── .gitignore
├── requirements.txt
├── docs/
│   ├── project-structure.md
│   └── workflow.md
├── data/
│   └── README.md
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
├── results/
│   └── README.md
└── out/
    └── generated figures (created at runtime)
```

## What is in the codebase

### Python scripts

The `scripts/` directory contains the main plotting and analysis scripts used to generate publication figures and tables.

- `fig1.py`: study design, DEG categories, cross-region overlap, core gene heatmap
- `fig2.py`: region- and sex-specific biological interpretation
- `fig3.py`: gene-set and regulatory analyses
- `fig4_submission.py`: final submission figure set
- `build_supp_figs_S1-S5.py`: supplementary figure generation
- `tables_T1_3.py` and `tables_T4.py`: summary tables
- `figdata.py` and `figstyle.py`: shared data and style helpers
- `fig3_enrichGO_recompute.R`: R-based enrichment recomputation used as part of the workflow

### R scripts

The `R_scripts/` directory contains figure generation code in R using `ggplot2` and associated packages.

- `theme_scm.R`: shared plotting theme, fonts, colors, and export utilities
- `Figure1.R` to `Figure4.R`: primary manuscript figures
- `Supplementary_Figure_S*.R`: supplementary figure panels

## Suggested project workflow

1. Prepare input data and metadata in `data/`.
2. Run the Python figure scripts in `scripts/` to generate figures and summary tables.
3. Use the R scripts in `R_scripts/` for publication-style figure panels and consistent formatting.
4. Save generated outputs under `results/` or `out/`.
5. Keep raw data separate from analyzed results to preserve reproducibility.

## Dependencies

### Python

Install the core packages used by the scripts:

```bash
python -m pip install numpy pandas matplotlib
```

Additional packages may be needed depending on local analysis inputs or if you extend the plotting pipeline.

### R

Install the plotting dependencies used in `R_scripts/theme_scm.R`:

```r
install.packages(c("ggplot2", "patchwork", "ggrepel", "scales", "ragg"))
```

## Data handling recommendation

This repository currently focuses on the analytical code and plotting logic. For a clean single-cell project workflow, use the following structure:

- `data/raw/`: original raw datasets and metadata
- `data/processed/`: cleaned matrices, annotations, and analysis tables
- `results/`: final figure exports and summary outputs
- `out/`: generated plots from scripts when run locally

The current repository does not contain the raw single-cell count matrices, which are typically kept outside the Git repository for size and privacy reasons.

## Notes about the biological analysis

The project centers on microglial states in Alzheimer’s disease, comparing disease-associated inflammatory states with homeostatic reference states across the HPC and PFC.

Key analyses in the code include:

- sex-shared vs sex-specific differential expression
- cross-region overlap of candidate genes
- functional enrichment of core gene programs
- transcription-factor network and activity analysis
- figure-ready summaries for manuscript and supplementary panels

## Repository maintenance

To keep this project organized over time:

- keep raw data out of version control
- store processed tables under `data/processed/`
- keep figure-generation scripts in their current folders
- save outputs in `results/` or `out/`
- document any new analysis modules in `docs/`

## License

This repository currently does not specify a license file. If publication or distribution is intended, add an explicit license before public release.

## Contribution guidance

When adding new analysis modules:

- place reusable functions in a dedicated utility file or script
- keep figure scripts modular and script-driven
- store generated plots outside the repo if they are large or generated outputs
- update this README when adding major analytical steps
