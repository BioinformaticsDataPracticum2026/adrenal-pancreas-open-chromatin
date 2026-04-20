# Adrenal-Pancreas Open Chromatin Comparative Project

Comparative epigenomics workflow for identifying conserved and species-specific open chromatin regions (OCRs) between human and mouse pancreas, then connecting those regions to biological processes, regulatory classes, and candidate transcription factors.

## Project Scope

This repository currently contains:

1. Task 2: Cross-species OCR mapping (HAL liftover + reciprocal best-hit pairing)
2. Task 3: Conserved and species-specific OCR partitioning from Task 2 outputs
3. Task 4: GO biological process enrichment using rGREAT
4. Task 5: Promoter vs enhancer classification (TSS +/- 2 kb rule)
5. Task 6: Motif analysis support scripts and result parsing

Phase 1 QC deliverables are stored in results/qc/ and related notes are in docs/.

## Repository Organization

```text
.
├── config/                      # Runtime configuration files
├── docs/                        # Methods, progress reports, project notes
├── scripts/                     # Executable scripts for Tasks 2-6
├── src/cross_species_ocr/       # Python package for Task 2 mapping pipeline
├── tests/                       # Unit tests
├── results/                     # Analysis outputs (tables, figures, mapping, GO, etc.)
├── data/                        # Local copy/cache area (ignored by git)
├── README.md
├── requirements.txt
└── .gitignore
```

## Environment and Dependencies

### Python setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### External tools

1. halLiftover for Task 2 mapping
2. bedtools for Task 3 processing and intersections
3. R + Bioconductor packages (rGREAT, GenomicRanges, rtracklayer, GenomeInfoDb) for Tasks 4-5
4. HOMER for Task 6 motif analysis

## Clone Repository

```bash
git clone https://github.com/BioinformaticsDataPracticum2026/adrenal-pancreas-open-chromatin
cd adrenal-pancreas-open-chromatin
```

## Quick Start

### Task 2: Cross-species mapping

1. Review config/task2_pancreas_mapping.yaml.
2. Discover inputs:

```bash
PYTHONPATH=src python3 scripts/run_task2_pancreas_mapping.py discover --config config/task2_pancreas_mapping.yaml
```

3. Run full mapping:

```bash
PYTHONPATH=src python3 scripts/run_task2_pancreas_mapping.py run --config config/task2_pancreas_mapping.yaml
```

4. Optional preprocessing-only run:

```bash
PYTHONPATH=src python3 scripts/run_task2_pancreas_mapping.py run --config config/task2_pancreas_mapping.yaml --skip-mapping
```

### Task 3: Conserved/species-specific OCRs

```bash
bash scripts/task_3_compare_ocrs_v3.sh
```

### Task 4: GO enrichment (rGREAT)

```bash
Rscript scripts/task4_GO.R
```

### Task 5: Promoter/enhancer classification

```bash
Rscript scripts/step5_promoter_enhancer.R
```

### Task 6: Motif analysis support

```bash
python3 scripts/task6_pre.py
python3 scripts/task6_analize.py
```

## Key Outputs

### Task 2 mapping outputs

1. results/mapping/orthologous_ocr_pairs.tsv
2. results/mapping/human_non_orthologous_ocr.tsv
3. results/mapping/mouse_non_orthologous_ocr.tsv
4. results/tables/task2_mapping_summary.tsv
5. results/figures/task2_mapping_counts.png
6. results/figures/task2_mapping_rates.png

### Task 3 outputs

1. results/mapping/conserved_human_in_mouse.bed
2. results/mapping/conserved_mouse_in_human.bed
3. results/mapping/human_specific.bed
4. results/mapping/mouse_specific.bed

### Task 4 outputs

Stored under results/task_4_go_analysis/.

### Task 5 outputs

Stored under results/task_5_enhancer_promoter/.

### Task 6 outputs

Stored under results/task6/.

## Documentation

1. docs/task2_cross_species_mapping.md
2. docs/task2_progress_report.md
3. docs/methods.md
4. docs/repo_conventions.md

## Reproducibility Notes

1. Use config-driven paths for Task 2 via config/task2_pancreas_mapping.yaml.
2. Keep large raw/intermediate data out of version control.
3. Write derived outputs to results/ with stable filenames.
4. Record software versions (HAL, bedtools, R packages, HOMER) in methods documentation.

## Known Limitations

1. Some Task 4/5/6 scripts currently include user-specific absolute paths and may need edits before rerun on another account.
2. Task 2 ortholog calling is conservative due to reciprocal best-hit strategy.
3. Temporary liftover artifacts can be large; keep results/mapping/tmp/ out of commits.

## Team Usage

1. Clone and create an environment.
2. Review and adjust config files as needed.
3. Run task scripts in order (Task 2 -> Task 3 -> Task 4/5 -> Task 6).
4. Keep method updates documented in docs/ and submit figures/tables from results/.
