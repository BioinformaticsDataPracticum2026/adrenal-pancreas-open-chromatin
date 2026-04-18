# Cross-Species OCR Mapping Pipeline for Pancreas

This repository now contains a reproducible Task 2 pipeline for building a cross-species open chromatin region (OCR) map between human and mouse pancreas ATAC-seq peak sets.

The implementation is intentionally scoped to Task 2 only. It prepares standardized OCR peak files, performs optional lightweight annotation, runs HAL-based cross-species liftover, classifies orthologous versus non-orthologous OCRs, and writes downstream-friendly tables and figures for later team tasks.

## Current Scope

Implemented here:
- Locate human and mouse pancreas peak files from the course data root
- Standardize peak files into clean OCR tables and BED outputs
- Optionally annotate OCRs with nearest gene and promoter/distal labels
- Lift OCRs between species using a HAL alignment via `halLiftover`
- Build orthologous and species-specific OCR outputs
- Write summary tables, figures, logs, and reproducible metadata

Not implemented here:
- Task 3 conserved-vs-species-specific biological interpretation
- Task 4 GO or pathway enrichment
- Task 5 promoter-vs-enhancer downstream analysis beyond simple labels
- Task 6 TF motif analysis

## Default Data Inputs

The default pancreas inputs are configured to use the reproducible IDR optimal peak calls:

- Human: `/ocean/projects/bio230007p/ikaplow/HumanAtac/Pancreas/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz`
- Mouse: `/ocean/projects/bio230007p/ikaplow/MouseAtac/Pancreas/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz`
- HAL alignment: `/ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal`

HAL species names in the provided alignment are configured as `Human` and `Mouse`.

## Repository Layout

```
config/
  task2_pancreas_mapping.yaml
docs/
  task2_cross_species_mapping.md
scripts/
  run_task2_pancreas_mapping.py
src/cross_species_ocr/
  __init__.py
  annotation.py
  cli.py
  config.py
  intervals.py
  logging_utils.py
  mapping.py
  peaks.py
  pipeline.py
  reporting.py
tests/
  test_standardize.py
results/
  mapping/
  tables/
  figures/
  logs/
```

## Environment

Python dependencies are listed in `requirements.txt`. The mapping step also requires an external HAL executable:

- `halLiftover`

If `halLiftover` is not already on `PATH`, load the appropriate module or supply an explicit binary path in the YAML config.

## How To Run

1. Create or activate your Python environment.
2. Adjust `config/task2_pancreas_mapping.yaml` only if you need non-default inputs or output paths.
3. Run input discovery:

```bash
PYTHONPATH=src python3 scripts/run_task2_pancreas_mapping.py discover --config config/task2_pancreas_mapping.yaml
```

4. Run the full Task 2 pipeline:

```bash
PYTHONPATH=src python3 scripts/run_task2_pancreas_mapping.py run --config config/task2_pancreas_mapping.yaml
```

5. If you want to validate preprocessing before HAL mapping is available:

```bash
PYTHONPATH=src python3 scripts/run_task2_pancreas_mapping.py run --config config/task2_pancreas_mapping.yaml --skip-mapping
```

## Main Outputs

- `results/mapping/human_pancreas_ocr.processed.bed`
- `results/mapping/mouse_pancreas_ocr.processed.bed`
- `results/mapping/human_pancreas_ocr.processed.tsv`
- `results/mapping/mouse_pancreas_ocr.processed.tsv`
- `results/mapping/orthologous_ocr_pairs.tsv`
- `results/mapping/human_non_orthologous_ocr.tsv`
- `results/mapping/mouse_non_orthologous_ocr.tsv`
- `results/tables/task2_mapping_summary.tsv`
- `results/figures/task2_mapping_counts.png`
- `results/figures/task2_mapping_rates.png`
- `results/logs/task2_pancreas_mapping.log`

## Design Notes

The pipeline keeps Task 2 outputs easy to reuse later by preserving stable OCR IDs, genomic coordinates, optional nearest-gene annotations, and explicit mapping-status labels.

## Task 4: GO biological process enrichment of species-specific and conserved open chromatin regions

### Overview
This task performs GO Biological Process enrichment analysis on species-specific and conserved open chromatin regions using `rGREAT`, a region-based enrichment framework designed for genomic intervals rather than pre-defined gene lists.

The goal is to identify biological processes associated with:
- mouse-specific OCRs
- human-specific OCRs
- conserved human OCRs mapped into mouse coordinates

### Analysis objective
We aimed to test whether different classes of open chromatin regions are preferentially associated with genes involved in specific biological processes.

Because these inputs are genomic regions rather than genes, we used `rGREAT`, which links regulatory regions to nearby or regulatory-domain-associated genes and evaluates GO enrichment in a region-centric manner.

---

## Input files

### Foreground peak sets
The foreground BED files are located in:

- `results/mapping/mouse_specific.bed`
- `results/mapping/human_specific.bed`
- `results/mapping/conserved_human_in_mouse.bed`

These represent:
- `mouse_specific.bed`: mouse-specific open chromatin regions
- `human_specific.bed`: human-specific open chromatin regions
- `conserved_human_in_mouse.bed`: conserved human open chromatin regions projected into mouse coordinates

### Background peak universes
The background OCR universes are located in:

- `results/mapping/mouse_pancreas_ocr.processed.bed`
- `results/mapping/human_pancreas_ocr.processed.bed`

These files represent the total experimentally observed OCR peak universes for each species and are used as species-matched backgrounds in `rGREAT`.

### Excluded file
The following file was not used in the final GO enrichment analysis:

- `results/mapping/conserved_mouse_in_human.bed`

This reciprocal conserved projection showed inconsistent promoter/enhancer composition after annotation and was therefore excluded from downstream GO enrichment to avoid introducing potentially unreliable signals.

---

## Background strategy

### Background used for each foreground set
We used species-matched OCR universes as background:

- `mouse_specific.bed` was tested against `results/mapping/mouse_pancreas_ocr.processed.bed`
- `human_specific.bed` was tested against `results/mapping/human_pancreas_ocr.processed.bed`
- `conserved_human_in_mouse.bed` was tested against `results/mapping/mouse_pancreas_ocr.processed.bed`

This design ensures that enrichment is evaluated relative to accessible regulatory regions detected in the corresponding species background.

---

## Software requirements

### R version
This analysis was designed to run in R 4.5.1

### Required R packages

#### Bioconductor packages
- `rGREAT`
- `rtracklayer`
- `GenomicRanges`
- `GenomeInfoDb`

#### CRAN packages
- `ggplot2`
- `dplyr`
- `stringr`

## Task 5: Promoter- and Enhancer-Associated OCR Classification

### Objective
To compare regulatory composition across species-specific and conserved open chromatin regions by classifying peaks into promoter-associated versus enhancer-associated categories.

### Input Peak Sets
- `results/mapping/mouse_specific.bed`
- `results/mapping/human_specific.bed`
- `results/mapping/conserved_human_in_mouse.bed`

### Excluded Set
- `results/mapping/conserved_mouse_in_human.bed` was excluded from final analysis due to inconsistent promoter/enhancer composition in prior QC.

### Reference TSS Annotations
- Human: `/ocean/projects/bio230007p/ikaplow/HumanGenomeInfo/gencode.v27.annotation.protTranscript.TSSsWithStrand_sorted.bed`
- Mouse: `/ocean/projects/bio230007p/ikaplow/MouseGenomeInfo/gencode.vM15.annotation.protTranscript.geneNames_TSSWithStrand_sorted.bed`

---

## Computational Method

### Core Genomic Framework
This analysis was implemented using the **Bioconductor GenomicRanges framework**.

- Peaks were imported as genomic intervals and represented as `GRanges` objects.
- TSS annotations were converted to stranded `GRanges`.
- Interval operations used GenomicRanges methods, including:
  - `keepStandardChromosomes()`
  - `sort()`
  - `promoters()`
  - `trim()`
  - `reduce()`
  - `overlapsAny()`

### Promoter and Enhancer Definition
- **Promoter-associated peaks**: peaks overlapping TSS ± 2 kb windows.
- **Candidate enhancer-associated peaks**: peaks not overlapping promoter windows.

This is an operational distance-based definition; enhancer labels here indicate non-promoter/distal OCRs rather than histone-mark-validated enhancer states.

### Species-Matched Assignment
- `mouse_specific.bed` and `conserved_human_in_mouse.bed` were annotated against mouse promoter windows.
- `human_specific.bed` was annotated against human promoter windows.

---

## Software
R packages used:
- `rtracklayer`
- `GenomicRanges`
- `GenomeInfoDb`
- `dplyr`
- `ggplot2`
- `scales`
- `patchwork`
- `tidyr`

---

## Main Script
- `scripts/task5_promoter_enhancer.R`

## Run Command (Bridges2)
```bash
cd /jet/home/wzhang37/step5_enhancer_promoter
module load r
Rscript task5_promoter_enhancer.R

## Outputs

- `results/task_5_enhancer_promoter/step5_peak_assignment.csv`
- `results/task_5_enhancer_promoter/step5_count_summary.csv`
- `results/task_5_enhancer_promoter/step5_proportion_summary.csv`
- `results/task_5_enhancer_promoter/proportion_count_panel.png`
- `results/task_5_enhancer_promoter/enhancer_promoter_ratio.png`
- `results/task_5_enhancer_promoter/all_three_panel.png`

## Summary

Using a GenomicRanges-based overlap framework, OCRs were classified into promoter-associated and candidate enhancer-associated groups via TSS-proximity rules. This enabled consistent comparison of promoter/distal composition across mouse-specific, human-specific, and conserved-human-in-mouse peak sets.