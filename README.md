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

## Task 4: GO biological process enrichment of species-specific open chromatin regions

### Goal
Identify biological processes associated with mouse-specific and human-specific open chromatin regions using a region-centric enrichment framework.

### Input peak sets
- `mouse_specific.bed`
- `human_specific.bed`

### Background peak universes
- `mouse_pancreas_ocr.processed.bed`
- `human_pancreas_ocr.processed.bed`

### Method
We performed GO Biological Process enrichment analysis using `rGREAT`, which is designed for genomic region-based input.

#### Foreground sets
- Mouse-specific open chromatin regions: `mouse_specific.bed`
- Human-specific open chromatin regions: `human_specific.bed`

#### Background strategy
We used species-matched OCR peak universes rather than the whole genome as background, because the goal was to evaluate functional enrichment within experimentally observed regulatory regions rather than across all genomic bases.

Specifically:
- Mouse-specific regions were tested against `mouse_pancreas_ocr.processed.bed`
- Human-specific regions were tested against `human_pancreas_ocr.processed.bed`

This background design provides a more biologically appropriate comparison framework by restricting enrichment to accessible regulatory regions detected in each species.

#### Note on conserved peak sets
Although conserved peak files were generated during earlier cross-species comparisons, they were not included in the final GO enrichment analysis. This decision was made because the conserved peak projections showed inconsistent regulatory composition after promoter/enhancer stratification, suggesting that they were not reliable for downstream functional enrichment.

### Main script
- `scripts/task4_GO.R`

### Output
Task 4 results are stored in:

- `results/task_4_go_analysis/`

This directory includes:
- GO BP enrichment result tables for `mouse_specific` and `human_specific`
- significant-term summary tables
- a combined significant-term table for species-specific peak sets
- dotplot visualization comparing mouse-specific and human-specific enrichment results

### Summary
The enrichment analysis identified biological processes associated with mouse-specific and human-specific open chromatin regions, highlighting distinct regulatory programs in the two species within pancreas OCR landscapes.


## Task 5: Compare candidate enhancers and candidate promoters

### Goal
Distinguish candidate promoter-associated and candidate enhancer-associated open chromatin regions, and compare their distributions across species-specific and conserved peak sets.

### Input peak sets
- `mouse_specific.bed`
- `human_specific.bed`
- `conserved_human_in_mouse.bed`
- `conserved_mouse_in_human.bed`

### Reference annotations
We used the instructor-provided TSS annotation resources:

- Human TSS annotation:  
  `/ocean/projects/bio230007p/ikaplow/HumanGenomeInfo/gencode.v27.annotation.protTranscript.TSSsWithStrand_sorted.bed`

- Mouse TSS annotation:  
  `/ocean/projects/bio230007p/ikaplow/MouseGenomeInfo/gencode.vM15.annotation.protTranscript.geneNames_TSSsWithStrand_sorted.bed`

### Method
To distinguish candidate promoters and candidate enhancers, we used annotated transcription start sites as the reference.

#### Promoter definition
Promoter-associated peaks were defined as peaks overlapping a **±2 kb window** around annotated TSSs.

#### Enhancer definition
Peaks that did **not** overlap promoter windows were classified as candidate enhancers.

#### Species-specific assignment
- `mouse_specific.bed` and `conserved_human_in_mouse.bed` were compared against **mouse TSS/promoter annotations**
- `human_specific.bed` and `conserved_mouse_in_human.bed` were compared against **human TSS/promoter annotations**

### Main script
- `scripts/step5_promoter_enhancer.R`

### Output
Task 5 results are stored in:

- `results/task_5_enhancer_promoter/`

This directory includes:
- `step5_peak_assignment.csv`
- `step5_count_summary.csv`
- `step5_proportion_summary.csv`
- `step5_promoter_enhancer_two_panel.png`
- `step5_promoter_enhancer_two_panel.pdf`

### Summary
Using annotated TSS windows, we classified open chromatin peaks into candidate promoter-associated and candidate enhancer-associated regions. This allowed us to compare promoter-versus-enhancer composition across mouse-specific, human-specific, and conserved peak sets.
