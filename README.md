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
- [Task 3 script](./scripts/task_3_compare_ocrs_v3.sh)
  

### Task 4: GO enrichment (rGREAT)

```bash
Rscript scripts/task4_GO.R
```
- [Task 4 script](./scripts/task4_GO.R)


### Task 5: Promoter/enhancer classification

```bash
Rscript scripts/task5_promoter_enhancer.R
```
- [Task 5 script](./scripts/task5_promoter_enhancer.R)

### Task 6: Motif analysis support

```bash
python3 scripts/task6_pre.py
python3 scripts/task6_analize.py
```
- [Task 6 script](./scripts/task6_analyze.py)
## Key Inputs & Outputs

### Task 2 mapping outputs

1. results/mapping/orthologous_ocr_pairs.tsv
2. results/mapping/human_non_orthologous_ocr.tsv
3. results/mapping/mouse_non_orthologous_ocr.tsv
4. results/tables/task2_mapping_summary.tsv
5. results/figures/task2_mapping_counts.png
6. results/figures/task2_mapping_rates.png

### Task 3 
#### Inputs:

1. results/mapping/human_pancreas_ocr.processed.bed
2. results/mapping/mouse_pancreas_ocr.processed.bed
3. results/mapping/tmp/human_to_mouse.lifted.bed
4. results/mapping/tmp/mouse_to_human.lifted.bed
5. results/mapping/orthologous_ocr_pairs.tsv
#### Outputs:

1. results/mapping/tmp/human_to_mouse.lifted.merged.bed
2. results/mapping/tmp/mouse_to_human.lifted.merged.bed
3. results/mapping/conserved_human_in_mouse.bed
4. results/mapping/conserved_mouse_in_human.bed
5. results/mapping/human_specific.bed
6. results/mapping/mouse_specific.bed

### Task 4
#### Inputs:

1. results/mapping/mouse_specific.bed
2. results/mapping/human_specific.bed
3. results/mapping/conserved_human_in_mouse.bed
4. results/mapping/mouse_pancreas_ocr.processed.bed
5. results/mapping/human_pancreas_ocr.processed.bed
#### Outputs:

1. results/task_4_go_analysis/mouse_specific_rGREAT_GO_BP.csv
2. results/task_4_go_analysis/mouse_specific_rGREAT_GO_BP_sig.csv
3. results/task_4_go_analysis/human_specific_rGREAT_GO_BP.csv
4. results/task_4_go_analysis/human_specific_rGREAT_GO_BP_sig.csv
5. results/task_4_go_analysis/conserved_human_in_mouse_rGREAT_GO_BP.csv
6. results/task_4_go_analysis/conserved_human_in_mouse_rGREAT_GO_BP_sig.csv
7. results/task_4_go_analysis/combined_mouse_human_conserved_sig_terms.csv
- [Task 4 outputs](./results/task_4_go_analysis/combined_top10_dotplots.png)
### Task 5
#### Inputs:
1. results/mapping/mouse_specific.bed
2. results/mapping/human_specific.bed
3. results/mapping/conserved_human_in_mouse.bed
4. results/mapping/human_tss.bed
5. results/mapping/mouse_tss.bed

#### Outputs:
1. results/task_5_enhancer_promoter/step5_peak_assignment.csv
2. results/task_5_enhancer_promoter/step5_count_summary.csv
3. results/task_5_enhancer_promoter/step5_proportion_summary.csv
- [Task 5 outputs](./results/task_5_enhancer_promoter/all_three_panel.png)

### Task 6
#### Inputs:
1. results/task6/human_ortho_unique.bed
2. results/task6/mouse_ortho_unique.bed
3. /ocean/projects/bio230007p/ikaplow/HumanGenomeInfo/hg38.fa
4. /ocean/projects/bio230007p/ikaplow/MouseGenomeInfo/mm10.fa

#### Outputs:
1. results/task6/homer_human_promoters/
2. results/task6/homer_human_enhancers/
3. results/task6/homer_mouse_promoters/
4. results/task6/homer_mouse_enhancers/
5. results/task6/task6_final_result.png
Stored under results/task6/.
- [Task 6 outputs](.results/task6)

## Reproducibility Notes

1. Use config-driven paths for Task 2 via config/task2_pancreas_mapping.yaml.
2. Keep large raw/intermediate data out of version control.
3. Write derived outputs to results/ with stable filenames.
4. Record software versions (HAL, bedtools, R packages, HOMER) in methods documentation.

## Known Limitations

1. Task 6 scripts currently include user-specific absolute paths and may need edits before rerun on another account.
2. Task 2 ortholog calling is conservative due to reciprocal best-hit strategy.
3. Temporary liftover artifacts can be large; keep results/mapping/tmp/ out of commits.

## Team Usage

1. Clone and create an environment.
2. Review and adjust config files as needed.
3. Run task scripts in order (Task 2 -> Task 3 -> Task 4/5 -> Task 6).
4. Keep method updates documented in docs/ and submit figures/tables from results/.

## Demo Videos

### Task 2
- Video: [Task 2 Demo (Cross-Species Mapping)](https://youtu.be/-IAszN7KlHI?si=o2VrlK41F1mDozUq)

### Task 3
- Video: [Task 3 Demo (OCR Comparison)](https://youtu.be/ZyHVatZ3sHM?si=M-0xn3bfwic-X8yj)

### Task 4
- Video: [Task 4 Demo (GO BP Enrichment)](https://youtu.be/CB8uWG67Mug)

### Task 5
- Video: [Task 5 Demo (Promoter vs Enhancer Classification)](https://youtu.be/V_zYp4Uemik)

### Task 6
- Video: [Task 6 Demo (ranscription Factor Motif Enrichment)](https://www.youtube.com/watch?v=pofi3isRvKM)

## Tools Used 

This pipeline uses [bedtools](https://bedtools.readthedocs.io/) for genomic interval 
operations, including merging lifted peak fragments (`bedtools groupby`) and 
identifying species-specific OCRs via reciprocal overlap filtering 
(`bedtools intersect`).

> Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing 
> genomic features. *Bioinformatics*. 2010;26(6):841-842. 
> [doi:10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033)


This pipeline uses [rGREAT](https://bioconductor.org/packages/rGREAT/) for 
region-based GO Biological Process enrichment of OCR sets. 

> McLean CY, Bristor D, Hiller M, Clarke SL, Schaar BT, Lowe CB, Wenger AM, 
> Bejerano G. GREAT improves functional interpretation of cis-regulatory regions. 
> *Nature Biotechnology*. 2010;28(5):495-501. 
> [doi:10.1038/nbt.1630](https://doi.org/10.1038/nbt.1630)

This pipeline uses [GenomicRanges](https://bioconductor.org/packages/GenomicRanges/) 
for genomic interval representation and overlap operations in R. If you use this 
pipeline, please cite:

> Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M, Gentleman R, Morgan MT, 
> Carey VJ. Software for computing and annotating genomic ranges. 
> *PLoS Computational Biology*. 2013;9(8):e1003118. 
> [doi:10.1371/journal.pcbi.1003118](https://doi.org/10.1371/journal.pcbi.1003118)

This pipeline uses [rtracklayer](https://bioconductor.org/packages/rtracklayer/) 
to import and export genome annotation files (e.g., BED) in R. 

> Lawrence M, Gentleman R, Carey V. rtracklayer: an R package for interfacing 
> with genome browsers. *Bioinformatics*. 2009;25(14):1841-1842. 
> [doi:10.1093/bioinformatics/btp328](https://doi.org/10.1093/bioinformatics/btp328)

This pipeline uses the [Gene Ontology](https://geneontology.org/) resource for 
Biological Process term definitions and annotation. 

> The Gene Ontology Consortium. The Gene Ontology resource: enriching a GOld mine. 
> *Nucleic Acids Research*. 2021;49(D1):D325-D334. 
> [doi:10.1093/nar/gkaa1113](https://doi.org/10.1093/nar/gkaa1113)

## Authors
- Aryan Sharan Guda: aryanshg@andrew.cmu.edu
- Joe Wang: joewang@andrew.cmu.edu
- Wenxin Zhang: zhangwez@andrew.cmu.edu
- Mrunmayee Wankhede: mwankhed@andrew.cmu.edu

## How to Cite

If you use this pipeline in your work, please cite it as:

> Guda AS, Wang J, Zhang W, Wankhede M. Adrenal-Pancreas Open Chromatin Comparative Project: A pipeline for cross-species OCR mapping, functional enrichment, and motif discovery in human and mouse pancreas. 2026. https://github.com/BioinformaticsDataPracticum2026/adrenal-pancreas-open-chromatin

Please also cite the underlying tools used by this pipeline (bedtools, halLiftover, rGREAT, GenomicRanges, rtracklayer, HOMER, Gene Ontology) as listed in the [Tools Used](#tools-used) section.
