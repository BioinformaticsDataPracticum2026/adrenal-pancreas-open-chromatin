# Comparative Epigenomics of Open Chromatin in Human and Mouse: Adrenal Gland and Pancreas

## Project Overview

This project investigates the conservation and divergence of regulatory DNA (open chromatin regions, OCRs) between human and mouse in two tissues: **adrenal gland** and **pancreas**. Using ATAC-seq–derived OCRs, we aim to:

- Evaluate dataset quality and select the best tissue for downstream analysis
- Map regulatory regions across species
- Identify conserved vs species-specific OCRs
- Interpret the biological functions regulated by these regions
- Compare promoter-like and enhancer-like OCRs
- Discover transcription factors (TFs) that may drive tissue-specific regulatory programs

The project is structured as a reproducible, multi-step computational pipeline, with clear division of labor and deliverables for each team member.

---

## Directory Structure

```
project_root/
│
├── config/                  # Configuration files (YAML, etc.)
│   └── project_config.example.yaml
│
├── data/                    # Raw and processed data
│   ├── Alignments/          # Multi-species alignments (HAL, etc.)
│   ├── CIS-BP_2.00/         # TF motif data
│   ├── external/            # External resources
│   ├── HumanAtac/           # Human ATAC-seq data
│   ├── HumanGenomeInfo/     # Human genome info
│   ├── interim/             # Intermediate files
│   ├── metadata/            # Metadata files
│   ├── MouseAtac/           # Mouse ATAC-seq data
│   ├── MouseGenomeInfo/     # Mouse genome info
│   ├── processed/           # Processed data
│   ├── raw/                 # Raw data
│
├── docs/                    # Project documentation
│   ├── methods.md
│   ├── project_overview.md
│   ├── repo_conventions.md
│   └── task1_qc_workflow.md
│
├── notebooks/               # Jupyter notebooks for analysis
│
├── results/                 # Output results
│   ├── annotations/
│   ├── figures/
│   ├── logs/
│   ├── mapping/
│   ├── qc/
│   ├── tables/
│
├── scripts/                 # Standalone scripts
│
├── src/                     # Source code (modules, pipeline)
│
├── tests/                   # Unit and integration tests
│
├── .github/                 # GitHub workflows, issue templates
├── .gitignore
├── main.py                  # Main entry point (if needed)
├── project.code-workspace   # VS Code workspace settings
├── README.md                # This file
└── requirements.txt         # Python dependencies
```

---

## Data Resources

- **Human ATAC-seq**: `/ocean/projects/bio230007p/ikaplow/HumanAtac`
- **Mouse ATAC-seq**: `/ocean/projects/bio230007p/ikaplow/MouseAtac`
- **Human genome info**: `/ocean/projects/bio230007p/ikaplow/HumanGenomeInfo`
- **Mouse genome info**: `/ocean/projects/bio230007p/ikaplow/MouseGenomeInfo`
- **Multi-species alignments**: `/ocean/projects/bio230007p/ikaplow/Alignments`
- **TF motif data**: `/ocean/projects/bio230007p/ikaplow/CIS-BP_2.00`

---

## Project Phases and Tasks

### Phase 1: Data Quality Evaluation

- **Goal:** Assess the quality of each ATAC-seq dataset (human/mouse, adrenal/pancreas) and select the best tissue for downstream analysis.
- **QC metrics:** Mapping rates, read length, nucleosomal patterning, number of peaks, TSS enrichment, read depth, reproducibility, IDR peaks, biological plausibility.
- **Deliverable:** Written report section with figures/tables and tissue selection rationale.

### Phase 2: Downstream Biological Analysis

#### Task 2: Build a Cross-Species OCR Map
- Map and annotate OCRs in human and mouse for the selected tissue
- Generate a unified cross-species map using tools like `halLiftover` and `HALPER`

#### Task 3: Compare OCRs Between Species
- Classify OCRs as conserved or species-specific
- Quantify conservation and divergence

#### Task 4: Biological Process Enrichment
- Link OCRs to genes
- Perform GO/pathway enrichment (preferably using GREAT/rGREAT logic)

#### Task 5: Promoter vs Enhancer Analysis
- Classify OCRs as promoter-like or enhancer-like (e.g., by distance to TSS)
- Compare conservation and gene associations

#### Task 6: Transcription Factor Analysis
- Motif enrichment/footprinting to identify TFs binding OCRs
- Integrate with conservation and regulatory class analyses

---

## Workflow Summary

1. **QC all four datasets** (human/mouse × adrenal/pancreas)
2. **Select best tissue** for downstream analysis
3. **Map OCRs across species** (build cross-species regulatory map)
4. **Classify OCRs** (conserved vs species-specific)
5. **Biological interpretation** (GO/pathway enrichment)
6. **Promoter vs enhancer comparison**
7. **TF motif/footprinting analysis**
8. **Integrate results into a final report**

---

## Best Practices

- Use version control (Git) for all code and documentation
- Document all scripts, parameters, and intermediate files
- Use config files for reproducibility
- Keep raw data immutable; all processing should be reproducible from scripts
- Organize results and figures for easy integration into the final report
- Write clear, modular code and notebooks
- Use environment management (e.g., `requirements.txt`, conda) to ensure reproducibility

---

## Getting Started

1. **Clone the repository:**
  ```bash
  git clone <repo_url>
  cd adrenal-pancreas-open-chromatin
  ```
2. **Set up the environment:**
  ```bash
  python3 -m venv env
  source env/bin/activate
  pip install -r requirements.txt
  ```
3. **Configure project settings:**
  - Copy `config/project_config.example.yaml` to `config/project_config.yaml` and edit as needed.
4. **Run analyses:**
  - Follow the documentation in `docs/` and notebooks in `notebooks/` for each project phase.

---

## Contributing

- Use feature branches and pull requests for all changes
- Write informative commit messages
- Add/modify tests as needed
- Keep documentation up to date

---

## Contact

For questions, contact the project team via the repository issues or your course communication channel.

---

## License

This project is for educational purposes as part of a comparative epigenomics course. Please cite appropriately if reusing code or analyses.
