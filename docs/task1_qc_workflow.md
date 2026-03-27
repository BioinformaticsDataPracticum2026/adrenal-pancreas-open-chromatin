# Phase 1 Task 1: Human Adrenal Gland QC Workflow

## Workflow Steps

1. **Data Inventory**
   - Scan configured directories for candidate human adrenal gland open chromatin files
   - Prefer alternative adrenal dataset if available
   - Write a TSV inventory of candidate files

2. **Candidate File Selection**
   - Use filename heuristics and metadata to identify likely human adrenal gland files
   - Mark TODOs where filenames need confirmation

3. **Validation**
   - Load candidate BED/peak files
   - Validate intervals (chromosome, coordinates, width, duplicates, malformed rows)

4. **QC Metrics**
   - Compute:
     - Total intervals
     - Chromosome distribution
     - Interval width summary
     - Duplicate count
     - Score column summary (if present)
     - Fraction near promoters (optional, if GTF provided)

5. **Replicate Comparison (if applicable)**
   - If multiple adrenal files/replicates are found:
     - Compare counts
     - Compute overlap statistics
     - Write summary table and plots
   - If not, skip gracefully

6. **Outputs**
   - Inventory, QC tables, plots, and logs under `results/`
   - Markdown QC report

## Reproducibility
- All paths and filenames are configurable
- No hard-coded filenames
- Outputs are versioned and organized

See `README.md` and `config/project_config.example.yaml` for setup and usage.
