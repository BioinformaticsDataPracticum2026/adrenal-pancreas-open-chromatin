# Methods: Human Adrenal Gland QC Workflow

- **Data Inventory:** Directory scan, filename filtering, and candidate selection
- **BED/Peak Validation:** Interval checks, malformed row detection, duplicate counting
- **QC Metrics:** Interval count, width stats, chromosome distribution, score summary
- **Replicate Comparison:** Overlap and count comparison if replicates exist
- **Plotting:** Matplotlib-based histograms and barplots
- **Config:** YAML-based, all paths configurable
- **Logging:** Standardized logging for reproducibility
