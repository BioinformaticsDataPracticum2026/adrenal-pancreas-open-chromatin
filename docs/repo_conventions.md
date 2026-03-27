# Repository Conventions

## Branch Naming
- Use `phase1-task1-qc` for this phase
- Use descriptive feature branches for new work

## Commit Style
- Use concise, imperative commit messages
- Reference task/issue if applicable

## Notebook Usage
- Place exploratory and reporting notebooks in `notebooks/`
- Keep reusable logic in `src/`
- Notebooks should not contain hard-coded paths

## Output Handling
- All outputs go under `results/`
- Do not commit raw data or large intermediates
- Use `results/logs/` for logs, `results/qc/` for QC summaries, `results/figures/` for plots, `results/tables/` for tables
