from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def identify_ortholog_column_schema(df: pd.DataFrame) -> tuple[list[str], list[str]]:
	"""Identify and extract column names from ortholog table.
	
	Supports two possible schemas:
	1. Query/target format (current pipeline):
	   query_chrom, query_start, query_end, query_peak_id,
	   target_chrom, target_start, target_end, target_peak_id
	2. Species format (legacy):
	   human_chr, human_start, human_end, human_peak_id,
	   mouse_chr, mouse_start, mouse_end, mouse_peak_id
	
	Args:
		df: DataFrame loaded from orthologous_ocr_pairs.tsv.
	
	Returns:
		Tuple of (human_column_names, mouse_column_names).
	
	Raises:
		ValueError: If schema does not match either expected format.
	"""
	# Try modern query/target column names first
	current_schema_cols = {"query_chrom", "query_start", "query_end", "query_peak_id",
						  "target_chrom", "target_start", "target_end", "target_peak_id"}
	if current_schema_cols.issubset(df.columns):
		human_cols = ["query_chrom", "query_start", "query_end", "query_peak_id"]
		mouse_cols = ["target_chrom", "target_start", "target_end", "target_peak_id"]
		return human_cols, mouse_cols

	# Try legacy species-specific column names
	legacy_schema_cols = {"human_chr", "human_start", "human_end", "human_peak_id",
					   "mouse_chr", "mouse_start", "mouse_end", "mouse_peak_id"}
	if legacy_schema_cols.issubset(df.columns):
		human_cols = ["human_chr", "human_start", "human_end", "human_peak_id"]
		mouse_cols = ["mouse_chr", "mouse_start", "mouse_end", "mouse_peak_id"]
		return human_cols, mouse_cols

	# No recognized schema found
	raise ValueError(
		f"Unrecognized ortholog table schema. Found columns: {sorted(df.columns)}"
	)


def main() -> int:
	"""Main entry point.
	
	Returns:
		0 on success, non-zero on error.
	"""
	parser = argparse.ArgumentParser(description="Prepare Task 6 ortholog BED files from Task 2 orthologous pairs.")
	parser.add_argument(
		"--ortholog-table",
		default="results/mapping/orthologous_ocr_pairs.tsv",
		help="Path to orthologous_ocr_pairs.tsv from Task 2.",
	)
	parser.add_argument(
		"--outdir",
		default="results/task6",
		help="Output directory for BED files.",
	)
	args = parser.parse_args()

	ortholog_table = Path(args.ortholog_table)
	outdir = Path(args.outdir)
	outdir.mkdir(parents=True, exist_ok=True)

	if not ortholog_table.exists():
		print(f"Error: Ortholog table not found: {ortholog_table}", file=__import__('sys').stderr)
		return 1

	try:
		df = pd.read_table(ortholog_table)
		human_cols, mouse_cols = identify_ortholog_column_schema(df)

		human_ortho = df[human_cols].drop_duplicates().copy()
		mouse_ortho = df[mouse_cols].drop_duplicates().copy()

		human_out = outdir / "human_ortho_unique.bed"
		mouse_out = outdir / "mouse_ortho_unique.bed"

		human_ortho.to_csv(human_out, sep="\t", index=False, header=False)
		mouse_ortho.to_csv(mouse_out, sep="\t", index=False, header=False)

		print(f"Wrote {len(human_ortho)} human conserved OCRs to {human_out}")
		print(f"Wrote {len(mouse_ortho)} mouse conserved OCRs to {mouse_out}")
		return 0
	except ValueError as e:
		print(f"Error: {e}", file=__import__('sys').stderr)
		return 1


if __name__ == "__main__":
	raise SystemExit(main())