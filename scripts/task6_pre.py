from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def pick_columns(df: pd.DataFrame) -> tuple[list[str], list[str]]:
	"""Support both legacy and current ortholog table schemas."""
	if {"query_chrom", "query_start", "query_end", "query_peak_id", "target_chrom", "target_start", "target_end", "target_peak_id"}.issubset(df.columns):
		human_cols = ["query_chrom", "query_start", "query_end", "query_peak_id"]
		mouse_cols = ["target_chrom", "target_start", "target_end", "target_peak_id"]
		return human_cols, mouse_cols

	if {"human_chr", "human_start", "human_end", "human_peak_id", "mouse_chr", "mouse_start", "mouse_end", "mouse_peak_id"}.issubset(df.columns):
		human_cols = ["human_chr", "human_start", "human_end", "human_peak_id"]
		mouse_cols = ["mouse_chr", "mouse_start", "mouse_end", "mouse_peak_id"]
		return human_cols, mouse_cols

	raise ValueError(
		"Unrecognized orthologous table schema. Expected query/target columns or human/mouse columns."
	)


def main() -> int:
	parser = argparse.ArgumentParser(description="Prepare Task 6 ortholog BED files from Task 2 orthologous pairs.")
	parser.add_argument(
		"--ortholog-table",
		default="results/mapping/orthologous_ocr_pairs.tsv",
		help="Path to orthologous pair table TSV.",
	)
	parser.add_argument(
		"--outdir",
		default="results/task6",
		help="Output directory for species-specific ortholog BED files.",
	)
	args = parser.parse_args()

	ortholog_table = Path(args.ortholog_table)
	outdir = Path(args.outdir)
	outdir.mkdir(parents=True, exist_ok=True)

	if not ortholog_table.exists():
		raise FileNotFoundError(f"Ortholog table not found: {ortholog_table}")

	df = pd.read_table(ortholog_table)
	human_cols, mouse_cols = pick_columns(df)

	human_ortho = df[human_cols].drop_duplicates().copy()
	mouse_ortho = df[mouse_cols].drop_duplicates().copy()

	human_out = outdir / "human_ortho_unique.bed"
	mouse_out = outdir / "mouse_ortho_unique.bed"

	human_ortho.to_csv(human_out, sep="\t", index=False, header=False)
	mouse_ortho.to_csv(mouse_out, sep="\t", index=False, header=False)

	print(f"Wrote {len(human_ortho)} human conserved OCRs to {human_out}")
	print(f"Wrote {len(mouse_ortho)} mouse conserved OCRs to {mouse_out}")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())