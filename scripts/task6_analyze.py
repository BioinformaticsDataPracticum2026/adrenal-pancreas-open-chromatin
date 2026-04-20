"""
Task 6: Motif enrichment analysis for promoters vs enhancers.

This script parses HOMER motif discovery results and generates summary
visualizations comparing enriched motifs across promoter/enhancer groups.
"""

import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def parse_homer_results(txt_file: str, category: str, top_n: int = 5) -> pd.DataFrame:
    """
    Parse HOMER knownResults.txt file and extract top motifs.
    
    HOMER output format has header lines to skip. Motif names may include
    family annotations (e.g., "NeuroD1(bHLH)/ATOH1") which we simplify.
    
    Args:
        txt_file: Path to HOMER knownResults.txt file.
        category: Label for this result set (e.g., 'Human Promoters').
        top_n: Number of top motifs to return (default 5).
    
    Returns:
        DataFrame with top motifs including columns:
        - Motif_Name: Simplified motif name (TF only)
        - Log_P_value: -log10(p-value) from HOMER
        - Category: Category label
        - Abs_Log_P: Absolute value of log p-value for plotting
    """
    # Read HOMER output, skipping header lines
    df = pd.read_csv(txt_file, sep='\t', skiprows=1)
    
    # Extract motif name (column 0) and log p-value (column 3)
    df.rename(columns={
        df.columns[0]: 'Motif_Name', 
        df.columns[3]: 'Log_P_value' 
    }, inplace=True)
    
    # Simplify motif names by removing family annotations
    # E.g., "NeuroD1(bHLH)/ATOH1" becomes "NeuroD1"
    df['Motif_Name'] = df['Motif_Name'].apply(lambda x: str(x).split('/')[0])
    
    # Keep top N and prepare for plotting
    top_motifs = df.head(top_n).copy()
    top_motifs['Category'] = category
    top_motifs['Abs_Log_P'] = abs(top_motifs['Log_P_value'])
    
    return top_motifs


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate motif enrichment summary plots from HOMER results."
    )
    parser.add_argument(
        "--base-path",
        default="results/task6",
        help="Base directory containing HOMER result subdirectories.",
    )
    parser.add_argument(
        "--output",
        default="results/task6/motif_top5_panel.png",
        help="Output figure path.",
    )
    return parser.parse_args()


def main() -> int:
    """
    Main entry point: load HOMER results, parse, and plot.
    
    Expected subdirectories under base-path:
    - homer_human_promoters/knownResults.txt
    - homer_human_enhancers/knownResults.txt
    - homer_mouse_promoters/knownResults.txt
    - homer_mouse_enhancers/knownResults.txt
    
    Returns:
        0 on success, non-zero on error.
    """
    args = parse_arguments()
    base_path = Path(args.base_path)

    # Define required HOMER result files
    required_files = [
        base_path / "homer_human_promoters" / "knownResults.txt",
        base_path / "homer_human_enhancers" / "knownResults.txt",
        base_path / "homer_mouse_promoters" / "knownResults.txt",
        base_path / "homer_mouse_enhancers" / "knownResults.txt",
    ]

    # Check all files exist
    missing_files = [str(p) for p in required_files if not p.exists()]
    if missing_files:
        print("Error: Missing HOMER knownResults.txt files:")
        for f in missing_files:
            print(f"  {f}")
        return 1

    # Parse each HOMER result file
    human_promoter_motifs = parse_homer_results(
        str(required_files[0]), 'Human Promoters'
    )
    human_enhancer_motifs = parse_homer_results(
        str(required_files[1]), 'Human Enhancers'
    )
    mouse_promoter_motifs = parse_homer_results(
        str(required_files[2]), 'Mouse Promoters'
    )
    mouse_enhancer_motifs = parse_homer_results(
        str(required_files[3]), 'Mouse Enhancers'
    )

    # Combine all results into single dataset
    all_motifs = pd.concat([
        human_promoter_motifs,
        human_enhancer_motifs,
        mouse_promoter_motifs,
        mouse_enhancer_motifs,
    ])

    # Create faceted barplot with seaborn
    g = sns.FacetGrid(
        all_motifs,
        col="Category",
        col_wrap=2,
        sharex=False,
        sharey=False,
        height=4,
        aspect=1.5
    )
    g.map(sns.barplot, "Abs_Log_P", "Motif_Name", palette="viridis")
    g.set_titles(col_template="{col_name}", size=14, fontweight='bold')
    g.set_axis_labels("Absolute Log P-value", "Transcription Factor")
    
    plt.suptitle(
        "Top Enriched Motifs: Promoters vs. Enhancers",
        y=1.05,
        fontsize=16,
        fontweight='bold'
    )
    plt.tight_layout()

    # Save output
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved motif summary plot to {output_path}")
    plt.close()
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
