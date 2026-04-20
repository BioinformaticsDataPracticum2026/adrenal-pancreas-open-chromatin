from __future__ import annotations

from bisect import bisect_left
from collections import defaultdict
from pathlib import Path

import pandas as pd


def load_tss_index(path: str | Path) -> dict[str, dict[str, list]]:
    """
    Load and index TSS positions from a BED-like file.
    
    Builds a per-chromosome index of TSS positions for fast nearest-neighbor lookup.
    Expected format: chrom, start, end, gene_name, strand (0-indexed, BED format).
    
    Args:
        path: Path to TSS BED file.
    
    Returns:
        Dictionary mapping chrom -> {positions, genes, strands} for fast lookup.
    """
    chrom_to_rows: dict[str, list[tuple[int, str, str]]] = defaultdict(list)
    with open(path, 'r', encoding='utf-8') as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 4:
                continue
            chrom = fields[0]
            start = int(fields[1])
            gene = fields[3]
            strand = fields[4] if len(fields) > 4 else '.'
            chrom_to_rows[chrom].append((start, gene, strand))

    index: dict[str, dict[str, list]] = {}
    for chrom, rows in chrom_to_rows.items():
        rows.sort(key=lambda item: item[0])
        index[chrom] = {
            'positions': [row[0] for row in rows],
            'genes': [row[1] for row in rows],
            'strands': [row[2] for row in rows],
        }
    return index


def annotate_with_nearest_tss(df: pd.DataFrame, tss_index: dict[str, dict[str, list]], promoter_window_bp: int) -> pd.DataFrame:
    """
    Annotate peaks with nearest TSS and classify as promoter or distal.
    
    For each peak, finds the nearest TSS on the same chromosome using binary search.
    Classifies peak as 'promoter' if within promoter_window_bp of TSS, else 'distal'.
    
    Args:
        df: Peak DataFrame with chrom, start, end columns.
        tss_index: TSS index from load_tss_index().
        promoter_window_bp: Window size around TSS (e.g., 2000 bp for ±2kb).
    
    Returns:
        DataFrame with added columns:
        - nearest_gene: Gene name of nearest TSS (or 'NA')
        - distance_to_nearest_tss: Distance in bp (or pd.NA)
        - nearest_gene_strand: Strand of nearest TSS (or 'NA')
        - regulatory_class: 'promoter' or 'distal' or 'unannotated'
    """
    nearest_gene = []
    nearest_distance = []
    nearest_strand = []
    regulatory_class = []

    for row in df.itertuples(index=False):
        peak_midpoint = (int(row.start) + int(row.end)) // 2
        chrom_index = tss_index.get(row.chrom)
        
        if chrom_index is None:
            # Chromosome not in index (e.g., chrM, scaffolds)
            nearest_gene.append('NA')
            nearest_distance.append(pd.NA)
            nearest_strand.append('NA')
            regulatory_class.append('unannotated')
            continue

        # Binary search for nearest TSS using bisect
        positions = chrom_index['positions']
        idx = bisect_left(positions, peak_midpoint)
        
        # Collect candidate indices (position at or just before/after midpoint)
        candidate_indices = []
        if idx < len(positions):
            candidate_indices.append(idx)
        if idx > 0:
            candidate_indices.append(idx - 1)

        # Find closest TSS
        best_idx = min(candidate_indices, key=lambda i: abs(positions[i] - peak_midpoint))
        distance_to_tss = abs(positions[best_idx] - peak_midpoint)
        
        nearest_gene.append(chrom_index['genes'][best_idx])
        nearest_distance.append(distance_to_tss)
        nearest_strand.append(chrom_index['strands'][best_idx])
        regulatory_class.append('promoter' if distance_to_tss <= promoter_window_bp else 'distal')

    annotated = df.copy()
    annotated['nearest_gene'] = nearest_gene
    annotated['distance_to_nearest_tss'] = nearest_distance
    annotated['nearest_gene_strand'] = nearest_strand
    annotated['regulatory_class'] = regulatory_class
    return annotated
