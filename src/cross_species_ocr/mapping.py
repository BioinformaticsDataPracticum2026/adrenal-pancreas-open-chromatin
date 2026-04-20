from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pandas as pd

from .intervals import find_best_overlaps


def hal_liftover_available(binary: str) -> bool:
    """
    Check if halLiftover binary is available on PATH or at specified location.
    
    Args:
        binary: Name or path of halLiftover binary.
    
    Returns:
        True if binary can be found, False otherwise.
    """
    return shutil.which(binary) is not None


def run_hal_liftover(binary: str, hal_file: str, src_genome: str, src_bed: str | Path, dst_genome: str, out_bed: str | Path) -> None:
    """
    Execute halLiftover to map genomic intervals from source to target genome.
    
    Args:
        binary: Path or name of halLiftover executable.
        hal_file: Path to HAL alignment file.
        src_genome: Source genome name (must match HAL file).
        src_bed: Input BED file with intervals in source genome coordinates.
        dst_genome: Destination genome name (must match HAL file).
        out_bed: Output BED file path for lifted intervals.
    
    Raises:
        CalledProcessError: If halLiftover exits with non-zero code.
    """
    command = [binary, hal_file, src_genome, str(src_bed), dst_genome, str(out_bed)]
    subprocess.run(command, check=True)


def read_liftover_bed(path: str | Path, source_prefix: str, target_species: str) -> pd.DataFrame:
    """
    Parse halLiftover output BED file.
    
    halLiftover fragments peaks during lifting. This function reads the raw
    lifted intervals and filters to peaks matching the source prefix.
    Fragmentation is handled by downstream grouping (see pipeline.py).
    
    Args:
        path: Path to lifted BED file from halLiftover.
        source_prefix: Prefix to filter lifted peaks (e.g., 'human_ocr_').
        target_species: Label for the target species.
    
    Returns:
        DataFrame with lifted interval information.
        Columns: peak_id, source_peak_id, chrom, start, end, target_species, width.
    """
    rows = []
    with open(path, 'r', encoding='utf-8') as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 4:
                continue
            peak_id = fields[3]
            rows.append({
                'peak_id': peak_id,
                'source_peak_id': peak_id,
                'chrom': fields[0],
                'start': int(fields[1]),
                'end': int(fields[2]),
                'target_species': target_species,
            })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    # Filter to peaks from the source species
    df = df[df['source_peak_id'].str.startswith(source_prefix)].copy()
    df['width'] = df['end'] - df['start']
    return df


def build_pair_tables(
    forward_lifted: pd.DataFrame,
    reverse_lifted: pd.DataFrame,
    target_df: pd.DataFrame,
    reverse_target_df: pd.DataFrame,
    min_reciprocal_overlap: float,
    source_label: str,
    target_label: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build forward and reverse mapping tables with reciprocal best-hit filtering.
    
    Identifies orthologs using reciprocal best-hit logic:
    - Forward direction: finds best overlap for each source peak in target
    - Reverse direction: finds best overlap for each target peak in source
    - Reciprocal best-hit: marked if pair is best match in both directions
    
    Args:
        forward_lifted: Lifted peaks from source to target genome.
        reverse_lifted: Lifted peaks from target back to source genome.
        target_df: Target genome peaks (for forward direction).
        reverse_target_df: Source genome peaks (for reverse direction).
        min_reciprocal_overlap: Minimum reciprocal overlap threshold (0.0-1.0).
        source_label: Label for source species (e.g., 'human').
        target_label: Label for target species (e.g., 'mouse').
    
    Returns:
        Tuple of (forward_pairs_df, reverse_pairs_df) with reciprocal best-hit flags.
    """
    forward_best = find_best_overlaps(forward_lifted, target_df, min_reciprocal_overlap)
    reverse_best = find_best_overlaps(reverse_lifted, reverse_target_df, min_reciprocal_overlap)

    if forward_best.empty:
        return forward_best, reverse_best

    # Build lookup of reverse pairs for reciprocal matching
    reverse_lookup = {
        (row.query_peak_id, row.target_peak_id)
        for row in reverse_best.itertuples(index=False)
    }
    # Mark reciprocal best-hits: forward pair is valid if reverse exists
    forward_best['reciprocal_best_hit'] = [
        (row.target_peak_id, row.query_peak_id) in reverse_lookup
        for row in forward_best.itertuples(index=False)
    ]
    forward_best['source_species'] = source_label
    forward_best['target_species'] = target_label
    return forward_best, reverse_best
