from __future__ import annotations

import pandas as pd


def reciprocal_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> tuple[int, float, float]:
    """
    Calculate reciprocal overlap between two intervals.
    
    Args:
        a_start, a_end: Start and end positions of interval A.
        b_start, b_end: Start and end positions of interval B.
    
    Returns:
        Tuple of (overlap_bp, overlap_fraction_in_a, overlap_fraction_in_b).
        overlap_fraction_in_a is the fraction of A covered by the overlap.
        overlap_fraction_in_b is the fraction of B covered by the overlap.
    """
    overlap_bp = max(0, min(a_end, b_end) - max(a_start, b_start))
    if overlap_bp == 0:
        return 0, 0.0, 0.0
    fraction_a = overlap_bp / max(1, a_end - a_start)
    fraction_b = overlap_bp / max(1, b_end - b_start)
    return overlap_bp, fraction_a, fraction_b


def split_by_chrom(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """
    Split a DataFrame of genomic intervals by chromosome.
    
    Args:
        df: DataFrame with a 'chrom' column.
    
    Returns:
        Dictionary mapping chromosome names to subsets of the DataFrame,
        sorted by start position within each chromosome.
    """
    grouped: dict[str, pd.DataFrame] = {}
    for chrom, subset in df.sort_values(['chrom', 'start', 'end']).groupby('chrom', sort=False):
        grouped[str(chrom)] = subset.reset_index(drop=True)
    return grouped


def find_best_overlaps(query_df: pd.DataFrame, target_df: pd.DataFrame, min_reciprocal_overlap: float) -> pd.DataFrame:
    """
    Find the best (highest overlap) target interval for each query interval.
    
    Only returns matches where the reciprocal overlap (minimum of the two overlap
    fractions) meets or exceeds min_reciprocal_overlap threshold.
    
    Args:
        query_df: DataFrame of query intervals (must have peak_id, chrom, start, end).
        target_df: DataFrame of target intervals (must have peak_id, chrom, start, end).
        min_reciprocal_overlap: Minimum reciprocal overlap fraction (0.0 to 1.0).
    
    Returns:
        DataFrame of best query-target pairs with overlap statistics.
        Returns empty DataFrame if no matches found.
    """
    result_rows: list[dict[str, object]] = []
    query_by_chrom = split_by_chrom(query_df)
    target_by_chrom = split_by_chrom(target_df)

    for chrom, query_subset in query_by_chrom.items():
        target_subset = target_by_chrom.get(chrom)
        if target_subset is None or target_subset.empty:
            continue

        target_records = list(target_subset.itertuples(index=False))
        # Track position in target list for efficiency
        left_idx = 0
        
        for query_record in query_subset.itertuples(index=False):
            # Skip targets that end before this query starts
            while left_idx < len(target_records) and int(target_records[left_idx].end) <= int(query_record.start):
                left_idx += 1

            best_match_row = None
            best_overlap_score = -1.0
            
            # Scan targets from current position while they start before query ends
            scan_idx = left_idx
            while scan_idx < len(target_records) and int(target_records[scan_idx].start) < int(query_record.end):
                target_record = target_records[scan_idx]
                overlap_bp, overlap_frac_query, overlap_frac_target = reciprocal_overlap(
                    int(query_record.start),
                    int(query_record.end),
                    int(target_record.start),
                    int(target_record.end),
                )
                # Check if reciprocal overlap meets threshold
                min_overlap_frac = min(overlap_frac_query, overlap_frac_target)
                if min_overlap_frac >= min_reciprocal_overlap:
                    if min_overlap_frac > best_overlap_score:
                        best_overlap_score = min_overlap_frac
                        best_match_row = {
                            'query_peak_id': query_record.peak_id,
                            'query_chrom': query_record.chrom,
                            'query_start': int(query_record.start),
                            'query_end': int(query_record.end),
                            'target_peak_id': target_record.peak_id,
                            'target_chrom': target_record.chrom,
                            'target_start': int(target_record.start),
                            'target_end': int(target_record.end),
                            'overlap_bp': overlap_bp,
                            'query_overlap_fraction': overlap_frac_query,
                            'target_overlap_fraction': overlap_frac_target,
                            'reciprocal_overlap_score': min_overlap_frac,
                        }
                scan_idx += 1

            if best_match_row is not None:
                result_rows.append(best_match_row)

    return pd.DataFrame(result_rows)
