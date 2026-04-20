from __future__ import annotations

from pathlib import Path

import pandas as pd


NARROWPEAK_COLUMNS = [
    'chrom',
    'start',
    'end',
    'name',
    'score',
    'strand',
    'signal_value',
    'p_value',
    'q_value',
    'summit',
]


def discover_peak_file(candidates: list[str]) -> str:
    """
    Discover the first available peak file from a list of candidates.
    
    Args:
        candidates: List of file paths to check in order.
    
    Returns:
        Path to the first existing candidate file.
    
    Raises:
        FileNotFoundError: If no candidate file exists.
    """
    for candidate in candidates:
        if Path(candidate).exists():
            return candidate
    joined = '\n'.join(candidates)
    raise FileNotFoundError(f'Could not locate any configured peak file. Checked:\n{joined}')


def read_peak_table(path: str | Path) -> pd.DataFrame:
    """
    Read a narrowPeak format file into a DataFrame.
    
    Args:
        path: Path to narrowPeak file (optionally gzip-compressed).
    
    Returns:
        DataFrame with standard narrowPeak columns.
    """
    compression = 'gzip' if str(path).endswith('.gz') else None
    return pd.read_csv(
        path,
        sep='\t',
        header=None,
        names=NARROWPEAK_COLUMNS,
        compression=compression,
    )


def standardize_peak_table(df: pd.DataFrame, species: str) -> pd.DataFrame:
    """
    Standardize and clean a raw narrowPeak DataFrame.
    
    Performs type coercion, removes invalid peaks, removes duplicates,
    assigns unique peak IDs, and reorders columns.
    
    Args:
        df: Raw narrowPeak DataFrame.
        species: Species label (e.g., 'human', 'mouse') for peak ID generation.
    
    Returns:
        Standardized DataFrame with columns ordered for downstream use.
    """
    result = df.copy()
    result['start'] = result['start'].astype(int)
    result['end'] = result['end'].astype(int)
    result['score'] = pd.to_numeric(result['score'], errors='coerce').fillna(0).astype(int)
    result['signal_value'] = pd.to_numeric(result['signal_value'], errors='coerce')
    result['p_value'] = pd.to_numeric(result['p_value'], errors='coerce')
    result['q_value'] = pd.to_numeric(result['q_value'], errors='coerce')
    result['summit'] = pd.to_numeric(result['summit'], errors='coerce').fillna(-1).astype(int)
    # Remove invalid peaks where start >= end
    result = result[result['end'] > result['start']].copy()
    # Remove exact duplicates
    result = result.drop_duplicates(subset=['chrom', 'start', 'end']).reset_index(drop=True)
    # Assign unique peak IDs
    result.insert(0, 'peak_id', [f'{species}_ocr_{i:07d}' for i in range(1, len(result) + 1)])
    result['species'] = species
    result['width'] = result['end'] - result['start']
    return result[[
        'peak_id', 'species', 'chrom', 'start', 'end', 'width', 'score', 'signal_value',
        'p_value', 'q_value', 'summit', 'name', 'strand'
    ]]


def write_bed(df: pd.DataFrame, path: str | Path) -> None:
    """
    Write a DataFrame to BED format.
    
    Args:
        df: DataFrame containing at least chrom, start, end, peak_id, score columns.
        path: Output BED file path.
    """
    bed = df[['chrom', 'start', 'end', 'peak_id', 'score']].copy()
    bed['strand'] = '.'
    bed.to_csv(path, sep='\t', header=False, index=False)


def write_table(df: pd.DataFrame, path: str | Path) -> None:
    """
    Write a DataFrame to tab-delimited TSV format.
    
    Args:
        df: DataFrame to write.
        path: Output TSV file path.
    """
    df.to_csv(path, sep='\t', index=False)


def write_manifest(entries: list[dict[str, str]], path: str | Path) -> None:
    """
    Write a manifest of file inputs and outputs.
    
    Args:
        entries: List of dictionaries with manifest information.
        path: Output manifest file path.
    """
    pd.DataFrame(entries).to_csv(path, sep='\t', index=False)
