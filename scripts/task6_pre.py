import pandas as pd

# Load reciprocal best-hit pairs from Step 3
# Ensure this file path is correct for your local setup
ortho_file = "orthologous_ocr_pairs.tsv"
df = pd.read_table(ortho_file)

# Extract Human unique coordinates for conserved peaks
human_ortho = df[['human_chr', 'human_start', 'human_end', 'human_peak_id']].copy()
human_ortho.to_csv("results/task6_motifs/human_ortho_unique.bed", sep='\t', index=False, header=False)

# Extract Mouse unique coordinates for conserved peaks
mouse_ortho = df[['mouse_chr', 'mouse_start', 'mouse_end', 'mouse_peak_id']].copy()
mouse_ortho.to_csv("results/task6_motifs/mouse_ortho_unique.bed", sep='\t', index=False, header=False)

print("BED files for conserved OCRs generated successfully.")