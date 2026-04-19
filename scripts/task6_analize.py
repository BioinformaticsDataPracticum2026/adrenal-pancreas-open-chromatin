import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def get_top_motifs_txt(txt_file, category):
    # Read the text file and SKIP the first line (the HOMER warning message)
    df = pd.read_csv(txt_file, sep='\t', skiprows=1)
    
    # FOOLPROOF TRICK: Rename columns based on their position (Index 0 and Index 3) 
    df.rename(columns={
        df.columns[0]: 'Motif_Name', 
        df.columns[3]: 'Log_P_value' 
    }, inplace=True)
    
    # Clean the motif name (e.g., "NeuroD1(bHLH)/..." becomes just "NeuroD1")
    df['Motif_Name'] = df['Motif_Name'].apply(lambda x: str(x).split('/')[0])
    
    # Keep top 5, make the Log P-value absolute for plotting
    top5 = df.head(5).copy()
    top5['Category'] = category
    top5['Abs_Log_P'] = abs(top5['Log_P_value'])
    return top5

# Your exact base path
base_path = r"C:\Users\joewa\files\cmu\bioinfodata"

# Parse all 4 TXT files
human_pro = get_top_motifs_txt(f"{base_path}\\homer_human_promoters\\knownResults.txt", 'Human Promoters')
human_enh = get_top_motifs_txt(f"{base_path}\\homer_human_enhancers\\knownResults.txt", 'Human Enhancers')
mouse_pro = get_top_motifs_txt(f"{base_path}\\homer_mouse_promoters\\knownResults.txt", 'Mouse Promoters')
mouse_enh = get_top_motifs_txt(f"{base_path}\\homer_mouse_enhancers\\knownResults.txt", 'Mouse Enhancers')

# Combine into one dataset
all_data = pd.concat([human_pro, human_enh, mouse_pro, mouse_enh])

# Plot using Seaborn
g = sns.FacetGrid(all_data, col="Category", col_wrap=2, sharex=False, sharey=False, height=4, aspect=1.5)
g.map(sns.barplot, "Abs_Log_P", "Motif_Name", palette="viridis")
g.set_titles(col_template="{col_name}", size=14, fontweight='bold')
g.set_axis_labels("Absolute Log P-value", "Transcription Factor")
plt.suptitle("Top Enriched Motifs: Promoters vs. Enhancers", y=1.05, fontsize=16, fontweight='bold')
plt.tight_layout()
plt.show()