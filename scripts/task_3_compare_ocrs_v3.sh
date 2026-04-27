#!/bin/bash
#step 3: Cross-species comparison of OCRs (human vs mouse)
#author: mrun :)
#date: march 30, 2026

#load required module
source /etc/profile.d/modules.sh
module load bedtools

#output directory
OUTDIR="results/mapping"
mkdir -p $OUTDIR

#merge fragmented lifted regions
#HALliftover splits peaks into fragments during lifting
#reconstruct each peak by grouping fragments by peak ID
#take the min start and max end per chromosome

#merge fragmented human to mouse lifted regions
awk 'NR>1' results/mapping/tmp/human_to_mouse.lifted.bed | \
  sort -k4,4 -k1,1 -k2,2n | \
  bedtools groupby -g 1,4 -c 2,3 -o min,max | \
  awk '{print $1"\t"$3"\t"$4"\t"$2"\t1000\t."}' | \
  sort -k1,1 -k2,2n > results/mapping/tmp/human_to_mouse.lifted.merged.bed
echo "Done: merged human to mouse lifted regions"

#for mouse to human lifted regions, some peaks map to multiple
#chromosomes, only the chromosome with the most
#fragments per peak (best lift) before merging
awk '{print $4"\t"$1}' results/mapping/tmp/mouse_to_human.lifted.bed | \
  sort | uniq -c | \
  sort -k2,2 -k1,1rn | \
  sort -k2,2 -u | \
  awk '{print $2"\t"$3}' > results/mapping/tmp/mouse_best_chrom.txt

#filter mouse lifted file to best chromosome per peak, then merge
awk 'NR==FNR {best[$1]=$2; next} ($4 in best) && best[$4]==$1' \
  results/mapping/tmp/mouse_best_chrom.txt \
  results/mapping/tmp/mouse_to_human.lifted.bed | \
  sort -k4,4 -k1,1 -k2,2n | \
  bedtools groupby -g 1,4 -c 2,3 -o min,max | \
  awk '{print $1"\t"$3"\t"$4"\t"$2"\t1000\t."}' | \
  sort -k1,1 -k2,2n > results/mapping/tmp/mouse_to_human.lifted.merged.bed
echo "Done: merged mouse to human lifted regions"

#extract conserved OCRs from orthologous pairs
#only peaks that are best matches in both directions are considered truly conserved

#conserved human OCRs (query coordinates)
awk 'NR>1 {print $2"\t"$3"\t"$4"\t"$1}' \
  results/mapping/orthologous_ocr_pairs.tsv \
  > $OUTDIR/conserved_human_in_mouse.bed
echo "Done: conserved human OCRs extracted from orthologous pairs"

#conserved mouse OCRs (target coordinates)
awk 'NR>1 {print $6"\t"$7"\t"$8"\t"$5}' \
  results/mapping/orthologous_ocr_pairs.tsv \
  > $OUTDIR/conserved_mouse_in_human.bed
echo "Done: conserved mouse OCRs extracted from orthologous pairs"

#SPECIES-SPECIFIC OCRs USING BEDTOOLS
#human/mouse OCRs (lifted) that DO NOT overlap with the
#other species' OCRs (>=30% reciprocal overlap, i set this threshold, can be changed)
#represent species-specific regulatory regions

#human-specific OCRs:
#human OCRs (lifted to mouse genome) that DO NOT
#overlap with mouse OCRs (>=30% reciprocal overlap)
#represents human-specific regulatory regions
bedtools intersect \
  -a results/mapping/tmp/human_to_mouse.lifted.merged.bed \
  -b results/mapping/mouse_pancreas_ocr.processed.bed \
  -f 0.3 -r \
  -v > $OUTDIR/human_specific.bed
echo "Done: human-specific OCRs identified"

#mouse-specific OCRs:
#mouse OCRs (lifted to human genome) that DO NOT
#overlap with human OCRs (>=30% reciprocal overlap)
#represents mouse-specific regulatory regions
bedtools intersect \
  -a results/mapping/tmp/mouse_to_human.lifted.merged.bed \
  -b results/mapping/human_pancreas_ocr.processed.bed \
  -f 0.3 -r \
  -v > $OUTDIR/mouse_specific.bed
echo "Done: mouse-specific OCRs identified"

#SUMMARY STATSSSSS

#count number of regions in each category
#wc -l counts number of lines (each line = one OCR)
echo "Summary Time CHATTTTTTT"
echo "Conserved human in mouse: $(wc -l < $OUTDIR/conserved_human_in_mouse.bed)"
echo "Human-specific: $(wc -l < $OUTDIR/human_specific.bed)"
echo "Conserved mouse in human: $(wc -l < $OUTDIR/conserved_mouse_in_human.bed)"
echo "Mouse-specific: $(wc -l < $OUTDIR/mouse_specific.bed)"


