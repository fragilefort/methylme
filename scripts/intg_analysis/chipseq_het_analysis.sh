#!/usr/bin/env bash
set -e

WGBS_TABLE="/vol/COMPEPIWS/groups/shared/WGBS/wgbs2/differential/differential_methylation_data/diffMethTable_region_cmp1_tiling.csv"
KIDNEY_SEG="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/segmentation/kidney_15_segments.bed"
LIVER_SEG="/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/segmentation/liver_15_segments.bed"

echo "===  Filtering DMRs ==="
# FDR < 0.05 and |diff| >= 0.2 (20% methylation difference)
awk -F',' 'NR>1 && $10 < 0.05 && ($7 >= 0.2 || $7 <= -0.2) {print $2"\t"$3"\t"$4}' "$WGBS_TABLE" \
  | sort -k1,1 -k2,2n > dmrs.bed
echo "Total DMRs: $(wc -l < dmrs.bed)"

echo "===  Filtering Heterochromatin States (Het_P and Het_S) ==="
grep -E "Het_P|Het_S" "$KIDNEY_SEG" | sort -k1,1 -k2,2n > kidney_Het.bed
grep -E "Het_P|Het_S" "$LIVER_SEG"  | sort -k1,1 -k2,2n > liver_Het.bed

echo "Kidney Het regions: $(wc -l < kidney_Het.bed)"
echo "Liver Het regions: $(wc -l < liver_Het.bed)"

echo "===  Comparing DMRs with Liver Het ==="
LIVER_OVERLAP=$(bedtools intersect -a dmrs.bed -b liver_Het.bed -u | wc -l)
echo "DMRs overlapping Liver Het state: $LIVER_OVERLAP"

echo "===  Comparing DMRs with Kidney Het ==="
KIDNEY_OVERLAP=$(bedtools intersect -a dmrs.bed -b kidney_Het.bed -u | wc -l)
echo "DMRs overlapping Kidney Het state: $KIDNEY_OVERLAP"

echo "===  Merging Het Regions Across Tissues ==="
cat kidney_Het.bed liver_Het.bed | sort -k1,1 -k2,2n | bedtools merge > combined_Het.bed
echo "Combined Het regions: $(wc -l < combined_Het.bed)"

echo "===  Comparing DMRs with Merged Het Regions ==="
COMBINED_OVERLAP=$(bedtools intersect -a dmrs.bed -b combined_Het.bed -u | wc -l)
echo "DMRs overlapping combined Het state (either tissue): $COMBINED_OVERLAP"
