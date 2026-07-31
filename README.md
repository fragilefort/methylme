---
title: Wgbs2
author: 
date:
theme: united
highlight: tango
output:
  html_document:
    toc: true
    toc_depth: 3
    number_sections: true
    toc_float:
      collapsed: true
      smooth_scroll: false
tags: wgbs
HackMD Tutorial: https://hackmd.io/c/tutorials/%2FUFLeoGd_SmGy8acWOTs6JA
---
## Day 1
### System arc
| Node | CPUs | Memory |
| :--- | :---: | :---: |
| **bibigrid-master-aorkgacx5hn960m** | 2 | 2000 |
| **bibigrid-worker-aorkgacx5hn960m-0** | 14 | 112000 |
| **bibigrid-worker-aorkgacx5hn960m-1** | 14 | 112000 |
| **bibigrid-worker-aorkgacx5hn960m-2** | 14 | 112000 |
| **bibigrid-worker-aorkgacx5hn960m-3** | 14 | 112000 |
| **bibigrid-worker-aorkgacx5hn960m-4** | 14 | 112000 |

There are 1 master node and 5 worker nodes, the table above has both the CPU and the RAM for each node, the volume size is 4.9T. The tmp folder is located at `/vol/COMPEPIWS/tmp/`. 

### Task 3 

<img width="931" height="118" alt="Screenshot 2026-07-30 113519" src="https://github.com/user-attachments/assets/4681fac7-4eee-4338-8e80-0f932b7f5622" />

### Task 4

<img width="1895" height="538" alt="Screenshot 2026-07-30 114751" src="https://github.com/user-attachments/assets/d62ed496-f21b-4eac-9161-7f3cfff14322" />

##number of lines equal: 
58239 

## number of exon entries:

24556

## number of exons longer than 1000bp :
1663

## number of exons where gene_id is "Sox17":
0

## Count exons in chr2 :
0  

## frequency of each feature type in column 3:

   2260 3UTR
   3363 5UTR
  21278 CDS
  24556 exon
   2146 start_codon
   2139 stop_codon
   2497 transcript
### The new sorted file:
chr18   refGene 3UTR    10064401        10066047        .       -       .       gene_id "Rock1"; transcript_id "NM_009071"; exon_number "1"; exon_id "NM_009071.1"; gene_name "Rock1";
chr18   refGene 3UTR    10560484        10562941        .       +       .       gene_id "Greb1l"; transcript_id "NM_001083628"; exon_number "33"; exon_id "NM_001083628.33"; gene_name "Greb1l";
chr18   refGene 3UTR    10566512        10567498        .       -       .       gene_id "Esco1"; transcript_id "NM_001081222"; exon_number "1"; 

# Task 6

done
Submitted batch job 23680
Submitted batch job 23681
Submitted batch job 23682
Submitted batch job 23683
Submitted batch job 23684
Submitted batch job 23685
Submitted batch job 23686
Submitted batch job 23687
Submitted batch job 23688
Submitted batch job 23689
Submitted batch job 23690
Submitted batch job 23691
Submitted batch job 23692
Submitted batch job 23693
Submitted batch job 23694

  JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
             23692   compepi  test_13 ibha0000  R       0:03      1 bibigrid-worker-aorkgacx5hn960m-1
             23693   compepi  test_14 ibha0000  R       0:03      1 bibigrid-worker-aorkgacx5hn960m-1
             23694   compepi  test_15 ibha0000  R       0:03      1 bibigrid-worker-aorkgacx5hn960m-2
             23680   compepi   test_1 ibha0000  R       0:04      1 bibigrid-worker-aorkgacx5hn960m-1
             23681   compepi   test_2 ibha0000  R       0:04      1 bibigrid-worker-aorkgacx5hn960m-1
             23682   compepi   test_3 ibha0000  R       0:04      1 bibigrid-worker-aorkgacx5hn960m-1
# Task 7

Q4

Table	Rows	Columns
Liver	1,362,728	6
Kidney	1,362,728	6

Q5 — Column names of each table

Both tables share the exact same column names:

text
"Chr"  "Start"  "End"  "H3K27me3"  "H3K36me3"  "H3K9me3"

Q6 — Genome length from each dataset

Both liver and kidney give the same total: 2,725,456,000 bp (~2.73 Gb), consistent with the mouse genome (mm10 build), calculated as the sum
 of (End - Start) across all 2kb bins.

Q8 — Dimension of the new (vertically concatenated) data frame

2,725,456 rows × 7 columns — exactly double the row count of one table (1,362,728 × 2), plus one extra column (cell_type) added on top of
 the original 6.

The density plots:
<img width="2000" height="1200" alt="image" src="https://github.com/user-attachments/assets/d9677325-ef09-44e2-86ca-d0f905ca2425" />
<img width="2000" height="1200" alt="image" src="https://github.com/user-attachments/assets/2a281e3c-8f76-4efd-af38-cf6133c613d1" />


Q16 -  scree plot and cumulative plot 

<img width="2000" height="1500" alt="image" src="https://github.com/user-attachments/assets/7bb9c0d6-34e5-4788-a7c4-85f52f4c9c23" />
<img width="2000" height="1500" alt="image" src="https://github.com/user-attachments/assets/01fc194c-3780-4e78-a9ab-143eb7f85246" />


Q17 - 

6 variables × 6 components = 36 individual loading values, organized as a 6×6 matrix rather than a flat list.

Q18 PC1 vs PC2 loadings interpretation

Variable	PC1	PC2
H3K27me3_liver	0.394	0.314
H3K36me3_liver	0.350	-0.624
H3K9me3_liver	0.391	0.241
H3K27me3_kidney	0.449	0.310
H3K36me3_kidney	0.386	-0.572
H3K9me3_kidney	0.468	0.175


PC1: All six variables load positively with similar magnitude (0.35–0.47), meaning PC1 captures overall signal intensity shared across all
 three histone marks in both tissues — it doesn't distinguish between marks or tissues, just represents a general "high vs low signal" axis.

PC2: Splits the marks by type, not tissue. H3K36me3 loads strongly negative in both liver (-0.62) and kidney (-0.57), 
while H3K27me3 and H3K9me3 both load positively (~0.18–0.31) in both tissues.
 This means PC2 separates H3K36me3 (a mark typically associated with active transcription) from the two repressive marks (H3K27me3, H3K9me3),
 regardless of which tissue the signal came from.


<img width="2000" height="1500" alt="image" src="https://github.com/user-attachments/assets/0f321366-0a04-4894-aaa5-704196de05e2" />

## Genomic ranges

Q2 — Total number of bases covered

Two different interpretations were computed:

Method	Liver	Kidney
Sum of all range widths (sum(width(gr)))	2,726,818,728	2,726,818,728
Sum of widths after merging overlaps (sum(width(reduce(gr))))	2,725,456,021	2,725,456,021


Q5 — Overlapping regions between liver and kidney

text
findOverlaps(liver_gr, kidney_gr) → 4,088,142 hits
pintersect(...) → GRanges with 4,088,142 overlapping regions

# Pipeline
## 1. Main steps of the pipelines

- Quality control
- Read trimming
- Read alignment
- Deduplicate Alignments
- Extract methylation calls
- Sample report
- Alignment QC
