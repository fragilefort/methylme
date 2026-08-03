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

There is one master node and five worker nodes. The table above shows the available CPUs and RAM for each node. The shared volume size is 4.9T, and the temporary directory is located at `/vol/COMPEPIWS/tmp/`.

### Task 3

<img width="931" height="118" alt="Screenshot 2026-07-30 113519" src="https://github.com/user-attachments/assets/4681fac7-4eee-4338-8e80-0f932b7f5622" />

### Task 4

<img width="1895" height="538" alt="Screenshot 2026-07-30 114751" src="https://github.com/user-attachments/assets/d62ed496-f21b-4eac-9161-7f3cfff14322" />

#### Number of lines

58239

#### Number of exon entries

24556

#### Number of exons longer than 1000 bp

1663

#### Number of exons where gene_id is "Sox17"

0

#### Count of exons in chr2

0

#### Frequency of each feature type in column 3

```
   2260 3UTR
   3363 5UTR
  21278 CDS
  24556 exon
   2146 start_codon
   2139 stop_codon
   2497 transcript
```

#### Example lines from the new sorted file

```
chr18   refGene 3UTR    10064401        10066047        .       -       .       gene_id "Rock1"; transcript_id "NM_009071"; exon_number "1"; exon_id "NM_009071.1"; gene_name "Rock1";
chr18   refGene 3UTR    10560484        10562941        .       +       .       gene_id "Greb1l"; transcript_id "NM_001083628"; exon_number "33"; exon_id "NM_001083628.33"; gene_name "Greb1l";
chr18   refGene 3UTR    10566512        10567498        .       -       .       gene_id "Esco1"; transcript_id "NM_001081222"; exon_number "1";
```

### Task 6

The jobs were submitted successfully:

```
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
```

Example running jobs:

```
JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
23692 compepi  test_13 ibha0000  R       0:03      1 bibigrid-worker-aorkgacx5hn960m-1
23693 compepi  test_14 ibha0000  R       0:03      1 bibigrid-worker-aorkgacx5hn960m-1
23694 compepi  test_15 ibha0000  R       0:03      1 bibigrid-worker-aorkgacx5hn960m-2
23680 compepi   test_1 ibha0000  R       0:04      1 bibigrid-worker-aorkgacx5hn960m-1
23681 compepi   test_2 ibha0000  R       0:04      1 bibigrid-worker-aorkgacx5hn960m-1
23682 compepi   test_3 ibha0000  R       0:04      1 bibigrid-worker-aorkgacx5hn960m-1
```

### Task 7

The liver and kidney bed-count tables were loaded into R and checked before continuing with the analysis.

#### Q4 — Table dimensions

| Table | Rows | Columns |
|---|---:|---:|
| Liver | 1,362,728 | 6 |
| Kidney | 1,362,728 | 6 |

Both tables have the same dimensions. This is expected because liver and kidney were measured using the same genome-wide 2 kb binning scheme.

#### Q5 — Column names

Both tables have the same column names:

```
"Chr"  "Start"  "End"  "H3K27me3"  "H3K36me3"  "H3K9me3"
```

Each row represents one genomic bin. The first three columns describe its genomic position, while the last three columns contain the signal values for the histone marks.

#### Q6 — Genome length

The genome length was calculated by summing `End - Start` across all bins. Both datasets gave **2,725,456,000 bp** (approximately **2.73 Gb**). The same result for liver and kidney is expected because both datasets cover the same mouse genome bins.

#### Q8 — Dimensions after vertical concatenation

After adding a `cell_type` column and combining the liver and kidney tables with `rbind()`, the resulting data frame had **2,725,456 rows and 7 columns**. The number of rows doubled because both datasets were stacked, and the extra column stores whether each row comes from liver or kidney.

#### Q9 — Dimensions after reshaping to long format

After reshaping the three histone-mark columns from wide to long format, the data frame had **8,176,368 rows and 3 columns**. Each original row now occurs three times: once for H3K27me3, once for H3K36me3, and once for H3K9me3.

#### Q10 and Q11 — Density plots

<img width="2000" height="1200" alt="Density plot, full x-axis range" src="https://github.com/user-attachments/assets/d9677325-ef09-44e2-86ca-d0f905ca2425" />

<img width="2000" height="1200" alt="Density plot, x-axis restricted to 0-100" src="https://github.com/user-attachments/assets/2a281e3c-8f76-4efd-af38-cf6133c613d1" />

The first density plot is strongly right-skewed: most genomic bins have very low signal values, while a small number of bins have high values and form a long right-hand tail. This is expected for genome-wide histone-mark data because these marks are enriched at specific genomic regions rather than uniformly distributed across the genome.

Restricting the x-axis to 0–100 makes the main part of the distributions easier to see. The three marks have different shapes, showing that they are distributed differently across genomic bins. The liver and kidney panels have broadly similar overall distributions, suggesting that the major differences between tissues are likely to be at particular genomic loci rather than a large global shift in all bins.

#### Q16 — Scree plot and cumulative scree plot

<img width="2000" height="1500" alt="Scree plot" src="https://github.com/user-attachments/assets/7bb9c0d6-34e5-4788-a7c4-85f52f4c9c23" />

<img width="2000" height="1500" alt="Cumulative scree plot" src="https://github.com/user-attachments/assets/01fc194c-3780-4e78-a9ab-143eb7f85246" />

Applying the elbow method to both graphs where we take into account the PCA components until the variance becomes stabilized ( no longer changes much between PC), we can take the first 3 PC as together they represent approximately  **90%** of the total variance. 

#### Q17 — Principal-component loadings

The PCA loading matrix has dimensions **6 x 6**. There are therefore **36 loading values** in total: each of the six original variables has one loading for each of the six principal components.

#### Q18 — PC1 and PC2 loadings

| Variable | PC1 | PC2 |
|---|---:|---:|
| H3K27me3_liver | 0.394 | 0.314 |
| H3K36me3_liver | 0.350 | -0.624 |
| H3K9me3_liver | 0.391 | 0.241 |
| H3K27me3_kidney | 0.449 | 0.310 |
| H3K36me3_kidney | 0.386 | -0.572 |
| H3K9me3_kidney | 0.468 | 0.175 |

<img width="2000" height="1555" alt="PC1 versus PC2 loadings" src="https://github.com/user-attachments/assets/83242506-4b03-4cbc-84f2-d2ca39b1bae8" />

All six variables have positive PC1 loadings of a similar size (approximately 0.35–0.47). PC1 therefore represents a general signal-intensity axis rather than separating liver from kidney or one histone mark from another.

PC2 separates H3K36me3 from the two other marks. H3K36me3 has strong negative PC2 loadings in both liver and kidney, whereas H3K27me3 and H3K9me3 have positive PC2 loadings. In this dataset, PC2 mainly distinguishes the H3K36me3 signal pattern from the patterns of H3K27me3 and H3K9me3. The liver and kidney versions of the same mark are located close to each other, indicating that mark type contributes more strongly than tissue type to this PCA separation.

### Genomic Ranges

#### Q2 — Total bases covered

| Method | Liver | Kidney |
|---|---:|---:|
| Sum of all range widths: `sum(width(gr))` | 2,726,818,728 | 2,726,818,728 |
| Non-redundant coverage after `reduce()`: `sum(width(reduce(gr)))` | 2,725,456,021 | 2,725,456,021 |

The raw total adds the width of every range. The reduced total first merges overlapping or adjacent ranges, then calculates the total width. The two values are very similar, indicating that the 2 kb genomic bins are largely non-overlapping. The small difference results from ranges that touch or overlap and are merged by `reduce()`.

#### Q5 — Overlaps between liver and kidney ranges

```
findOverlaps(liver_gr, kidney_gr) -> 4,088,142 hits
pintersect(...) -> GRanges with 4,088,142 overlapping regions
```

The number of hits is larger than the number of bins in one dataset because `findOverlaps()` considers ranges that share a boundary position as overlapping. Since adjacent bins use an end coordinate equal to the next bin's start coordinate, a liver bin can match its corresponding kidney bin and neighbouring bins. This reflects the interval definition used by the function, rather than a biological difference between liver and kidney.

# Day 2

# Pipeline
## 1. Main steps of the pipelines

- Quality control
- Read trimming
- Read alignment
- Deduplicate Alignments
- Extract methylation calls
- Sample report
- Alignment QC

## 2 
- a) The data is single end 
- b) the mapper used is bismark 

# Day 3
## 1.3 Quality control
### fastqc report
- Which read length was used in the sequencing? 
Read length is about 150 bp
- What does Per-base sequence quality tell you? How many samples fail this check? If any, why? 
It tells us the number of reads with average quality scores. Shows if a subset of reads has poor quality. Zero samples fails this test.
- What does Per-base sequence content tell you? How many samples fail this check? If any, why? 
It tells us the proportion of each base position for which each of the four normal DNA bases has been called. All samples actually fail this check due to the bias towards the high propertion of thymine at each position. This makes sense because all the unmethylated cytosines are read as thymines after bisulfite treatment.

### samtools stats for deduplicated reads
- Run samtools flagstats for each BAM file and report how many reads were aligned.
```
liver_14.5.1_trimmed_bismark_bt2.deduplicated.bam : 48411840 + 0 mapped 
liver_14.5.2_trimmed_bismark_bt2.deduplicated.bam : 49586255 + 0 mapped
liver_15.5.1_trimmed_bismark_bt2.deduplicated.bam : 43031256 + 0 mapped 
liver_15.5.2_trimmed_bismark_bt2.deduplicated.bam : 38484289 + 0 mapped 

kidney_14.5.1_trimmed_bismark_bt2.deduplicated.bam : 48836159 + 0 mapped 
kidney_14.5.2_trimmed_bismark_bt2.deduplicated.bam : 44830803 + 0 mapped 
kidney_15.5.1_trimmed_bismark_bt2.deduplicated.bam : 39130857 + 0 mapped 
kidney_15.5.2_trimmed_bismark_bt2.deduplicated.bam : 43080895 + 0 mapped 
```
### Mthylation calls
- Count how many cytosines were called in each sample.

```
methylation_coverage/kidney_14.5.1_trimmed_bismark_bt2.deduplicated.bismark.cov.gz: 2459863
methylation_coverage/kidney_14.5.2_trimmed_bismark_bt2.deduplicated.bismark.cov.gz: 2456188
methylation_coverage/kidney_15.5.1_trimmed_bismark_bt2.deduplicated.bismark.cov.gz: 2480893
methylation_coverage/kidney_15.5.2_trimmed_bismark_bt2.deduplicated.bismark.cov.gz: 2489813

methylation_coverage/liver_14.5.1_trimmed_bismark_bt2.deduplicated.bismark.cov.gz: 2460789
methylation_coverage/liver_14.5.2_trimmed_bismark_bt2.deduplicated.bismark.cov.gz: 2461990
methylation_coverage/liver_15.5.1_trimmed_bismark_bt2.deduplicated.bismark.cov.gz: 2453361
methylation_coverage/liver_15.5.2_trimmed_bismark_bt2.deduplicated.bismark.cov.gz: 2451935
```

### Multiqc report
- Check the mapping efficiencies for all samples. Are the efficiencies reasonable? Is there any substantial difference between the samples? 
Yes the efficinecy are reasonable, almost all reads are aligned uniquely, with little no. of reads are abmiguous or not aligned. This is consitient across all samples in liver and kidney.

- Which sample has the highest read duplication rate? 
`liver_14.5.2` with a precentage of duplicate reads equal 16.2%

- Does kidney or liver tissue have a globally higher CpG methylation rate? 
Kidney tissue samples has more overall methylatio rate at the CpG context with close to 80% compare to liver samples with approx 60% CpG methylation rate.

- Which sample has the lowest mean coverage? 
`kidney_14.5.2` with 28.7X mean coverage
