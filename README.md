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

# General notes 
- The script mentioned in here are suppose to run in the core environment on the elixir nodes
- All the paths in the script are absolute path so it can be run from anywhere as long as you are on elexir servers where the data are hosted
- We will not attach any code snippets here, code are hosted here: [https://github.com/fragilefort/methylme](https://github.com/fragilefort/methylme). Most scripts are in the scripts directory but you can use any fuzzy findig tool (simply grep) to find any script mentioned in here
# Day 1

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

```
Number of lines

58239

Number of exon entries

24556

Number of exons longer than 1000 bp

1663

Number of exons where gene_id is "Sox17"

0

Count of exons in chr2

0

```

```
# Frequency of each feature type in column 3
   2260 3UTR
   3363 5UTR
  21278 CDS
  24556 exon
   2146 start_codon
   2139 stop_codon
   2497 transcript
```


```
# Example lines from the new sorted file

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

#### Q4  Table dimensions

| Table | Rows | Columns |
|---|---:|---:|
| Liver | 1,362,728 | 6 |
| Kidney | 1,362,728 | 6 |

Both tables have the same dimensions. This is expected because liver and kidney were measured using the same genome-wide 2 kb binning scheme.

#### Q5  Column names

Both tables have the same column names:

```
"Chr"  "Start"  "End"  "H3K27me3"  "H3K36me3"  "H3K9me3"
```

Each row represents one genomic bin. The first three columns describe its genomic position, while the last three columns contain the signal values for the histone marks.

#### Q6 Genome length

The genome length was calculated by summing `End - Start` across all bins. Both datasets gave **2,725,456,000 bp** (approximately **2.73 Gb**). The same result for liver and kidney is expected because both datasets cover the same mouse genome bins.

#### Q8 Dimensions after vertical concatenation

After adding a `cell_type` column and combining the liver and kidney tables with `rbind()`, the resulting data frame had **2,725,456 rows and 7 columns**. The number of rows doubled because both datasets were stacked, and the extra column stores whether each row comes from liver or kidney.

#### Q9 Dimensions after reshaping to long format

After reshaping the three histone-mark columns from wide to long format, the data frame had **8,176,368 rows and 3 columns**. Each original row now occurs three times: once for H3K27me3, once for H3K36me3, and once for H3K9me3.

#### Q10 and Q11 Density plots

<img width="2000" height="1200" alt="Density plot, full x-axis range" src="https://github.com/user-attachments/assets/d9677325-ef09-44e2-86ca-d0f905ca2425" />

<img width="2000" height="1200" alt="Density plot, x-axis restricted to 0-100" src="https://github.com/user-attachments/assets/2a281e3c-8f76-4efd-af38-cf6133c613d1" />

The first density plot is strongly right-skewed: most genomic bins have very low signal values, while a small number of bins have high values and form a long right-hand tail. This is expected for genome-wide histone-mark data because these marks are enriched at specific genomic regions rather than uniformly distributed across the genome.

Restricting the x-axis to 0–100 makes the main part of the distributions easier to see. The three marks have different shapes, showing that they are distributed differently across genomic bins. The liver and kidney panels have broadly similar overall distributions, suggesting that the major differences between tissues are likely to be at particular genomic loci rather than a large global shift in all bins.

#### Q16 Scree plot and cumulative scree plot

<img width="2000" height="1500" alt="Scree plot" src="https://github.com/user-attachments/assets/7bb9c0d6-34e5-4788-a7c4-85f52f4c9c23" />

<img width="2000" height="1500" alt="Cumulative scree plot" src="https://github.com/user-attachments/assets/01fc194c-3780-4e78-a9ab-143eb7f85246" />

Applying the elbow method to both graphs where we take into account the PCA components until the variance becomes stabilized ( no longer changes much between PC), we can take the first 3 PC as together they represent approximately  **90%** of the total variance. 

#### Q17 Principal-component loadings

The PCA loading matrix has dimensions **6 x 6**. There are therefore **36 loading values** in total: each of the six original variables has one loading for each of the six principal components.

#### Q18 PC1 and PC2 loadings

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

#### Q2 Total bases covered

| Method | Liver | Kidney |
|---|---:|---:|
| Sum of all range widths: `sum(width(gr))` | 2,726,818,728 | 2,726,818,728 |
| Non-redundant coverage after `reduce()`: `sum(width(reduce(gr)))` | 2,725,456,021 | 2,725,456,021 |

The raw total adds the width of every range. The reduced total first merges overlapping or adjacent ranges, then calculates the total width. The two values are very similar, indicating that the 2 kb genomic bins are largely non-overlapping. The small difference results from ranges that touch or overlap and are merged by `reduce()`.

#### Q5 Overlaps between liver and kidney ranges

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

## 2. Config file
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
### Exploratory Analysis using IGV
- Q4)
  - No the only chromosomes covered are chromosomes 18 and 19
  - The heatmap for unzoomed look:
  <img width="1781" height="788" alt="Screenshot 2026-08-03 122645" src="https://github.com/user-attachments/assets/f2e159e1-3f6a-498f-908e-073195004638" />
  - The general state of The genome-wide IGV view shows that the DNA methylation landscape is generally highly methylated across the chromosome in both liver and kidney samples. Most bedGraph signals are in the upper part of the 0–100% methylation scale, which is expected for CpG methylation in differentiated mammalian tissue. However, methylation is not uniform: there are many local regions with lower signal, visible as valleys in the tracks. These low-methylation regions may correspond to regulatory regions such as CpG islands and promoters. The overall patterns appear broadly similar among replicates however from the heatmap point of view, we can see that generally the kidney have higher methylation levels than the levels of the liver a little bit also specific local differences need to be examined at higher zoom levels.
  - As we can see in those 2 screenshots we found both a hypomethylated and hypermethylated promoter respectively:
  <img width="1626" height="818" alt="Screenshot 2026-08-03 123058" src="https://github.com/user-attachments/assets/e37c85a3-f0c5-40d4-b85d-c84e8ca600e6" />
  <img width="1787" height="751" alt="Screenshot 2026-08-03 124946" src="https://github.com/user-attachments/assets/6b8b6f21-c476-471c-9c78-d67f950c4e84" />
### Running the pipeline
- How does RnBeads compute differentially methylated sites and regions?
  - RnBeads computes differentially methylated sites (DMCs) by comparing methylation levels at each CpG between sample groups defined in the sample annotation—for example, tissue groups. For each CpG, it calculates the mean methylation in each group, the difference between group means, a log2 methylation quotient, and a p-value for the group comparison. By default, RnBeads uses a **two-sided Welch t-test**; it can alternatively use limma linear modelling, including adjustment covariates. It corrects site-level p-values for multiple testing using false-discovery-rate (FDR) adjustment.
  - RnBeads computes differentially methylated regions (DMRs) from CpG-level results. First, it assigns CpGs to predefined genomic annotations, such as genes, promoters, CpG islands, or tiling regions. For each region, it summarizes the group methylation levels across the CpGs in the region and reports the mean methylation difference and log2 quotient.To evaluate regional significance, RnBeads combines the p-values of the CpGs in each region using a **generalized Fisher method**, then calculates an FDR-adjusted regional p-value. The output also includes the number of CpGs contributing to that region and, for sequencing datasets, coverage-related summaries. Therefore, robust DMRs are regions with a meaningful methylation difference, low regional FDR, adequate coverage, and support from multiple CpGs rather than one isolated CpG.
  - RnBeads uses Welch’s t-test at individual CpGs because it compares mean methylation between two groups without assuming equal variances. This is important because methylation variability can differ between groups and between CpG sites. For DMRs, RnBeads does not run Fisher’s test to assume equal variance. Instead, it combines the unadjusted CpG-level p-values within each genomic region using a generalized Fisher method to produce one regional p-value. This increases power when several CpGs in the same region show concordant evidence of differential methylation, while accounting for the expected correlation among CpGs. The method also separately summarizes regional methylation effect sizes, such as the mean difference between groups.

### RnBeads 
#### Data import
- Which genome assembly was used? mm10 
- How many CpG sites with methylation calls were annotated in the dataset? 1248150
#### Quality control 
- Look at the coverage information. Were there any samples that had signifciantly fewer CpGs covered or the read depth per CpG was significantly lower? No, none of the samples show significantly lower CpG coverage or read depth

#### Preprocessing 
- What flitering steps were applied? 
    * SNP Removal: Removed 85,372 sites that overlap with single-nucleotide polymorphisms (SNPs).
    * High Coverage Outliers: Checked for sites exceeding 50 times the 0.95-quantile of coverage values (0 sites removed).
    * Low Coverage Masking: Masked 394,214 sites with a coverage depth of less than 10 by setting them to NA.
    * Sex Chromosomes: Checked for sites located on sex chromosomes (0 sites removed).
    * Missing Value Filtering: Removed 117,698 sites containing any missing values (more than 0% missing across samples).
- How many CpGs survived the filtering? Out of 1,248,150 total sites, 203,070 sites were removed, leaving 1,045,080 sites
Yes, the removed CpGs have a massive spike at zero methylation (0% methylated), whereas the retained CpGs are mostly highly methylated (around 90% methylated).

#### Differential methylation 
- Inspect difefrential methylation between kidney and liver. Which tissue generally exhibits higher methylation levels? Kidney generally exhibits higher methylation levels.
- In which type of genomic regions do you see the largest differences? Promoters
- Inspect the table of difefrentially methylated promoters. What does the combinedRank column mean? The combinedRank column evaluates overall differential methylation by assigning a promoter its worst (highest) rank across absolute difference, relative fold-change, and statistical significance, ensuring top-ranked regions is doing well across all three metrics.

- List the names of the top 5 differentially methylated promoters (according to combinedRank). Are they hypomethylated or hypermethylated in kidney compared to liver? 

1. 5830416P10Rik (ENSMUSG00000097636)
2. Olfr1425 (ENSMUSG00000067526)
3. Anxa1 (ENSMUSG00000024659)
4. Ms4a6c (ENSMUSG00000079419)
5. Cyp26c1 (ENSMUSG00000062432)

Four of the top five promoters (5830416P10Rik, Olfr1425, Anxa1, and Ms4a6c) are hypermethylated in kidney compared to liver, while one promoter (Cyp26c1) is hypomethylated in kidney compared to liver.


# Day 4
## Differential methylation
- Sort the resulting table by the combinedRank column and note the coordinates of these regions.

```
--- Top 10 Ranked Mystery Regions ---
     Chromosome    Start      End Strand combinedRank mean.mean.diff
5236      chr19 40536439 40537317      +           80      0.6731252
5920      chr19 53464400 53465046      +           96      0.7303905
3239      chr18 82569663 82569913      +          100      0.8010435
4937      chr19 32487643 32488124      +          115      0.7295919
4307      chr19 17014480 17015015      +          118      0.6360558
4899      chr19 32200198 32200547      +          122      0.7336288
5028      chr19 34844841 34845136      +          132      0.7288215
2220      chr18 61789026 61789440      +          135      0.6272329
4474      chr19 21897242 21897659      +          142      0.6230632
4012      chr19  9938564  9939500      +          149      0.6196579
     comb.p.adj.fdr
5236   0.0003747531
5920   0.0004823422
3239   0.0001302144
4937   0.0003747531
4307   0.0005185470
4899   0.0005254040
5028   0.0003747531
2220   0.0005294051
4474   0.0005185470
4012   0.0004823422
```


# 3.3 Exploratory Analysis

## PCA and scree plot:

<img width="2700" height="1800" alt="PCA_scree_plot" src="https://github.com/user-attachments/assets/ba089df7-8a83-441c-a693-caa13828612d" />

<img width="3000" height="2100" alt="PCA_top_10000_variable_CpGs" src="https://github.com/user-attachments/assets/1060628f-b1a6-493c-8c8b-658795ea5392" />

## Density plot comparing methylation levels in the different sample groups:


<img width="3000" height="2100" alt="Density_methylation_by_tissue" src="https://github.com/user-attachments/assets/02ecc4a1-748a-43c8-9b43-83f0f4d5fc49" />


<img width="3300" height="2100" alt="Density_methylation_tissue_by_timepoint" src="https://github.com/user-attachments/assets/bca8612d-0dc9-4566-9408-afedc0a34a94" />

- Yes the distribution is different as we can see that the Kidney has a higher methylation than the liver and also the difference in time-points didn't create any major difference in the distribution.

## Density plot comparing methylation levels in the different GC-contexts:


<img width="3000" height="2100" alt="Density_methylation_by_CpG_context" src="https://github.com/user-attachments/assets/d3dd599b-a219-4a65-b337-3803ed0f899d" />

- Yes we can see distribution is different as we can see that at 0 beta methylation value the highest density is for the islands while the shores and shelf has low density and as we increase the beta methylation value the density for islands decreases until we have at 1 beta methylation value, the shelf then shores then islands with highest distribution respectively.

## Heatmap for the 1,000 most variable sites:


<img width="3000" height="3600" alt="Heatmap_top_1000_variable_CpGs" src="https://github.com/user-attachments/assets/2f80d8b1-bda0-4010-a262-6b75aa65812d" />

## Heatmap for the 100 most variable of these mystery regions:


<img width="3000" height="3000" alt="Heatmap_top_100_variable_mystery_regions" src="https://github.com/user-attachments/assets/38b78ccf-444b-44de-86dd-4a0e8c1f41f2" />

## Boxplots to compare the samples to each other:

<img width="3000" height="2100" alt="Boxplot_global_methylation_by_group" src="https://github.com/user-attachments/assets/b56b1f6b-2e10-49de-b607-9d8484393984" />

<img width="4200" height="2100" alt="Boxplot_global_methylation_by_sample" src="https://github.com/user-attachments/assets/95814f41-8aac-4972-b96a-bd2e625db756" />

- Yes we can see that also like density plot, we have a higher methylation level for kidneys with respect to the liver.

## 9

## 10

## 11

## Comparing the segmentation to the raw methylation signal:


<img width="1803" height="697" alt="task 3 3 12" src="https://github.com/user-attachments/assets/970ad10a-1b1f-4080-87b2-01c2e0dd5016" />

<img width="1777" height="771" alt="task 3 3 12 (2)" src="https://github.com/user-attachments/assets/8bfe6762-292c-4b92-9693-5d9ff4675ccc" />

- As we can see also here it confirms that Kidney has slightly higher methylation than liver and the raw signal corresponds to raw signal files especially with UMR and HMR



# 3.4 Exploring differential methylation

## Explore differential methylation in IGV. Explore the top 5 regions you identified

- Mystery region 1
  - a) The region appears to be hypomethylated in kidney compared with liver. The liver track shows a clear methylation signal over the region, whereas the kidney track has little or no signal.
  - b) The region is intergenic, because no RefSeq gene overlaps the displayed region.
  - c) No annotated RefSeq gene is directly overlapping the region in this view.
  - d) In liver, the region overlaps an LMR (low-methylated region). In kidney, it lies within a PMD (partially methylated domain).

- Mystery region 2
  - a) The region appears to be more methylated in kidney than in liver overall, although the liver signal changes across the displayed interval.
  - b) The region is located in the gene body of Mirt1.
  - c) The nearby and overlapping gene is Mirt1.
  - d) In liver, the mystery region overlaps an LMR. In kidney, the whole displayed locus lies in a PMD.

- Mystery region 3
  - a) The region appears to be hypomethylated in kidney compared with liver.
  - b) The region is located in the gene body of Mbp. The gene is on the reverse strand, as indicated by the left-facing arrows, but the mystery region still falls within its genomic span.
  - c) The nearby and overlapping gene is Mbp.
  - d) In liver, the region overlaps an HMR (high-methylated region). In kidney, it falls in a PMD.

- Mystery region 4
  - a) The region appears to be hypomethylated in kidney compared with liver. Liver has a stronger methylation signal across the locus, while kidney shows little signal.
  - b) The region is in the gene body of Minpp1.
  - c) The nearby and overlapping gene is Minpp1.
  - d) In liver, the region overlaps an HMR. It is close to the transition from a UMR to an HMR, but the blue mystery interval itself is mainly within the HMR. In kidney, it lies in a PMD.

- Mystery region 5
  - a) The region appears to be hypomethylated in kidney compared with liver.
  - b) The region is in the gene body of Prune2.
  - c) The nearby and overlapping gene is Prune2.
  - d) In liver, the region overlaps an LMR. In kidney, it lies within a PMD.
- Overall interpretation:
  - The five mystery regions are tissue-specific DMRs. Region 1 is an intergenic LMR and is therefore a candidate distal enhancer-like DMR. Regions 2 and 5 are intragenic LMRs and may represent intragenic enhancer-like elements. Regions 3 and 4 are gene-body DMRs overlapping HMRs and are less likely to be enhancer-associated. None of the regions clearly overlaps a promoter, since none lies at an annotated transcription start site or within a UMR.

# Summary questions:

## 1. Is normalization performed for WGBS data? Explain your answer.

WGBS data are not normalized in the same way as RNA-seq count data. RNA-seq commonly uses library-size normalization because the main measurement is the number of reads assigned to a gene. In WGBS, methylation is usually measured at each CpG as a proportion of methylated reads:

\[
\text{Methylation level} = \frac{\text{methylated reads}}{\text{methylated reads} + \text{unmethylated reads}}
\]

For example, if 8 out of 10 reads at a CpG are methylated, the estimated methylation level is 80%. This calculation already accounts for the local number of reads at that CpG.

However, WGBS data still require quality control and processing. CpGs with very low coverage are normally filtered because a methylation estimate based on one or two reads is unreliable. During differential methylation analysis, the coverage at each CpG should also be considered: 50% methylation based on 2 reads is much less reliable than 50% methylation based on 100 reads. It may also be necessary to correct for technical factors such as batch effects or incomplete bisulfite conversion.

Therefore, WGBS does not usually use a simple global read-depth normalization such as CPM. Instead, it uses methylation proportions, coverage filtering, and statistical models that account for methylated and unmethylated read counts.

## 2. What are the mystery regions?

The mystery regions are differentially methylated regions (DMRs) between kidney and liver. Based on their location relative to genes and their methylation-segmentation category, they can be interpreted as follows:


Region 1 is the strongest candidate for a distal enhancer-like DMR because it is intergenic and overlaps an LMR in liver. Regions 2 and 5 are located within gene bodies but also overlap liver LMRs. Since enhancers can occur inside genes, these are possible intragenic enhancer-like DMRs.

Regions 3 and 4 are located in gene bodies and overlap liver HMRs. They are better described as gene-body DMRs than as enhancers. None of the mystery regions clearly overlaps an annotated transcription start site or a UMR, so none is a strong candidate for a classical promoter DMR.

The enhancer-like labels are predictions based on methylation patterns and genomic position. Confirmation would require additional evidence, such as chromatin accessibility, enhancer-associated histone marks, gene-expression correlation, or chromatin-contact data.

## 3. Explain what UMR, LMR, HMR, and PMD are.

### UMR: Unmethylated Region

A UMR is a genomic region with very low or nearly absent DNA methylation. UMRs are often CpG-rich and commonly occur at CpG islands and active gene promoters.

### LMR: Low-Methylated Region

An LMR has low, but not completely absent, DNA methylation. LMRs are often associated with distal regulatory elements, including transcription-factor binding sites and enhancer-like regions. However, an LMR alone is not sufficient to prove that a region is an enhancer.

### HMR: Highly Methylated Region

An HMR is a region with high DNA methylation. These regions make up much of the genomic background outside unmethylated promoters, CpG islands, and some regulatory elements. HMRs can occur in gene bodies or intergenic regions.

### PMD: Partially Methylated Domain

A PMD is a large genomic domain with intermediate or variable methylation rather than consistently high methylation. PMDs are often associated with heterochromatin, low gene density, low gene expression, and cellular proliferation history. A DMR inside a PMD should be interpreted carefully because it may reflect a broad domain-level methylation difference rather than a local regulatory event.


# Integrative analysis using IGV


## Question 1: General methylation states of promoters and gene bodies

In general, active gene promoters are lowly methylated or unmethylated, whereas the bodies of actively transcribed genes are often highly methylated.

At active promoters, low DNA methylation allows transcription factors and RNA polymerase II to bind DNA. These promoters commonly overlap unmethylated regions (UMRs), CpG islands, active promoter chromatin states, and H3K4me3 peaks. In contrast, promoter hypermethylation is generally associated with reduced accessibility and transcriptional silencing.

Gene bodies frequently show high DNA methylation even when the gene is actively expressed. This is not contradictory: gene-body methylation can suppress inappropriate internal transcription start sites and is often found together with the transcription-elongation mark H3K36me3.


| Genomic region | Typical methylation state in an active gene | Main interpretation |
|---|---|---|
| Promoter / transcription start site | Low or unmethylated | Supports transcription-factor and RNA polymerase binding |
| Gene body | Often highly methylated | Compatible with active transcription and may suppress cryptic internal initiation |
| Silenced promoter | Often highly methylated | Usually associated with closed chromatin and reduced transcription |

<img width="1536" height="710" alt="q1) liver_specific_diff_meth" src="https://github.com/user-attachments/assets/9fd48a0d-93a2-49fd-b7b2-a2f38997baae" />

<img width="1536" height="710" alt="q1) liver_specific_diff_meth_gene_2" src="https://github.com/user-attachments/assets/431121f2-d51b-4b19-8e3d-1e9deade99e4" />

<img width="1536" height="710" alt="q1) liver_specific_diff_meth_gene_3" src="https://github.com/user-attachments/assets/1c18dd16-3026-4c91-980d-47bf456a3ba0" />


---

## Question 2: Overlap between MethylSeekR and ChromHMM states

The DNA-methylation segmentation produced by MethylSeekR is consistent with the chromatin-state segmentation generated by ChromHMM.

| MethylSeekR state | Typical overlapping ChromHMM state | Interpretation |
|---|---|---|
| UMR: unmethylated region | Active promoter (`7_Pr_A`) or flanking promoter (`6_Pr_F`) | Active or potentially active promoter, often near a transcription start site |
| LMR: low-methylated region | Active enhancer (`8_Enh_A`) or weak enhancer (`9_Enh_W`) | Candidate distal or intragenic regulatory element |
| HMR: highly methylated region | Strong transcription (`12_Tx_S`), weak transcription (`11_Tx_W`), or transcription transition (`10_Tx_T`) | Often found in transcribed gene bodies |
| PMD: partially methylated domain | No signal / quiescent (`15_NS`) or weak heterochromatin (`14_Het_W`) | Broad inactive or weakly active chromatin domain |

UMRs were observed at promoter regions, LMRs at enhancer-like regions, HMRs across transcribed gene bodies, and PMDs in broad inactive or heterochromatic regions. This agreement supports the biological consistency of the segmentation tracks.

---

## Question 3: ATAC-seq peaks and chromatin states

ATAC-seq peaks preferentially overlap open and regulatory chromatin states. The clearest overlap is with active promoters and enhancer states.

- Strong ATAC-seq peaks were observed at active and flanking promoters (`7_Pr_A` and `6_Pr_F`). These regions are accessible because transcription factors and RNA polymerase need access to the DNA.
- ATAC-seq peaks also overlap active and weak enhancer states (`8_Enh_A` and `9_Enh_W`). These peaks often occur at tissue-specific regulatory regions.
- In contrast, ATAC-seq signal is weak or absent in transcribed gene-body states (`10_Tx_T`, `11_Tx_W`, and `12_Tx_S`) and in repressed or quiescent chromatin (`14_Het_W` and `15_NS`).

Therefore, chromatin accessibility is concentrated at promoters and enhancers rather than being distributed uniformly across gene bodies.

---

## Question 4: Examples of unexpressed genes

Several examples of unexpressed genes were identified. These genes show little or no RNA-seq signal and are associated with repressive chromatin features.

<img width="1536" height="710" alt="q4) gene1_example" src="https://github.com/user-attachments/assets/0e86b468-0750-4a71-84ee-0de5db477caf" />

<img width="1536" height="710" alt="q4) gene1_example_rest_of_png" src="https://github.com/user-attachments/assets/96bdb286-2be4-4a08-8740-b674764335c6" />

<img width="1536" height="710" alt="q4) gene2_example" src="https://github.com/user-attachments/assets/0e5fae02-e369-4426-b64f-7fb18a8521c0" />

<img width="1536" height="710" alt="q4) gene2_example_rest_of_png" src="https://github.com/user-attachments/assets/58d57539-09a0-43b3-a1dc-295084921d8c" />

<img width="1536" height="710" alt="q4) gene3_example" src="https://github.com/user-attachments/assets/bfca6821-89eb-4c41-b2f6-97b8c2b5d95a" />

<img width="1536" height="710" alt="q4) gene3_example_rest_of_png" src="https://github.com/user-attachments/assets/abd4d34a-856a-495f-af1e-e91ac01d42e9" />


### Example 1: *Mro*

*Mro* has no detectable RNA-seq signal in the inspected liver tracks. Its promoter and gene body are highly methylated. The locus is enriched for repressive histone marks, including H3K27me3 and H3K9me3, while active marks such as H3K27ac and H3K36me3 are absent or near background. The promoter is classified as a bivalent / poised promoter state (`4_Pr_B`), indicating the coexistence of H3K4me3 with repressive H3K27me3. This state is compatible with a gene that is poised but currently silent.

### Example 2: *Mapk4*

*Mapk4* is also not expressed in the inspected liver sample. Its promoter is unmethylated, but the locus is classified as bivalent / poised (`4_Pr_B`) and is enriched for H3K27me3. H3K27ac is absent, so the promoter lacks the active acetylation mark needed for active transcription. This shows that unmethylated promoter DNA alone is not sufficient for expression when repressive chromatin marks are present.

### Example 3: *Gm32661* and *Gm32801*

These loci are unexpressed and occur in broad repressive chromatin. They show enrichment of H3K9me3 and H3K27me3, little or no active-mark signal, low or absent ATAC-seq accessibility, and PMD / quiescent or weak-heterochromatin segmentation. These features are consistent with heterochromatic silencing.

| Gene / locus | RNA-seq expression | Promoter / body methylation | Major chromatin features | Interpretation |
|---|---|---|---|---|
| *Mro* | Not detected | Promoter and body highly methylated | H3K27me3 and H3K9me3; bivalent / poised promoter | Repressed or poised silent gene |
| *Mapk4* | Not detected | Promoter unmethylated; body highly methylated | H3K4me3 plus H3K27me3; no H3K27ac; bivalent promoter | Poised but inactive gene |
| *Gm32661* / *Gm32801* | Not detected | Broad PMD-associated methylation pattern | H3K9me3 and H3K27me3; no active marks; no ATAC peaks | Heterochromatin-associated silencing |

---

## Question 5: Associations between epigenetic layers and expression

The data show clear associations between DNA methylation, chromatin accessibility, histone modifications, chromatin state, and gene expression. These associations depend strongly on genomic context.

### Promoters

At promoters, DNA methylation is generally negatively associated with expression. Active promoters are usually unmethylated, accessible by ATAC-seq, enriched for H3K4me3 and H3K27ac, and classified as active promoter or flanking-promoter chromatin. Inactive promoters often show hypermethylation and/or repressive histone marks such as H3K27me3 and H3K9me3.

### Gene bodies

Gene-body methylation has a different relationship with expression. Actively transcribed genes can have highly methylated bodies and strong H3K36me3 signal. Thus, high methylation within a gene body does not necessarily imply that the gene is inactive. ATAC-seq accessibility is generally lower across gene bodies than at promoters and enhancers.

### Enhancers

Enhancers tend to show low or intermediate methylation, ATAC-seq accessibility, and enhancer-associated histone marks such as H3K4me1 and H3K27ac. These features are positively associated with expression of their target genes, although target-gene assignment requires additional evidence.

| Region type | DNA methylation | ATAC-seq accessibility | Typical histone marks | Expected expression state |
|---|---|---|---|---|
| Active promoter | Low / unmethylated | High | H3K4me3 and H3K27ac | High |
| Silenced promoter | Often high, or repressed despite low methylation | Low | H3K27me3 and/or H3K9me3 | Low or absent |
| Active gene body | Often high | Usually low | H3K36me3 | High |
| Enhancer | Low or intermediate | High | H3K4me1, often H3K27ac | Associated with active target-gene regulation |
| Heterochromatin / PMD | Variable or partially methylated | Low | H3K9me3 and/or H3K27me3 | Low or absent |

A useful example is tissue-specific regulation of *Tnfaip8*: the liver locus has promoter / enhancer accessibility and active chromatin features, whereas the kidney locus is relatively closed and has lower expression.

---

## Question 6: Example of multi-omic regulation

The **Tcf4 locus** on chromosome 18 provides an example of tissue-specific regulation supported by differential expression, DNA methylation, chromatin accessibility, and chromatin segmentation.

The `DEGs_sig_01FDR_2FC.bed` track marks **Tcf4** as a differentially expressed gene. Around the displayed locus, the WGBS tracks show different methylation profiles in kidney and liver. In particular, the region around approximately 69.794–69.795 Mb has a much stronger kidney methylation signal than liver, and it overlaps one of the DMRs in the `mystery_diffmeth_sorted.bed` track. Thus, this is a differentially methylated regulatory region rather than simply a methylation difference elsewhere in the genome.

The chromatin data also differ strongly between tissues. In liver, the complete interval is classified as `15_NS` (no signal / quiescent) in the 15-state ChromHMM segmentation, and there are no liver-specific ATAC-seq peaks in this view. In kidney, however, a kidney-specific ATAC-seq peak overlaps the DMR. The same area is classified as a flanking-promoter state (`6_Pr_F`) and adjacent weak-enhancer states (`9_Enh_W`) in the kidney 15-state segmentation.

Together, these observations indicate that the region has a kidney-specific regulatory chromatin configuration: it is accessible in kidney and is annotated as promoter-flanking / enhancer-like chromatin, whereas the corresponding region is quiescent in liver. The overlapping DMR suggests that DNA methylation differs between the tissues at the same regulatory locus. Since Tcf4 is also differentially expressed, this provides a candidate example of coordinated regulation by DNA methylation and chromatin accessibility.


<img width="1536" height="710" alt="igv_6_gene1" src="https://github.com/user-attachments/assets/85046b31-6a62-4aed-a1cd-6eaf2e4c521a" />

---


# Summarize signals over regions of interest

- a) Mind the gene direction and use 1500bp downstream and 500bp upstream of the TSS. Explain in
1–2 sentences how your script handles coordinate calculation for genes on the minus (−) strand
versus the plus (+) strand to prevent 0-length or reversed interval errors

For plus-strand genes, the transcription start site (TSS) is located at the gene start coordinate (`$2`). The promoter interval is therefore defined as 500 bp upstream and 1,500 bp downstream of the TSS:

```text
start = $2 - 500
end   = $2 + 1500
```

For minus-strand genes, the TSS is located at the gene end coordinate (`$3`). Since transcription proceeds toward decreasing genomic coordinates on the minus strand, the promoter interval is calculated as:

```text
start = $3 - 1500
end   = $3 + 500
```

The script replaces negative start coordinates with `0` and writes an interval only if `end > start`. This prevents invalid genomic coordinates, zero-length intervals, and reversed intervals.

  b) Combine computeMatrix with plotHeatmap/plotProfile from deepTools22 to generate plots. e.g.
histone mark/accessibility/DNA methylation signals across all/some of the above mentioned ROI. What
do you observe? Explain your findings from a gene-regulation perspective

<img width="1450" height="693" alt="image" src="https://github.com/user-attachments/assets/f9b587f8-0d64-4295-aa85-e8fdd920796c" />

<img width="1422" height="693" alt="image" src="https://github.com/user-attachments/assets/52814ddd-533e-4242-952e-445b1bcd0497" />

<img width="1422" height="693" alt="image" src="https://github.com/user-attachments/assets/d49a70e7-d83d-4a27-9661-29931d1fc621" />


---

## 1. Promoter profiles in liver and kidney

The promoter heatmaps and average profiles show a clear signal centered at the TSS. ATAC-seq has a narrow, sharp accessibility peak at or very close to the TSS. This represents a nucleosome-free region, where DNA is accessible to transcription factors and RNA polymerase II.

H3K4me3 is enriched around the promoter, with a broader peak shifted slightly downstream of the TSS. This pattern is consistent with H3K4me3 on the +1 nucleosome immediately downstream of an active promoter. H3K4me3 is therefore associated with transcription initiation and active promoter chromatin.

The WGBS profile shows a methylation minimum at the TSS. This local reduction in methylation corresponds to an unmethylated region (UMR), which is expected at active promoters because methylation can interfere with transcription-factor binding and transcription initiation. Moving downstream from the TSS into the gene body, the methylation level rises again. This is consistent with gene-body hypermethylation, which can coexist with active transcription.

The same basic promoter architecture is visible in both tissues, but the liver profile has stronger mean ATAC-seq and H3K4me3 signal than kidney in these plots. The stronger liver signals are consistent with more widespread or stronger promoter activity in the profiled fetal liver samples. .

### Promoter interpretation

| Signal | Pattern around the TSS | Gene-regulation interpretation |
|---|---|---|
| ATAC-seq | Sharp peak at the TSS | Open, nucleosome-depleted promoter that permits transcription-factor and polymerase binding |
| H3K4me3 | Broad enrichment around the TSS, often slightly downstream | Active promoter and transcription-initiation-associated chromatin |
| WGBS methylation | Valley at the TSS, increasing into the gene body | Unmethylated promoter; methylated gene body is compatible with transcription |

So we can conclude by saying active promoters are characterized by the convergence of low DNA methylation, high accessibility, and H3K4me3 enrichment at the TSS.

---

## 2. Kidney enhancer profiles

The kidney enhancer heatmaps show the expected enhancer-associated combination of broad H3K4me1 enrichment, localized H3K27ac signal, reduced methylation, and ATAC-seq accessibility.

H3K4me1 forms a broad signal across the enhancer region. This mark identifies enhancer-associated chromatin and can be present at both primed and active enhancers. H3K27ac is more concentrated near the enhancer center and distinguishes active enhancers from merely poised ones. Therefore, regions with both H3K4me1 and H3K27ac are likely active kidney enhancers.

The WGBS signal shows a local decrease in DNA methylation at enhancer centers, consistent with low-methylated regions (LMRs). This may occur because transcription-factor binding and open chromatin are associated with localized loss or exclusion of DNA methylation.

The ATAC-seq profile is centered on the enhancer. The profile may show accessible shoulders surrounding a lower central signal. This shape is compatible with transcription-factor footprinting: a bound factor can protect its exact binding site from Tn5 insertion, while nearby DNA remains accessible.

### Kidney enhancer

| Signal | Observed pattern | Gene-regulation interpretation |
|---|---|---|
| H3K4me1 | Broad enrichment across enhancer regions | Marks primed or active enhancer chromatin |
| H3K27ac | More localized enrichment near the center | Indicates active enhancer activity |
| WGBS methylation | Local methylation decrease / LMR | Compatible with regulatory DNA and transcription-factor binding |
| ATAC-seq | Central accessibility with possible footprint pattern | Open enhancer chromatin; possible factor occupancy at the center |

So we can conclude by saying kidney enhancer ROIs have a regulatory signature of H3K4me1, H3K27ac, chromatin accessibility, and local hypomethylation. These features are consistent with enhancer activity in kidney.

---

## 3. Liver enhancer profiles

The liver enhancer plots show the same general enhancer architecture, but with very strong H3K27ac signal in the highest-intensity rows of the heatmap. Such a strong enrichment is consistent with a subset of highly active liver enhancers and may include super-enhancer-like regions. These regions are candidates for regulating genes important for liver identity and fetal liver function.

As in kidney, H3K4me1 is broad across the enhancer domain, whereas H3K27ac is more focused at the active enhancer core. WGBS shows lower methylation at enhancer centers, consistent with LMRs. The ATAC-seq profile shows central accessibility and can contain a central dip with accessible flanking shoulders, again consistent with transcription-factor occupancy at the central motif.

This pattern is biologically plausible for liver-specific regulatory elements. Liver transcription factors, such as HNF4A, FOXA-family factors, or CEBP-family factors, could bind enhancer DNA, promote accessible chromatin, and recruit co-activators that deposit H3K27ac.

### Liver enhancer

| Signal | Observed pattern | Gene-regulation interpretation |
|---|---|---|
| H3K4me1 | Broad enhancer-associated profile | Primed or active enhancer domain |
| H3K27ac | Strong, focused signal; especially intense in the top heatmap rows | Active enhancer core; strongest regions may be super-enhancer-like |
| WGBS methylation | Local methylation valley / LMR | Regulatory DNA with low methylation, consistent with factor binding |
| ATAC-seq | Accessible flanks and possible central footprint | Open chromatin and possible transcription-factor binding |

So we can conclude by saying liver enhancer ROIs are enriched for features of active enhancers. The very strong H3K27ac signal at a subset of loci suggests particularly active liver regulatory elements.

---

## 4. Gene-regulation perspective

The three deepTools analyses support a coherent model of gene regulation:

```text
Active promoter:
low DNA methylation + ATAC-seq accessibility + H3K4me3
→ transcription initiation

Active enhancer:
LMR / local hypomethylation + ATAC-seq accessibility + H3K4me1 + H3K27ac
→ tissue-specific activation of target genes

Gene body:
increased DNA methylation downstream of the TSS, often with H3K36me3
→ compatible with transcription elongation
```

Promoters are generally shared regulatory units required for transcription initiation, so both tissues show the same basic promoter profile. Enhancers provide more tissue specificity: enhancer activity differs according to which transcription factors, chromatin marks, and accessible regions are present in each tissue.

The inverse relationship between WGBS methylation and ATAC-seq accessibility is especially clear at regulatory elements. Promoters and enhancers that are accessible tend to be locally hypomethylated, while methylated chromatin tends to be less accessible. This relationship is an association rather than proof that methylation alone causes accessibility changes; both can be shaped by transcription-factor binding and chromatin-remodeling activity.

---
