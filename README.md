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
