#!/bin/bash
# =============================================================================
# RNA-seq processing: Sp8 KO/WT, Sp6 KO/WT, DKO/WT (mouse ectoderm, mm10).
#
# Expected in input/:
#   <prefix>_1/2.fastq.gz for all 18 samples below (Sp8 SRA accessions,
#   Sp6 and DKO local prefixes)
#   GRCm38_mouse (HISAT2 index, files alongside this prefix)
#   gencode.vM25.basic.annotation.gtf
#
# Tools: fastp, HISAT2, samtools, deepTools (bamCoverage), StringTie,
# prepDE.py3 (https://ccb.jhu.edu/software/stringtie/dl/prepDE.py3 — put it
# on PATH or next to this script).
#
# Samples (raw prefix -> label):
#   SRR11015074 -> RNASEQ_ect_sp8_WT_rep1   SRR11015077 -> RNASEQ_ect_sp8_KO_rep1
#   SRR11015075 -> RNASEQ_ect_sp8_WT_rep2   SRR11015078 -> RNASEQ_ect_sp8_KO_rep2
#   SRR11015076 -> RNASEQ_ect_sp8_WT_rep3   SRR11015079 -> RNASEQ_ect_sp8_KO_rep3
#   rnasq_sp6_wt_4 -> RNASEQ_ect_sp6_WT_rep1   rnasq_sp6_mt_4 -> RNASEQ_ect_sp6_KO_rep1
#   rnasq_sp6_wt_5 -> RNASEQ_ect_sp6_WT_rep2   rnasq_sp6_mt_5 -> RNASEQ_ect_sp6_KO_rep2
#   rnasq_sp6_wt_6 -> RNASEQ_ect_sp6_WT_rep3   rnasq_sp6_mt_6 -> RNASEQ_ect_sp6_KO_rep3
#   RNA_WT_32 -> RNASEQ_ect_DKOexp_WT_rep1     RNA_DKO_17 -> RNASEQ_ect_DKO_rep1
#   RNA_WT_53 -> RNASEQ_ect_DKOexp_WT_rep2     RNA_DKO_40 -> RNASEQ_ect_DKO_rep2
#   RNA_WT_56 -> RNASEQ_ect_DKOexp_WT_rep3     RNA_DKO_41 -> RNASEQ_ect_DKO_rep3
# =============================================================================

set -euo pipefail

mkdir -p output/trimmed output/align output/bigwig output/counts

# =============================================================================
# 1. Quality trimming (fastp) — paired-end, low-complexity filter (-y)
# =============================================================================

fastp -i input/SRR11015074_1.fastq.gz -I input/SRR11015074_2.fastq.gz -o output/trimmed/SRR11015074_1.fastq.gz -O output/trimmed/SRR11015074_2.fastq.gz --detect_adapter_for_pe -y -p -h output/trimmed/RNASEQ_ect_sp8_WT_rep1_fastp.html --thread 16
fastp -i input/SRR11015075_1.fastq.gz -I input/SRR11015075_2.fastq.gz -o output/trimmed/SRR11015075_1.fastq.gz -O output/trimmed/SRR11015075_2.fastq.gz --detect_adapter_for_pe -y -p -h output/trimmed/RNASEQ_ect_sp8_WT_rep2_fastp.html --thread 16
fastp -i input/SRR11015076_1.fastq.gz -I input/SRR11015076_2.fastq.gz -o output/trimmed/SRR11015076_1.fastq.gz -O output/trimmed/SRR11015076_2.fastq.gz --detect_adapter_for_pe -y -p -h output/trimmed/RNASEQ_ect_sp8_WT_rep3_fastp.html --thread 16
fastp -i input/SRR11015077_1.fastq.gz -I input/SRR11015077_2.fastq.gz -o output/trimmed/SRR11015077_1.fastq.gz -O output/trimmed/SRR11015077_2.fastq.gz --detect_adapter_for_pe -y -p -h output/trimmed/RNASEQ_ect_sp8_KO_rep1_fastp.html --thread 16
fastp -i input/SRR11015078_1.fastq.gz -I input/SRR11015078_2.fastq.gz -o output/trimmed/SRR11015078_1.fastq.gz -O output/trimmed/SRR11015078_2.fastq.gz --detect_adapter_for_pe -y -p -h output/trimmed/RNASEQ_ect_sp8_KO_rep2_fastp.html --thread 16
fastp -i input/SRR11015079_1.fastq.gz -I input/SRR11015079_2.fastq.gz -o output/trimmed/SRR11015079_1.fastq.gz -O output/trimmed/SRR11015079_2.fastq.gz --detect_adapter_for_pe -y -p -h output/trimmed/RNASEQ_ect_sp8_KO_rep3_fastp.html --thread 16

fastp -i input/rnasq_sp6_wt_4_1.fastq.gz -I input/rnasq_sp6_wt_4_2.fastq.gz -o output/trimmed/rnasq_sp6_wt_4_1.fastq.gz -O output/trimmed/rnasq_sp6_wt_4_2.fastq.gz --detect_adapter_for_pe -y -p -h output/trimmed/RNASEQ_ect_sp6_WT_rep1_fastp.html --thread 16
fastp -i input/rnasq_sp6_wt_5_1.fastq.gz -I input/rnasq_sp6_wt_5_2.fastq.gz -o output/trimmed/rnasq_sp6_wt_5_1.fastq.gz -O output/trimmed/rnasq_sp6_wt_5_2.fastq.gz --detect_adapter_for_pe -y -p -h output/trimmed/RNASEQ_ect_sp6_WT_rep2_fastp.html --thread 16
fastp -i input/rnasq_sp6_wt_6_1.fastq.gz -I input/rnasq_sp6_wt_6_2.fastq.gz -o output/trimmed/rnasq_sp6_wt_6_1.fastq.gz -O output/trimmed/rnasq_sp6_wt_6_2.fastq.gz --detect_adapter_for_pe -y -p -h output/trimmed/RNASEQ_ect_sp6_WT_rep3_fastp.html --thread 16
fastp -i input/rnasq_sp6_mt_4_1.fastq.gz -I input/rnasq_sp6_mt_4_2.fastq.gz -o output/trimmed/rnasq_sp6_mt_4_1.fastq.gz -O output/trimmed/rnasq_sp6_mt_4_2.fastq.gz --detect_adapter_for_pe -y -p -h output/trimmed/RNASEQ_ect_sp6_KO_rep1_fastp.html --thread 16
fastp -i input/rnasq_sp6_mt_5_1.fastq.gz -I input/rnasq_sp6_mt_5_2.fastq.gz -o output/trimmed/rnasq_sp6_mt_5_1.fastq.gz -O output/trimmed/rnasq_sp6_mt_5_2.fastq.gz --detect_adapter_for_pe -y -p -h output/trimmed/RNASEQ_ect_sp6_KO_rep2_fastp.html --thread 16
fastp -i input/rnasq_sp6_mt_6_1.fastq.gz -I input/rnasq_sp6_mt_6_2.fastq.gz -o output/trimmed/rnasq_sp6_mt_6_1.fastq.gz -O output/trimmed/rnasq_sp6_mt_6_2.fastq.gz --detect_adapter_for_pe -y -p -h output/trimmed/RNASEQ_ect_sp6_KO_rep3_fastp.html --thread 16

fastp -i input/RNA_WT_32_1.fastq.gz -I input/RNA_WT_32_2.fastq.gz -o output/trimmed/RNA_WT_32_1.fastq.gz -O output/trimmed/RNA_WT_32_2.fastq.gz --detect_adapter_for_pe -y -p -h output/trimmed/RNASEQ_ect_DKOexp_WT_rep1_fastp.html --thread 16
fastp -i input/RNA_WT_53_1.fastq.gz -I input/RNA_WT_53_2.fastq.gz -o output/trimmed/RNA_WT_53_1.fastq.gz -O output/trimmed/RNA_WT_53_2.fastq.gz --detect_adapter_for_pe -y -p -h output/trimmed/RNASEQ_ect_DKOexp_WT_rep2_fastp.html --thread 16
fastp -i input/RNA_WT_56_1.fastq.gz -I input/RNA_WT_56_2.fastq.gz -o output/trimmed/RNA_WT_56_1.fastq.gz -O output/trimmed/RNA_WT_56_2.fastq.gz --detect_adapter_for_pe -y -p -h output/trimmed/RNASEQ_ect_DKOexp_WT_rep3_fastp.html --thread 16
fastp -i input/RNA_DKO_17_1.fastq.gz -I input/RNA_DKO_17_2.fastq.gz -o output/trimmed/RNA_DKO_17_1.fastq.gz -O output/trimmed/RNA_DKO_17_2.fastq.gz --detect_adapter_for_pe -y -p -h output/trimmed/RNASEQ_ect_DKO_rep1_fastp.html --thread 16
fastp -i input/RNA_DKO_40_1.fastq.gz -I input/RNA_DKO_40_2.fastq.gz -o output/trimmed/RNA_DKO_40_1.fastq.gz -O output/trimmed/RNA_DKO_40_2.fastq.gz --detect_adapter_for_pe -y -p -h output/trimmed/RNASEQ_ect_DKO_rep2_fastp.html --thread 16
fastp -i input/RNA_DKO_41_1.fastq.gz -I input/RNA_DKO_41_2.fastq.gz -o output/trimmed/RNA_DKO_41_1.fastq.gz -O output/trimmed/RNA_DKO_41_2.fastq.gz --detect_adapter_for_pe -y -p -h output/trimmed/RNASEQ_ect_DKO_rep3_fastp.html --thread 16

# =============================================================================
# 2. Alignment (HISAT2, --dta-cufflinks) + drop unmapped + sort + index
# =============================================================================

hisat2 -t -p 16 --dta-cufflinks -x input/GRCm38_mouse -1 output/trimmed/SRR11015074_1.fastq.gz -2 output/trimmed/SRR11015074_2.fastq.gz --summary-file output/align/RNASEQ_ect_sp8_WT_rep1_hisat2.txt | samtools view -@ 16 -ShuF 4 - | samtools sort -@ 16 -o output/align/RNASEQ_ect_sp8_WT_rep1.bam - && samtools index output/align/RNASEQ_ect_sp8_WT_rep1.bam -@ 16
hisat2 -t -p 16 --dta-cufflinks -x input/GRCm38_mouse -1 output/trimmed/SRR11015075_1.fastq.gz -2 output/trimmed/SRR11015075_2.fastq.gz --summary-file output/align/RNASEQ_ect_sp8_WT_rep2_hisat2.txt | samtools view -@ 16 -ShuF 4 - | samtools sort -@ 16 -o output/align/RNASEQ_ect_sp8_WT_rep2.bam - && samtools index output/align/RNASEQ_ect_sp8_WT_rep2.bam -@ 16
hisat2 -t -p 16 --dta-cufflinks -x input/GRCm38_mouse -1 output/trimmed/SRR11015076_1.fastq.gz -2 output/trimmed/SRR11015076_2.fastq.gz --summary-file output/align/RNASEQ_ect_sp8_WT_rep3_hisat2.txt | samtools view -@ 16 -ShuF 4 - | samtools sort -@ 16 -o output/align/RNASEQ_ect_sp8_WT_rep3.bam - && samtools index output/align/RNASEQ_ect_sp8_WT_rep3.bam -@ 16
hisat2 -t -p 16 --dta-cufflinks -x input/GRCm38_mouse -1 output/trimmed/SRR11015077_1.fastq.gz -2 output/trimmed/SRR11015077_2.fastq.gz --summary-file output/align/RNASEQ_ect_sp8_KO_rep1_hisat2.txt | samtools view -@ 16 -ShuF 4 - | samtools sort -@ 16 -o output/align/RNASEQ_ect_sp8_KO_rep1.bam - && samtools index output/align/RNASEQ_ect_sp8_KO_rep1.bam -@ 16
hisat2 -t -p 16 --dta-cufflinks -x input/GRCm38_mouse -1 output/trimmed/SRR11015078_1.fastq.gz -2 output/trimmed/SRR11015078_2.fastq.gz --summary-file output/align/RNASEQ_ect_sp8_KO_rep2_hisat2.txt | samtools view -@ 16 -ShuF 4 - | samtools sort -@ 16 -o output/align/RNASEQ_ect_sp8_KO_rep2.bam - && samtools index output/align/RNASEQ_ect_sp8_KO_rep2.bam -@ 16
hisat2 -t -p 16 --dta-cufflinks -x input/GRCm38_mouse -1 output/trimmed/SRR11015079_1.fastq.gz -2 output/trimmed/SRR11015079_2.fastq.gz --summary-file output/align/RNASEQ_ect_sp8_KO_rep3_hisat2.txt | samtools view -@ 16 -ShuF 4 - | samtools sort -@ 16 -o output/align/RNASEQ_ect_sp8_KO_rep3.bam - && samtools index output/align/RNASEQ_ect_sp8_KO_rep3.bam -@ 16

hisat2 -t -p 16 --dta-cufflinks -x input/GRCm38_mouse -1 output/trimmed/rnasq_sp6_wt_4_1.fastq.gz -2 output/trimmed/rnasq_sp6_wt_4_2.fastq.gz --summary-file output/align/RNASEQ_ect_sp6_WT_rep1_hisat2.txt | samtools view -@ 16 -ShuF 4 - | samtools sort -@ 16 -o output/align/RNASEQ_ect_sp6_WT_rep1.bam - && samtools index output/align/RNASEQ_ect_sp6_WT_rep1.bam -@ 16
hisat2 -t -p 16 --dta-cufflinks -x input/GRCm38_mouse -1 output/trimmed/rnasq_sp6_wt_5_1.fastq.gz -2 output/trimmed/rnasq_sp6_wt_5_2.fastq.gz --summary-file output/align/RNASEQ_ect_sp6_WT_rep2_hisat2.txt | samtools view -@ 16 -ShuF 4 - | samtools sort -@ 16 -o output/align/RNASEQ_ect_sp6_WT_rep2.bam - && samtools index output/align/RNASEQ_ect_sp6_WT_rep2.bam -@ 16
hisat2 -t -p 16 --dta-cufflinks -x input/GRCm38_mouse -1 output/trimmed/rnasq_sp6_wt_6_1.fastq.gz -2 output/trimmed/rnasq_sp6_wt_6_2.fastq.gz --summary-file output/align/RNASEQ_ect_sp6_WT_rep3_hisat2.txt | samtools view -@ 16 -ShuF 4 - | samtools sort -@ 16 -o output/align/RNASEQ_ect_sp6_WT_rep3.bam - && samtools index output/align/RNASEQ_ect_sp6_WT_rep3.bam -@ 16
hisat2 -t -p 16 --dta-cufflinks -x input/GRCm38_mouse -1 output/trimmed/rnasq_sp6_mt_4_1.fastq.gz -2 output/trimmed/rnasq_sp6_mt_4_2.fastq.gz --summary-file output/align/RNASEQ_ect_sp6_KO_rep1_hisat2.txt | samtools view -@ 16 -ShuF 4 - | samtools sort -@ 16 -o output/align/RNASEQ_ect_sp6_KO_rep1.bam - && samtools index output/align/RNASEQ_ect_sp6_KO_rep1.bam -@ 16
hisat2 -t -p 16 --dta-cufflinks -x input/GRCm38_mouse -1 output/trimmed/rnasq_sp6_mt_5_1.fastq.gz -2 output/trimmed/rnasq_sp6_mt_5_2.fastq.gz --summary-file output/align/RNASEQ_ect_sp6_KO_rep2_hisat2.txt | samtools view -@ 16 -ShuF 4 - | samtools sort -@ 16 -o output/align/RNASEQ_ect_sp6_KO_rep2.bam - && samtools index output/align/RNASEQ_ect_sp6_KO_rep2.bam -@ 16
hisat2 -t -p 16 --dta-cufflinks -x input/GRCm38_mouse -1 output/trimmed/rnasq_sp6_mt_6_1.fastq.gz -2 output/trimmed/rnasq_sp6_mt_6_2.fastq.gz --summary-file output/align/RNASEQ_ect_sp6_KO_rep3_hisat2.txt | samtools view -@ 16 -ShuF 4 - | samtools sort -@ 16 -o output/align/RNASEQ_ect_sp6_KO_rep3.bam - && samtools index output/align/RNASEQ_ect_sp6_KO_rep3.bam -@ 16

hisat2 -t -p 16 --dta-cufflinks -x input/GRCm38_mouse -1 output/trimmed/RNA_WT_32_1.fastq.gz -2 output/trimmed/RNA_WT_32_2.fastq.gz --summary-file output/align/RNASEQ_ect_DKOexp_WT_rep1_hisat2.txt | samtools view -@ 16 -ShuF 4 - | samtools sort -@ 16 -o output/align/RNASEQ_ect_DKOexp_WT_rep1.bam - && samtools index output/align/RNASEQ_ect_DKOexp_WT_rep1.bam -@ 16
hisat2 -t -p 16 --dta-cufflinks -x input/GRCm38_mouse -1 output/trimmed/RNA_WT_53_1.fastq.gz -2 output/trimmed/RNA_WT_53_2.fastq.gz --summary-file output/align/RNASEQ_ect_DKOexp_WT_rep2_hisat2.txt | samtools view -@ 16 -ShuF 4 - | samtools sort -@ 16 -o output/align/RNASEQ_ect_DKOexp_WT_rep2.bam - && samtools index output/align/RNASEQ_ect_DKOexp_WT_rep2.bam -@ 16
hisat2 -t -p 16 --dta-cufflinks -x input/GRCm38_mouse -1 output/trimmed/RNA_WT_56_1.fastq.gz -2 output/trimmed/RNA_WT_56_2.fastq.gz --summary-file output/align/RNASEQ_ect_DKOexp_WT_rep3_hisat2.txt | samtools view -@ 16 -ShuF 4 - | samtools sort -@ 16 -o output/align/RNASEQ_ect_DKOexp_WT_rep3.bam - && samtools index output/align/RNASEQ_ect_DKOexp_WT_rep3.bam -@ 16
hisat2 -t -p 16 --dta-cufflinks -x input/GRCm38_mouse -1 output/trimmed/RNA_DKO_17_1.fastq.gz -2 output/trimmed/RNA_DKO_17_2.fastq.gz --summary-file output/align/RNASEQ_ect_DKO_rep1_hisat2.txt | samtools view -@ 16 -ShuF 4 - | samtools sort -@ 16 -o output/align/RNASEQ_ect_DKO_rep1.bam - && samtools index output/align/RNASEQ_ect_DKO_rep1.bam -@ 16
hisat2 -t -p 16 --dta-cufflinks -x input/GRCm38_mouse -1 output/trimmed/RNA_DKO_40_1.fastq.gz -2 output/trimmed/RNA_DKO_40_2.fastq.gz --summary-file output/align/RNASEQ_ect_DKO_rep2_hisat2.txt | samtools view -@ 16 -ShuF 4 - | samtools sort -@ 16 -o output/align/RNASEQ_ect_DKO_rep2.bam - && samtools index output/align/RNASEQ_ect_DKO_rep2.bam -@ 16
hisat2 -t -p 16 --dta-cufflinks -x input/GRCm38_mouse -1 output/trimmed/RNA_DKO_41_1.fastq.gz -2 output/trimmed/RNA_DKO_41_2.fastq.gz --summary-file output/align/RNASEQ_ect_DKO_rep3_hisat2.txt | samtools view -@ 16 -ShuF 4 - | samtools sort -@ 16 -o output/align/RNASEQ_ect_DKO_rep3.bam - && samtools index output/align/RNASEQ_ect_DKO_rep3.bam -@ 16

# =============================================================================
# 3. BigWig generation (deepTools bamCoverage, RPGC normalisation)
# =============================================================================

bamCoverage -b output/align/RNASEQ_ect_sp8_WT_rep1.bam -o output/bigwig/RNASEQ_ect_sp8_WT_rep1.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 -p 16
bamCoverage -b output/align/RNASEQ_ect_sp8_WT_rep2.bam -o output/bigwig/RNASEQ_ect_sp8_WT_rep2.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 -p 16
bamCoverage -b output/align/RNASEQ_ect_sp8_WT_rep3.bam -o output/bigwig/RNASEQ_ect_sp8_WT_rep3.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 -p 16
bamCoverage -b output/align/RNASEQ_ect_sp8_KO_rep1.bam -o output/bigwig/RNASEQ_ect_sp8_KO_rep1.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 -p 16
bamCoverage -b output/align/RNASEQ_ect_sp8_KO_rep2.bam -o output/bigwig/RNASEQ_ect_sp8_KO_rep2.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 -p 16
bamCoverage -b output/align/RNASEQ_ect_sp8_KO_rep3.bam -o output/bigwig/RNASEQ_ect_sp8_KO_rep3.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 -p 16

bamCoverage -b output/align/RNASEQ_ect_sp6_WT_rep1.bam -o output/bigwig/RNASEQ_ect_sp6_WT_rep1.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 -p 16
bamCoverage -b output/align/RNASEQ_ect_sp6_WT_rep2.bam -o output/bigwig/RNASEQ_ect_sp6_WT_rep2.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 -p 16
bamCoverage -b output/align/RNASEQ_ect_sp6_WT_rep3.bam -o output/bigwig/RNASEQ_ect_sp6_WT_rep3.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 -p 16
bamCoverage -b output/align/RNASEQ_ect_sp6_KO_rep1.bam -o output/bigwig/RNASEQ_ect_sp6_KO_rep1.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 -p 16
bamCoverage -b output/align/RNASEQ_ect_sp6_KO_rep2.bam -o output/bigwig/RNASEQ_ect_sp6_KO_rep2.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 -p 16
bamCoverage -b output/align/RNASEQ_ect_sp6_KO_rep3.bam -o output/bigwig/RNASEQ_ect_sp6_KO_rep3.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 -p 16

bamCoverage -b output/align/RNASEQ_ect_DKOexp_WT_rep1.bam -o output/bigwig/RNASEQ_ect_DKOexp_WT_rep1.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 -p 16
bamCoverage -b output/align/RNASEQ_ect_DKOexp_WT_rep2.bam -o output/bigwig/RNASEQ_ect_DKOexp_WT_rep2.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 -p 16
bamCoverage -b output/align/RNASEQ_ect_DKOexp_WT_rep3.bam -o output/bigwig/RNASEQ_ect_DKOexp_WT_rep3.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 -p 16
bamCoverage -b output/align/RNASEQ_ect_DKO_rep1.bam -o output/bigwig/RNASEQ_ect_DKO_rep1.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 -p 16
bamCoverage -b output/align/RNASEQ_ect_DKO_rep2.bam -o output/bigwig/RNASEQ_ect_DKO_rep2.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 -p 16
bamCoverage -b output/align/RNASEQ_ect_DKO_rep3.bam -o output/bigwig/RNASEQ_ect_DKO_rep3.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 -p 16

# =============================================================================
# 4. Transcript quantification (StringTie, expression-only mode -e)
# =============================================================================

stringtie -o output/counts/RNASEQ_ect_sp8_WT_rep1.gtf -e -A output/counts/RNASEQ_ect_sp8_WT_rep1.tab -v -p 12 -G input/gencode.vM25.basic.annotation.gtf output/align/RNASEQ_ect_sp8_WT_rep1.bam
stringtie -o output/counts/RNASEQ_ect_sp8_WT_rep2.gtf -e -A output/counts/RNASEQ_ect_sp8_WT_rep2.tab -v -p 12 -G input/gencode.vM25.basic.annotation.gtf output/align/RNASEQ_ect_sp8_WT_rep2.bam
stringtie -o output/counts/RNASEQ_ect_sp8_WT_rep3.gtf -e -A output/counts/RNASEQ_ect_sp8_WT_rep3.tab -v -p 12 -G input/gencode.vM25.basic.annotation.gtf output/align/RNASEQ_ect_sp8_WT_rep3.bam
stringtie -o output/counts/RNASEQ_ect_sp8_KO_rep1.gtf -e -A output/counts/RNASEQ_ect_sp8_KO_rep1.tab -v -p 12 -G input/gencode.vM25.basic.annotation.gtf output/align/RNASEQ_ect_sp8_KO_rep1.bam
stringtie -o output/counts/RNASEQ_ect_sp8_KO_rep2.gtf -e -A output/counts/RNASEQ_ect_sp8_KO_rep2.tab -v -p 12 -G input/gencode.vM25.basic.annotation.gtf output/align/RNASEQ_ect_sp8_KO_rep2.bam
stringtie -o output/counts/RNASEQ_ect_sp8_KO_rep3.gtf -e -A output/counts/RNASEQ_ect_sp8_KO_rep3.tab -v -p 12 -G input/gencode.vM25.basic.annotation.gtf output/align/RNASEQ_ect_sp8_KO_rep3.bam

stringtie -o output/counts/RNASEQ_ect_sp6_WT_rep1.gtf -e -A output/counts/RNASEQ_ect_sp6_WT_rep1.tab -v -p 12 -G input/gencode.vM25.basic.annotation.gtf output/align/RNASEQ_ect_sp6_WT_rep1.bam
stringtie -o output/counts/RNASEQ_ect_sp6_WT_rep2.gtf -e -A output/counts/RNASEQ_ect_sp6_WT_rep2.tab -v -p 12 -G input/gencode.vM25.basic.annotation.gtf output/align/RNASEQ_ect_sp6_WT_rep2.bam
stringtie -o output/counts/RNASEQ_ect_sp6_WT_rep3.gtf -e -A output/counts/RNASEQ_ect_sp6_WT_rep3.tab -v -p 12 -G input/gencode.vM25.basic.annotation.gtf output/align/RNASEQ_ect_sp6_WT_rep3.bam
stringtie -o output/counts/RNASEQ_ect_sp6_KO_rep1.gtf -e -A output/counts/RNASEQ_ect_sp6_KO_rep1.tab -v -p 12 -G input/gencode.vM25.basic.annotation.gtf output/align/RNASEQ_ect_sp6_KO_rep1.bam
stringtie -o output/counts/RNASEQ_ect_sp6_KO_rep2.gtf -e -A output/counts/RNASEQ_ect_sp6_KO_rep2.tab -v -p 12 -G input/gencode.vM25.basic.annotation.gtf output/align/RNASEQ_ect_sp6_KO_rep2.bam
stringtie -o output/counts/RNASEQ_ect_sp6_KO_rep3.gtf -e -A output/counts/RNASEQ_ect_sp6_KO_rep3.tab -v -p 12 -G input/gencode.vM25.basic.annotation.gtf output/align/RNASEQ_ect_sp6_KO_rep3.bam

stringtie -o output/counts/RNASEQ_ect_DKOexp_WT_rep1.gtf -e -A output/counts/RNASEQ_ect_DKOexp_WT_rep1.tab -v -p 12 -G input/gencode.vM25.basic.annotation.gtf output/align/RNASEQ_ect_DKOexp_WT_rep1.bam
stringtie -o output/counts/RNASEQ_ect_DKOexp_WT_rep2.gtf -e -A output/counts/RNASEQ_ect_DKOexp_WT_rep2.tab -v -p 12 -G input/gencode.vM25.basic.annotation.gtf output/align/RNASEQ_ect_DKOexp_WT_rep2.bam
stringtie -o output/counts/RNASEQ_ect_DKOexp_WT_rep3.gtf -e -A output/counts/RNASEQ_ect_DKOexp_WT_rep3.tab -v -p 12 -G input/gencode.vM25.basic.annotation.gtf output/align/RNASEQ_ect_DKOexp_WT_rep3.bam
stringtie -o output/counts/RNASEQ_ect_DKO_rep1.gtf -e -A output/counts/RNASEQ_ect_DKO_rep1.tab -v -p 12 -G input/gencode.vM25.basic.annotation.gtf output/align/RNASEQ_ect_DKO_rep1.bam
stringtie -o output/counts/RNASEQ_ect_DKO_rep2.gtf -e -A output/counts/RNASEQ_ect_DKO_rep2.tab -v -p 12 -G input/gencode.vM25.basic.annotation.gtf output/align/RNASEQ_ect_DKO_rep2.bam
stringtie -o output/counts/RNASEQ_ect_DKO_rep3.gtf -e -A output/counts/RNASEQ_ect_DKO_rep3.tab -v -p 12 -G input/gencode.vM25.basic.annotation.gtf output/align/RNASEQ_ect_DKO_rep3.bam

# =============================================================================
# 5. Count matrix generation (prepDE.py3) — one sample list + one call per experiment
# =============================================================================

cat > output/counts/sample_list_sp8.txt <<'EOF'
RNASEQ_ect_sp8_WT_rep1 output/counts/RNASEQ_ect_sp8_WT_rep1.gtf
RNASEQ_ect_sp8_WT_rep2 output/counts/RNASEQ_ect_sp8_WT_rep2.gtf
RNASEQ_ect_sp8_WT_rep3 output/counts/RNASEQ_ect_sp8_WT_rep3.gtf
RNASEQ_ect_sp8_KO_rep1 output/counts/RNASEQ_ect_sp8_KO_rep1.gtf
RNASEQ_ect_sp8_KO_rep2 output/counts/RNASEQ_ect_sp8_KO_rep2.gtf
RNASEQ_ect_sp8_KO_rep3 output/counts/RNASEQ_ect_sp8_KO_rep3.gtf
EOF
python3 prepDE.py3 -i output/counts/sample_list_sp8.txt -g output/counts/gene_count_matrix_sp8.csv -t output/counts/transcript_count_matrix_sp8.csv

cat > output/counts/sample_list_sp6.txt <<'EOF'
RNASEQ_ect_sp6_WT_rep1 output/counts/RNASEQ_ect_sp6_WT_rep1.gtf
RNASEQ_ect_sp6_WT_rep2 output/counts/RNASEQ_ect_sp6_WT_rep2.gtf
RNASEQ_ect_sp6_WT_rep3 output/counts/RNASEQ_ect_sp6_WT_rep3.gtf
RNASEQ_ect_sp6_KO_rep1 output/counts/RNASEQ_ect_sp6_KO_rep1.gtf
RNASEQ_ect_sp6_KO_rep2 output/counts/RNASEQ_ect_sp6_KO_rep2.gtf
RNASEQ_ect_sp6_KO_rep3 output/counts/RNASEQ_ect_sp6_KO_rep3.gtf
EOF
python3 prepDE.py3 -i output/counts/sample_list_sp6.txt -g output/counts/gene_count_matrix_sp6.csv -t output/counts/transcript_count_matrix_sp6.csv

cat > output/counts/sample_list_dko.txt <<'EOF'
RNASEQ_ect_DKOexp_WT_rep1 output/counts/RNASEQ_ect_DKOexp_WT_rep1.gtf
RNASEQ_ect_DKOexp_WT_rep2 output/counts/RNASEQ_ect_DKOexp_WT_rep2.gtf
RNASEQ_ect_DKOexp_WT_rep3 output/counts/RNASEQ_ect_DKOexp_WT_rep3.gtf
RNASEQ_ect_DKO_rep1 output/counts/RNASEQ_ect_DKO_rep1.gtf
RNASEQ_ect_DKO_rep2 output/counts/RNASEQ_ect_DKO_rep2.gtf
RNASEQ_ect_DKO_rep3 output/counts/RNASEQ_ect_DKO_rep3.gtf
EOF
python3 prepDE.py3 -i output/counts/sample_list_dko.txt -g output/counts/gene_count_matrix_DKO.csv -t output/counts/transcript_count_matrix_DKO.csv

