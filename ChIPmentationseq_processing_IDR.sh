#!/bin/bash
# =============================================================================
# ChIPmentation processing + IDR peak calling for Sp6 (PE) and Sp8 (SE), mm10.
#
# Expected in input/:
#   chipd_sp6_rep1_1/2.fastq.gz, chipd_sp6_rep2_1/2.fastq.gz,
#   input_sp6_rep1_1/2.fastq.gz, input_sp6_rep2_1/2.fastq.gz   (Sp6, PE)
#   rochp_sp8_rep1/2.fastq.gz, roinp_sp8_rep1/2.fastq.gz       (Sp8, SE)
#   GRCm38.primary_assembly.genome.fa (+ bwa index alongside it)
#   mm10-blacklist.v2.bed, tss.sorted.bed, atac_ectoderm.bw
#
# Tools: fastp, bwa, samtools, sambamba, bedtools, deepTools, MACS2 (env
# macs2_python_37), IDR (env IDR_python3). Edit env names below if different.
# =============================================================================

set -euo pipefail

mkdir -p output/trimmed output/align output/bigwig output/heatmaps output/matrix output/annot \
         output/idr/sp6/CONS_PEAKS output/idr/sp6/OPT_PEAKS output/idr/sp6/TRACKS \
         output/idr/sp8/CONS_PEAKS output/idr/sp8/OPT_PEAKS output/idr/sp8/TRACKS

# =============================================================================
# 1. Quality trimming (fastp)
# =============================================================================

fastp -i input/chipd_sp6_rep1_1.fastq.gz -I input/chipd_sp6_rep1_2.fastq.gz -o output/trimmed/chipd_sp6_rep1_1.fastq.gz -O output/trimmed/chipd_sp6_rep1_2.fastq.gz --detect_adapter_for_pe -p -h output/trimmed/chipd_sp6_rep1_fastp.html --thread 16
fastp -i input/chipd_sp6_rep2_1.fastq.gz -I input/chipd_sp6_rep2_2.fastq.gz -o output/trimmed/chipd_sp6_rep2_1.fastq.gz -O output/trimmed/chipd_sp6_rep2_2.fastq.gz --detect_adapter_for_pe -p -h output/trimmed/chipd_sp6_rep2_fastp.html --thread 16
fastp -i input/input_sp6_rep1_1.fastq.gz -I input/input_sp6_rep1_2.fastq.gz -o output/trimmed/input_sp6_rep1_1.fastq.gz -O output/trimmed/input_sp6_rep1_2.fastq.gz --detect_adapter_for_pe -p -h output/trimmed/input_sp6_rep1_fastp.html --thread 16
fastp -i input/input_sp6_rep2_1.fastq.gz -I input/input_sp6_rep2_2.fastq.gz -o output/trimmed/input_sp6_rep2_1.fastq.gz -O output/trimmed/input_sp6_rep2_2.fastq.gz --detect_adapter_for_pe -p -h output/trimmed/input_sp6_rep2_fastp.html --thread 16

fastp -i input/rochp_sp8_rep1.fastq.gz -o output/trimmed/rochp_sp8_rep1.fastq.gz -p -h output/trimmed/rochp_sp8_rep1_fastp.html --thread 16
fastp -i input/rochp_sp8_rep2.fastq.gz -o output/trimmed/rochp_sp8_rep2.fastq.gz -p -h output/trimmed/rochp_sp8_rep2_fastp.html --thread 16
fastp -i input/roinp_sp8_rep1.fastq.gz -o output/trimmed/roinp_sp8_rep1.fastq.gz -p -h output/trimmed/roinp_sp8_rep1_fastp.html --thread 16
fastp -i input/roinp_sp8_rep2.fastq.gz -o output/trimmed/roinp_sp8_rep2.fastq.gz -p -h output/trimmed/roinp_sp8_rep2_fastp.html --thread 16

# =============================================================================
# 2. Alignment (bwa) + coordinate sort (samtools)
#    Sp6 paired-end -> bwa mem. Sp8 single-end -> bwa aln | bwa samse.
# =============================================================================

bwa mem -t 16 input/GRCm38.primary_assembly.genome.fa output/trimmed/chipd_sp6_rep1_1.fastq.gz output/trimmed/chipd_sp6_rep1_2.fastq.gz | samtools sort -@ 16 -o output/align/chipd_sp6_rep1_sorted.bam -
bwa mem -t 16 input/GRCm38.primary_assembly.genome.fa output/trimmed/chipd_sp6_rep2_1.fastq.gz output/trimmed/chipd_sp6_rep2_2.fastq.gz | samtools sort -@ 16 -o output/align/chipd_sp6_rep2_sorted.bam -
bwa mem -t 16 input/GRCm38.primary_assembly.genome.fa output/trimmed/input_sp6_rep1_1.fastq.gz output/trimmed/input_sp6_rep1_2.fastq.gz | samtools sort -@ 16 -o output/align/input_sp6_rep1_sorted.bam -
bwa mem -t 16 input/GRCm38.primary_assembly.genome.fa output/trimmed/input_sp6_rep2_1.fastq.gz output/trimmed/input_sp6_rep2_2.fastq.gz | samtools sort -@ 16 -o output/align/input_sp6_rep2_sorted.bam -

bwa aln -t 16 input/GRCm38.primary_assembly.genome.fa output/trimmed/rochp_sp8_rep1.fastq.gz | bwa samse input/GRCm38.primary_assembly.genome.fa - output/trimmed/rochp_sp8_rep1.fastq.gz | samtools sort -@ 16 -o output/align/rochp_sp8_rep1_sorted.bam -
bwa aln -t 16 input/GRCm38.primary_assembly.genome.fa output/trimmed/rochp_sp8_rep2.fastq.gz | bwa samse input/GRCm38.primary_assembly.genome.fa - output/trimmed/rochp_sp8_rep2.fastq.gz | samtools sort -@ 16 -o output/align/rochp_sp8_rep2_sorted.bam -
bwa aln -t 16 input/GRCm38.primary_assembly.genome.fa output/trimmed/roinp_sp8_rep1.fastq.gz | bwa samse input/GRCm38.primary_assembly.genome.fa - output/trimmed/roinp_sp8_rep1.fastq.gz | samtools sort -@ 16 -o output/align/roinp_sp8_rep1_sorted.bam -
bwa aln -t 16 input/GRCm38.primary_assembly.genome.fa output/trimmed/roinp_sp8_rep2.fastq.gz | bwa samse input/GRCm38.primary_assembly.genome.fa - output/trimmed/roinp_sp8_rep2.fastq.gz | samtools sort -@ 16 -o output/align/roinp_sp8_rep2_sorted.bam -

# =============================================================================
# 3. Filter to uniquely-mapping, non-duplicate reads (sambamba)
# =============================================================================

sambamba view -h -t 2 -f bam -F "[XA] == null and not unmapped and not duplicate" output/align/chipd_sp6_rep1_sorted.bam > output/align/chipd_sp6_rep1_filt.bam
sambamba view -h -t 2 -f bam -F "[XA] == null and not unmapped and not duplicate" output/align/chipd_sp6_rep2_sorted.bam > output/align/chipd_sp6_rep2_filt.bam
sambamba view -h -t 2 -f bam -F "[XA] == null and not unmapped and not duplicate" output/align/input_sp6_rep1_sorted.bam > output/align/input_sp6_rep1_filt.bam
sambamba view -h -t 2 -f bam -F "[XA] == null and not unmapped and not duplicate" output/align/input_sp6_rep2_sorted.bam > output/align/input_sp6_rep2_filt.bam
sambamba view -h -t 2 -f bam -F "[XA] == null and not unmapped and not duplicate" output/align/rochp_sp8_rep1_sorted.bam > output/align/rochp_sp8_rep1_filt.bam
sambamba view -h -t 2 -f bam -F "[XA] == null and not unmapped and not duplicate" output/align/rochp_sp8_rep2_sorted.bam > output/align/rochp_sp8_rep2_filt.bam
sambamba view -h -t 2 -f bam -F "[XA] == null and not unmapped and not duplicate" output/align/roinp_sp8_rep1_sorted.bam > output/align/roinp_sp8_rep1_filt.bam
sambamba view -h -t 2 -f bam -F "[XA] == null and not unmapped and not duplicate" output/align/roinp_sp8_rep2_sorted.bam > output/align/roinp_sp8_rep2_filt.bam

# =============================================================================
# 4. Restrict to canonical chromosomes and index (samtools)
# =============================================================================

samtools view -b output/align/chipd_sp6_rep1_filt.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY -@ 16 > output/align/chipd_sp6_rep1_filt_chr.bam && samtools index output/align/chipd_sp6_rep1_filt_chr.bam -@ 16
samtools view -b output/align/chipd_sp6_rep2_filt.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY -@ 16 > output/align/chipd_sp6_rep2_filt_chr.bam && samtools index output/align/chipd_sp6_rep2_filt_chr.bam -@ 16
samtools view -b output/align/input_sp6_rep1_filt.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY -@ 16 > output/align/input_sp6_rep1_filt_chr.bam && samtools index output/align/input_sp6_rep1_filt_chr.bam -@ 16
samtools view -b output/align/input_sp6_rep2_filt.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY -@ 16 > output/align/input_sp6_rep2_filt_chr.bam && samtools index output/align/input_sp6_rep2_filt_chr.bam -@ 16
samtools view -b output/align/rochp_sp8_rep1_filt.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY -@ 16 > output/align/rochp_sp8_rep1_filt_chr.bam && samtools index output/align/rochp_sp8_rep1_filt_chr.bam -@ 16
samtools view -b output/align/rochp_sp8_rep2_filt.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY -@ 16 > output/align/rochp_sp8_rep2_filt_chr.bam && samtools index output/align/rochp_sp8_rep2_filt_chr.bam -@ 16
samtools view -b output/align/roinp_sp8_rep1_filt.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY -@ 16 > output/align/roinp_sp8_rep1_filt_chr.bam && samtools index output/align/roinp_sp8_rep1_filt_chr.bam -@ 16
samtools view -b output/align/roinp_sp8_rep2_filt.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY -@ 16 > output/align/roinp_sp8_rep2_filt_chr.bam && samtools index output/align/roinp_sp8_rep2_filt_chr.bam -@ 16

# =============================================================================
# 5. BigWig generation, RPGC-normalised (deepTools bamCoverage)
# =============================================================================

bamCoverage -b output/align/chipd_sp6_rep1_filt_chr.bam -o output/bigwig/chipd_sp6_rep1.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 --extendReads 150 --centerReads -p 16
bamCoverage -b output/align/chipd_sp6_rep2_filt_chr.bam -o output/bigwig/chipd_sp6_rep2.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 --extendReads 150 --centerReads -p 16
bamCoverage -b output/align/input_sp6_rep1_filt_chr.bam -o output/bigwig/input_sp6_rep1.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 --extendReads 150 --centerReads -p 16
bamCoverage -b output/align/input_sp6_rep2_filt_chr.bam -o output/bigwig/input_sp6_rep2.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 --extendReads 150 --centerReads -p 16
bamCoverage -b output/align/rochp_sp8_rep1_filt_chr.bam -o output/bigwig/rochp_sp8_rep1.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 --extendReads 150 --centerReads -p 16
bamCoverage -b output/align/rochp_sp8_rep2_filt_chr.bam -o output/bigwig/rochp_sp8_rep2.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 --extendReads 150 --centerReads -p 16
bamCoverage -b output/align/roinp_sp8_rep1_filt_chr.bam -o output/bigwig/roinp_sp8_rep1.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 --extendReads 150 --centerReads -p 16
bamCoverage -b output/align/roinp_sp8_rep2_filt_chr.bam -o output/bigwig/roinp_sp8_rep2.bw --binSize 20 --effectiveGenomeSize 1870000000 --normalizeUsing RPGC --smoothLength 60 --extendReads 150 --centerReads -p 16

# =============================================================================
# 6. Replicate QC: Pearson / Spearman correlation (deepTools)
# =============================================================================

multiBamSummary bins --bamfiles output/align/rochp_sp8_rep1_filt_chr.bam output/align/rochp_sp8_rep2_filt_chr.bam output/align/roinp_sp8_rep1_filt_chr.bam output/align/roinp_sp8_rep2_filt_chr.bam output/align/chipd_sp6_rep1_filt_chr.bam output/align/chipd_sp6_rep2_filt_chr.bam output/align/input_sp6_rep1_filt_chr.bam output/align/input_sp6_rep2_filt_chr.bam --labels Sp8_rep1 Sp8_rep2 Sp8_inp_rep1 Sp8_inp_rep2 Sp6_rep1 Sp6_rep2 Sp6_inp_rep1 Sp6_inp_rep2 -o output/align/sp6_sp8_multibamsummary.npz -p 16 -bl input/mm10-blacklist.v2.bed

plotCorrelation --corData output/align/sp6_sp8_multibamsummary.npz --corMethod pearson --removeOutliers --whatToPlot heatmap --plotNumbers -o output/heatmaps/corr_pearson_heatmap.pdf
plotCorrelation --corData output/align/sp6_sp8_multibamsummary.npz --corMethod spearman --removeOutliers --whatToPlot heatmap --plotNumbers -o output/heatmaps/corr_spearman_heatmap.pdf

# =============================================================================
# 7. BAM -> BEDPE (samtools/bedtools), for MACS2 BEDPE-mode peak calling
#    Sp6 (PE): name-sort -> fixmate -> coord-sort -> bamtobed -bedpe
#    Sp8 (SE): degenerate BEDPE (read spans itself)
# =============================================================================

samtools sort -n -@ 16 output/align/chipd_sp6_rep1_filt_chr.bam | samtools fixmate -m - - | samtools sort -@ 16 - | bedtools bamtobed -i stdin -bedpe | sort -k1,1 -k2,2n > output/idr/sp6/chipd_sp6_rep1.bedpe
samtools sort -n -@ 16 output/align/chipd_sp6_rep2_filt_chr.bam | samtools fixmate -m - - | samtools sort -@ 16 - | bedtools bamtobed -i stdin -bedpe | sort -k1,1 -k2,2n > output/idr/sp6/chipd_sp6_rep2.bedpe
samtools sort -n -@ 16 output/align/input_sp6_rep1_filt_chr.bam | samtools fixmate -m - - | samtools sort -@ 16 - | bedtools bamtobed -i stdin -bedpe | sort -k1,1 -k2,2n > output/idr/sp6/input_sp6_rep1.bedpe
samtools sort -n -@ 16 output/align/input_sp6_rep2_filt_chr.bam | samtools fixmate -m - - | samtools sort -@ 16 - | bedtools bamtobed -i stdin -bedpe | sort -k1,1 -k2,2n > output/idr/sp6/input_sp6_rep2.bedpe

bedtools bamtobed -i output/align/rochp_sp8_rep1_filt_chr.bam | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$1,$2,$3,$4,$5,$6}' | sort -k1,1 -k2,2n > output/idr/sp8/rochp_sp8_rep1.bedpe
bedtools bamtobed -i output/align/rochp_sp8_rep2_filt_chr.bam | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$1,$2,$3,$4,$5,$6}' | sort -k1,1 -k2,2n > output/idr/sp8/rochp_sp8_rep2.bedpe
bedtools bamtobed -i output/align/roinp_sp8_rep1_filt_chr.bam | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$1,$2,$3,$4,$5,$6}' | sort -k1,1 -k2,2n > output/idr/sp8/roinp_sp8_rep1.bedpe
bedtools bamtobed -i output/align/roinp_sp8_rep2_filt_chr.bam | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$1,$2,$3,$4,$5,$6}' | sort -k1,1 -k2,2n > output/idr/sp8/roinp_sp8_rep2.bedpe

# =============================================================================
# 8. Peak calling (MACS2, p=0.02, BEDPE mode) + IDR (soft threshold 0.1)
#    ENCODE IDR framework: conservative (true replicates vs pool) and
#    optimal (per-pseudoreplicate) peak sets, filtered at global IDR >= 1.
# =============================================================================

# ---- Sp6 ----

cat output/idr/sp6/input_sp6_rep1.bedpe output/idr/sp6/input_sp6_rep2.bedpe > output/idr/sp6/input_pooled.bedpe
cat output/idr/sp6/chipd_sp6_rep1.bedpe output/idr/sp6/chipd_sp6_rep2.bedpe > output/idr/sp6/chipd_sp6_pooled.bedpe

micromamba activate macs2_python_37

macs2 callpeak -f BEDPE -t output/idr/sp6/chipd_sp6_rep1.bedpe -c output/idr/sp6/input_sp6_rep1.bedpe -p 0.02 --keep-dup auto -g mm -n sp6_rep1 --outdir output/idr/sp6/CONS_PEAKS
macs2 callpeak -f BEDPE -t output/idr/sp6/chipd_sp6_rep2.bedpe -c output/idr/sp6/input_sp6_rep2.bedpe -p 0.02 --keep-dup auto -g mm -n sp6_rep2 --outdir output/idr/sp6/CONS_PEAKS
macs2 callpeak -f BEDPE -t output/idr/sp6/chipd_sp6_rep1.bedpe output/idr/sp6/chipd_sp6_rep2.bedpe -c output/idr/sp6/input_pooled.bedpe -p 0.02 --keep-dup auto -g mm -n sp6_pooled --outdir output/idr/sp6/CONS_PEAKS

sort -k8,8nr output/idr/sp6/CONS_PEAKS/sp6_rep1_peaks.narrowPeak | head -n 500000 > output/idr/sp6/CONS_PEAKS/sp6_rep1_top.bed
sort -k8,8nr output/idr/sp6/CONS_PEAKS/sp6_rep2_peaks.narrowPeak | head -n 500000 > output/idr/sp6/CONS_PEAKS/sp6_rep2_top.bed
sort -k8,8nr output/idr/sp6/CONS_PEAKS/sp6_pooled_peaks.narrowPeak | head -n 500000 > output/idr/sp6/CONS_PEAKS/sp6_pooled_top.bed

shuf output/idr/sp6/chipd_sp6_rep1.bedpe | split -d -l $(( ($(wc -l < output/idr/sp6/chipd_sp6_rep1.bedpe) + 1) / 2 )) - output/idr/sp6/OPT_PEAKS/sp6_rep1_ps
shuf output/idr/sp6/chipd_sp6_rep2.bedpe | split -d -l $(( ($(wc -l < output/idr/sp6/chipd_sp6_rep2.bedpe) + 1) / 2 )) - output/idr/sp6/OPT_PEAKS/sp6_rep2_ps
shuf output/idr/sp6/chipd_sp6_pooled.bedpe | split -d -l $(( ($(wc -l < output/idr/sp6/chipd_sp6_pooled.bedpe) + 1) / 2 )) - output/idr/sp6/OPT_PEAKS/sp6_pool_ps

macs2 callpeak -f BEDPE -t output/idr/sp6/OPT_PEAKS/sp6_rep1_ps00 -c output/idr/sp6/input_sp6_rep1.bedpe -p 0.02 --keep-dup auto -g mm -n sp6_rep1_ps00 --outdir output/idr/sp6/OPT_PEAKS
macs2 callpeak -f BEDPE -t output/idr/sp6/OPT_PEAKS/sp6_rep1_ps01 -c output/idr/sp6/input_sp6_rep1.bedpe -p 0.02 --keep-dup auto -g mm -n sp6_rep1_ps01 --outdir output/idr/sp6/OPT_PEAKS
macs2 callpeak -f BEDPE -t output/idr/sp6/OPT_PEAKS/sp6_rep2_ps00 -c output/idr/sp6/input_sp6_rep2.bedpe -p 0.02 --keep-dup auto -g mm -n sp6_rep2_ps00 --outdir output/idr/sp6/OPT_PEAKS
macs2 callpeak -f BEDPE -t output/idr/sp6/OPT_PEAKS/sp6_rep2_ps01 -c output/idr/sp6/input_sp6_rep2.bedpe -p 0.02 --keep-dup auto -g mm -n sp6_rep2_ps01 --outdir output/idr/sp6/OPT_PEAKS
macs2 callpeak -f BEDPE -t output/idr/sp6/OPT_PEAKS/sp6_pool_ps00 -c output/idr/sp6/input_pooled.bedpe -p 0.02 --keep-dup auto -g mm -n sp6_pool_ps00 --outdir output/idr/sp6/OPT_PEAKS
macs2 callpeak -f BEDPE -t output/idr/sp6/OPT_PEAKS/sp6_pool_ps01 -c output/idr/sp6/input_pooled.bedpe -p 0.02 --keep-dup auto -g mm -n sp6_pool_ps01 --outdir output/idr/sp6/OPT_PEAKS

sort -k8,8nr output/idr/sp6/OPT_PEAKS/sp6_rep1_ps00_peaks.narrowPeak | head -n 500000 > output/idr/sp6/OPT_PEAKS/sp6_rep1_ps00_top.bed
sort -k8,8nr output/idr/sp6/OPT_PEAKS/sp6_rep1_ps01_peaks.narrowPeak | head -n 500000 > output/idr/sp6/OPT_PEAKS/sp6_rep1_ps01_top.bed
sort -k8,8nr output/idr/sp6/OPT_PEAKS/sp6_rep2_ps00_peaks.narrowPeak | head -n 500000 > output/idr/sp6/OPT_PEAKS/sp6_rep2_ps00_top.bed
sort -k8,8nr output/idr/sp6/OPT_PEAKS/sp6_rep2_ps01_peaks.narrowPeak | head -n 500000 > output/idr/sp6/OPT_PEAKS/sp6_rep2_ps01_top.bed
sort -k8,8nr output/idr/sp6/OPT_PEAKS/sp6_pool_ps00_peaks.narrowPeak | head -n 500000 > output/idr/sp6/OPT_PEAKS/sp6_pool_ps00_top.bed
sort -k8,8nr output/idr/sp6/OPT_PEAKS/sp6_pool_ps01_peaks.narrowPeak | head -n 500000 > output/idr/sp6/OPT_PEAKS/sp6_pool_ps01_top.bed

micromamba deactivate

bedtools intersect -v -a output/idr/sp6/CONS_PEAKS/sp6_rep1_top.bed -b input/mm10-blacklist.v2.bed > output/idr/sp6/CONS_PEAKS/sp6_rep1_top_bl.bed
bedtools intersect -v -a output/idr/sp6/CONS_PEAKS/sp6_rep2_top.bed -b input/mm10-blacklist.v2.bed > output/idr/sp6/CONS_PEAKS/sp6_rep2_top_bl.bed
bedtools intersect -v -a output/idr/sp6/CONS_PEAKS/sp6_pooled_top.bed -b input/mm10-blacklist.v2.bed > output/idr/sp6/CONS_PEAKS/sp6_pooled_top_bl.bed
bedtools intersect -v -a output/idr/sp6/OPT_PEAKS/sp6_rep1_ps00_top.bed -b input/mm10-blacklist.v2.bed > output/idr/sp6/OPT_PEAKS/sp6_rep1_ps00_top_bl.bed
bedtools intersect -v -a output/idr/sp6/OPT_PEAKS/sp6_rep1_ps01_top.bed -b input/mm10-blacklist.v2.bed > output/idr/sp6/OPT_PEAKS/sp6_rep1_ps01_top_bl.bed
bedtools intersect -v -a output/idr/sp6/OPT_PEAKS/sp6_rep2_ps00_top.bed -b input/mm10-blacklist.v2.bed > output/idr/sp6/OPT_PEAKS/sp6_rep2_ps00_top_bl.bed
bedtools intersect -v -a output/idr/sp6/OPT_PEAKS/sp6_rep2_ps01_top.bed -b input/mm10-blacklist.v2.bed > output/idr/sp6/OPT_PEAKS/sp6_rep2_ps01_top_bl.bed
bedtools intersect -v -a output/idr/sp6/OPT_PEAKS/sp6_pool_ps00_top.bed -b input/mm10-blacklist.v2.bed > output/idr/sp6/OPT_PEAKS/sp6_pool_ps00_top_bl.bed
bedtools intersect -v -a output/idr/sp6/OPT_PEAKS/sp6_pool_ps01_top.bed -b input/mm10-blacklist.v2.bed > output/idr/sp6/OPT_PEAKS/sp6_pool_ps01_top_bl.bed

micromamba activate IDR_python3

idr --input-file-type narrowPeak --rank p.value --soft-idr-threshold 0.1 -s output/idr/sp6/CONS_PEAKS/sp6_rep1_top_bl.bed output/idr/sp6/CONS_PEAKS/sp6_rep2_top_bl.bed -p output/idr/sp6/CONS_PEAKS/sp6_pooled_top_bl.bed --verbose -l output/idr/sp6/CONS_PEAKS/idr_conservative.log -o output/idr/sp6/CONS_PEAKS/sp6_rep1_rep2_idr.txt --plot
idr --input-file-type narrowPeak --rank p.value --soft-idr-threshold 0.1 -s output/idr/sp6/OPT_PEAKS/sp6_rep1_ps00_top_bl.bed output/idr/sp6/OPT_PEAKS/sp6_rep1_ps01_top_bl.bed -p output/idr/sp6/CONS_PEAKS/sp6_rep1_top_bl.bed --verbose -l output/idr/sp6/OPT_PEAKS/idr_opt_rep1.log -o output/idr/sp6/OPT_PEAKS/sp6_rep1_ps_idr.txt --plot
idr --input-file-type narrowPeak --rank p.value --soft-idr-threshold 0.1 -s output/idr/sp6/OPT_PEAKS/sp6_rep2_ps00_top_bl.bed output/idr/sp6/OPT_PEAKS/sp6_rep2_ps01_top_bl.bed -p output/idr/sp6/CONS_PEAKS/sp6_rep2_top_bl.bed --verbose -l output/idr/sp6/OPT_PEAKS/idr_opt_rep2.log -o output/idr/sp6/OPT_PEAKS/sp6_rep2_ps_idr.txt --plot
idr --input-file-type narrowPeak --rank p.value --soft-idr-threshold 0.1 -s output/idr/sp6/OPT_PEAKS/sp6_pool_ps00_top_bl.bed output/idr/sp6/OPT_PEAKS/sp6_pool_ps01_top_bl.bed -p output/idr/sp6/CONS_PEAKS/sp6_pooled_top_bl.bed --verbose -l output/idr/sp6/OPT_PEAKS/idr_opt_pool.log -o output/idr/sp6/OPT_PEAKS/sp6_pool_ps_idr.txt --plot

micromamba deactivate

awk 'BEGIN{OFS="\t"} $12>=1 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' output/idr/sp6/CONS_PEAKS/sp6_rep1_rep2_idr.txt | sort -k7n,7n > output/idr/sp6/TRACKS/sp6_idrConsPeaks.bed
awk 'BEGIN{OFS="\t"} $12>=1 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' output/idr/sp6/OPT_PEAKS/sp6_rep1_ps_idr.txt | sort -k7n,7n > output/idr/sp6/OPT_PEAKS/sp6_idrOptPeaks_rep1_ps.bed
awk 'BEGIN{OFS="\t"} $12>=1 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' output/idr/sp6/OPT_PEAKS/sp6_rep2_ps_idr.txt | sort -k7n,7n > output/idr/sp6/OPT_PEAKS/sp6_idrOptPeaks_rep2_ps.bed
awk 'BEGIN{OFS="\t"} $12>=1 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' output/idr/sp6/OPT_PEAKS/sp6_pool_ps_idr.txt | sort -k7n,7n > output/idr/sp6/OPT_PEAKS/sp6_idrOptPeaks_pool_ps.bed

N_T=$(wc -l < output/idr/sp6/TRACKS/sp6_idrConsPeaks.bed)
N_1=$(wc -l < output/idr/sp6/OPT_PEAKS/sp6_idrOptPeaks_rep1_ps.bed)
N_2=$(wc -l < output/idr/sp6/OPT_PEAKS/sp6_idrOptPeaks_rep2_ps.bed)
N_P=$(wc -l < output/idr/sp6/OPT_PEAKS/sp6_idrOptPeaks_pool_ps.bed)
echo "Sp6 IDR peaks -- conservative: ${N_T}, optimal(pool): ${N_P}, optimal(rep1): ${N_1}, optimal(rep2): ${N_2}"
awk -v a="${N_T}" -v b="${N_P}" 'BEGIN{m=(a>b?a:b); n=(a<b?a:b); printf "Sp6 rescue ratio: %.3f\n", m/n}'
awk -v a="${N_1}" -v b="${N_2}" 'BEGIN{m=(a>b?a:b); n=(a<b?a:b); printf "Sp6 self-consistency ratio: %.3f\n", m/n}'

# ---- Sp8 ----

cat output/idr/sp8/roinp_sp8_rep1.bedpe output/idr/sp8/roinp_sp8_rep2.bedpe > output/idr/sp8/input_pooled.bedpe
cat output/idr/sp8/rochp_sp8_rep1.bedpe output/idr/sp8/rochp_sp8_rep2.bedpe > output/idr/sp8/rochp_sp8_pooled.bedpe

micromamba activate macs2_python_37

macs2 callpeak -f BEDPE -t output/idr/sp8/rochp_sp8_rep1.bedpe -c output/idr/sp8/roinp_sp8_rep1.bedpe -p 0.02 --keep-dup auto -g mm -n sp8_rep1 --outdir output/idr/sp8/CONS_PEAKS
macs2 callpeak -f BEDPE -t output/idr/sp8/rochp_sp8_rep2.bedpe -c output/idr/sp8/roinp_sp8_rep2.bedpe -p 0.02 --keep-dup auto -g mm -n sp8_rep2 --outdir output/idr/sp8/CONS_PEAKS
macs2 callpeak -f BEDPE -t output/idr/sp8/rochp_sp8_rep1.bedpe output/idr/sp8/rochp_sp8_rep2.bedpe -c output/idr/sp8/input_pooled.bedpe -p 0.02 --keep-dup auto -g mm -n sp8_pooled --outdir output/idr/sp8/CONS_PEAKS

sort -k8,8nr output/idr/sp8/CONS_PEAKS/sp8_rep1_peaks.narrowPeak | head -n 500000 > output/idr/sp8/CONS_PEAKS/sp8_rep1_top.bed
sort -k8,8nr output/idr/sp8/CONS_PEAKS/sp8_rep2_peaks.narrowPeak | head -n 500000 > output/idr/sp8/CONS_PEAKS/sp8_rep2_top.bed
sort -k8,8nr output/idr/sp8/CONS_PEAKS/sp8_pooled_peaks.narrowPeak | head -n 500000 > output/idr/sp8/CONS_PEAKS/sp8_pooled_top.bed

shuf output/idr/sp8/rochp_sp8_rep1.bedpe | split -d -l $(( ($(wc -l < output/idr/sp8/rochp_sp8_rep1.bedpe) + 1) / 2 )) - output/idr/sp8/OPT_PEAKS/sp8_rep1_ps
shuf output/idr/sp8/rochp_sp8_rep2.bedpe | split -d -l $(( ($(wc -l < output/idr/sp8/rochp_sp8_rep2.bedpe) + 1) / 2 )) - output/idr/sp8/OPT_PEAKS/sp8_rep2_ps
shuf output/idr/sp8/rochp_sp8_pooled.bedpe | split -d -l $(( ($(wc -l < output/idr/sp8/rochp_sp8_pooled.bedpe) + 1) / 2 )) - output/idr/sp8/OPT_PEAKS/sp8_pool_ps

macs2 callpeak -f BEDPE -t output/idr/sp8/OPT_PEAKS/sp8_rep1_ps00 -c output/idr/sp8/roinp_sp8_rep1.bedpe -p 0.02 --keep-dup auto -g mm -n sp8_rep1_ps00 --outdir output/idr/sp8/OPT_PEAKS
macs2 callpeak -f BEDPE -t output/idr/sp8/OPT_PEAKS/sp8_rep1_ps01 -c output/idr/sp8/roinp_sp8_rep1.bedpe -p 0.02 --keep-dup auto -g mm -n sp8_rep1_ps01 --outdir output/idr/sp8/OPT_PEAKS
macs2 callpeak -f BEDPE -t output/idr/sp8/OPT_PEAKS/sp8_rep2_ps00 -c output/idr/sp8/roinp_sp8_rep2.bedpe -p 0.02 --keep-dup auto -g mm -n sp8_rep2_ps00 --outdir output/idr/sp8/OPT_PEAKS
macs2 callpeak -f BEDPE -t output/idr/sp8/OPT_PEAKS/sp8_rep2_ps01 -c output/idr/sp8/roinp_sp8_rep2.bedpe -p 0.02 --keep-dup auto -g mm -n sp8_rep2_ps01 --outdir output/idr/sp8/OPT_PEAKS
macs2 callpeak -f BEDPE -t output/idr/sp8/OPT_PEAKS/sp8_pool_ps00 -c output/idr/sp8/input_pooled.bedpe -p 0.02 --keep-dup auto -g mm -n sp8_pool_ps00 --outdir output/idr/sp8/OPT_PEAKS
macs2 callpeak -f BEDPE -t output/idr/sp8/OPT_PEAKS/sp8_pool_ps01 -c output/idr/sp8/input_pooled.bedpe -p 0.02 --keep-dup auto -g mm -n sp8_pool_ps01 --outdir output/idr/sp8/OPT_PEAKS

sort -k8,8nr output/idr/sp8/OPT_PEAKS/sp8_rep1_ps00_peaks.narrowPeak | head -n 500000 > output/idr/sp8/OPT_PEAKS/sp8_rep1_ps00_top.bed
sort -k8,8nr output/idr/sp8/OPT_PEAKS/sp8_rep1_ps01_peaks.narrowPeak | head -n 500000 > output/idr/sp8/OPT_PEAKS/sp8_rep1_ps01_top.bed
sort -k8,8nr output/idr/sp8/OPT_PEAKS/sp8_rep2_ps00_peaks.narrowPeak | head -n 500000 > output/idr/sp8/OPT_PEAKS/sp8_rep2_ps00_top.bed
sort -k8,8nr output/idr/sp8/OPT_PEAKS/sp8_rep2_ps01_peaks.narrowPeak | head -n 500000 > output/idr/sp8/OPT_PEAKS/sp8_rep2_ps01_top.bed
sort -k8,8nr output/idr/sp8/OPT_PEAKS/sp8_pool_ps00_peaks.narrowPeak | head -n 500000 > output/idr/sp8/OPT_PEAKS/sp8_pool_ps00_top.bed
sort -k8,8nr output/idr/sp8/OPT_PEAKS/sp8_pool_ps01_peaks.narrowPeak | head -n 500000 > output/idr/sp8/OPT_PEAKS/sp8_pool_ps01_top.bed

micromamba deactivate

bedtools intersect -v -a output/idr/sp8/CONS_PEAKS/sp8_rep1_top.bed -b input/mm10-blacklist.v2.bed > output/idr/sp8/CONS_PEAKS/sp8_rep1_top_bl.bed
bedtools intersect -v -a output/idr/sp8/CONS_PEAKS/sp8_rep2_top.bed -b input/mm10-blacklist.v2.bed > output/idr/sp8/CONS_PEAKS/sp8_rep2_top_bl.bed
bedtools intersect -v -a output/idr/sp8/CONS_PEAKS/sp8_pooled_top.bed -b input/mm10-blacklist.v2.bed > output/idr/sp8/CONS_PEAKS/sp8_pooled_top_bl.bed
bedtools intersect -v -a output/idr/sp8/OPT_PEAKS/sp8_rep1_ps00_top.bed -b input/mm10-blacklist.v2.bed > output/idr/sp8/OPT_PEAKS/sp8_rep1_ps00_top_bl.bed
bedtools intersect -v -a output/idr/sp8/OPT_PEAKS/sp8_rep1_ps01_top.bed -b input/mm10-blacklist.v2.bed > output/idr/sp8/OPT_PEAKS/sp8_rep1_ps01_top_bl.bed
bedtools intersect -v -a output/idr/sp8/OPT_PEAKS/sp8_rep2_ps00_top.bed -b input/mm10-blacklist.v2.bed > output/idr/sp8/OPT_PEAKS/sp8_rep2_ps00_top_bl.bed
bedtools intersect -v -a output/idr/sp8/OPT_PEAKS/sp8_rep2_ps01_top.bed -b input/mm10-blacklist.v2.bed > output/idr/sp8/OPT_PEAKS/sp8_rep2_ps01_top_bl.bed
bedtools intersect -v -a output/idr/sp8/OPT_PEAKS/sp8_pool_ps00_top.bed -b input/mm10-blacklist.v2.bed > output/idr/sp8/OPT_PEAKS/sp8_pool_ps00_top_bl.bed
bedtools intersect -v -a output/idr/sp8/OPT_PEAKS/sp8_pool_ps01_top.bed -b input/mm10-blacklist.v2.bed > output/idr/sp8/OPT_PEAKS/sp8_pool_ps01_top_bl.bed

micromamba activate IDR_python3

idr --input-file-type narrowPeak --rank p.value --soft-idr-threshold 0.1 -s output/idr/sp8/CONS_PEAKS/sp8_rep1_top_bl.bed output/idr/sp8/CONS_PEAKS/sp8_rep2_top_bl.bed -p output/idr/sp8/CONS_PEAKS/sp8_pooled_top_bl.bed --verbose -l output/idr/sp8/CONS_PEAKS/idr_conservative.log -o output/idr/sp8/CONS_PEAKS/sp8_rep1_rep2_idr.txt --plot
idr --input-file-type narrowPeak --rank p.value --soft-idr-threshold 0.1 -s output/idr/sp8/OPT_PEAKS/sp8_rep1_ps00_top_bl.bed output/idr/sp8/OPT_PEAKS/sp8_rep1_ps01_top_bl.bed -p output/idr/sp8/CONS_PEAKS/sp8_rep1_top_bl.bed --verbose -l output/idr/sp8/OPT_PEAKS/idr_opt_rep1.log -o output/idr/sp8/OPT_PEAKS/sp8_rep1_ps_idr.txt --plot
idr --input-file-type narrowPeak --rank p.value --soft-idr-threshold 0.1 -s output/idr/sp8/OPT_PEAKS/sp8_rep2_ps00_top_bl.bed output/idr/sp8/OPT_PEAKS/sp8_rep2_ps01_top_bl.bed -p output/idr/sp8/CONS_PEAKS/sp8_rep2_top_bl.bed --verbose -l output/idr/sp8/OPT_PEAKS/idr_opt_rep2.log -o output/idr/sp8/OPT_PEAKS/sp8_rep2_ps_idr.txt --plot
idr --input-file-type narrowPeak --rank p.value --soft-idr-threshold 0.1 -s output/idr/sp8/OPT_PEAKS/sp8_pool_ps00_top_bl.bed output/idr/sp8/OPT_PEAKS/sp8_pool_ps01_top_bl.bed -p output/idr/sp8/CONS_PEAKS/sp8_pooled_top_bl.bed --verbose -l output/idr/sp8/OPT_PEAKS/idr_opt_pool.log -o output/idr/sp8/OPT_PEAKS/sp8_pool_ps_idr.txt --plot

micromamba deactivate

awk 'BEGIN{OFS="\t"} $12>=1 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' output/idr/sp8/CONS_PEAKS/sp8_rep1_rep2_idr.txt | sort -k7n,7n > output/idr/sp8/TRACKS/sp8_idrConsPeaks.bed
awk 'BEGIN{OFS="\t"} $12>=1 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' output/idr/sp8/OPT_PEAKS/sp8_rep1_ps_idr.txt | sort -k7n,7n > output/idr/sp8/OPT_PEAKS/sp8_idrOptPeaks_rep1_ps.bed
awk 'BEGIN{OFS="\t"} $12>=1 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' output/idr/sp8/OPT_PEAKS/sp8_rep2_ps_idr.txt | sort -k7n,7n > output/idr/sp8/OPT_PEAKS/sp8_idrOptPeaks_rep2_ps.bed
awk 'BEGIN{OFS="\t"} $12>=1 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' output/idr/sp8/OPT_PEAKS/sp8_pool_ps_idr.txt | sort -k7n,7n > output/idr/sp8/OPT_PEAKS/sp8_idrOptPeaks_pool_ps.bed

N_T=$(wc -l < output/idr/sp8/TRACKS/sp8_idrConsPeaks.bed)
N_1=$(wc -l < output/idr/sp8/OPT_PEAKS/sp8_idrOptPeaks_rep1_ps.bed)
N_2=$(wc -l < output/idr/sp8/OPT_PEAKS/sp8_idrOptPeaks_rep2_ps.bed)
N_P=$(wc -l < output/idr/sp8/OPT_PEAKS/sp8_idrOptPeaks_pool_ps.bed)
echo "Sp8 IDR peaks -- conservative: ${N_T}, optimal(pool): ${N_P}, optimal(rep1): ${N_1}, optimal(rep2): ${N_2}"
awk -v a="${N_T}" -v b="${N_P}" 'BEGIN{m=(a>b?a:b); n=(a<b?a:b); printf "Sp8 rescue ratio: %.3f\n", m/n}'
awk -v a="${N_1}" -v b="${N_2}" 'BEGIN{m=(a>b?a:b); n=(a<b?a:b); printf "Sp8 self-consistency ratio: %.3f\n", m/n}'

# =============================================================================
# 9. Heatmaps at IDR optimal (pooled-pseudoreplicate) peaks (deepTools)
#    Sp8 signal, Sp6 signal, ATAC-seq, each ±3 kb from peak center
# =============================================================================

awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3}' output/idr/sp6/OPT_PEAKS/sp6_idrOptPeaks_pool_ps.bed | sort -k1,1 -k2,2n > output/idr/sp6/sp6_peaks_4c_sorted.bed
awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3}' output/idr/sp8/OPT_PEAKS/sp8_idrOptPeaks_pool_ps.bed | sort -k1,1 -k2,2n > output/idr/sp8/sp8_peaks_4c_sorted.bed

computeMatrix reference-point -S output/bigwig/rochp_sp8_rep1.bw output/bigwig/rochp_sp8_rep2.bw -R output/idr/sp8/sp8_peaks_4c_sorted.bed -a 3000 -b 3000 -p 16 --referencePoint center -o output/matrix/sp8_peaks_in_sp8.gz
plotHeatmap -m output/matrix/sp8_peaks_in_sp8.gz -o output/heatmaps/sp8_peaks_in_sp8.pdf --startLabel center --colorMap bwr --legendLocation none

computeMatrix reference-point -S output/bigwig/chipd_sp6_rep1.bw output/bigwig/chipd_sp6_rep2.bw -R output/idr/sp6/sp6_peaks_4c_sorted.bed -a 3000 -b 3000 -p 16 --referencePoint center -o output/matrix/sp6_peaks_in_sp6.gz
plotHeatmap -m output/matrix/sp6_peaks_in_sp6.gz -o output/heatmaps/sp6_peaks_in_sp6.pdf --startLabel center --colorMap bwr --legendLocation none

computeMatrix reference-point -S input/atac_ectoderm.bw -R output/idr/sp8/sp8_peaks_4c_sorted.bed -a 3000 -b 3000 -p 16 --referencePoint center -o output/matrix/sp8_peaks_in_atac.gz
plotHeatmap -m output/matrix/sp8_peaks_in_atac.gz -o output/heatmaps/sp8_peaks_in_atac.pdf --startLabel center --colorList blue,white,red --legendLocation none --zMax 1

computeMatrix reference-point -S input/atac_ectoderm.bw -R output/idr/sp6/sp6_peaks_4c_sorted.bed -a 3000 -b 3000 -p 16 --referencePoint center -o output/matrix/sp6_peaks_in_atac.gz
plotHeatmap -m output/matrix/sp6_peaks_in_atac.gz -o output/heatmaps/sp6_peaks_in_atac.pdf --startLabel center --colorList blue,white,red --legendLocation none --zMax 1

# =============================================================================
# 10. Proximal / distal annotation of IDR optimal peaks (bedtools closest,
#     nearest TSS, signed distance in column 9, threshold = +/- 2 kb)
# =============================================================================

bedtools closest -D ref -t first -a output/idr/sp6/sp6_peaks_4c_sorted.bed -b input/tss.sorted.bed > output/annot/sp6_distance_to_tss.bed
awk '{ if ($9>=-2000 && $9<=2000) print > "output/annot/sp6_tss_proximal_2kb.bed"; else print > "output/annot/sp6_tss_distal_2kb.bed" }' output/annot/sp6_distance_to_tss.bed

bedtools closest -D ref -t first -a output/idr/sp8/sp8_peaks_4c_sorted.bed -b input/tss.sorted.bed > output/annot/sp8_distance_to_tss.bed
awk '{ if ($9>=-2000 && $9<=2000) print > "output/annot/sp8_tss_proximal_2kb.bed"; else print > "output/annot/sp8_tss_distal_2kb.bed" }' output/annot/sp8_distance_to_tss.bed

echo "Done. IDR optimal peaks: output/idr/{sp6,sp8}/OPT_PEAKS/{sp6,sp8}_idrOptPeaks_pool_ps.bed"
echo "Proximal/distal BEDs:    output/annot/{sp6,sp8}_tss_{proximal,distal}_2kb.bed"
