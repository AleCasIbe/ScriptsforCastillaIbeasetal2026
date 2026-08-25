# Software, databases and reference data

Companion record for the Sp6/Sp8 limb ectoderm manuscript.

`sessionInfo.txt` (written by `figures_and_source_data.R`) records the **R**
environment only. This file records everything else: command-line tools, web
services, databases and reference assemblies. Together they are the complete
computational provenance for the Code Availability statement.

---

## 1. R environment

Full detail in [`output/sessionInfo.txt`](output/sessionInfo.txt), regenerated on
every run of `figures_and_source_data.R`. Key packages as installed in the
environment that produced the figures and Source Data:

| Package | Version | Used for |
|---|---|---|
| R | 4.4.2 | Figure and Source Data generation |
| Bioconductor | 3.20 | — |
| DESeq2 | 1.46.0 | Differential expression |
| SPIA | 2.58.0 | Pathway impact analysis |
| KEGGREST | 1.46.0 | KEGG pathway membership |
| org.Mm.eg.db | 3.20.0 | Symbol → Entrez mapping |
| GenomicRanges | 1.58.0 | Peak / VISTA interval overlaps |
| Seurat | 5.2.1 | scRNA-seq classification and plots |
| ComplexHeatmap | 2.25.2 | L2FC heatmaps |
| circlize | 0.4.16 | Heatmap colour ramps |
| ggplot2 | 4.0.1 | All ggplot figures |
| eulerr | 7.0.2 | Euler / Venn diagrams |
| ggpubr | 0.6.2 | Table panels |
| pheatmap | 1.0.13 | QC heatmaps |
| monarchr | 2.1 | Monarch knowledge-graph queries |
| homologene | 1.4.68.19.3.27 | Human → mouse orthologs |
| biomaRt | 2.62.1 | Annotation lookups |
| data.table | 1.17.8 | MGI / HPO table handling |
| xlsx / writexl / readxl | 0.6.5 / 1.5.4 / 1.4.5 | Excel I/O |
| sessioninfo | 1.2.3 | Environment capture |

---

## 2. Command-line tools

| Tool | Version | Step |
|---|---|---|
| fastp | 0.23.4 | Adapter/quality trimming (ChIPmentation and RNA-seq) |
| BWA | 0.7.17 | `bwa mem` (Sp6, paired-end); `bwa aln`/`samse` (Sp8, single-end) |
| SAMtools | 1.19 | Coordinate sorting |
| sambamba | 1.0 | Unique, non-duplicate read filtering |
| MACS2 | 2.2.7.1 | Narrow peak calling (`-p 0.02 --keep-dup auto`) |
| IDR | 2.0.4 | Reproducible peak sets (threshold 0.1) |
| bedtools | 2.31.1   | `closest`, `merge`, `intersect` |
| deepTools | 3.5.6 | `multiBamSummary`, `plotCorrelation`, `computeMatrix`, `plotHeatmap`, `bamCoverage` |
| seqMINER | 1.3.4 | K-means ranked clustering (10 clusters, seed 74676790) |
| MEME Suite | 5.0.1 | MEME-ChIP, FIMO, AME |
| FastQC | 0.12.1 | RNA-seq read QC |
| HISAT2 | 2.2.1 | RNA-seq alignment (`--dta-cufflinks`) |
| StringTie | 3.0.3 | Transcript assembly and quantification |
| Cell Ranger | 3.0.1 | scRNA-seq preprocessing (10x Genomics) |

---

## 3. Web services


| Service | Version  | Used for |
|---|---|---|---|
| GREAT | 4.0.4 | Gene assignment (basal plus extension) and GO enrichment |
| DAVID | 6.8 | GO Biological Process enrichment |
| Cistrome Analysis Pipeline | Conservation Plot 1.0.0 | PhastCons conservation profiles |
| MEME-ChIP (web) | 5.0.1 | De novo motif discovery |


---

## 4. Databases and annotation

| Resource | Version / release | Used for |
|---|---|---|---|
| GRCm38 / mm10 primary assembly | GRCm38  | ChIP and RNA-seq alignment |
| Ensembl annotation | GRCm38.93 | Cell Ranger reference |
| ENCODE blacklist (mm10) | v2 | Peak filtering |
| HOCOMOCO | H12 CORE (human + mouse) | Known-motif matching |
| KEGG (live REST) | current at run time  | Pathway gene membership for heatmaps |
| MGI | `MGI_DiseaseGeneModel.rpt` file dated 2026-06-25 | Mouse limb / ectodermal phenotypes |
| Human Phenotype Ontology | release v2025-09-01 | Human limb / ectodermal phenotypes |
| Monarch Initiative | via monarchr 2.1 | uPheno gene–phenotype associations |
| VISTA Enhancer Browser | `experiments.tsv` file dated 2025-02-12 | Curated enhancer overlap |

---

## 5. Public datasets

| Accession | Description | Used for |
|---|---|---|
| GSE145657 | Single-cell ATAC-seq of E11.5 forelimb bud, WT and *Hox13* mutant (10x Chromium); Desanlis et al. 2020 | Ectodermal accessibility track — sample GSM4323465, file `GSM4323465_WT.8_treat_pileup.bw` |
| GSE149368 | E10.5 mouse limb bud scRNA-seq; Allou et al. 2021 | Ectodermal expression classification |

---

## 6. Analysis scripts

| Script | Purpose |
|---|---|
| `figures_and_source_data.R` | Master script: figure panels, Supplementary Data 1–8, Source Data workbook, `sessionInfo.txt` |
| `1_ChIP_processing.R` | ChIP-seq processing |

---



