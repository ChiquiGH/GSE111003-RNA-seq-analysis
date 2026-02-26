# GSE111003 – RNA-seq analysis (β-glucan, early response)

This repository contains the scripts used to analyze bulk RNA-seq data from GEO (GSE111003),
focusing on early transcriptional responses to β-glucan stimulation (T0, 4h, 24h) with time-matched RPMI controls.

## Repository structure
- `scripts/`: analysis pipeline steps (R scripts)
- `SRC/`: helper functions
- `report/`: R Markdown used to generate the final PDF report
- `results/`: output tables and figures

## Main methods
- Filtering: edgeR `filterByExpr()`
- Normalization: TMM (edgeR)
- Multivariate analysis: PCA (logCPM TMM, scaled)
- Differential expression: DESeq2 (`~ Donor + Condition`)

## Data availability
Raw mmseq files were downloaded from GEO (GSE111003) and are not included in this repository due to size.