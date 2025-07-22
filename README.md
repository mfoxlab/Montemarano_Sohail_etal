# Project Title

A short description of your project goes here.

## Table of Contents

-   [R and Package Versions](#1-versions)
-   [Raw Files](#2-raw)
-   [R Scripts](#3-scripts)

## 1. Versions 

This project was run using R 4.4.1. This version is recommended for full reproducibility and compatibility with packages.

All scripts will have the `library()` commands, but not the `install::packages()` ones. To install all proper packages necessary for this project, please follow the instructions below.

How to use:

1.  Download the [`renv.lock`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/renv.lock) file and place into your working directory.
2.  In R or RStudio, set your working directory to your project folder.
3.  In your console, run:

```         
install.packages("renv") # just once, if not already installed
renv::restore()
```

This reads the [`renv.lock`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/renv.lock) file and downloads/installs all the required packages that have been used (and the exact versions used) into a local renv library. This does not conflict with your global R packages.

4.  Further R scripts may be run in this environment.

*Note: Due to the complexity of installing BiocManager packages, `Go.db` (v 3.19.1) was not added to the `renv.lock` file. Please manually install this package.*

## 2. Raw Files

Next, please download the following two files and place them into your working directory:

[`counts.Rdata`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/Raw%20Files/counts.Rdata)

-   Contains raw counts of 56,980 genes (after filtering in Galaxy?) for numerous group conditions.

[`meta.Rdata`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/Raw%20Files/meta.Rdata)

-   Contains metadata for variables involved in the RNASeq analysis: sex, drug, VTA target, etc.

## 3. R Scripts

Download and run the following files in your working directory:

*Note: Where there is a line named* `workdir = ""` *or something similar, please write your file path to set your own working directory. Skip the* `load(file = paste0(workdir,"/FILENAME_workspace.Rdata"))` *line in your first run.*

1\. [`WGCNA.R`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/R%20Scripts/WGCNA.R)

-   Perform quality control to remove lowly-expressed reads
-   Sets up differential expression based on comparisons between sex, drug (fentanyl vs. saline), and VTA projection target (PFC, NAC, or INPUT (total VTA)).
-   Produces a `DEoutputp05` file that ???
-   Produces three `DEgenesummary` files for VTA-PFC, VTA-NAC, and VTA-INPUT genes.
-   *Note: Gene information for VTA to amygdala neurons (DEgene_summary_AMY) can be ignored. Although this RNASeq data contains VTA-AMY neuron information, this project does not consider AMY as a variable to create WGCNA networks.*

2.  `logCPM_PFC.R`, `logCPM_NAC.R`, and `logCPM_INPUT.R`

-   Builds WGCNA networks with only differentially expressed genes as a function of drug, projection target and sex..
-   Makes `netREGION_8.Rdata` to create gene networks, modules, hub genes, and `REGION_gene_summary.txt` files for each projection. Three moduleEigengengesDE .txt files will be further used to create circos data and visualizations.

3.  `circosplotPFCdata.R`, `circosplotNAC.Rdata`, and `circosplotINPUTdata.R`

-   Organizes data to be used in circos plots. Merges log fold-change and P-Value data from sex-combined data with female and male data.

4.  `circosplotPFC.R`, `circosplotNAC.R`, and `circosplotINPUT.R`

-   *Note: Please download and place `colors.txt` into your working directory before running.*

-   Generates four PDF plots:

    1.  circos plot with significant modules (unlabeled)
    2.  circos plot with significant modules (labeled)
    3.  circos plot with all modules (unlabeled)
    4.  circos plot with all modules (labeled)

-   Run this entire script at once (as objects in this environment are overwritten to create each new plot).

5.  `go_PFC.R`, `go_NAC.R`, `go_INPUT.R`

-   Performs gene ontology and KEGG analyses via gprofiler2. Outputs 2 .txt files.
-   Visualizes the top enriched GO terms in dot plots.
