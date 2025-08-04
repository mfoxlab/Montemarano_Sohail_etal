# Project Title

A short description of your project goes here.

## Table of Contents

-   [R and Package Versions](#1-versions)
-   [Raw Files](#2-raw-files)
-   [R Scripts](#3-r-scripts)

## 1. R and Package Versions

This project was run using R 4.4.1. This version is recommended for full reproducibility and compatibility with packages.

All scripts will have the `library()` commands, but not the `install::packages()` ones. To install all proper packages necessary for this project, please download and run [`installpackages.R`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/R%20Scripts/installpackages.R). Exact package versions are noted, however, running this script as is will also reproduce the data.

## 2. Raw Files

Next, please download the following files and place them into your working directory:

[`counts.Rdata`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/Raw%20Files/counts.Rdata)

-   Contains raw counts of 56,980 genes (after filtering in Galaxy?) for numerous group conditions (Used for Step 1)

[`meta.Rdata`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/Raw%20Files/meta.Rdata)

-   Contains metadata for variables involved in the RNASeq analysis: sex, drug, VTA target, etc. (Step 1)

[`modulecolors.txt`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/Raw%20Files/modulecolors.txt)

-   Color key for gene modules and creating plots (Steps 4 and 5).

[`manual_go_terms.txt`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/Raw%20Files/manual_go_terms.txt)

-   For use in gene ontology analysis (Step 5).

## 3. R Scripts

Download and run the following files in your working directory:

*Note: Where there is a line named* `workdir = ""` *or something similar, please write your file path to set your own working directory. Skip the* `load(file = paste0(workdir,"/FILENAME_workspace.Rdata"))` *line in your first run.*

1.  [`DElimma.R`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/R%20Scripts/DElimma.R)

-   Perform quality control to remove lowly-expressed reads
-   Sets up differential expression based on comparisons between sex, drug (fentanyl vs. saline), and VTA projection target (PFC, NAC, or INPUT (total VTA)).
-   Produces a `DEoutputp05` file that ???
-   Produces three `DEgenesummary` files for VTA-PFC, VTA-NAC, and VTA-INPUT genes.
-   *Note: Gene information for VTA to amygdala neurons (DEgene_summary_AMY) can be ignored. Although this RNASeq data contains VTA-AMY neuron information, this project does not consider AMY as a variable to create WGCNA networks.*

2.  [`WGCNA_PFC.R`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/R%20Scripts/WGCNA_PFC.R), [`WGCNA_NAC.R`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/R%20Scripts/WGCNA_NAC.R), and [`WGCNA_INPUT.R`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/R%20Scripts/WGCNA_INPUT.R)

-   *Note: Mac users are recommended to download the* `netREGION_8.Rdata` *file respective to each VTA projection. The script will have two commented out lines of code that need to be run. Run those two lines and skip the rest of that section.*

-   Builds WGCNA networks with only differentially expressed genes as a function of drug, projection target and sex.
-   Makes `netREGION_8.Rdata` to create gene networks, modules, hub genes, and `REGION_gene_summary.txt` files for each projection. Three moduleEigengengesDE .txt files will be further used to create circos data and visualizations.

3.  [`circosplotPFCdata.R`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/R%20Scripts/circosplotPFCdata.R), [`circosplotNAC.Rdata`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/R%20Scripts/circosplotNACdata.R), and [`circosplotINPUTdata.R`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/R%20Scripts/circosplotINPUTdata.R)

-   Organizes data to be used in circos plots. Merges log fold-change and P-Value data from sex-combined data with female and male data.

4.  [`circosplotPFC.R`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/R%20Scripts/circosplotPFC.R), [`circosplotNAC.R`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/R%20Scripts/circosplotNAC.R), and [`circosplotINPUT.R`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/R%20Scripts/circosplotINPUT.R)

-   Generates four PDF plots:

    1.  circos plot with significant modules (unlabeled)
    2.  circos plot with significant modules (labeled)
    3.  circos plot with all modules (unlabeled)
    4.  circos plot with all modules (labeled)

-   Run this entire script at once (as objects in this environment are overwritten to create each new plot).

5.  [`go_PFC.R`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/R%20Scripts/go_PFC.R), [`go_NAC.R`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/R%20Scripts/go_NAC.R), [`go_INPUT.R`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/R%20Scripts/go_INPUT.R)

-   *Note: Rerunning these GOst queries may be subject to change over time. For exact results, download the* `REGION_go_all_modules_cleaned_colorkey_7_30_2025.Rdata` *file respective to each VTA projection. The script will have a line of code where it will load it into your workspace.*

-   Performs gene ontology and KEGG analyses via gprofiler2. Outputs 2 .txt files.

-   Visualizes the top enriched GO terms in dot plots.

6.  [`RRHO.R`](https://github.com/mfoxlab/Montemarano_Sohail_etal/blob/main/R%20Scripts/RRHO.R)

-   Uses rank-rank hypergeometric overlap between two ranked gene lists to plot heatmaps. All combinations of projection, sex, and drug are compared.
