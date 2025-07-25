# Install all libraries

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")           # used v3.60.6

BiocManager::install("edgeR")           # used v4.2.2

install.packages("data.table")          # used v1.17.6

install.packages("dplyr")               # used v1.1.4

BiocManager::install("biomaRt")         # used v2.60.1

install.packages("WGCNA")               # used v1.73
 
install.packages("pheatmap")            # used v1.0.13

BiocManager::install("GO.db")           # used v3.19.1

BiocManager::install("impute")          # used v1.78.0

BiocManager::install("preprocessCore")  # used v1.66.0

install.packages("circlize")            # used v0.4.16

install.packages("readxl")              # used v1.4.3

BiocManager::install("ComplexHeatmap")  # used v2.20.0

install.packages("gprofiler2")          # used v0.2.3
# GO VERSION "e112_eg59_p19_25aa4782"


