# Build WGCNA gene networks using the data with only DEGs as function of drug or projection target
# Hajra Sohail
# 2025-06-13


# SETUP ------------------------------------------------------------------------------------------
workdir = ""
workdirINPUT = ""
setwd(workdirINPUT)

load(file = paste0(workdirINPUT, "/logCPM_INPUT_workspace.RData"))

library(WGCNA)
library(limma)
library(dplyr)
library(biomaRt) 
library(pheatmap) 


allowWGCNAThreads()

# FINDING THE CORRECT SOFT THRESHOLD ------------------------------------------------------------------------------------------

#define the power range 
powers <- c(1:30, seq(from = 35, to = 50, by = 5)) 

#transpose the data so that genes and samples are correctly arrangedd
logcpm_INPUT_sig<- t(logcpm_INPUT_sig)

#perform the soft-thresholding analysis for VTA->INPUT
powerTable <- pickSoftThreshold(logcpm_INPUT_sig, powerVector = powers, verbose = 5)

#plot the results for this set
par(mfrow = c(1, 2)) 
plot(powerTable$fitIndices[, 1], -sign(powerTable$fitIndices[, 3]) * powerTable$fitIndices[, 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", 
     main = paste("Scale Independence for INPUT"))
text(powerTable$fitIndices[, 1], -sign(powerTable$fitIndices[, 3]) * powerTable$fitIndices[, 2], labels = powers, col = "red") 
plot(powerTable$fitIndices[, 1], powerTable$fitIndices[, 5], 
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = paste("Mean Connectivity for INPUT")) 
text(powerTable$fitIndices[, 1], powerTable$fitIndices[, 5], labels = powers, col = "red")


#find the modules using blockwise
powerthres = 8
netINPUT_8 = blockwiseModules(logcpm_INPUT_sig, corType = "bicor", maxPOutliers = 0.1, 
                            power = powerthres, networkType = "signed", 
                            minModuleSize = 30, reassignThreshold = 0,
                            mergeCutHeight = 0.25, numericLabels = TRUE, minMEtoStay = 0, pamRespectsDendro = FALSE, 
                            saveTOMs = TRUE, saveTOMFileBase = "8TOM.INPUT",
                            loadTOM = T, verbose = 5)



save(netINPUT_8, file = "netINPUT_8.Rdata")


# CREATE MODULES AND GENE SUMMARIES ------------------------------------------------------------------------------------------

#relate modules to traits
sample_names <- rownames(logcpm_INPUT_sig)
traits <- ifelse(grepl("FENT", sample_names), 1, 0)
trait_data <- data.frame(Sample = sample_names, Trait = traits)
rownames(trait_data) <- sample_names

#define module colors
INPUT_module_colors<- netINPUT_8$colors
#look at expression as function of fentanyl treatment
design <- model.matrix(~ trait_data$Trait)
expr_data<- t(logcpm_INPUT_sig)
results <- list()
module_list <- unique(INPUT_module_colors) 
results <- list()

for(module in module_list) { 
  module_genes <- names(INPUT_module_colors)[INPUT_module_colors == module] 
  module_expr_data <- expr_data[rownames(expr_data) %in% module_genes, ] 
  gene_index <- which(rownames(expr_data) %in% module_genes) 
  roast_result <- roast(expr_data, index = gene_index, design = design)
  results[[as.character(module)]] <- roast_result
}



#identify MEs and make adjacency heatmaps

#identify module eigengenes for INPUT
MEsINPUT = moduleEigengenes(
  logcpm_INPUT_sig,
  INPUT_module_colors,
  impute = TRUE,
  nPC = 1,
  align = "along average",
  excludeGrey = FALSE,
  grey = if (is.numeric(INPUT_module_colors))
    0
  else
    "grey",
  subHubs = TRUE,
  softPower = 8,
  scale = TRUE,
  verbose = 5,
  indent = 0
)

save(MEsINPUT, file = "MEsINPUT.RData")

eigengenes <- MEsINPUT$eigengenes
adjacency_matrix <- cor(eigengenes, use = "pairwise.complete.obs")


#plot the adjacency matrix using pheatmap
pdf("input_adjacency_matrix_heatmap.pdf")
pheatmap(adjacency_matrix,
         main = "Adjacency Matrix of Module Eigengenes_INPUT",
         color = colorRampPalette(c("blue", "white", "red"))(50), # Customize the color palette
         cluster_rows = TRUE, # Cluster the rows
         cluster_cols = TRUE, # Cluster the columns
         display_numbers = FALSE)

dev.off()


#save/ export the file with module membership for all the genes


#load in the DE tables to append to the gene summary 
DEsINPUT <- read.table(paste0(workdir, "/DEgene_summary_INPUT.txt"), sep = "\t", header=TRUE)


geneModuleMembership = list()
MEsINPUTAve = MEsINPUT$averageExpr
nSamples = nrow(logcpm_INPUT_sig)
ModuleMembership = as.data.frame(bicor(logcpm_INPUT_sig, MEsINPUTAve, use = "p"))
MMPValue = as.data.frame(corPvalueStudent(as.matrix(ModuleMembership), nSamples))
names(ModuleMembership) = gsub("AE", paste0("INPUT", ".MM"), names(ModuleMembership))
names(MMPValue) = gsub("AE", "p.MM", names(ModuleMembership))
geneModuleMembership = list(membership=as.matrix(ModuleMembership), pvals=as.matrix(MMPValue))

save(geneModuleMembership, file ="INPUT_gene_module_membership.Rdata") 

netGenes = data.frame(ensembl_gene_id = colnames(logcpm_INPUT_sig), module.INPUT = netINPUT_8$colors)


#convert membership matrix to data frame for easy joining 
membership_df <- as.data.frame(geneModuleMembership$membership, stringsAsFactors = FALSE)
membership_df$ensembl_gene_id <- rownames(membership_df)


outp <- DEsINPUT %>%
  full_join(netGenes, by = "ensembl_gene_id") %>%
  full_join(membership_df, by = "ensembl_gene_id") 


#export the data 
write.table(outp, "INPUT_gene_summary.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#find hub genes 
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
topHubGenesList <- list()
moduleNames <- colnames(ModuleMembership)

for (module in moduleNames) {
  moduleMembershipValues <- ModuleMembership[[module]]
  
  rankedGenes <- order(moduleMembershipValues, decreasing = TRUE)
  
  
  top10Genes <- rankedGenes[1:10]
  
  
  top10GeneIDs <- rownames(ModuleMembership)[top10Genes]
  top10kME <- moduleMembershipValues[top10Genes]
  
  geneModuleMembership<-as.data.frame(geneModuleMembership)
  annotations <- getBM(
    attributes = c("ensembl_gene_id", "mgi_symbol"),
    filters = "ensembl_gene_id",
    values = rownames(geneModuleMembership),
    mart = ensembl
  )
  
  
  topHubGenesDF <- data.frame(
    ensembl_gene_id = top10GeneIDs,
    kME = top10kME
  )
  topHubGenesDF <- merge(topHubGenesDF, annotations, by = "ensembl_gene_id", all.x = TRUE)
  
  topHubGenesList[[module]] <- topHubGenesDF
}
write.table(topHubGenesList, "tophubgeneslist_input.txt", sep="\t", row.names=FALSE)




###
#alternative strategy to determine fentanyl associated modules 


combined_df <- cbind(MEsINPUT$eigengenes, Condition=trait_data$Trait)
combined_df$Condition <- factor(combined_df$Condition, levels = c(0, 1), labels = c('SAL', 'FENT'))
kME <- signedKME(logcpm_INPUT_sig, MEsINPUT$eigengenes)
hub_genes <- apply(kME, 2, function(x) names(x)[which.max(x)])
names(hub_genes) <- paste0("ME", sub("kME", "", names(hub_genes)))


#create design matrix

design <- model.matrix(~ Condition, data = combined_df)
eigengenes_df <- t(MEsINPUT$eigengenes)

fit <- lmFit(eigengenes_df, design)
fit <- eBayes(fit)

# Get the results table
results <- topTable(fit, coef="ConditionFENT", number=nrow(eigengenes_df)) 
results$ModuleEigengene <- rownames(results)
results$ensembl_gene_id <- sapply(rownames(results), function(x) hub_genes[x]) 

#annotate with mgi_symbols
results_annotated <- merge(results, annotations, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = TRUE)


write.table(results_annotated, "INPUTmoduleEigengenesDE.txt", sep = "\t", row.names = FALSE, quote = FALSE)
