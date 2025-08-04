# Build WGCNA gene networks using the data with only DEGs as function of drug or projection target
# Hajra Sohail
# 2025-06-13


# SETUP ------------------------------------------------------------------------------------------
workdir = "/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/Annalisa/New WGCNA"
workdirPFC = "/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/Annalisa/New PFC"
setwd(workdirPFC)

load(file = paste0(workdirPFC, "/WGCNA_PFC_workspace.RData"))

#load(paste0(workdir, "/logCPM_PFC_sig_Annalisa.Rdata"))

library(WGCNA)
library(limma)
library(dplyr)
library(biomaRt) 
library(pheatmap) 


allowWGCNAThreads()

# FINDING THE CORRECT SOFT THRESHOLD ------------------------------------------------------------------------------------------
# MAC USERS: Please run the two lines below. Skip the rest of the code in this section. Move onto 'CREATE MODULES AND GENE SUMMARIES'.
#load(paste0(workdirPFC, "/netPFC_8.Rdata"))
#logcpm_PFC_sig <-t (logcpm_PFC_sig) 

#define the power range 
powers <- c(1:30, seq(from = 35, to = 50, by = 5)) 

#transpose the data so that genes and samples are correctly arranged
logcpm_PFC_sig<- t(logcpm_PFC_sig)

#perform the soft-thresholding analysis for VTA->PFC
powerTable <- pickSoftThreshold(logcpm_PFC_sig, powerVector = powers, verbose = 5)

#plot the results for this set
par(mfrow = c(1, 2)) 
plot(powerTable$fitIndices[, 1], -sign(powerTable$fitIndices[, 3]) * powerTable$fitIndices[, 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", 
     main = paste("Scale Independence for PFC"))
text(powerTable$fitIndices[, 1], -sign(powerTable$fitIndices[, 3]) * powerTable$fitIndices[, 2], labels = powers, col = "red") 
plot(powerTable$fitIndices[, 1], powerTable$fitIndices[, 5], 
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = paste("Mean Connectivity for PFC")) 
text(powerTable$fitIndices[, 1], powerTable$fitIndices[, 5], labels = powers, col = "red")


#find the modules using blockwise
powerthres = 8
netPFC_8 = blockwiseModules(logcpm_PFC_sig, corType = "bicor", maxPOutliers = 0.1, 
                            power = powerthres, networkType = "signed", 
                            minModuleSize = 30, reassignThreshold = 0,
                            mergeCutHeight = 0.25, numericLabels = TRUE, minMEtoStay = 0, pamRespectsDendro = FALSE, 
                            saveTOMs = TRUE, saveTOMFileBase = "8TOM.PFC",
                            loadTOM = T, verbose = 5)



save(netPFC_8, file = "netPFC_8.Rdata")


# CREATE MODULES AND GENE SUMMARIES ------------------------------------------------------------------------------------------

#relate modules to traits
sample_names <- rownames(logcpm_PFC_sig)
traits <- ifelse(grepl("FENT", sample_names), 1, 0)
trait_data <- data.frame(Sample = sample_names, Trait = traits)
rownames(trait_data) <- sample_names

#define module colors
PFC_module_colors<- netPFC_8$colors
#look at expression as function of fentanyl treatment
design <- model.matrix(~ trait_data$Trait)
expr_data<- t(logcpm_PFC_sig)
results <- list()
module_list <- unique(PFC_module_colors) 
results <- list()

for(module in module_list) { 
  module_genes <- names(PFC_module_colors)[PFC_module_colors == module] 
  module_expr_data <- expr_data[rownames(expr_data) %in% module_genes, ] 
  gene_index <- which(rownames(expr_data) %in% module_genes) 
  roast_result <- roast(expr_data, index = gene_index, design = design)
  results[[as.character(module)]] <- roast_result
}
#add sex metadata for later circos comparisons
sexes<-ifelse(grepl("_M", sample_names),1,0)
sex_data<-data.frame(Sample=sample_names, Sex=sexes)
rownames(sex_data)<-sample_names


#identify MEs and make adjacency heatmaps

#identify module eigengenes for PFC
MEsPFC = moduleEigengenes(
  logcpm_PFC_sig,
  PFC_module_colors,
  impute = TRUE,
  nPC = 1,
  align = "along average",
  excludeGrey = FALSE,
  grey = if (is.numeric(PFC_module_colors))
    0
  else
    "grey",
  subHubs = TRUE,
  softPower = 8,
  scale = TRUE,
  verbose = 5,
  indent = 0
)

save(MEsPFC, file = "MEsPFC.RData")

eigengenes <- MEsPFC$eigengenes
adjacency_matrix <- cor(eigengenes, use = "pairwise.complete.obs")


#plot the adjacency matrix using pheatmap
pdf("pfc_adjacency_matrix_heatmap.pdf")
pheatmap(adjacency_matrix,
         main = "Adjacency Matrix of Module Eigengenes_PFC",
         color = colorRampPalette(c("blue", "white", "red"))(50), # Customize the color palette
         cluster_rows = TRUE, # Cluster the rows
         cluster_cols = TRUE, # Cluster the columns
         display_numbers = FALSE)

dev.off()


#save/ export the file with module membership for all the genes


#load in the DE tables to append to the gene summary 
DEsPFC <- read.table(paste0(workdir, "/DEgene_summary_PFC.txt"), sep = "\t", header=TRUE)


geneModuleMembership = list()
MEsPFCAve = MEsPFC$averageExpr
nSamples = nrow(logcpm_PFC_sig)
ModuleMembership = as.data.frame(bicor(logcpm_PFC_sig, MEsPFCAve, use = "p"))
MMPValue = as.data.frame(corPvalueStudent(as.matrix(ModuleMembership), nSamples))
names(ModuleMembership) = gsub("AE", paste0("PFC", ".MM"), names(ModuleMembership))
names(MMPValue) = gsub("AE", "p.MM", names(ModuleMembership))
geneModuleMembership = list(membership=as.matrix(ModuleMembership), pvals=as.matrix(MMPValue))

save(geneModuleMembership, file ="PFC_gene_module_membership.Rdata") 

netGenes = data.frame(ensembl_gene_id = colnames(logcpm_PFC_sig), module.PFC = netPFC_8$colors)


#convert membership matrix to data frame for easy joining 
membership_df <- as.data.frame(geneModuleMembership$membership, stringsAsFactors = FALSE)
membership_df$ensembl_gene_id <- rownames(membership_df)


outp <- DEsPFC %>%
  full_join(netGenes, by = "ensembl_gene_id") %>%
  full_join(membership_df, by = "ensembl_gene_id") 


#export the data 
write.table(outp, "PFC_gene_summary.txt", sep = "\t", row.names = FALSE, quote = FALSE)


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
write.table(topHubGenesList, "tophubgeneslist_pfc.txt", sep="\t", row.names=FALSE)




###
#alternative strategy to determine fentanyl associated modules 


combined_df <- cbind(MEsPFC$eigengenes, Condition=trait_data$Trait)
combined_df$Condition <- factor(combined_df$Condition, levels = c(0, 1), labels = c('SAL', 'FENT'))
kME <- signedKME(logcpm_PFC_sig, MEsPFC$eigengenes)
hub_genes <- apply(kME, 2, function(x) names(x)[which.max(x)])
names(hub_genes) <- paste0("ME", sub("kME", "", names(hub_genes)))


#sex related 
sex_combined_df <- cbind(MEsPFC$eigengenes, Sex = sex_data$Sex)
sex_combined_df$Sex <- factor(sex_combined_df$Sex, levels = c(0, 1), labels = c('FEMALE', 'MALE'))
sex_condition_combined <- interaction(sex_combined_df$Sex, combined_df$Condition, drop = TRUE)
sex_condition_combined <- factor(sex_condition_combined, levels = c("FEMALE.SAL", "FEMALE.FENT", "MALE.SAL", "MALE.FENT"))
sex_drug_combined_df <- cbind(MEsPFC$eigengenes, SexCondition = sex_condition_combined)
sex_combined_df$SexCondition <- interaction(sex_combined_df$Sex, combined_df$Condition, drop = TRUE)
sex_combined_df$SexCondition <- factor(sex_combined_df$SexCondition, 
                                       levels = c("FEMALE.SAL", "FEMALE.FENT", "MALE.SAL", "MALE.FENT"))
# Step 4: Subset for female and male data separately
female_data <- sex_combined_df[sex_combined_df$Sex == "FEMALE", ]
male_data <- sex_combined_df[sex_combined_df$Sex == "MALE", ]

design <- model.matrix(~ Condition, data = combined_df)
eigengenes_df <- t(MEsPFC$eigengenes)

fit <- lmFit(eigengenes_df, design)
fit <- eBayes(fit,trend=TRUE)

# Get the results table
results <- topTable(fit, coef="ConditionFENT", number=nrow(eigengenes_df)) 
results$ModuleEigengene <- rownames(results)
results$ensembl_gene_id <- sapply(rownames(results), function(x) hub_genes[x]) 

#annotate with mgi_symbols
results_annotated <- merge(results, annotations, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = TRUE)


female_data$SexCondition <- factor(female_data$SexCondition, 
                                   levels = c("FEMALE.SAL", "FEMALE.FENT"))
male_data$SexCondition <- factor(male_data$SexCondition, 
                                 levels = c("MALE.SAL", "MALE.FENT"))

# For females, relevel the factor so FEMALE.SAL is the reference
female_data$SexCondition <- relevel(female_data$SexCondition, ref = "FEMALE.SAL")

female_design<-model.matrix(~0 +SexCondition, data = female_data)
female_eigengenes <- t(female_data[, 1:ncol(MEsPFC$eigengenes)])
female_lm_fit <- lmFit(female_eigengenes, female_design)

male_design <- model.matrix(~0 +SexCondition, data = male_data)
male_eigengenes <- t(MEsPFC$eigengenes[rownames(male_data), ])
male_lm_fit <- lmFit(male_eigengenes, male_design)


# Define contrasts (FENT vs SAL) for females
female_contrasts <- makeContrasts(
  FEM_FENT_vs_SAL = SexConditionFEMALE.FENT - SexConditionFEMALE.SAL, 
  levels = female_design
)
male_contrasts <- makeContrasts(
  MALE_FENT_vs_SAL = SexConditionMALE.FENT - SexConditionMALE.SAL, 
  levels = male_design
)

female_fit <- contrasts.fit(female_lm_fit, female_contrasts)
female_fit <- eBayes(female_fit,trend=TRUE)
male_fit <- contrasts.fit(male_lm_fit, male_contrasts)
male_fit <- eBayes(male_fit,trend=TRUE)

female_results<-topTable(female_fit, coef="FEM_FENT_vs_SAL",number=nrow(female_eigengenes))
male_results<-topTable(male_fit, coef="MALE_FENT_vs_SAL", number=nrow(male_eigengenes))
female_results$ModuleEigengene<-rownames(female_results)
female_results$ensembl_gene_id<-sapply(rownames(female_results), function(x) hub_genes[x])

male_results$ModuleEigengene<-rownames(male_results)
male_results$ensembl_gene_id<-sapply(rownames(male_results),function(x) hub_genes[x])

write.table(results_annotated, "PFC_sexcomb_moduleEigengenesDE.txt", sep = "\t", row.names = FALSE, quote = FALSE) #module eigengene data for sexes combined
write.table(female_results, "PFC_F_moduleEigengenesDE.txt",sep="\t", row.names=FALSE, quote=FALSE) #module eigengene data for females
write.table(male_results,"PFC_M_moduleEigengenesDE.txt", sep="\t", row.names=FALSE, quote=FALSE) #module eigengene data for males

#now with sexes separated 


save.image(file = "WGCNA_PFC_workspace.RData")

