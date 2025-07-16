# Creating data tables for circos plots
# Hajra Sohail
# 2025-06-13

# PATH DIRECTORY, LOAD FILES ----------------------------------------------

workdir = ""
workdirNAC = ""
setwd(workdir)

#load the workspace
load(file = paste0(workdir,"/circosplotNACdataworkspace.Rdata"))

library(dplyr)

#load each of the three data categories
NAC_sexcomb_moduleEigengenesDE <- read.table(paste0(workdirNAC, "/NAC_sexcomb_moduleEigengenesDE.txt"), sep = '\t', header = TRUE, quote = '')
NAC_M_moduleEgigngenesDE <- read.table(paste0(workdirNAC, "/NAC_M_moduleEigengenesDE.txt"), sep = '\t', header = TRUE, quote = '')
NAC_F_moduleEgigngenesDE <- read.table(paste0(workdirNAC, "/NAC_F_moduleEigengenesDE.txt"), sep = '\t', header = TRUE, quote = '')


# CREATE A DATA TABLE FOR ALL MODULES IN THE NAC ----------------------------------------------

all_modules <- NAC_sexcomb_moduleEigengenesDE %>%
  dplyr::select('ModuleEigengene', logFC, P.Value)

#append the category of data to the columns so you know what it's referring to
colnames(all_modules)[c(2:3)] <- paste(colnames(all_modules)[c(2:3)], 'NAC_FENT', sep = '_')

#prep female data
female_data <- NAC_F_moduleEgigngenesDE %>%
  dplyr::select('ModuleEigengene', logFC, P.Value)
colnames(female_data)[c(2:3)] <- paste(colnames(female_data)[c(2:3)], 'NAC_F_FENT', sep = '_')

#join female data
all_modules <- full_join(all_modules, female_data, by = "ModuleEigengene")

#prep male data
male_data <- NAC_M_moduleEgigngenesDE %>%
  dplyr::select('ModuleEigengene', logFC, P.Value)
colnames(male_data)[c(2:3)] <- paste(colnames(male_data)[c(2:3)], 'NAC_M_FENT', sep = '_')

#join male data (final version)
all_modules <- full_join(all_modules, male_data, by = "ModuleEigengene")

#now, rearrange 
all_modules <- all_modules %>%
  select(c(1, 2, 4, 6, 3, 5, 7))


write.table(all_modules, file = "NAC_circosdata_all_MEs.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)

# CREATE A DATA TABLE FOR ALL SIGNIFICANT MODULES IN THE NAC ----------------------------------------------

significant_modules <- NAC_sexcomb_moduleEigengenesDE %>%
  filter(P.Value < 0.05) %>%
  dplyr::select('ModuleEigengene', logFC, P.Value)
  
#remove the negligible ME0 as it is not really a module and should not be considered significant
significant_modules <- subset(significant_modules, ModuleEigengene != 'ME0') 

#append the category of data to the columns so you know what it's referring to
colnames(significant_modules)[c(2:3)] <- paste(colnames(significant_modules)[c(2:3)], 'NAC_FENT', sep = '_')

#join female data
significant_modules <- inner_join(significant_modules, female_data, by = "ModuleEigengene")

#join male data (final version)
significant_modules <- inner_join(significant_modules, male_data, by = "ModuleEigengene")

#now, rearrange 
significant_modules <- significant_modules %>%
  select(c(1, 2, 4, 6, 3, 5, 7))


write.table(significant_modules, file = "NAC_circosdata.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)



# Save the workspace
save.image(file = paste0(workdir, "/circosplotNACdataworkspace.Rdata"))

