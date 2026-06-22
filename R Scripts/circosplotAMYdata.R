# Creating data tables for circos plots
# Meghan Hnilo copied 06-22-26

# PATH DIRECTORY, LOAD FILES ----------------------------------------------
workdir = ""
workdirAMY = ""
setwd(workdir)

#load the workspace
load(file = paste0(workdir,"/circosplotAMYdataworkspace.Rdata"))

library(dplyr)

#load each of the three data categories
AMY_sexcomb_moduleEigengenesDE <- read.table(paste0(workdirAMY, "/AMY_sexcomb_moduleEigengenesDE.txt"), sep = '\t', header = TRUE, quote = '')
AMY_M_moduleEgigngenesDE <- read.table(paste0(workdirAMY, "/AMY_M_moduleEigengenesDE.txt"), sep = '\t', header = TRUE, quote = '')
AMY_F_moduleEgigngenesDE <- read.table(paste0(workdirAMY, "/AMY_F_moduleEigengenesDE.txt"), sep = '\t', header = TRUE, quote = '')


# CREATE A DATA TABLE FOR ALL MODULES IN THE AMY ----------------------------------------------

all_modules <- AMY_sexcomb_moduleEigengenesDE %>%
  dplyr::select('ModuleEigengene', logFC, P.Value)

#append the category of data to the columns so you know what it's referring to
colnames(all_modules)[c(2:3)] <- paste(colnames(all_modules)[c(2:3)], 'AMY_FENT', sep = '_')

#prep female data
female_data <- AMY_F_moduleEgigngenesDE %>%
  dplyr::select('ModuleEigengene', logFC, P.Value)
colnames(female_data)[c(2:3)] <- paste(colnames(female_data)[c(2:3)], 'AMY_F_FENT', sep = '_')

#join female data
all_modules <- full_join(all_modules, female_data, by = "ModuleEigengene")

#prep male data
male_data <- AMY_M_moduleEgigngenesDE %>%
  dplyr::select('ModuleEigengene', logFC, P.Value)
colnames(male_data)[c(2:3)] <- paste(colnames(male_data)[c(2:3)], 'AMY_M_FENT', sep = '_')

#join male data (final version)
all_modules <- full_join(all_modules, male_data, by = "ModuleEigengene")

#now, rearrange 
all_modules <- all_modules %>%
  dplyr::select(c(1, 2, 4, 6, 3, 5, 7))


write.table(all_modules, file = "AMY_circosdata_all_MEs.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)

# CREATE A DATA TABLE FOR ALL SIGNIFICANT MODULES IN THE AMY ----------------------------------------------

significant_modules <- AMY_sexcomb_moduleEigengenesDE %>%
  filter(P.Value < 0.05) %>%
  dplyr::select('ModuleEigengene', logFC, P.Value)

#remove the negligible ME0 as it is not really a module and should not be considered significant
significant_modules <- subset(significant_modules, ModuleEigengene != 'ME0') 

#append the category of data to the columns so you know what it's referring to
colnames(significant_modules)[c(2:3)] <- paste(colnames(significant_modules)[c(2:3)], 'AMY_FENT', sep = '_')

#join female data
significant_modules <- inner_join(significant_modules, female_data, by = "ModuleEigengene")

#join male data (final version)
significant_modules <- inner_join(significant_modules, male_data, by = "ModuleEigengene")

#now, rearrange 
significant_modules <- significant_modules %>%
  dplyr::select(c(1, 2, 4, 6, 3, 5, 7))


write.table(significant_modules, file = "AMY_circosdata.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)



# Save the workspace
save.image(file = paste0(workdir, "/circosplotAMYdataworkspace.Rdata"))
