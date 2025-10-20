# Creating data tables for circos plots
# Hajra Sohail
# 2025-06-13
# 2025-10-20-meg edit to select also sig modules regardless of sex

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
  dplyr::select('ModuleEigengene', logFC, P.Value) %>%
  dplyr::rename(logFC_NAC_FENT = logFC, P.Value_NAC_FENT = P.Value)


#prep female data
female_data <- NAC_F_moduleEgigngenesDE %>%
  dplyr::select('ModuleEigengene', logFC, P.Value)%>%
  dplyr::rename(logFC_NAC_F_FENT = logFC, P.Value_NAC_F_FENT = P.Value)



#prep male data
male_data <- NAC_M_moduleEgigngenesDE %>%
  dplyr::select('ModuleEigengene', logFC, P.Value)%>%
  dplyr::rename(logFC_NAC_M_FENT = logFC, P.Value_NAC_M_FENT = P.Value)

#join all
all_modules <- all_modules %>%
  full_join(female_data, by = "ModuleEigengene") %>%
  full_join(male_data, by = "ModuleEigengene")

write.table(all_modules, file = "NAC_circosdata_all_MEs.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)

# CREATE A DATA TABLE FOR ALL SIGNIFICANT MODULES IN THE NAC ----------------------------------------------

significant_modules <- all_modules %>%
  filter((P.Value_NAC_FENT < 0.05) |
           (P.Value_NAC_F_FENT < 0.05) |
           (P.Value_NAC_M_FENT < 0.05)) %>%
  filter(ModuleEigengene != "ME0")

significant_modules <- significant_modules %>%
  dplyr::select(ModuleEigengene,
                logFC_NAC_FENT, logFC_NAC_F_FENT, logFC_NAC_M_FENT,
                P.Value_NAC_FENT, P.Value_NAC_F_FENT, P.Value_NAC_M_FENT)

# Export significant modules
write.table(significant_modules, 
            file = "NAC_circosdata.txt",
            row.names = FALSE, sep = "\t", col.names = TRUE, quote = FALSE)

write.table(significant_modules, file = "NAC_circosdata.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# Save the workspace
save.image(file = paste0(workdir, "/circosplotNACdataworkspace.Rdata"))


