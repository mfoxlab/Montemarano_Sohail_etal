# Creating data tables for circos plots
# Hajra Sohail
# 2025-06-13
#2025-10-20 meg update to include sig modules in single sexes

# PATH DIRECTORY, LOAD FILES ----------------------------------------------
workdir = ""
workdirPFC = ""
setwd(workdir)

#load the workspace
load(file = paste0(workdir,"/circosplotPFCdataworkspace.Rdata"))

library(dplyr)

#load each of the three data categories
PFC_sexcomb_moduleEigengenesDE <- read.table(paste0(workdirPFC, "/PFC_sexcomb_moduleEigengenesDE.txt"), sep = '\t', header = TRUE, quote = '')
PFC_M_moduleEgigngenesDE <- read.table(paste0(workdirPFC, "/PFC_M_moduleEigengenesDE.txt"), sep = '\t', header = TRUE, quote = '')
PFC_F_moduleEgigngenesDE <- read.table(paste0(workdirPFC, "/PFC_F_moduleEigengenesDE.txt"), sep = '\t', header = TRUE, quote = '')


# CREATE A DATA TABLE FOR ALL MODULES IN THE PFC ----------------------------------------------

all_modules <- PFC_sexcomb_moduleEigengenesDE %>%
  dplyr::select('ModuleEigengene', logFC, P.Value)%>%
  dplyr::rename(logFC_PFC_FENT = logFC, P.Value_PFC_FENT = P.Value)


#prep female data
female_data <- PFC_F_moduleEgigngenesDE %>%
  dplyr::select('ModuleEigengene', logFC, P.Value)%>%
  dplyr::rename(logFC_PFC_F_FENT =logFC,  P.Value_PFC_F_FENT = P.Value)



#prep male data
male_data <- PFC_M_moduleEgigngenesDE %>%
  dplyr::select('ModuleEigengene', logFC, P.Value)%>%
  dplyr::rename(logFC_PFC_M_FENT = logFC, P.Value_PFC_M_FENT = P.Value)



  #join all
  all_modules <- all_modules %>%
    full_join(female_data, by = "ModuleEigengene") %>%
    full_join(male_data, by = "ModuleEigengene")


write.table(all_modules, file = "PFC_circosdata_all_MEs.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)

# CREATE A DATA TABLE FOR ALL SIGNIFICANT MODULES IN THE PFC ----------------------------------------------
significant_modules <- all_modules %>%
  filter((P.Value_PFC_FENT < 0.05) |
           (P.Value_PFC_F_FENT < 0.05) |
           (P.Value_PFC_M_FENT < 0.05)) %>%
  filter(ModuleEigengene != "ME0")

significant_modules <- significant_modules %>%
  dplyr::select(ModuleEigengene,
                logFC_PFC_FENT, logFC_PFC_F_FENT, logFC_PFC_M_FENT,
                P.Value_PFC_FENT, P.Value_PFC_F_FENT, P.Value_PFC_M_FENT)




write.table(significant_modules, file = "PFC_circosdata.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)


# Save the workspace
save.image(file = paste0(workdir, "/circosplotPFCdataworkspace.Rdata"))


