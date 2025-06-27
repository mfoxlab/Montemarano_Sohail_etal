# Creating data tables for circos plots
# Hajra Sohail
# 2025-06-26

# PATH DIRECTORY, LOAD FILES ----------------------------------------------

workdir = "/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/Annalisa/New INPUT"
setwd(workdir)
load(file = paste0(workdir,"/go_INPUT_workspace.Rdata"))


# Load libraries, data setup
library(gprofiler2)
library(dplyr)

INPUT_gene_summary <- read.table(paste0(workdir, "/INPUT_gene_summary.txt"), sep = '\t', header = TRUE, quote = '')
INPUT_sig_modules <- c(24, 22, 16)


'''
Single module





test_module <- INPUT_gene_summary %>%
  filter(module.INPUT == "16") %>%
  dplyr::select(mgi_symbol) %>%
  distinct() %>%
  pull()
  
# Convert to a format g:Profiler recognizes
converted <- gconvert(
  query = test_module,
  organism = "mmusculus",
  target = "ENSG",  # or try "ENSMUSG" if it accepts it
  mthreshold = Inf
)

head(converted)
recognized_ids <- unique(converted$target)


# Run GO
gost_result  <- gost(
  query = recognized_ids,
  organism = "mmusculus",   # mouse genome
  sources = c("GO:BP", "GO:MF", "GO:CC"),
  correction_method = "fdr", # similar to BH correction
  significant = TRUE,       # only return significant results
  evcodes = TRUE,
  user_threshold = 0.05     # adjust if needed
)



gost_result_clean <- gost_result$result



gost_result_clean <- gost_result$result %>%
  mutate(Genes = sapply(intersection, paste, collapse = "|"))

'''





# 1. RUNNING BP, CC, MF GENE ONTOLOGY WITH GOST ----------------------------------------------

INPUT_each_go_result <- list()

for (module in INPUT_sig_modules) {
  
  current_module <- INPUT_gene_summary %>%
    filter(module.INPUT == module) %>%
    dplyr::select(mgi_symbol) %>% 
    distinct() %>%
    pull()
  
  # Run GO
  gost_result  <- gost(
    query = current_module,
    organism = "mmusculus",   
    sources = c("GO:BP", "GO:MF", "GO:CC"),
    correction_method = "fdr", 
    significant = TRUE, 
    evcodes = TRUE,
    user_threshold = 0.05     
  )
  
  # Adding data to the table gost_result_df
  if (!is.null(gost_result) && !is.null(gost_result$result)) {
    gost_result_df <- gost_result$result
    gost_result_df$module <- module  # Add module info
    INPUT_each_go_result[[as.character(module)]] <- gost_result_df
    
  } else {
    message(paste("Module", module, "outputted a NULL result")
    )
    
  }
  # Compile all the results together 
  INPUT_go_all_modules <- bind_rows(INPUT_each_go_result)
  
}


INPUT_go_all_modules_cleaned <- INPUT_go_all_modules %>%
  select(-query) %>%
  select(module, everything()) %>%
  mutate(across(where(is.list), ~ sapply(., paste, collapse = "|")))




# VISUALIZING RESULTS
INPUT_go_overview <- INPUT_go_all_modules_cleaned %>%
  group_by(module) %>%
  count(source)

library(ggplot2)
ggplot(INPUT_go_overview, aes(x = factor(module), y = n, fill = source)) +
  geom_bar(stat = "identity", 
           position = position_dodge(preserve = "single"), 
           width = 0.7) +
  labs(title = "Number of GO Terms by Module and Ontology: INPUT",
       x = "Module",
       y = "Count of Enriched Terms",
       fill = "GO Category") +
  theme_minimal(base_size = 14)


# GRAPH 2.
# Optional: select top 5 terms per module
top_terms <- INPUT_go_all_modules_cleaned %>%
  group_by(module) %>%
  slice_min(order_by = p_value, n = 5) %>%
  ungroup()

theme_minimal(base_size = 12)

top_terms <- INPUT_go_all_modules %>%
  group_by(module) %>%
  slice_min(order_by = p_value, n = 5) %>%
  ungroup()

ggplot(top_terms, aes(x = as.factor(module), 
                      y = reorder(term_name, p_value))) +
  geom_point(aes(size = intersection_size, color = -log10(p_value))) +
  scale_color_viridis_c() +
  labs(title = "GO Enrichment Dot Plot: INPUT",
       x = "Module",
       y = "GO Term",
       size = "# Genes in Term",
       color = "-log10(p-value)") +
  theme_minimal(base_size = 12)


write.table(INPUT_go_all_modules_cleaned, file = "INPUT_go_masterlist_concat.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)

# 2. RUNNING KEGG WITH GOST ----------------------------------------------

INPUT_each_kegg_result <- list()

'''
# Single module 

test_module <- INPUT_gene_summary %>%
    filter(module.INPUT == 16) %>%
    dplyr::select(mgi_symbol) %>% 
    distinct() %>%
    pull()


# Perform KEGG enrichment
results <- gost(
  query = test_module,
  organism = "mmusculus",        
  sources = c("KEGG")           
)

testresults <- results$result

'''


for (module in INPUT_sig_modules) {
  
  current_module <- INPUT_gene_summary %>%
    filter(module.INPUT == module) %>%
    dplyr::select(mgi_symbol) %>% 
    distinct() %>%
    pull()
  
  # Perform KEGG enrichment
  kegg_result <- gost(
    query = current_module,
    organism = "mmusculus",        
    sources = c("KEGG"),
    evcodes = TRUE
  )
  
  # Adding data to the table kegg_result_df
  if (!is.null(kegg_result) && !is.null(kegg_result$result)) {
    kegg_result_df <- kegg_result$result
    kegg_result_df$module <- module  # Add module info
    INPUT_each_kegg_result[[as.character(module)]] <- kegg_result_df
    
  } else {
    message(paste("Module", module, "outputted a NULL result")
    )
    
  }
  # Compile all the results together 
  INPUT_kegg_all_modules <- bind_rows(INPUT_each_kegg_result)
  
}


INPUT_kegg_all_modules_cleaned <- INPUT_kegg_all_modules %>%
  select(-query) %>%
  select(module, everything()) %>%
  mutate(across(where(is.list), ~ sapply(., paste, collapse = "|")))


# VISUALIZING RESULTS
INPUT_kegg_overview <- INPUT_kegg_all_modules_cleaned %>%
  count(term_name)



write.table(INPUT_kegg_all_modules_cleaned, file = "INPUT_kegg_masterlist_concat.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)


# Save the workspace --------------------------------------------------------------------------------
save.image(file = paste0(workdir, "/go_INPUT_workspace.Rdata"))
