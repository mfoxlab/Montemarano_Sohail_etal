# Gene Ontology Analysis
# Hajra Sohail
# 2025-06-26

# PATH DIRECTORY, LOAD FILES ----------------------------------------------
workdir = ""
setwd(workdir)
load(file = paste0(workdir,"/go_NAC_workspace.Rdata"))


# Load libraries, data setup
library(gprofiler2)
library(dplyr)
library(ggplot2)


NAC_gene_summary <- read.table(paste0(workdir, "/NAC_gene_summary.txt"), sep = '\t', header = TRUE, quote = '')
NAC_sig_modules <- c(39, 38, 57, 64, 18, 52)

# Import the colors.txt to get color keys for the graphs
colordir = "/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/Annalisa/New WGCNA"
color_data <- read.delim(file.path(colordir, "modulecolors.txt"), header = TRUE, stringsAsFactors = FALSE, quote = '', fill = TRUE)

color_key <- color_data %>%
  rename_at('Module', ~'module') 

color_key$module <- gsub("ME","", as.character(color_key$module))

color_key <- color_key %>%
  mutate(module = as.numeric(module))

# 1. RUNNING BP, CC, MF GENE ONTOLOGY WITH GOST ----------------------------------------------

'''
Single module





test_module <- NAC_gene_summary %>%
  filter(module.NAC == "16") %>%
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



NAC_each_go_result <- list()

for (module in NAC_sig_modules) {
  
  current_module <- NAC_gene_summary %>%
    filter(module.NAC == module) %>%
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
    NAC_each_go_result[[as.character(module)]] <- gost_result_df
    
  } else {
    message(paste("Module", module, "outputted a NULL result")
    )
    
  }
  # Compile all the results together 
  NAC_go_all_modules <- bind_rows(NAC_each_go_result)
  
}

NAC_go_all_modules_cleaned <- NAC_go_all_modules %>%
  select(-query) %>%
  select(module, everything(), -significant, significant) %>%
  mutate(across(where(is.list), ~ sapply(., paste, collapse = "|")))

NAC_go_all_modules_cleaned_colorkey <- inner_join(NAC_go_all_modules_cleaned, color_key, by = "module") %>%
  select(module, Hex.Color, Color.Name, everything())

# GRAPH 1: VISUALIZING OVERALL RESULTS ----------------------------------------------
'''
NAC_go_overview <- NAC_go_all_modules_cleaned %>%
  group_by(module) %>%
  count(source)

ggplot(NAC_go_overview, aes(x = factor(module), y = n, fill = source)) +
  geom_bar(stat = "identity", 
           position = position_dodge(preserve = "single"), 
           width = 0.7) +
  labs(title = "Number of GO Terms by Module and Ontology: NAC",
       x = "Module",
       y = "Count of Enriched Terms",
       fill = "GO Category") +
  theme_minimal(base_size = 14)
'''

# GRAPH 2: ENRICHMENT DOT PLOT, NUMBERED MODULE LABELS ----------------------------------------------
# Select top 3 terms per module

top_terms <- NAC_go_all_modules_cleaned_colorkey %>%
  group_by(module) %>%
  #filter(intersection_size >= 5) %>%
  filter(source == "GO:BP") %>%
  slice_min(order_by = p_value, n = 3) %>%
  ungroup()

theme_minimal(base_size = 12)
#top_terms$intersection_size_capped <- pmin(top_terms$intersection_size, 100)

'''
# Numbered modules version

ggplot(top_terms, aes(x = as.factor(module), 
                      y = reorder(term_name, -module))) +
  geom_point(aes(size = intersection_size, color = -log10(p_value))) +
  scale_color_viridis_c() +
  scale_size_continuous(
    range = c(5, 10),  # adjust dot size appearance; tweak as needed
    breaks = c(2, 4),
    labels = c("2", "4")
  ) +
  labs(title = "GO:BP Enrichment Dot Plot: NAC",
       x = "Module",
       y = "GO Term",
       size = "# Genes in Term",
       color = "-log10(p-value)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

'''
# GRAPH 3: ENRICHMENT DOT PLOT, COLORED MODULE LABELS ----------------------------------------------
# Ensure x-axis factor levels are consistent
top_terms$Color.Name <- factor(top_terms$Color.Name, levels = unique(top_terms$Color.Name))

# Create a module label dataset (1 row per module)
module_labels <- top_terms %>%
  distinct(Color.Name, Hex.Color) %>%
  mutate(y = -1)  # Position below the lowest term

# Plot
NAC_go_plot_colorkey <- ggplot(top_terms, aes(x = Color.Name, y = reorder(term_name, -module))) +
  
  # Colored rectangles under x-axis
  geom_tile(data = module_labels,
            aes(x = Color.Name, y = y, fill = Hex.Color),
            width = 0.9, height = 0.5,
            inherit.aes = FALSE) +
  
  # Dot plot
  geom_point(aes(size = intersection_size, color = -log10(p_value))) +
  
  # Scales
  scale_color_viridis_c() +
  scale_fill_identity() +
  scale_size_continuous(
    range = c(5, 10),
    breaks = c(2, 4),
    labels = c("2", "4")
  ) +
  
  # Axes & labels
  labs(title = "GO:BP Enrichment Dot Plot: NAC",
       x = "Module",
       y = "GO Term",
       size = "# Genes in Term",
       color = "-log10(p-value)") +
  
  # Theme
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)#,
   # plot.margin = margin(b = 30, t = 30, r = 10, l = 10)  
  ) +
  
  # Expand y-axis to make room for bottom color tiles
  scale_y_discrete(expand = expansion(add = 1.5))

NAC_go_plot_colorkey

# Saving and exporting
pdf("NAC GO BP Enrichment Dot Plot Colored.pdf", width = 11, height = 8.5)
print(NAC_go_plot_colorkey)
dev.off()

write.table(NAC_go_all_modules_cleaned_colorkey, file = "NAC_go_masterlist_concat.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)

# GRAPH 4: ENRICHMENT DOT PLOT, MANUALLY SELECTED GO TERMS WITH COLORED LABELS ----------------------------------------------
# Import the file from the same directory as the color key

manual_go_terms <- read.delim(file.path(colordir, "manual_go_terms.txt"), header = TRUE, stringsAsFactors = FALSE, quote = '', fill = TRUE)
manual_go_terms <- manual_go_terms %>%
  rename_at('GO_ID', ~'term_id')

# inner_join by term_id
NAC_go_manual_terms <- inner_join(NAC_go_all_modules_cleaned_colorkey, manual_go_terms, by = "term_id")     # Shows no results!

# Visualize
top_terms_manual <- NAC_go_manual_terms %>%
  group_by(module) %>%
  #filter(intersection_size >= 5) %>%
  #filter(source.y == "GO:BP") %>%
  #slice_min(order_by = p_value, n = 3) %>%           # No filtering for this graph!
  ungroup()



theme_minimal(base_size = 12)
top_terms_manual$intersection_size_capped <- pmin(top_terms_manual$intersection_size, 75)

# Ensure x-axis factor levels are consistent
top_terms_manual$Color.Name <- factor(top_terms_manual$Color.Name, levels = unique(top_terms_manual$Color.Name))

# Create a module label dataset (1 row per module)
module_labels_manual <- top_terms_manual %>%
  distinct(Color.Name, Hex.Color) %>%
  mutate(y = -1)  # Position below the lowest term

# Plot
NAC_go_plot_colorkey_manual <- ggplot(top_terms_manual, aes(x = Color.Name, y = reorder(term_name.y, -module))) +
  
  # Colored rectangles under x-axis
  geom_tile(data = module_labels_manual,
            aes(x = Color.Name, y = y, fill = Hex.Color),
            width = 0.9, height = 0.5,
            inherit.aes = FALSE) +
  
  # Dot plot
  geom_point(aes(size = intersection_size_capped, color = -log10(p_value))) +
  
  # Scales
  scale_color_viridis_c() +
  scale_fill_identity() +
  scale_size_continuous(
    range = c(2, 8),
    breaks = c(3, 25, 50, 75),
    labels = c("3", "25", "50", "75+")
  ) +
  
  # Axes & labels
  labs(title = "GO Enrichment Dot Plot: NAC (Manually Selected Terms)",
       x = "Module",
       y = "GO Term",
       size = "# Genes in Term",
       color = "-log10(p-value)") +
  
  # Theme
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(b = 30, t = 30, r = 10, l = 10)  
  ) +
  
  # Expand y-axis to make room for bottom color tiles
  scale_y_discrete(expand = expansion(add = 1.5))

NAC_go_plot_colorkey_manual

# Saving and exporting
pdf("NAC GO Enrichment Dot Plot Colored Manually Selected Terms.pdf", width = 11, height = 8.5)
print(NAC_go_plot_colorkey_manual)
dev.off()


# 2. RUNNING KEGG WITH GOST ----------------------------------------------

NAC_each_kegg_result <- list()

'''
# Single module 

test_module <- NAC_gene_summary %>%
    filter(module.NAC == 11) %>%
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


for (module in NAC_sig_modules) {
  
  current_module <- NAC_gene_summary %>%
    filter(module.NAC == module) %>%
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
    NAC_each_kegg_result[[as.character(module)]] <- kegg_result_df
    
  } else {
    message(paste("Module", module, "outputted a NULL result")
    )
    
  }
  # Compile all the results together 
  NAC_kegg_all_modules <- bind_rows(NAC_each_kegg_result)
  
}


NAC_kegg_all_modules_cleaned <- NAC_kegg_all_modules %>%
  select(-query) %>%
  select(module, everything()) %>%
  mutate(across(where(is.list), ~ sapply(., paste, collapse = "|")))


# VISUALIZING RESULTS
NAC_kegg_overview <- NAC_kegg_all_modules_cleaned %>%
  count(term_name)




write.table(NAC_kegg_all_modules_cleaned, file = "NAC_kegg_masterlist_concat.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)


# Save the workspace --------------------------------------------------------------------------------
save.image(file = paste0(workdir, "/go_NAC_workspace.Rdata"))
