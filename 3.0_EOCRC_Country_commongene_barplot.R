
Gene_list_country_EO <- merge(
  setDT(clin.EOCRC@clinical.data[, .(Country, Panel)]), 
  setDT(assays_genes),   by = "Panel",   allow.cartesian = TRUE) %>% 
  distinct(Country, Hugo_Symbol) %>%  
  group_by(Country) %>% 
  summarize(Gene_list = list(unique(Hugo_Symbol)), .groups = "drop")

gene_lists <- lapply(Gene_list_country_EO$Gene_list, unlist)
common_genes <- Reduce(intersect, gene_lists)

mutation_table <- data.table(
  gene = common_genes,
  Canada = NA_real_,
  China = NA_real_,
  US = NA_real_,
  Netherlands = NA_real_,
  Nigeria = NA_real_,
  France = NA_real_,
  Spain = NA_real_,
  Korea = NA_real_
)

generate_country_frequency <- function(countries, maf,Gene_common){
  
  Alter_Gene = maf@data[Hugo_Symbol == Gene_common,] %>% 
    filter(!duplicated(Tumor_Sample_Barcode))
  
  mut_samples <- left_join(Alter_Gene[, c('Tumor_Sample_Barcode', 'Hugo_Symbol')],
                           maf@clinical.data[, c('Tumor_Sample_Barcode', "Country")],
                           by = 'Tumor_Sample_Barcode')
  
  Gene_panels_samples <- maf@clinical.data %>% 
    filter(Panel %in% gene_panel_map[[Gene_common]])
  
  
  for (guojia in countries) {
    
    mut_country_samples <- mut_samples %>% 
      filter(Country == guojia)  %>% 
      summarise(unique_samples_count = n_distinct(Tumor_Sample_Barcode)) %>% 
      pull(unique_samples_count)
    
    panels_country_samples <- Gene_panels_samples %>% 
      filter(Country == guojia) %>% 
      summarise(unique_samples_count = n_distinct(Tumor_Sample_Barcode)) %>% 
      pull(unique_samples_count)
    
    mutation_table[gene == Gene_common, (guojia) := mut_country_samples / panels_country_samples]
    
  }
  mutation_table  <<- mutation_table
}

for (i in common_genes) {
  generate_country_frequency(countries = unique(clin.EOCRC@clinical.data$Country), 
                             maf = clin.EOCRC, Gene_common = i)
}

mutation_table_EO <- mutation_table %>%
  mutate(Age_class = 'EOCRC') %>% 
  filter(rowSums(across(c("Canada", "China", "US", 
                          "Netherlands", "Nigeria", 
                          "France", "Spain", "Korea"), ~ .x > 0)) >= 8)
write.csv(mutation_table_EO,
          file = paste0(dir,'table/4.0_EOCRC_Country_commongenes_mutation_table_EO.csv'))
