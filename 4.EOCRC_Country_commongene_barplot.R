library(ggsci)
library(ggplot2)
library(tidyverse)
library(treeio)
library(ape)
library(ggnewscale)
library(MetBrewer)
library(tidyr)
library(ggthemes)

countries_EO <- unique(rt@clinical.data$Country)

Gene_list_country_EO <- left_join(clin.EOCRC@data[,.(Hugo_Symbol,Tumor_Sample_Barcode)],
                     clin.EOCRC@clinical.data[,.(Tumor_Sample_Barcode,Country)],
          by = 'Tumor_Sample_Barcode') %>% 
  group_by(Country) %>% 
  reframe(
    Gene_list = list(Hugo_Symbol)
  )

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

generate_country_frequency <- function(countries, maf,jiyin){
 
  Alter_Gene = maf@data[Hugo_Symbol == jiyin,] %>% 
    filter(!duplicated(Tumor_Sample_Barcode))
  
 
  mut_samples <- merge(Alter_Gene[, c('Tumor_Sample_Barcode', 'Hugo_Symbol')],
                                           maf@clinical.data[, c('Tumor_Sample_Barcode', "Country")],
                                           by = 'Tumor_Sample_Barcode')

  Gene_panels <- unique(panel_hugo_symbol[Hugo_Symbol == jiyin]$Panel)
  
  
  for (guojia in countries) {

    mut_country_samples <- mut_samples %>% 
      filter(Country == guojia)  %>% 
      summarise(unique_samples_count = n_distinct(Tumor_Sample_Barcode)) %>% 
      pull(unique_samples_count)

    panels_country_samples <- maf@clinical.data %>% 
      filter(Panel %in% Gene_panels,
             Country == guojia) %>% 
      summarise(unique_samples_count = n_distinct(Tumor_Sample_Barcode)) %>% 
      pull(unique_samples_count)

    mutation_table[gene == jiyin, (guojia) := mut_country_samples / panels_country_samples]
    
  }
  mutation_table  <<- mutation_table
}


for (i in common_genes) {
    generate_country_frequency(countries = countries_EO, 
                             maf = clin.EOCRC, jiyin = i)
}

  mutation_table_EO <- mutation_table %>%
    mutate(Age_class = 'EOCRC')

