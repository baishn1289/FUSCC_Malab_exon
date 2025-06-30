
uniquegenes <- function(m1, m2, m1Name = NULL, m2Name = NULL,
                        minMut = 5){
  
  m1.gs <- getGeneSummary(x = m1)
  m2.gs <- getGeneSummary(x = m2)
  
  m1.genes = as.character(m1.gs[AlteredSamples >= minMut, 
                                Hugo_Symbol])
  m2.genes = as.character(m2.gs[AlteredSamples >= minMut, 
                                Hugo_Symbol])
  uniquegenes =  intersect(m1.genes, m2.genes)
  
  return(uniquegenes)
  
}


logstic_country_mafcompare <- function(m1, m2, Clinical.data, Gene_logistic,Country_logistic,TMB_status) {
  
  Alter_m1_Gene <- as.data.table(m1@data)[Hugo_Symbol == Gene_logistic]
  Alter_m1_Gene <- unique(Alter_m1_Gene, by = "Tumor_Sample_Barcode")
  
  Alter_m2_Gene <- as.data.table(m2@data)[Hugo_Symbol == Gene_logistic]
  Alter_m2_Gene <- unique(Alter_m2_Gene, by = "Tumor_Sample_Barcode")
  
  Gene_panels <- Clinical.data%>% 
    filter(Panel %in% gene_panel_map[[Gene_logistic]]) %>% 
    pull(Panel)
  
  mut_samples_m1 <- length(unique(Alter_m1_Gene$Tumor_Sample_Barcode))
  panels_samples_m1 <- length(unique(m1@clinical.data[m1@clinical.data$Panel %in% Gene_panels, ]$Tumor_Sample_Barcode))
  m1_freq <- mut_samples_m1 / panels_samples_m1
  
  mut_samples_m2 <- length(unique(Alter_m2_Gene$Tumor_Sample_Barcode))
  panels_samples_m2 <- length(unique(m2@clinical.data[m2@clinical.data$Panel %in% Gene_panels, ]$Tumor_Sample_Barcode))
  m2_freq <- mut_samples_m2 / panels_samples_m2
  
  if ( m1_freq < 0.01 | m2_freq<0.01 ) {
    print(paste(Gene_logistic,Country_logistic,TMB_status,
                "Low mutation frequencies, returning from function.",
                sep = '____'))
    return()
  }
  
  imp_list_Gene_clinical_info <-lapply(imp_list, function(df) {
    dt <- as.data.table(df)
    dt[Country == Country_logistic]
    dt[, Gene_status := as.integer(Tumor_Sample_Barcode %in% c(Alter_m1_Gene$Tumor_Sample_Barcode, Alter_m2_Gene$Tumor_Sample_Barcode))]
    dt[, Gene_status := factor(Gene_status, levels = c(0, 1))]
    dt[Panel %in% Gene_panels]
  })
  
  
  each_gene_mutation_counts = gene_mutation_counts %>% 
    dplyr::select(all_of(Gene_logistic), Tumor_Sample_Barcode,Age_class, TMB) %>%
    filter(Tumor_Sample_Barcode %in% imp_list_Gene_clinical_info[[1]]$Tumor_Sample_Barcode)  %>% 
    mutate(mutation_rate = .data[[Gene_logistic]] / TMB)
  
  each_gene_mutation_EO = each_gene_mutation_counts %>% 
    filter(Age_class == 'EOCRC') %>% pull(mutation_rate)
  each_gene_mutation_LO = each_gene_mutation_counts%>% 
    filter(Age_class == 'LOCRC')%>% pull(mutation_rate)
  
  EOCRC_rate = mean(each_gene_mutation_EO)
  LOCRC_rate = mean(each_gene_mutation_LO)
  
  
  mutation_rate_pvalue <- wilcox.test(each_gene_mutation_EO, each_gene_mutation_LO)$p.value
  
  valid_vars <- "Age_class"
  covariates <- c("Race", "Panel", "Sex", "Tumor_site", "Sample_type", "Histology_type")
  
  for (cov in covariates) {
    
    cov_exists <- all(sapply(imp_list_Gene_clinical_info, function(df) cov %in% colnames(df)))
    if (!cov_exists) next
    
    
    sufficient_levels <- all(sapply(imp_list_Gene_clinical_info, function(df) {
      length(unique(df[[cov]])) >= 2
    }))
    
    if (sufficient_levels) {
      valid_vars <- c(valid_vars, cov)
    }
  }
  
  
  if (length(valid_vars) == 0) {
    formula_base <- "Gene_status ~ 1" 
  } else {
    formula_base <- paste("Gene_status ~", paste(valid_vars, collapse = " + "))
  }
  
  
  logistic_models <- lapply(imp_list_Gene_clinical_info, function(df) {
    tryCatch({
      glm(as.formula(formula_base), family = binomial(link = "logit"), data = df)
    }, error = function(e) {
      
      return(NULL)
    })
  })
  
  
  logistic_models <- Filter(Negate(is.null), logistic_models)
  
  
  if (length(logistic_models) == 0) {
    return()
  }
  
  
  summary_pooled <- summary(pool(logistic_models))
  
  locrc_eocrc_results <- summary_pooled %>% 
    filter(term == "Age_classEOCRC") %>% 
    mutate(conf.low = estimate - 1.959964 * std.error,
           conf.high = estimate + 1.959964 * std.error)%>% 
    mutate(OR = exp(estimate),
           `95% lower limit` = exp(conf.low),
           `95% upper limit` = exp(conf.high))
  
  df_Gene_test <- data.table::data.table(
    Hugo_Symbol = Gene_logistic,
    EOCRC_total_case= panels_samples_m1,
    LOCRC_total_case= panels_samples_m2,
    EOCRC = mut_samples_m1, 
    LOCRC = mut_samples_m2, 
    EOCRC_freq = m1_freq,
    LOCRC_freq = m2_freq,
    pval = locrc_eocrc_results$p.value, 
    OR = locrc_eocrc_results$OR,
    `95% upper limit` = locrc_eocrc_results$`95% upper limit`,
    `95% lower limit` = locrc_eocrc_results$`95% lower limit`,
    EOCRC_mutation_rate = EOCRC_rate,
    LOCRC_mutation_rate = LOCRC_rate,
    mutation_rate_pval = mutation_rate_pvalue,
    Country = Country_logistic,
    TMB_Status = TMB_status,
    Formula = formula_base
    
  )
  
  return(list(df_Gene_test))
  
}



generate_bubble_table <- function(guojia,Minmut){
  
  
  country_tsb<- rt@clinical.data[Country == guojia]$Tumor_Sample_Barcode
  
  
  if (length(intersect(clin.EOCRC@clinical.data$Tumor_Sample_Barcode, country_tsb)) > 0 &
      length(intersect(clin.LOCRC@clinical.data$Tumor_Sample_Barcode, country_tsb)) > 0){
    
    
    
    rt_country_EO <- subsetMaf(clin.EOCRC,tsb=country_tsb)
    rt_country_LO <- subsetMaf(clin.LOCRC,tsb=country_tsb)
    
    bubble_genes <- uniquegenes(m1 = rt_country_EO, m2 = rt_country_LO, minMut = Minmut)
    
    logistic_country_overall <- list()
    cl <- makeCluster(10)
    registerDoParallel(cl)
    
    logistic_country_overall <- foreach(gene = bubble_genes, .combine =  c, 
                                        .packages = c("data.table", "tidyverse","mice","stats"),
                                        .export = c("logstic_country_mafcompare",'imp_list',
                                                    'gene_panel_map','gene_mutation_counts')) %dopar% {
                                                      
                                                      logstic_country_mafcompare(m1 = rt_country_EO, m2 = rt_country_LO,
                                                                                 Clinical.data = rbind(rt_country_EO@clinical.data,
                                                                                                       rt_country_LO@clinical.data),  
                                                                                 Gene_logistic =gene,Country_logistic = guojia,TMB_status = 'Overall')
                                                      
                                                    }
    
    logistic_country_overall <- rbindlist(logistic_country_overall)
    stopCluster(cl) 
    
    
  }
  
  if (length(intersect(clin.EOCRC_hyper@clinical.data$Tumor_Sample_Barcode, country_tsb)) > 0 &
      length(intersect(clin.LOCRC_hyper@clinical.data$Tumor_Sample_Barcode, country_tsb)) > 0) {
    
    rt_country_EO_hyper <- subsetMaf(clin.EOCRC_hyper,tsb=country_tsb)
    rt_country_LO_hyper <- subsetMaf(clin.LOCRC_hyper,tsb=country_tsb)
    
    bubble_genes_hyper <- uniquegenes(m1 = rt_country_EO_hyper, m2 = rt_country_LO_hyper,
                                      minMut = Minmut)
    
    logistic_country_hyper <- list()
    cl <- makeCluster(10)
    registerDoParallel(cl)
    
    logistic_country_hyper <- foreach(gene = bubble_genes_hyper, .combine =  c, 
                                      .packages = c("data.table", "tidyverse","mice","stats"),
                                      .export = c("logstic_country_mafcompare",'imp_list',
                                                  'gene_panel_map','gene_mutation_counts')) %dopar% {
                                                    
                                                    logstic_country_mafcompare(m1 = rt_country_EO_hyper, m2 = rt_country_LO_hyper,
                                                                               Clinical.data = rbind(rt_country_EO_hyper@clinical.data,
                                                                                                     rt_country_LO_hyper@clinical.data), 
                                                                               Gene_logistic =gene,Country_logistic = guojia,TMB_status = 'Hyper')
                                                    
                                                  }
    
    logistic_country_hyper <- rbindlist(logistic_country_hyper)
    stopCluster(cl) 
    
    
  } 
  
  if (length(intersect(clin.EOCRC_nonhyper@clinical.data$Tumor_Sample_Barcode, country_tsb)) > 0 &
      length(intersect(clin.LOCRC_nonhyper@clinical.data$Tumor_Sample_Barcode, country_tsb)) > 0) {
    
    rt_country_EO_nonhyper <- subsetMaf(clin.EOCRC_nonhyper,tsb=country_tsb)
    rt_country_LO_nonhyper <- subsetMaf(clin.LOCRC_nonhyper,tsb=country_tsb)
    
    bubble_genes_nonhyper <- uniquegenes(m1 = rt_country_EO_nonhyper, m2 = rt_country_LO_nonhyper,
                                         minMut = Minmut)
    
    
    logistic_country_nonhyper <- list()
    cl <- makeCluster(10)
    registerDoParallel(cl)
    
    logistic_country_nonhyper <- foreach(gene = bubble_genes_nonhyper, .combine =  c, 
                                         .packages = c("data.table", "tidyverse","mice","stats"),
                                         .export = c("logstic_country_mafcompare",'imp_list',
                                                     'gene_panel_map','gene_mutation_counts')) %dopar% {
                                                       
                                                       logstic_country_mafcompare(m1 = rt_country_EO_nonhyper, m2 = rt_country_LO_nonhyper,
                                                                                  Clinical.data = rbind(rt_country_EO_nonhyper@clinical.data,
                                                                                                        rt_country_LO_nonhyper@clinical.data), 
                                                                                  Gene_logistic =gene,Country_logistic = guojia,TMB_status = 'Non-hyper')
                                                       
                                                     }
    
    logistic_country_nonhyper <- rbindlist(logistic_country_nonhyper)
    stopCluster(cl) 
    
  }
  
  return(rbind(logistic_country_overall,logistic_country_hyper,logistic_country_nonhyper))
}

France_data <-  generate_bubble_table('France',Minmut = 2)
Nigeria_data <-  generate_bubble_table('Nigeria',Minmut = 2)
Netherlands_data <-  generate_bubble_table('Netherlands',Minmut = 2)
Canada_data <-  generate_bubble_table('Canada',Minmut = 2)
Spain_data <-  generate_bubble_table('Spain',Minmut = 2)
China_data <-   generate_bubble_table('China',Minmut = 2)
US_data <-  generate_bubble_table('US',Minmut = 2)

total_bubble_table <- rbind(France_data,Nigeria_data,Netherlands_data,
                            Canada_data,Spain_data,China_data,US_data) %>% 
  rename(Gene = Hugo_Symbol )


filter_gene_country <- total_bubble_table %>%
  filter(TMB_Status == "Overall" & EOCRC_freq < 0.05)  %>%
  dplyr::select(Gene, Country,EOCRC_freq)

total_bubble_table_filter <- total_bubble_table %>%
  anti_join(filter_gene_country, by = c("Gene", "Country"))


{
  
  gene_country_counts <- total_bubble_table_filter %>%
    group_by(Gene) %>%
    summarise(Country_count = n_distinct(Country))
  
  genes_in_6_countries <- gene_country_counts %>%
    filter(Country_count >= 6) %>%
    dplyr::select(Gene)
  }

bubble_table = total_bubble_table_filter %>%
  filter(Gene %in%  genes_in_6_countries$Gene)%>% 
  mutate(adjPval = p.adjust(p = pval, method = "BH")) 

bubble_table$p_sig <- as.character(symnum(bubble_table$adjPval, 
                                          cutpoints = c(0, 0.01, 0.05, 0.25, 1), 
                                          symbols = c("***", "**", "*", "")))

p1 <- ggplot(data = bubble_table, 
             aes(x = Gene, y = TMB_Status)) +
  geom_point(aes(color = factor(case_when(
    OR > 1 ~ 'OR > 1',
    OR < 1 ~ 'OR < 1',
    OR == 1 ~ 'No significance'
  ), levels = c('OR > 1', 'OR < 1', 'No significance'))), 
  size = 6, alpha = 0.7, show.legend = c(size = FALSE))  +
  
  geom_text(aes(label = p_sig), vjust = 0.8, hjust = 0.5, size = 3, color = "black"
  ) +
  facet_grid(Country ~ ., scales = "free_y", space = "free") +  
  scale_y_discrete(position = "left") +  
  scale_color_manual(values = c('OR > 1' = '#d53e4f', 'OR < 1' = '#4393c3', 'No significance' = 'gray')) +
  theme_few() +
  labs(y ='TMB Category') +
  theme(
    strip.background = element_blank(), 
    strip.text.y.right = element_text(size = 13, angle = 0),
    panel.spacing = unit(0, "lines"),   
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    axis.text.y = element_text(color = "black", size = 10, angle = 0),
    axis.text.x = element_text(color = "black", size = 10, angle = 45, vjust = 0.5, hjust = 0.5),
    legend.text = element_text(color = "black", size = 10),
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

p1

ggsave(paste(dir,'picture/8_intersect_gene_among_countries.pdf',sep = ''),
       width = 12,height = 6,dpi = 300)

write.csv(total_bubble_table,file = paste0(dir,'table/8_total_bubble_data.csv'),
          quote = FALSE, row.names = FALSE)

