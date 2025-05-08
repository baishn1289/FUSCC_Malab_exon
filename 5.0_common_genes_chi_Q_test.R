#EO_LO_Country_common_genes_Chi_Q_test

logstic_country_commongenes <- function(m1, m2, Clinical.data,Gene_logistic_common,Country_common) {
  
  Alter_m1_Gene = m1@data[Hugo_Symbol == Gene_logistic_common,] %>% 
    left_join(m1@clinical.data[,.(Tumor_Sample_Barcode,Country)],by='Tumor_Sample_Barcode') %>% 
    filter(Country == Country_common) %>% 
    filter(!duplicated(Tumor_Sample_Barcode))
  
  Alter_m2_Gene = m2@data[Hugo_Symbol == Gene_logistic_common,] %>% 
    left_join(m2@clinical.data[,.(Tumor_Sample_Barcode,Country)],by='Tumor_Sample_Barcode') %>% 
    filter(Country == Country_common) %>% 
    filter(!duplicated(Tumor_Sample_Barcode))
  
  Gene_panels <- Clinical.data %>% 
    filter(Country == Country_common,Panel %in% gene_panel_map[[Gene_logistic_common]]) %>% 
    pull(Panel)
  
  mut_samples_m1 <- length(unique(Alter_m1_Gene$Tumor_Sample_Barcode))
  
  panels_samples_m1 <- length(unique(m1@clinical.data %>% 
                                       filter(Country == Country_common,Panel %in% Gene_panels) %>% 
                                       pull(Tumor_Sample_Barcode)))
  
  m1_freq <- mut_samples_m1 / panels_samples_m1
  
  mut_samples_m2 <- length(unique(Alter_m2_Gene$Tumor_Sample_Barcode))
  
  
  panels_samples_m2 <- length(unique(m2@clinical.data %>% 
                                       filter(Country == Country_common,Panel %in% Gene_panels) %>% 
                                       pull(Tumor_Sample_Barcode)))
  
  
  m2_freq <- mut_samples_m2 / panels_samples_m2
  
  imp_list_Gene_clinical_info <-lapply(imp_list, function(df) {
    dt <- as.data.table(df)
    dt <- dt[Tumor_Sample_Barcode %in% Clinical.data$Tumor_Sample_Barcode]
    dt <- dt[Country == Country_common]
    dt[, Gene_status := as.integer(Tumor_Sample_Barcode %in% c(Alter_m1_Gene$Tumor_Sample_Barcode, Alter_m2_Gene$Tumor_Sample_Barcode))]
    dt[, Gene_status := factor(Gene_status, levels = c(0, 1))]
    dt <- dt[Panel %in% Gene_panels]
    return(dt)
  })
  
  valid_vars <- "Age_class"
  
  covariates <- c( "Panel", "Sex", "Tumor_site", "Sample_type", "Histology_type")
  
  
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
    Hugo_Symbol = Gene_logistic_common,
    EOCRC_total_case= panels_samples_m1,
    LOCRC_total_case= panels_samples_m2,
    EOCRC = mut_samples_m1, 
    LOCRC = mut_samples_m2, 
    EOCRC_freq = m1_freq,
    LOCRC_freq = m2_freq,
    pval = locrc_eocrc_results$p.value, 
    OR = locrc_eocrc_results$OR,
    se = locrc_eocrc_results$std.error,
    beta = locrc_eocrc_results$estimate,
    `95% upper limit` = locrc_eocrc_results$`95% upper limit`,
    `95% lower limit` = locrc_eocrc_results$`95% lower limit`,
    Formula = formula_base,
    Country = Country_common
  )
  
  return(list(df_Gene_test))
}


result_list <- list()

for (country in unique(rt@clinical.data$Country)[1:7]) {
  for (gene in common_EO_LO_gene) {
    res <- logstic_country_commongenes(m1 = clin.EOCRC,
                                       m2 = clin.LOCRC,
                                       Clinical.data = rbind(clin.EOCRC@clinical.data, clin.LOCRC@clinical.data), 
                                       Gene_logistic_common = gene,
                                       Country_common = country)
    if (!is.null(res)) {
      result_list <- append(result_list, res)
    }
  }
}

result_df_Country <- data.table::rbindlist(result_list, fill = TRUE)

write_csv(result_df_Country,file = paste0(dir,'table/5.1_common_gene_across_country.csv'))                                  


#cochran Q test
{
  result_df_Country_Qtest <- result_df_Country %>%
    group_by(Hugo_Symbol) %>%
    filter(n_distinct(Country) == 7) %>%
    do(calculate_heterogeneity(.)) %>%
    ungroup()
  
}
write.csv(result_df_Country_Qtest,file = paste0(dir,'table/5.1_common_gene_across_country_Qtest.csv'))


#chi-sequare test
{
  result_Country_chi <-   result_df_Country %>%
    mutate(EOCRC_not_mut_case = EOCRC_total_case- EOCRC) %>% 
    dplyr::select(Hugo_Symbol, EOCRC_not_mut_case,EOCRC, Country) %>%   
    group_by(Hugo_Symbol) %>%
    filter(n_distinct(Country) == 7) %>%
    ungroup()  %>%
    split(.$Hugo_Symbol) %>%
    map_dfr(~ {
      contingency_table <- .x %>% dplyr::select(EOCRC_not_mut_case, EOCRC)
      
      expected_counts <- chisq.test(contingency_table)$expected
      
      
      if (any(expected_counts < 5)) {
        print(expected_counts)
        print(paste("Using Fisher test for:", unique(.x$Hugo_Symbol)))
        result_df_chi <- fisher.test(contingency_table, workspace = 2e7)
        method_used <- "fisher"
      } else {
        result_df_chi <- chisq.test(contingency_table)
        method_used <- "chi-square"
      }
      
      data.frame(Hugo_Symbol = unique(.x$Hugo_Symbol),
                 p.value = result_df_chi$p.value,
                 method = method_used)
    })
  
  write.csv(result_Country_chi,file = paste0(dir,'table/5.1_common_gene_across_country_chi.csv'))
  
}


result_df_Country_all <- result_df_Country %>%
  dplyr::select(Hugo_Symbol, EOCRC_freq, LOCRC_freq, Country) %>% 
  group_by(Hugo_Symbol) %>% 
  pivot_wider(names_from = 'Country', values_from = c('EOCRC_freq', 'LOCRC_freq')) %>%
  rename_with(~ gsub("^(.*)_freq_(.*)$", "\\2_\\1_freq", .), starts_with("EOCRC")) %>%
  rename_with(~ gsub("^(.*)_freq_(.*)$", "\\2_\\1_freq", .), starts_with("LOCRC")) %>% 
  left_join(result_Country_chi,by='Hugo_Symbol') %>% 
  left_join(result_df_Country_Qtest,by = 'Hugo_Symbol') %>% 
  dplyr::select(-df,-method) %>% 
  mutate(Q_test_adjpvalue = p.adjust(p = Q_test_p_value, method = "BH")%>% 
           format_Pval(),
         Pvalue = p.adjust(p = p.value, method = "BH")%>% 
           format_Pval()
  ) %>% 
  arrange(Q_test_p_value)

result_df_Country_all_publish <- result_df_Country_all[result_df_Country_all$Q_test_adjpvalue<0.05,] %>% 
  dplyr::select(-p.value,-Q_test_p_value)

write_csv(result_df_Country_all,file = paste0(dir,'table/5.1_common_gene_across_country_all.csv'))
write_csv(result_df_Country_all_publish,file = paste0(dir,'table/5.1_common_gene_across_country_all_publish.csv'))
