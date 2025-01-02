#EO_LO_Country_common_genes_Chi_Q_test

logstic_country_commongenes <- function(m1, m2, Clinical.data,Gene_test,Guojia) {
  
  Alter_m1_Gene = m1@data[Hugo_Symbol == Gene_test,] %>% 
                    left_join(m1@clinical.data[,.(Tumor_Sample_Barcode,Country)],by='Tumor_Sample_Barcode') %>% 
                    filter(Country == Guojia) %>% 
                    filter(!duplicated(Tumor_Sample_Barcode))
  
  Alter_m2_Gene = m2@data[Hugo_Symbol == Gene_test,] %>% 
                    left_join(m2@clinical.data[,.(Tumor_Sample_Barcode,Country)],by='Tumor_Sample_Barcode') %>% 
                    filter(Country == Guojia) %>% 
                    filter(!duplicated(Tumor_Sample_Barcode))
  
  Gene_panels = Clinical.data[, c('Tumor_Sample_Barcode', 'Age_class', 'Sex', 'Country','Panel', 'Race')] %>% 
                  filter(Country == Guojia,Panel %in% gene_panel_map[[Gene_test]]) %>% 
                  pull(Panel)
  
  mut_samples_m1 = length(unique(Alter_m1_Gene$Tumor_Sample_Barcode))
  
  panels_samples_m1 = length(unique(m1@clinical.data %>% 
                                       filter(Country == Guojia,Panel %in% Gene_panels) %>% 
                                       pull(Tumor_Sample_Barcode)))
  
  m1_freq <- mut_samples_m1 / panels_samples_m1

  mut_samples_m2 <- length(unique(Alter_m2_Gene$Tumor_Sample_Barcode))
  
  panels_samples_m2 <- length(unique(m2@clinical.data %>% 
                                       filter(Country == Guojia,Panel %in% Gene_panels) %>% 
                                       pull(Tumor_Sample_Barcode)))
  
  m2_freq <- mut_samples_m2 / panels_samples_m2
  
  Gene_clinical_info  <-  Clinical.data[, Gene_status := ifelse(Tumor_Sample_Barcode %in% c(Alter_m1_Gene$Tumor_Sample_Barcode,
                                                                                            Alter_m2_Gene$Tumor_Sample_Barcode), 1, 0)]%>% 
                            mutate(Age_class = factor(Age_class,levels = c('LOCRC','EOCRC'))) %>% 
                            filter(Country == Guojia,Panel %in% Gene_panels)
  
  formula_base <- "Gene_status ~ Age_class"
  
  if (length(unique(Gene_clinical_info$Race)) > 1) {
    formula_base <- paste(formula_base, "+ Race")
  }
  if (length(unique(Gene_clinical_info$Panel)) > 1) {
    formula_base <- paste(formula_base, "+ Panel")
  }
  if (length(unique(Gene_clinical_info$Country)) > 1) {
    formula_base <- paste(formula_base, "+ Country")
  }
  if (length(unique(Gene_clinical_info$Sex)) > 1) {
    formula_base <- paste(formula_base, "+ Sex")
  }
  if (length(unique(Gene_clinical_info$Tumor_site)) > 1) {
    formula_base <- paste(formula_base, "+ Tumor_site")
  }
  if (length(unique(Gene_clinical_info$Sample_type)) > 1) {
    formula_base <- paste(formula_base, "+ Sample_type")
  }
  if (length(unique(Gene_clinical_info$Histology_type)) > 1) {
    formula_base <- paste(formula_base, "+ Histology_type")
  }

  fit <- glm(as.formula(formula_base), 
             family = binomial(link = "logit"), 
             data = Gene_clinical_info)
  
  OR_age_class <- exp(coef(fit)["Age_classEOCRC"])
  pval_age_class <- summary(fit)$coefficients["Age_classEOCRC", "Pr(>|z|)"]
  coef_age_class <- coef(fit)["Age_classEOCRC"]
  se_age_class <- summary(fit)$coefficients["Age_classEOCRC", "Std. Error"]
  ci_low <- exp(coef_age_class - 1.96 * se_age_class)
  ci_up  <- exp(coef_age_class + 1.96 * se_age_class)
  
  df_Gene_test <- data.table::data.table(
                    Hugo_Symbol = Gene_test,
                    EOCRC_total_case= panels_samples_m1,
                    LOCRC_total_case= panels_samples_m2,
                    EOCRC = mut_samples_m1, 
                    LOCRC = mut_samples_m2, 
                    EOCRC_freq = m1_freq,
                    LOCRC_freq = m2_freq,
                    pval = pval_age_class, 
                    OR = OR_age_class,
                    ci.up = ci_up,
                    ci.low = ci_low,
                    se = se_age_class,
                    beta = coef_age_class,
                    formula = formula_base,
                    Country = Guojia
                   )

  result_df <<- rbind(result_df, df_Gene_test)  
}

result_df <- data.table::data.table(
              Hugo_Symbol = character(),
              EOCRC_total_case= numeric(),
              LOCRC_total_case= numeric(),
              EOCRC = numeric(), 
              LOCRC = numeric(), 
              EOCRC_freq = numeric(),
              LOCRC_freq = numeric(),
              pval = numeric(), 
              OR = numeric(),
              ci.up =  numeric(),
              ci.low =  numeric(),
              se = numeric(),
              beta = numeric(),
              formula =character(),
              Country = character()
            )

for (country in unique(rt@clinical.data$Country)[1:7]){
  
  for (genes in common_EO_LO_gene) {
    
    logstic_country_commongenes(m1 = clin.EOCRC, m2 = clin.LOCRC,
                               Clinical.data = rbind(clin.EOCRC@clinical.data,
                                                     clin.LOCRC@clinical.data), 
                               Gene_test = genes, Guojia = country)
  }
}

  result_df_Country <- result_df

  write_csv(result_df_Country,file = './table/5.1_common_gene_across_country.csv')



{  #cochran Q test
  
  result_df_Country_Qtest <- result_df_Country %>%
    group_by(Hugo_Symbol) %>%
    filter(n_distinct(Country) == 7) %>%
    do(calculate_heterogeneity(.)) %>%
    ungroup()
  
}

  write.csv(result_df_Country_Qtest,file = './table/5.1_common_gene_across_country_Qtest.csv')



{  #chi-sequare test
  
  result_Country_chi <-   result_df_Country %>%
                            mutate(EOCRC_not_mut_case = EOCRC_total_case- EOCRC) %>% 
                            dplyr::select(Hugo_Symbol, EOCRC_not_mut_case,EOCRC, Country) %>%   
                            group_by(Hugo_Symbol) %>%
                            filter(n_distinct(Country) == 7) %>%
                            ungroup()  %>%
                            split(.$Hugo_Symbol) %>%
                            map_dfr(~ {
                              contingency_table <- .x %>% dplyr::select(EOCRC_not_mut_case, EOCRC)
                              
                              tryCatch({
                                result_df_chi <- chisq.test(contingency_table)
                                data.frame(Hugo_Symbol = unique(.x$Hugo_Symbol),
                                           p.value = result_df_chi$p.value,
                                           method = "chi-square")
                                
                              }, warning = function(w) {
                                print(paste("warning:", unique(.x$Hugo_Symbol)))
                                result_df_chi <- fisher.test(contingency_table, workspace = 2e7)
                                data.frame(Hugo_Symbol = unique(.x$Hugo_Symbol),
                                           p.value = result_df_chi$p.value,
                                           method = "fisher")
                              })
                            })
  
  write.csv(result_Country_chi,file = './table/5.1_common_gene_across_country_chi.csv')
  
}

  result_df_Country_all <- result_df_Country %>%
                            dplyr::select(Hugo_Symbol, EOCRC_freq, LOCRC_freq, Country) %>% 
                            group_by(Hugo_Symbol) %>% 
                            pivot_wider(names_from = 'Country', values_from = c('EOCRC_freq', 'LOCRC_freq')) %>%
                            rename_with(~ gsub("^(.*)_freq_(.*)$", "\\2_\\1_freq", .), starts_with("EOCRC")) %>%
                            rename_with(~ gsub("^(.*)_freq_(.*)$", "\\2_\\1_freq", .), starts_with("LOCRC")) %>% 
                            left_join(result_Country_chi,by='Hugo_Symbol') %>% 
                            left_join(result_df_Country_Qtest,by = 'Hugo_Symbol') %>% 
                            dplyr::select(-df,-method)

  write_csv(result_df_Country_all,file = './table/5.1_common_gene_across_country_all.csv')
