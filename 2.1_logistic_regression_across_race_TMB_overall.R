
logstic_mafcompare_race <- function(m1, m2, Clinical.data,Gene_logistic,tmb_status,Race_logistic) {
  
  Alter_m1_Gene <- as.data.table(m1@data)[Hugo_Symbol == Gene_logistic]
  Alter_m1_Gene <- unique(Alter_m1_Gene, by = "Tumor_Sample_Barcode")
  
  Alter_m2_Gene <- as.data.table(m2@data)[Hugo_Symbol == Gene_logistic]
  Alter_m2_Gene <- unique(Alter_m2_Gene, by = "Tumor_Sample_Barcode")
  
  Gene_panels <- Clinical.data %>% 
    filter(Panel %in% gene_panel_map[[Gene_logistic]]) %>% 
    pull(Panel)
  
  mut_samples_m1 <- length(unique(Alter_m1_Gene$Tumor_Sample_Barcode))
  panels_samples_m1 <- length(unique(m1@clinical.data[m1@clinical.data$Panel %in% Gene_panels, ]$Tumor_Sample_Barcode))
  m1_freq <- mut_samples_m1 / panels_samples_m1
  
  mut_samples_m2 <- length(unique(Alter_m2_Gene$Tumor_Sample_Barcode))
  panels_samples_m2 <- length(unique(m2@clinical.data[m2@clinical.data$Panel %in% Gene_panels, ]$Tumor_Sample_Barcode))
  m2_freq <- mut_samples_m2 / panels_samples_m2
  
  if ( m1_freq < 0.01 | m2_freq<0.01 ) {
    print(paste(Gene_logistic,
                "Low mutation frequencies, returning from function.",
                sep = '____'))
    return()
  }
  
  
  imp_list_Gene_clinical_info <-lapply(imp_list, function(df) {
    dt <- as.data.table(df)
    dt <- dt[Tumor_Sample_Barcode %in% Clinical.data$Tumor_Sample_Barcode]
    dt[, Gene_status := as.integer(Tumor_Sample_Barcode %in% c(Alter_m1_Gene$Tumor_Sample_Barcode, Alter_m2_Gene$Tumor_Sample_Barcode))]
    dt[, Gene_status := factor(Gene_status, levels = c(0, 1))]
    dt <- dt[Panel %in% Gene_panels]
    return(dt)
  })
  
  
  if (any(sapply(imp_list_Gene_clinical_info[[1]][, c("Sex", "Tumor_site", "Panel",
                                                      'Histology_type','Sample_type')], 
                 function(x) length(unique(x)) <= 1))) {
    return()
  }
  
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
  
  logistic_models <- lapply(imp_list_Gene_clinical_info, function(df) {
    tryCatch({
      glm(as.formula("Gene_status ~ Age_class + Sex + Tumor_site + Panel+ Histology_type+Sample_type"), family = binomial(link = "logit"),
          data = df)
    }, error = function(e) {
      return(NULL)
    })
  })
  
  
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
    se = locrc_eocrc_results$std.error,
    beta = locrc_eocrc_results$estimate,
    `95% upper limit` = locrc_eocrc_results$`95% upper limit`,
    `95% lower limit` = locrc_eocrc_results$`95% lower limit`,
    EOCRC_mutation_rate = EOCRC_rate,
    LOCRC_mutation_rate = LOCRC_rate,
    mutation_rate_pval = mutation_rate_pvalue,
    TMB_Status = tmb_status,
    Race = Race_logistic
  )
  
  return(list(df_Gene_test))
  
}



{
  logstic_race_model <- function(race,source_EO,source_LO,tmb_status,Mutation_number =5){
    
    race_EO = subsetMaf(source_EO,tsb = source_EO@clinical.data[Race == race]$Tumor_Sample_Barcode)
    race_LO = subsetMaf(source_LO,tsb = source_LO@clinical.data[Race == race]$Tumor_Sample_Barcode)
    
    race_clinical_data = rbind(race_EO@clinical.data,race_LO@clinical.data)
    
    compare_genes_list <- uniquegenes(m1=race_EO,m2 = race_LO,minMut = Mutation_number)
    
    result_df <- list()
    
    cl <- makeCluster(12)
    registerDoParallel(cl)
    
    result_df <- foreach(gene = compare_genes_list, .combine =  c, 
                         .packages = c("data.table", "tidyverse","mice","stats"),
                         .export = c("logstic_mafcompare_race",'imp_list',
                                     'gene_panel_map','gene_mutation_counts')) %dopar% {
                                       
                                       logstic_mafcompare_race(m1 = race_EO, m2 = race_LO, 
                                                               Clinical.data = race_clinical_data, 
                                                               tmb_status = tmb_status, Race_logistic = race, 
                                                               Gene_logistic =gene)
                                       
                                     }
    
    stopCluster(cl) 
    
    
    result_df_modified <- rbindlist(result_df) %>% 
      mutate(adjPval = p.adjust(p = pval, method = "BH"),
             mutation_rate_adjPval =  p.adjust(p = mutation_rate_pval, method = "BH")) %>% 
      mutate(across(c(OR, `95% upper limit`, `95% lower limit`), ~ round(., digits = 2)),
             across(c(pval, adjPval, mutation_rate_pval,mutation_rate_adjPval), ~ format_Pval(.))
      )
    return(result_df_modified)
  }
  
  race_compare_results <- lapply(races_names, function(r) {
    logstic_race_model(
      race = r,
      source_EO = clin.EOCRC,
      source_LO = clin.LOCRC,
      tmb_status = 'Overall',
      Mutation_number = 5
    )
  })
  
  result_df_race <- rbindlist(race_compare_results)
  
  write.csv(result_df_race,
            file = paste0(dir,'table/3.1_race_overall_original_data.csv'))
  
}

#cochran Q test
{
  
  calculate_heterogeneity <- function(data) {
    
    betas <- data$beta
    ses <- data$se
    
    weights <- 1 / (ses^2)
    
    weighted_mean <- sum(weights * betas) / sum(weights)
    
    Q <- sum(weights * (betas - weighted_mean)^2)
    
    df <- length(betas) - 1
    
    p_value <-pchisq(Q, df, lower.tail = FALSE) 
    
    
    return(data.frame(Q_stat = Q, df = df, Q_test_p_value = p_value))
  }
  
  result_df_race_Qtest <- result_df_race %>%
    group_by(Hugo_Symbol) %>%
    filter(n_distinct(Race) == 3) %>%
    filter(any(adjPval < 0.05)) %>%
    do(calculate_heterogeneity(.)) %>%
    ungroup()
  
}

write.csv(result_df_race_Qtest,
          file = paste0(dir,'table/3.1_result_df_race_Qtest.csv'))

#chi-sequare test
{
  
  result_race_chi <-  result_df_race %>%
    mutate(EOCRC_not_mut_case = EOCRC_total_case- EOCRC) %>% 
    dplyr::select(Hugo_Symbol, EOCRC_not_mut_case,EOCRC, Race) %>%   
    group_by(Hugo_Symbol) %>%
    filter(n_distinct(Race) == 3) %>%
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
  
  write.csv(result_race_chi,
            file = paste0(dir,'table/3.1_result_race_chi.csv'))
  
}



{
  
  Asia_data <- result_df_race[Race == "Asian Pacific lslander"] %>%
    dplyr::select(Hugo_Symbol, EOCRC_freq, LOCRC_freq, 
                  adjPval, OR, `95% upper limit`, `95% lower limit`,
                  EOCRC_mutation_rate,LOCRC_mutation_rate,
                  mutation_rate_adjPval) %>%
    rename_with(~ paste("Asia", ., sep = "_"), -Hugo_Symbol)
  
  Black_data <- result_df_race[Race == "Black"]%>% 
    dplyr::select(Hugo_Symbol, EOCRC_freq, LOCRC_freq, 
                  adjPval, OR, `95% upper limit`, `95% lower limit`,
                  EOCRC_mutation_rate,LOCRC_mutation_rate,
                  mutation_rate_adjPval)%>%
    rename_with(~ paste("Black", ., sep = "_"), -Hugo_Symbol)
  
  White_data <- result_df_race[Race == "White"]%>% 
    dplyr::select(Hugo_Symbol, EOCRC_freq, LOCRC_freq, 
                  adjPval, OR, `95% upper limit`, `95% lower limit`,
                  EOCRC_mutation_rate,LOCRC_mutation_rate,
                  mutation_rate_adjPval)%>%
    rename_with(~ paste("White", ., sep = "_"), -Hugo_Symbol)
  
  
  Race_data_all <- full_join(
    Asia_data,Black_data,by='Hugo_Symbol'  ) %>% 
    full_join(White_data,by='Hugo_Symbol') %>% 
    left_join(result_race_chi,by='Hugo_Symbol') %>%
    left_join(result_df_race_Qtest,by='Hugo_Symbol')%>%
    dplyr::select(-df) %>% 
    mutate(Q_test_adjpvalue = p.adjust(p = Q_test_p_value, method = "BH")%>% 
             format_Pval(),
           Pvalue = p.adjust(p = p.value, method = "BH")%>% 
             format_Pval()
           
    ) %>% 
    arrange(Q_test_p_value)
  
  Race_data_all_publish <- Race_data_all[Q_test_adjpvalue<0.05] %>% 
    select(-p.value,-Q_test_p_value,-method)
  
  write_csv(Race_data_all,
            file = paste0(dir,'table/3.1_Race_data_all.csv'))
  write_csv(Race_data_all_publish,
            file = paste0(dir,'table/3.1_Race_data_all_publish.csv'))
}



