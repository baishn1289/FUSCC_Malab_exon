#EO_LO_Country_common_genes_Chi_Q_test




#相当于是把2.2代码的种族替换成国家进行比较即可


logstic_mafcompare_country <- function(m1, m2, Clinical.data,Gene_test,Guojia) {
  
  
  Alter_m1_Gene = m1@data[Hugo_Symbol == Gene_test,] %>% 
    filter(!duplicated(Tumor_Sample_Barcode))
  
  Alter_m2_Gene = m2@data[Hugo_Symbol == Gene_test,] %>% 
    filter(!duplicated(Tumor_Sample_Barcode))
  
  
  Gene_panels <- panel_hugo_symbol[Hugo_Symbol == Gene_test] %>%
    left_join(Clinical.data[, c('Tumor_Sample_Barcode', 'Age_class', 'Sex', 'Country', 'Race')], 
              by = "Tumor_Sample_Barcode") %>%
    filter(!duplicated(Tumor_Sample_Barcode),
           !is.na(Age_class)) %>%pull(Panel)
    # group_by(Panel) %>%
    # reframe(
    #   Age_class_count = n_distinct(Age_class)
    # ) %>%
    # filter(Age_class_count ==2 ) %>% 
    # pull(Panel)
  
  
  
  # if (length(Gene_panels) == 0) {
  #   print(paste(Guojia,Gene_test, '没有一个中心有两种Age_class', sep = '_'))
  #   
  #   singlesex_panel_record_log <<- rbind(singlesex_panel_record_log, 
  #                                        data.frame(Gene_test = Gene_test,
  #                                                   stringsAsFactors = FALSE))
  #   
  #   return()
  # }
  
  
  mut_samples_m1 <- length(unique(Alter_m1_Gene$Tumor_Sample_Barcode))
  
  
  panels_samples_m1 <- length(unique(m1@clinical.data[m1@clinical.data$Panel %in% Gene_panels, ]$Tumor_Sample_Barcode))
  
  
  m1_freq <- mut_samples_m1 / panels_samples_m1
  
  # ------------
  mut_samples_m2 <- length(unique(Alter_m2_Gene$Tumor_Sample_Barcode))
  
  
  panels_samples_m2 <- length(unique(m2@clinical.data[m2@clinical.data$Panel %in% Gene_panels, ]$Tumor_Sample_Barcode))
  
  
  m2_freq <- mut_samples_m2 / panels_samples_m2
  
  
  Gene_clinical_info  <-  Clinical.data[, Gene_status := ifelse(Tumor_Sample_Barcode %in% c(Alter_m1_Gene$Tumor_Sample_Barcode,
                                                                                            Alter_m2_Gene$Tumor_Sample_Barcode), 1, 0)]%>% 
    left_join(TMB_mutload[,.(Tumor_Sample_Barcode,TMB)],by = 'Tumor_Sample_Barcode') %>% 
    mutate(Age_class = factor(Age_class,levels = c('LOCRC','EOCRC'))) %>% 
    filter(Panel %in% Gene_panels) 
  
  formula_base <- "Gene_status ~ Age_class+ TMB"
  
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
  # if (length(unique(Gene_clinical_info$Histology_type)) > 1) {
  #   formula_base <- paste(formula_base, "+ Histology_type")
  # }
  
  
  fit <- glm(as.formula(formula_base), family = binomial(link = "logit"), data = Gene_clinical_info)
  
  
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



logstic_country_model <- function(country,source_EO,source_LO){
  Country_EO = subsetMaf(source_EO,tsb = source_EO@clinical.data[Country == country]$Tumor_Sample_Barcode)
  Country_LO = subsetMaf(source_LO,tsb = source_LO@clinical.data[Country == country]$Tumor_Sample_Barcode)
  
  
  # Country_EO = subsetMaf(clin.EOCRC,tsb = clin.EOCRC@clinical.data[Country == "Nigeria"]$Tumor_Sample_Barcode)
  # Country_LO = subsetMaf(clin.LOCRC,tsb = clin.LOCRC@clinical.data[Country == "Nigeria"]$Tumor_Sample_Barcode)
  
  Country_clinical_data = rbind(Country_EO@clinical.data,Country_LO@clinical.data)
  
 
  
  # race_hugo_symbol <- panel_hugo_symbol[Hugo_Symbol %in% compare_genes_list] %>%
  #   merge(rt@clinical.data[, c('Tumor_Sample_Barcode', 'Age_class','Country',"Sex","Race")], 
  #         by = "Tumor_Sample_Barcode") %>%
  #   group_by(Hugo_Symbol) %>%
  #   reframe(
  #     Panel_count = n_distinct(Panel),
  #     Country_count = n_distinct(Country),
  #     Race_count = n_distinct(Race),
  #     Sex_count = n_distinct(Sex)
  #   ) %>%
  #   pull(Hugo_Symbol)
  
  
  for (i in common_EO_LO_gene) {
    
      logstic_mafcompare_country(m1 = Country_EO, m2 = Country_LO,
                              Clinical.data = Country_clinical_data, Gene_test = i,
                              Guojia = country)
    
    # logstic_mafcompare_country(m1 = Country_EO, m2 = Country_LO,
    #                            Clinical.data = rt@clinical.data[Country == "Nigeria"], Gene_test = i,
    #                            Guojia = "Nigeria")
  
    }
  
}



singlesex_panel_record_log <- data.frame(Gene_test = character(),
                                         stringsAsFactors = FALSE)

result_df_Country <- data.table::data.table(
  Hugo_Symbol = NA,
  EOCRC_total_case= NA,
  LOCRC_total_case= NA,
  EOCRC = NA, 
  LOCRC = NA, 
  EOCRC_freq = NA,
  LOCRC_freq = NA,
  pval = NA, 
  OR = NA,
  ci.up =  NA,
  ci.low =  NA,
  se = NA,
  beta = NA,
  formula =NA,
  Country = NA,
  adjPval = NA
)


for (i in unique(rt@clinical.data$Country)[1:7]){
  
  problematic_genes <- c() 
  
  result_df <- data.table::data.table(
    Hugo_Symbol = NA,
    EOCRC_total_case= NA,
    LOCRC_total_case= NA,
    EOCRC = NA, 
    LOCRC = NA, 
    EOCRC_freq = NA,
    LOCRC_freq = NA,
    pval = NA, 
    OR = NA,
    ci.up =  NA,
    ci.low =  NA,
    se = NA,
    beta = NA,
    formula =NA,
    Country = NA
  )
  
  logstic_country_model(country = i,
                     source_EO=clin.EOCRC,
                     source_LO=clin.LOCRC)
  
  result_df= na.omit(result_df)
  
  
  write_csv(result_df,file = paste('./table/5.1_common_genes',i,'original.csv',sep = '_'))
  
  #此处变量需再斟酌(荷兰只有结肠癌，因此部位只有1个，是否考虑第一个代码中保留左右结肠的分类，AACR数据的colon单独做一个分类
  #像许多国家也只有1个assay；无法纳入逻辑回归)
  result_df_modified <- result_df %>% 
    # filter(!Hugo_Symbol %in% problematic_genes) %>%
    # filter(formula == 'Gene_status ~ Age_class+ TMB + Sex') %>% 
    mutate(adjPval = p.adjust(p = pval, method = "fdr")) 
  
  write_csv(result_df_modified,
            file = paste('./table/5.1_common_genes',i,'modified.csv',sep = '_'))
  
  result_df_Country <-  rbind(result_df_Country,result_df_modified)
  
  
}

result_df_Country <- na.omit(result_df_Country)

write_csv(result_df_Country,file = './table/5.1_common_gene_across_country.csv')



{#cochran Q test
  
  #先筛出在三个种族中都有的基因，然后再筛出在至少1个种族中显著的基因
  #对这个基因收集各个种族的beta和se，做Q检验
  result_df_Country_Qtest <- result_df_Country %>%
    group_by(Hugo_Symbol) %>%
    filter(n_distinct(Country) == 7) %>%
    # filter(any(adjPval < 0.05)) %>%
    do(calculate_heterogeneity(.)) %>%
    ungroup()
  
}

write.csv(result_df_Country_Qtest,file = './table/5.1_common_gene_across_country_Qtest.csv')



{#chi-sequare test
  
  result_Country_chi <-   result_df_Country %>%
    mutate(EOCRC_not_mut_case = EOCRC_total_case- EOCRC) %>% 
    select(Hugo_Symbol, EOCRC_not_mut_case,EOCRC, Country) %>%   
    group_by(Hugo_Symbol) %>%
    filter(n_distinct(Country) == 7) %>%
    ungroup()  %>%
    split(.$Hugo_Symbol) %>%
    map_dfr(~ {
      contingency_table <- .x %>% select(EOCRC_not_mut_case, EOCRC)
      
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
  
  write.csv(result_race_chi,file = './table/5.1_common_gene_across_country_chi.csv')
  
}





