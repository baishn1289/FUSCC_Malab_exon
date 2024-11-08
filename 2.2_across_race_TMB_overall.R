
# library(gt)
library(tidyverse)
library(devtools)

logstic_mafcompare_race <- function(m1, m2, Clinical.data,Gene_test,tmb_status,zhongzu) {
  
  
  Alter_m1_Gene = m1@data[Hugo_Symbol == Gene_test,] %>% 
    filter(!duplicated(Tumor_Sample_Barcode))
  
  Alter_m2_Gene = m2@data[Hugo_Symbol == Gene_test,] %>% 
    filter(!duplicated(Tumor_Sample_Barcode))
  
  
  Gene_panels <- panel_hugo_symbol[Hugo_Symbol == Gene_test] %>%
    left_join(Clinical.data[, c('Tumor_Sample_Barcode', 'Age_class', 'Sex', 'Country', 'Race')], 
              by = "Tumor_Sample_Barcode") %>%
    filter(!duplicated(Tumor_Sample_Barcode),
           !is.na(Age_class)) %>%
    group_by(Panel) %>%
    reframe(
      Age_class_count = n_distinct(Age_class)
    ) %>%
    filter(Age_class_count ==2 ) %>% 
    pull(Panel)
  
  
  
  if (length(Gene_panels) == 0) {
    print(paste(tmb_status,zhongzu,Gene_test, '没有一个中心有两种性别', sep = '_'))
    
    singlesex_panel_record_log <<- rbind(singlesex_panel_record_log, 
                                         data.frame(Gene_test = Gene_test,
                                                    stringsAsFactors = FALSE))
    
    return()
  }
  
  
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
    TMB_Status = tmb_status,
    Race = zhongzu
  )
  
  
  result_df <<- rbind(result_df, df_Gene_test)  
  
  
}



{
  #logstic_regression
  
  
  logstic_race_model <- function(race,source_EO,source_LO,tmb_status,Mutation_number =5){
    race_EO = subsetMaf(source_EO,tsb = source_EO@clinical.data[Race == race]$Tumor_Sample_Barcode)
    race_LO = subsetMaf(source_LO,tsb = source_LO@clinical.data[Race == race]$Tumor_Sample_Barcode)
    
    race_clinical_data = rbind(race_EO@clinical.data,race_LO@clinical.data)
    
    compare_genes_list <- uniquegenes(m1=race_EO,m2 = race_LO,minMut = Mutation_number)
    
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
    
    
    for (i in compare_genes_list) {
      result <- tryCatch({
        logstic_mafcompare_race(m1 = race_EO, m2 = race_LO, Clinical.data = race_clinical_data, Gene_test = i,
                                tmb_status = tmb_status, zhongzu = race)
        NULL  # 如果没有错误，返回NULL
      }, warning = function(w) {
        problematic_genes <<- c(problematic_genes, i)  # 记录产生警告的基因
        return(NULL)  # 返回NULL以继续循环
      })
    }
    
  }
  
  
  
  
  singlesex_panel_record_log <- data.frame(Gene_test = character(),
                                           stringsAsFactors = FALSE)
  
  result_df_race <- data.table::data.table(
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
    TMB_Status = NA,
    Race = NA,
    adjPval = NA
  )
  
  formula_base <- "Gene_status ~ Age_class+ TMB"
  
  for (i in c("Asian Pacific lslander","Black","White")){
    
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
      TMB_Status = NA,
      Race = NA
    )
    
    logstic_race_model(race = i,
                       source_EO=clin.EOCRC,
                       source_LO=clin.LOCRC,
                       tmb_status = 'Overall')
    
    result_df= na.omit(result_df)
    
    write_csv(result_df,file = paste('./table/2.2_logsitic',i,'original_TMBoverall.csv',sep = '_'))
    
    
    result_df_modified <- result_df %>% 
      filter(!Hugo_Symbol %in% problematic_genes) %>% 
      filter(formula == 'Gene_status ~ Age_class+ TMB + Panel + Country + Sex + Tumor_site + Sample_type') %>% 
      mutate(adjPval = p.adjust(p = pval, method = "fdr")) 
    
    write_csv(result_df_modified,file = paste('./table/2.2_logsitic',i,'modified_TMBoverall.csv',sep = '_'))
         
    result_df_race <-  rbind(result_df_race,result_df_modified)

    
  }
  
  
  
  result_df_race <- na.omit(result_df_race)
  
  write_csv(result_df_race,file = './table/2.2_logsitic_result_df_race.csv')
  
}
  
  


{#cochran Q test
    
    calculate_heterogeneity <- function(data) {
      
      betas <- data$beta
      ses <- data$se
      
      weights <- 1 / (ses^2)
      
      weighted_mean <- sum(weights * betas) / sum(weights)
      
      Q <- sum(weights * (betas - weighted_mean)^2)
      
      df <- length(betas) - 1
      
      p_value <- pchisq(Q, df, lower.tail = FALSE)
      
      return(data.frame(Q_stat = Q, df = df, Q_test_p_value = p_value))
    }
    
    
    #先筛出在三个种族中都有的基因，然后再筛出在至少1个种族中显著的基因
    #对这个基因收集各个种族的beta和se，做Q检验
    result_df_race_Qtest <- result_df_race %>%
      group_by(Hugo_Symbol) %>%
      filter(n_distinct(Race) == 3) %>%
      filter(any(adjPval < 0.05)) %>%
      do(calculate_heterogeneity(.)) %>%
      ungroup()
    
  }

  write.csv(result_df_race_Qtest,file = './table/2.2_result_df_race_Qtest.csv')

  {#chi-sequare test
  
  result_race_chi <-  result_df_race %>%
      mutate(EOCRC_not_mut_case = EOCRC_total_case- EOCRC) %>% 
      select(Hugo_Symbol, EOCRC_not_mut_case,EOCRC, Race) %>%   
      group_by(Hugo_Symbol) %>%
      filter(n_distinct(Race) == 3) %>%
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
        result_df_chi <- fisher.test(contingency_table)
        data.frame(Hugo_Symbol = unique(.x$Hugo_Symbol),
                   p.value = result_df_chi$p.value,
                   method = "fisher")
      })
    })
    
    write.csv(result_race_chi,file = './table/2.2_result_race_chi.csv')
    
  }
  


{#数据整合
  
  Asia_data <- read.csv('./table/2.2_logsitic_Asian Pacific lslander_modified_TMBoverall.csv') %>%
    select(Hugo_Symbol, EOCRC_freq, LOCRC_freq, adjPval, OR, ci.up, ci.low) %>%
    rename_with(~ paste("Asia", ., sep = "_"), -Hugo_Symbol)
  
  Black_data <- read.csv('./table/2.2_logsitic_Black_modified_TMBoverall.csv')%>% 
    select(Hugo_Symbol,EOCRC_freq,LOCRC_freq,adjPval,OR,ci.up,ci.low)%>%
    rename_with(~ paste("Blcak", ., sep = "_"), -Hugo_Symbol)
  
  White_data <- read.csv('./table/2.2_logsitic_White_modified_TMBoverall.csv')%>% 
    select(Hugo_Symbol,EOCRC_freq,LOCRC_freq,adjPval,OR,ci.up,ci.low)%>%
    rename_with(~ paste("White", ., sep = "_"), -Hugo_Symbol)
  
  
  Race_data_all <- full_join(
    Asia_data,Black_data,by='Hugo_Symbol'  ) %>% 
    full_join(White_data,by='Hugo_Symbol') %>% 
    left_join(result_race_chi,by='Hugo_Symbol') %>%
    left_join(result_df_race_Qtest,by='Hugo_Symbol')%>%
    select(-df) %>% 
    arrange(Q_test_p_value)
    
  write_csv(Race_data_all,file = './table/2.2_Race_data_all.csv')
}
    


