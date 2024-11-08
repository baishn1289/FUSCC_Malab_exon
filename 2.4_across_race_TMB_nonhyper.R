#across_race_nonhypermutated_CRC


singlesex_panel_record_log <- data.frame(Gene_test = character(),
                                         stringsAsFactors = FALSE)

result_df_race_nonhyper <- data.table::data.table(
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
                     source_EO=clin.EOCRC_nonhyper,
                     source_LO=clin.LOCRC_nonhyper,
                     tmb_status = 'Nonhyper',
                     Mutation_number = 2
  )
  
  result_df= na.omit(result_df)
  
  write_csv(result_df,file = paste('./table/2.2_logsitic',i,'original_TMBnonhyper.csv',sep = '_'))
  
  
  result_df_modified <- result_df %>% 
    filter(!Hugo_Symbol %in% problematic_genes) %>% 
    filter(formula == 'Gene_status ~ Age_class+ TMB + Panel + Country + Sex + Tumor_site + Sample_type') %>% 
    mutate(adjPval = p.adjust(p = pval, method = "fdr")) 
  
  write_csv(result_df_modified,file = paste('./table/2.2_logsitic',i,'modified_TMBnonhyper.csv',sep = '_'))
  
  result_df_race_nonhyper <-  rbind(result_df_race_nonhyper,result_df_modified)
}

result_df_race_nonhyper <- na.omit(result_df_race_nonhyper)

write_csv(result_df_race_nonhyper,file = './table/2.2_logsitic_result_df_race_TMB_nonhyper.csv')


#####
result_df_race_Qtest_nonhyper <- result_df_race_nonhyper %>%
  group_by(Hugo_Symbol) %>%
  filter(n_distinct(Race) == 3) %>%
  filter(any(adjPval < 0.05)) %>%
  do(calculate_heterogeneity(.)) %>%
  ungroup()


write.csv(result_df_race_Qtest,file = './table/2.2_result_df_race_Qtest_nonhyper.csv')

####
result_race_chi_nonhyper <-  result_df_race_nonhyper %>%
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

write.csv(result_race_chi,file = './table/2.2_result_race_chi_nonhyper.csv')






{#数据整合
  
  Asia_data_nonhyper <- read.csv('./table/2.2_logsitic_Asian Pacific lslander_modified_TMBnonhyper.csv') %>%
    select(Hugo_Symbol, EOCRC_freq, LOCRC_freq, adjPval, OR, ci.up, ci.low) %>%
    rename_with(~ paste("Asia", ., sep = "_"), -Hugo_Symbol)
  
  Black_data_nonhyper <- read.csv('./table/2.2_logsitic_Black_modified_TMBnonhyper.csv')%>% 
    select(Hugo_Symbol,EOCRC_freq,LOCRC_freq,adjPval,OR,ci.up,ci.low)%>%
    rename_with(~ paste("Blcak", ., sep = "_"), -Hugo_Symbol)
  
  White_data_nonhyper <- read.csv('./table/2.2_logsitic_White_modified_TMBnonhyper.csv')%>% 
    select(Hugo_Symbol,EOCRC_freq,LOCRC_freq,adjPval,OR,ci.up,ci.low)%>%
    rename_with(~ paste("White", ., sep = "_"), -Hugo_Symbol)
  
  
  Race_data_all_nonhyper <- full_join(
    Asia_data_nonhyper,Black_data_nonhyper,by='Hugo_Symbol'  ) %>% 
    full_join(White_data_nonhyper,by='Hugo_Symbol') %>% 
    left_join(result_race_chi_nonhyper,by='Hugo_Symbol') %>%
    left_join(result_df_race_Qtest_nonhyper,by='Hugo_Symbol')%>%
    select(-df) %>% 
    arrange(Q_test_p_value)
  
  write_csv(Race_data_all_nonhyper,file = './table/2.2_Race_data_all_nonhyper.csv')
}


