
race_compare_results_nonhyper <- lapply(races_names, function(r) {
  logstic_race_model(
    race = r,
    source_EO = clin.EOCRC_nonhyper,
    source_LO = clin.LOCRC_nonhyper,
    tmb_status = 'Non-hyper',
    Mutation_number = 5
  )
})

result_df_race_nonhyper <- rbindlist(race_compare_results_nonhyper)

write.csv(result_df_race_nonhyper,
          file = paste0(dir,'table/3.3_race_nonhyper_original_data.csv'))

result_df_race_Qtest_nonhyper <- result_df_race_nonhyper %>%
  group_by(Hugo_Symbol) %>%
  filter(n_distinct(Race) == 3) %>%
  filter(any(adjPval < 0.05)) %>%
  do(calculate_heterogeneity(.)) %>%
  ungroup()

write.csv(result_df_race_Qtest_nonhyper,
          file = paste0(dir,'table/3.3_result_df_race_Qtest_nonhyper.csv'))

result_race_chi_nonhyper <-  result_df_race_nonhyper %>%
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

write.csv(result_race_chi_nonhyper,
          file = paste0(dir,'table/3.3_result_race_chi_nonhyper.csv'))



{
  
  Asia_data_nonhyper <- result_df_race_nonhyper[Race == "Asian Pacific lslander"] %>%
    dplyr::select(Hugo_Symbol, EOCRC_freq, LOCRC_freq, 
                  adjPval, OR, `95% upper limit`, `95% lower limit`,
                  EOCRC_mutation_rate,LOCRC_mutation_rate,
                  mutation_rate_adjPval) %>%
    rename_with(~ paste("Asia", ., sep = "_"), -Hugo_Symbol)
  
  Black_data_nonhyper <- result_df_race_nonhyper[Race == "Black"] %>% 
    dplyr::select(Hugo_Symbol, EOCRC_freq, LOCRC_freq, 
                  adjPval, OR, `95% upper limit`, `95% lower limit`,
                  EOCRC_mutation_rate,LOCRC_mutation_rate,
                  mutation_rate_adjPval) %>%
    rename_with(~ paste("Black", ., sep = "_"), -Hugo_Symbol)
  
  White_data_nonhyper <- result_df_race_nonhyper[Race == "White"] %>% 
    dplyr::select(Hugo_Symbol, EOCRC_freq, LOCRC_freq, 
                  adjPval, OR, `95% upper limit`, `95% lower limit`,
                  EOCRC_mutation_rate,LOCRC_mutation_rate,
                  mutation_rate_adjPval) %>%
    rename_with(~ paste("White", ., sep = "_"), -Hugo_Symbol)
  
  
  Race_data_all_nonhyper <- full_join(
    Asia_data_nonhyper,Black_data_nonhyper,by='Hugo_Symbol'  ) %>% 
    full_join(White_data_nonhyper,by='Hugo_Symbol') %>% 
    left_join(result_race_chi_nonhyper,by='Hugo_Symbol') %>%
    left_join(result_df_race_Qtest_nonhyper,by='Hugo_Symbol')%>%
    dplyr::select(-df) %>% 
    mutate(Q_test_adjpvalue = p.adjust(p = Q_test_p_value, method = "BH")%>% 
             format_Pval(),
           Pvalue = p.adjust(p = p.value, method = "BH")%>% 
             format_Pval()
    ) %>% 
    arrange(Q_test_p_value) 
  
  Race_data_all_nonhyper_publish <- Race_data_all_nonhyper[Q_test_adjpvalue < 0.05]%>% 
    select(-p.value,-Q_test_p_value,-method)
  
  write_csv(Race_data_all_nonhyper,
            file = paste0(dir,'table/3.3_Race_data_all_nonhyper.csv'))
  write_csv(Race_data_all_nonhyper_publish,
            file = paste0(dir,'table/3.3_Race_data_all_nonhyper_publish.csv'))
  
}


