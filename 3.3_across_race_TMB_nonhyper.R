#across_race_nonhypermutated_CRC

  result_df <- result_df_race_template
  
  formula_base <- "Gene_status ~ Age_class"
  
  logstic_race_model(
    race = "Asian Pacific lslander", source_EO=clin.EOCRC_nonhyper, source_LO=clin.LOCRC_nonhyper,
    tmb_status = 'Nonhyper',Mutation_number = 5
  )
  
  logstic_race_model(
    race = "White", source_EO=clin.EOCRC_nonhyper, source_LO=clin.LOCRC_nonhyper,
    tmb_status = 'Nonhyper',Mutation_number =5
  )
  
  logstic_race_model(
    race = "Black", source_EO=clin.EOCRC_nonhyper, source_LO=clin.LOCRC_nonhyper,
    tmb_status = 'Nonhyper',Mutation_number =5
  )


  result_df_race_nonhyper <- rbind(result_df_modified_Black_Nonhyper,
                                  result_df_modified_White_Nonhyper,
                                  `result_df_modified_Asian Pacific lslander_Nonhyper`)

####Q-test
  result_df_race_Qtest_nonhyper <- result_df_race_nonhyper %>%
                                    group_by(Hugo_Symbol) %>%
                                    filter(n_distinct(Race) == 3) %>%
                                    filter(any(adjPval < 0.05)) %>%
                                    do(calculate_heterogeneity(.)) %>%
                                    ungroup()
  
  write.csv(result_df_race_Qtest_nonhyper,file = './table/3.3_result_df_race_Qtest_nonhyper.csv')

####chi-seq_test
  result_race_chi_nonhyper <-  result_df_race_nonhyper %>%
                                mutate(EOCRC_not_mut_case = EOCRC_total_case- EOCRC) %>% 
                                dplyr::select(Hugo_Symbol, EOCRC_not_mut_case,EOCRC, Race) %>%   
                                group_by(Hugo_Symbol) %>%
                                filter(n_distinct(Race) == 3) %>%
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
                                              result_df_chi <- fisher.test(contingency_table)
                                              data.frame(Hugo_Symbol = unique(.x$Hugo_Symbol),
                                                         p.value = result_df_chi$p.value,
                                                         method = "fisher")
                                            })
                                          })

  write.csv(result_race_chi_nonhyper,file = './table/3.3_result_race_chi_nonhyper.csv')

{
  Asia_data_nonhyper <- `result_df_modified_Asian Pacific lslander_Nonhyper` %>%
                        dplyr::select(Hugo_Symbol, EOCRC_freq, LOCRC_freq, 
                                      adjPval, OR, ci.up, ci.low,
                                      EOCRC_mutation_rate,LOCRC_mutation_rate,
                                      mutation_rate_pval) %>%
                        rename_with(~ paste("Asia", ., sep = "_"), -Hugo_Symbol)
  
  Black_data_nonhyper <- result_df_modified_Black_Nonhyper %>% 
                          dplyr::select(Hugo_Symbol, EOCRC_freq, LOCRC_freq, 
                                        adjPval, OR, ci.up, ci.low,
                                        EOCRC_mutation_rate,LOCRC_mutation_rate,
                                        mutation_rate_pval)%>%
                          rename_with(~ paste("Blcak", ., sep = "_"), -Hugo_Symbol)
                        
  White_data_nonhyper <- result_df_modified_White_Nonhyper %>% 
                          dplyr::select(Hugo_Symbol, EOCRC_freq, LOCRC_freq, 
                                        adjPval, OR, ci.up, ci.low,
                                        EOCRC_mutation_rate,LOCRC_mutation_rate,
                                        mutation_rate_pval)%>%
                          rename_with(~ paste("White", ., sep = "_"), -Hugo_Symbol)

  Race_data_all_nonhyper <- full_join(  Asia_data_nonhyper,Black_data_nonhyper,by='Hugo_Symbol'  ) %>% 
                                        full_join(White_data_nonhyper,by='Hugo_Symbol') %>% 
                                        left_join(result_race_chi_nonhyper,by='Hugo_Symbol') %>%
                                        left_join(result_df_race_Qtest_nonhyper,by='Hugo_Symbol')%>%
                                        dplyr::select(-df) %>% 
                                        arrange(Q_test_p_value)
  
  write_csv(Race_data_all_nonhyper,file = './table/3.3_Race_data_all_nonhyper.csv')
}