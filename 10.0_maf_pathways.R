
pws_EO = pathways(maf = clin.EOCRC,plotType = 'bar', col = "#f39c12",
                  pathdb='smgbp')
pws_LO= pathways(maf=clin.LOCRC,plotType = 'bar', col = "#f39c12",
                 pathdb='smgbp')


CairoPDF(file =paste0(dir,'picture/11_maf_EO_pathways.pdf'),
         height =5,width=10)
pathways(maf = clin.EOCRC,plotType = 'bar', col = "#f39c12",
         pathdb='smgbp')
dev.off()

CairoPDF(file = paste0(dir,'picture/11_maf_LO_pathways.pdf'),
         height =5,width=10)
pathways(maf=clin.LOCRC,plotType = 'bar', col = "#f39c12",
         pathdb='smgbp')
dev.off()


compare_pathway_fisher <- function(pws1, pws2) {
  
  results <- list()
  
  for (pathway in intersect(pws1$Pathway, pws2$Pathway)) {
    
    mutated1 <- pws1$Mutated_samples[pws1$Pathway == pathway]
    total1 <- nrow(clin.EOCRC@clinical.data)
    
    mutated2 <- pws2$Mutated_samples[pws2$Pathway == pathway]
    total2 <- nrow(clin.LOCRC@clinical.data)
    
    
    contingency_table <- matrix(c(mutated1, total1 - mutated1, 
                                  mutated2, total2 - mutated2), 
                                nrow = 2, 
                                byrow = TRUE)
    
    
    expected_counts <- chisq.test(contingency_table)$expected
    
    
    if (any(expected_counts < 5)) {
      result_df_chi <- fisher.test(contingency_table, workspace = 2e7)
      method_used <- "fisher"
    } else {
      result_df_chi <- chisq.test(contingency_table)
      method_used <- "chi-square"
    }
    
    
    cal_result <- data.frame(
      Pathway_name = pathway,
      p.value = result_df_chi$p.value,
      method = method_used,
      EO_pathway_mut = mutated1,
      EO_pathway_all = total1,
      LO_pathway_mut = mutated2,
      LO_pathway_all = total2,
      EO_freq = mutated1 / total1,
      LO_freq = mutated2 / total2
    )
    
    
    results[[length(results) + 1]] <- cal_result
  }
  
  
  results_df <- do.call(rbind, results)
  return(results_df)
}


pathways_compare_results <- compare_pathway_fisher(pws_EO, pws_LO) %>% 
  mutate(adjPval = p.adjust(p = p.value, method = "BH"))

pathways_compare_results_publish <- pathways_compare_results %>% 
  mutate(across(c(p.value, adjPval), ~ format_Pval(.)))

write.csv(pathways_compare_results,
          file=paste0(dir,"table/11.pathways_compare_results.csv"))
write.csv(pathways_compare_results_publish,
          file=paste0(dir,"table/11.publish_pathways_compare_results.csv"))
