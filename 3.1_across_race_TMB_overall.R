library(tidyverse)
library(patchwork)

logstic_mafcompare_race <- function(m1, m2, Clinical.data,Gene_test,tmb_status,zhongzu) {
  
  Alter_m1_Gene = m1@data[Hugo_Symbol == Gene_test,] %>% 
                  filter(!duplicated(Tumor_Sample_Barcode))
  
  Alter_m2_Gene = m2@data[Hugo_Symbol == Gene_test,] %>% 
                  filter(!duplicated(Tumor_Sample_Barcode))
  
  Gene_panels <- Clinical.data[, c('Tumor_Sample_Barcode', 'Age_class', 'Sex', 'Country','Panel', 'Race')] %>% 
                  filter(Panel %in% gene_panel_map[[Gene_test]]) %>% 
                  pull(Panel)
  
  mut_samples_m1 <- length(unique(Alter_m1_Gene$Tumor_Sample_Barcode))

  panels_samples_m1 <- length(unique(m1@clinical.data[m1@clinical.data$Panel %in% Gene_panels, ]$Tumor_Sample_Barcode))
  
  m1_freq <- mut_samples_m1 / panels_samples_m1
  
  # ------------
  mut_samples_m2 <- length(unique(Alter_m2_Gene$Tumor_Sample_Barcode))

  panels_samples_m2 <- length(unique(m2@clinical.data[m2@clinical.data$Panel %in% Gene_panels, ]$Tumor_Sample_Barcode))

  m2_freq <- mut_samples_m2 / panels_samples_m2
  
  if ( m1_freq < 0.01 | m2_freq<0.01 ) {
    print(paste(Gene_test,tmb_status,zhongzu,
                "Low mutation frequencies, returning from function.",
                sep = '____'))
    return()
  }
  
  Gene_clinical_info  <-  Clinical.data[, Gene_status := ifelse(Tumor_Sample_Barcode %in% c(Alter_m1_Gene$Tumor_Sample_Barcode,
                                                                                            Alter_m2_Gene$Tumor_Sample_Barcode), 1, 0)]%>% 
                          mutate(Age_class = factor(Age_class,levels = c('LOCRC','EOCRC'))) %>% 
                          filter(Panel %in% Gene_panels) 

  each_gene_mutation_counts = gene_mutation_counts %>% 
                              dplyr::select(all_of(Gene_test), Tumor_Sample_Barcode,Age_class, TMB,Race) %>%
                              filter(Tumor_Sample_Barcode %in% Gene_clinical_info$Tumor_Sample_Barcode)  %>% 
                              mutate(mutation_rate = .data[[Gene_test]] / TMB)
  
  each_gene_mutation_EO = each_gene_mutation_counts %>% 
                           filter(Age_class == 'EOCRC') %>% pull(mutation_rate)
  each_gene_mutation_LO = each_gene_mutation_counts%>% 
                           filter(Age_class == 'LOCRC')%>% pull(mutation_rate)
  
  EOCRC_rate = mean(each_gene_mutation_EO)
  LOCRC_rate = mean(each_gene_mutation_LO)

  mutation_rate_pvalue <- wilcox.test(each_gene_mutation_EO, each_gene_mutation_LO)$p.value
  
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
                  EOCRC_mutation_rate = EOCRC_rate,
                  LOCRC_mutation_rate = LOCRC_rate,
                  mutation_rate_pval = mutation_rate_pvalue,
                  formula = formula_base,
                  TMB_Status = tmb_status,
                  Race = zhongzu  )

  result_df <<- rbind(result_df, df_Gene_test)  

}

{    #logstic_regression

logstic_race_model <- function(race,source_EO,source_LO,tmb_status,Mutation_number =5){
    
    race_EO = subsetMaf(source_EO,tsb = source_EO@clinical.data[Race == race]$Tumor_Sample_Barcode)
    race_LO = subsetMaf(source_LO,tsb = source_LO@clinical.data[Race == race]$Tumor_Sample_Barcode)
    
    race_clinical_data = rbind(race_EO@clinical.data,race_LO@clinical.data)
    
    compare_genes_list <- uniquegenes(m1=race_EO,m2 = race_LO,minMut = Mutation_number)
    
    problematic_genes <- c()
    
    for (i in compare_genes_list) {
      tryCatch({
        logstic_mafcompare_race(m1 = race_EO, m2 = race_LO, 
                                Clinical.data = race_clinical_data, Gene_test = i,
                                tmb_status = tmb_status, zhongzu = race)
        NULL 
      }, warning = function(w) {
        problematic_genes <- c(problematic_genes, i)  
        return(NULL)
      })
      
    }

    result_df_modified <- result_df[Race == race] %>% 
                          filter(!Hugo_Symbol %in% problematic_genes) %>% 
                          filter(formula == 'Gene_status ~ Age_class + Panel + Country + Sex + Tumor_site + Sample_type + Histology_type') %>% 
                          mutate(adjPval = p.adjust(p = pval, method = "fdr"),
                                 mutation_rate_adjPval =  p.adjust(p = mutation_rate_pval, method = "fdr"))  
    
    assign(paste('result_df_modified', race, tmb_status,sep = '_'), 
           result_df_modified, envir = .GlobalEnv)
    
    write_csv(result_df_modified,file = paste('./table/3.1_logsitic_modified',
                                              race,tmb_status,
                                              '.csv',sep = '_'))
}
 
  result_df_race_template <- data.table::data.table(
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
                              EOCRC_mutation_rate = numeric(),
                              LOCRC_mutation_rate = numeric(),
                              mutation_rate_pval = numeric(),
                              formula =numeric(),
                              TMB_Status = character(),
                              Race = character() )

  result_df <- result_df_race_template
  
  formula_base <- "Gene_status ~ Age_class"
  
  logstic_race_model(race = "Asian Pacific lslander", source_EO=clin.EOCRC, source_LO=clin.LOCRC,
                      tmb_status = 'Overall',Mutation_number =5  )
  
  logstic_race_model(race = "White", source_EO=clin.EOCRC, source_LO=clin.LOCRC,
                      tmb_status = 'Overall',Mutation_number =5  )
     
  logstic_race_model(race = "Black", source_EO=clin.EOCRC, source_LO=clin.LOCRC,
                      tmb_status = 'Overall',Mutation_number =5 )
   
  result_df_race <- rbind(result_df_modified_Black_Overall,
                          result_df_modified_White_Overall,
                          `result_df_modified_Asian Pacific lslander_Overall`)
}
  
{    #cochran Q test
    
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
    
  result_df_race_Qtest <- result_df_race %>%
                          group_by(Hugo_Symbol) %>%
                          filter(n_distinct(Race) == 3) %>%
                          filter(any(adjPval < 0.05)) %>%
                          do(calculate_heterogeneity(.)) %>%
                          ungroup()
    
}

  write.csv(result_df_race_Qtest,file = './table/3.1_result_df_race_Qtest.csv')

{    #chi-sequare test
  
    result_race_chi <-  result_df_race %>%
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
    
    write.csv(result_race_chi,file = './table/3.1_result_race_chi.csv')
 
}
  
{
  Asia_data <- `result_df_modified_Asian Pacific lslander_Overall` %>%
                dplyr::select(Hugo_Symbol, EOCRC_freq, LOCRC_freq, 
                              adjPval, OR, ci.up, ci.low,
                              EOCRC_mutation_rate,LOCRC_mutation_rate,
                              mutation_rate_adjPval) %>%
                rename_with(~ paste("Asia", ., sep = "_"), -Hugo_Symbol)
  
  Black_data <- result_df_modified_Black_Overall%>% 
                dplyr::select(Hugo_Symbol, EOCRC_freq, LOCRC_freq, 
                              adjPval, OR, ci.up, ci.low,
                              EOCRC_mutation_rate,LOCRC_mutation_rate,
                              mutation_rate_adjPval)%>%
                rename_with(~ paste("Blcak", ., sep = "_"), -Hugo_Symbol)
  
  White_data <- result_df_modified_White_Overall%>% 
                dplyr::select(Hugo_Symbol, EOCRC_freq, LOCRC_freq, 
                              adjPval, OR, ci.up, ci.low,
                              EOCRC_mutation_rate,LOCRC_mutation_rate,
                              mutation_rate_adjPval)%>%
                rename_with(~ paste("White", ., sep = "_"), -Hugo_Symbol)

  Race_data_all <- full_join( Asia_data,Black_data,by='Hugo_Symbol'  ) %>% 
                              full_join(White_data,by='Hugo_Symbol') %>% 
                              left_join(result_race_chi,by='Hugo_Symbol') %>%
                              left_join(result_df_race_Qtest,by='Hugo_Symbol')%>%
                              dplyr::select(-df) %>% 
                              arrange(Q_test_p_value)
    
  write_csv(Race_data_all,file = './table/3.1_Race_data_all.csv')
}
