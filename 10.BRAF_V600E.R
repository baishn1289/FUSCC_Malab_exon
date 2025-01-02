library(tidyverse)
library(patchwork)
library(foreach)
library(data.table)

logstic_mafcompare_hotspot <- function(m1, m2, Clinical.data, Gene_test, aachange, TMB_Status) {
  
  m1_AA_changed_samples = m1@data[Hugo_Symbol == Gene_test & HGVSp_Short == aachange,]  %>% 
                          filter(!duplicated(Tumor_Sample_Barcode))
   
  m2_AA_changed_samples = m2@data[Hugo_Symbol == Gene_test & HGVSp_Short == aachange,] %>% 
                          filter(!duplicated(Tumor_Sample_Barcode))
  
  Gene_panels <- Clinical.data[, c('Tumor_Sample_Barcode', 'Age_class', 'Sex', 'Country','Panel', 'Race')] %>% 
                  filter(Panel %in% gene_panel_map[[Gene_test]]) %>% 
                  pull(Panel)
  
  panels_samples_m1 <- length(unique(m1@clinical.data[m1@clinical.data$Panel %in% Gene_panels, ]$Tumor_Sample_Barcode)) 
  m1_freq <- nrow(m1_AA_changed_samples) / panels_samples_m1
  
  # ------------

  panels_samples_m2 <- length(unique(m2@clinical.data[m2@clinical.data$Panel %in% Gene_panels, ]$Tumor_Sample_Barcode))
  m2_freq <- nrow(m2_AA_changed_samples) / panels_samples_m2
 
  Gene_clinical_info  <-  Clinical.data[, Gene_status := ifelse(Tumor_Sample_Barcode %in% c(m1_AA_changed_samples$Tumor_Sample_Barcode,
                                                                                            m2_AA_changed_samples$Tumor_Sample_Barcode), 1, 0)]%>% 
                            mutate(Age_class = factor(Age_class,levels = c('LOCRC','EOCRC'))) %>% 
                            filter(Panel %in% Gene_panels)
 
  each_gene_mutation_counts = gene_mutation_counts %>% 
                              dplyr::select(Panel,Tumor_Sample_Barcode,Age_class, TMB) %>% 
                              filter(Tumor_Sample_Barcode %in% Gene_clinical_info$Tumor_Sample_Barcode) %>%
                              mutate(hotspot = ifelse(Tumor_Sample_Barcode %in% c(m1_AA_changed_samples$Tumor_Sample_Barcode,
                                                                                  m2_AA_changed_samples$Tumor_Sample_Barcode), 1, 0)) %>% 
                              mutate(mutation_rate = hotspot / TMB)
  
  each_gene_mutation_EO = each_gene_mutation_counts %>% 
                           filter(Age_class == 'EOCRC') %>% 
                           pull(mutation_rate)
  each_gene_mutation_LO = each_gene_mutation_counts%>% 
                           filter(Age_class == 'LOCRC')%>% 
                            pull(mutation_rate)
  
  EOCRC_rate = mean(each_gene_mutation_EO)
  LOCRC_rate = mean(each_gene_mutation_LO)
  
  mutation_rate_pvalue <- wilcox.test(each_gene_mutation_EO, each_gene_mutation_LO)$p.value
  
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
                    Hotspot = aachange,
                    EOCRC_total_case= panels_samples_m1,
                    LOCRC_total_case= panels_samples_m2,
                    EOCRC = nrow(m1_AA_changed_samples), 
                    LOCRC = nrow(m2_AA_changed_samples), 
                    EOCRC_freq = m1_freq,
                    LOCRC_freq = m2_freq,
                    pval = pval_age_class, 
                    or = OR_age_class,
                    ci.up = ci_up,
                    ci.low = ci_low,
                    EOCRC_mutation_rate = EOCRC_rate,
                    LOCRC_mutation_rate = LOCRC_rate,
                    mutation_rate_pval = mutation_rate_pvalue,
                    formula = formula_base,
                    TMB_status = TMB_Status
                    )
  
  result_df <<- rbind(result_df, df_Gene_test)  
}

  result_df <- data.table::data.table(
                Hugo_Symbol = character(),
                Hotspot = character(),
                EOCRC_total_case= numeric(),
                LOCRC_total_case= numeric(),
                EOCRC = numeric(), 
                LOCRC = numeric(), 
                EOCRC_freq = numeric(),
                LOCRC_freq = numeric(),
                pval = numeric(), 
                or = numeric(),
                ci.up = numeric(),
                ci.low = numeric(),
                EOCRC_mutation_rate = numeric(),
                LOCRC_mutation_rate = numeric(),
                mutation_rate_pval = numeric(),
                formula = character(),
                TMB_status = character()
                 )

  formula_base <- "Gene_status ~ Age_class"
  
  logstic_mafcompare_hotspot(m1 = clin.EOCRC, m2 = clin.LOCRC, Clinical.data = rt@clinical.data,
                             Gene_test = 'BRAF', aachange = 'p.V600E', TMB_Status = 'Overall')
  
  logstic_mafcompare_hotspot(m1 = clin.EOCRC_hyper, m2 = clin.LOCRC_hyper, 
                             Clinical.data = rbind(clin.EOCRC_hyper@clinical.data,clin.LOCRC_hyper@clinical.data),
                             Gene_test = 'BRAF', aachange = 'p.V600E', TMB_Status = 'Hypermutated')
  
  logstic_mafcompare_hotspot(m1 = clin.EOCRC_nonhyper, m2 = clin.LOCRC_nonhyper,
                             Clinical.data = rbind(clin.EOCRC_nonhyper@clinical.data,clin.LOCRC_nonhyper@clinical.data),
                             Gene_test = 'BRAF', aachange = 'p.V600E', TMB_Status = 'Non-hypermutated')
  
  
  
  logstic_mafcompare_hotspot(m1 = clin.EOCRC, m2 = clin.LOCRC, Clinical.data = rt@clinical.data,
                             Gene_test = 'KRAS', aachange = 'p.G12D', TMB_Status = 'Overall')

  logstic_mafcompare_hotspot(m1 = clin.EOCRC_hyper, m2 = clin.LOCRC_hyper, 
                              Clinical.data = rbind(clin.EOCRC_hyper@clinical.data,clin.LOCRC_hyper@clinical.data),
                              Gene_test = 'KRAS', aachange = 'p.G12D', TMB_Status = 'Hypermutated')

  logstic_mafcompare_hotspot(m1 = clin.EOCRC_nonhyper, m2 = clin.LOCRC_nonhyper, 
                           Clinical.data = rbind(clin.EOCRC_nonhyper@clinical.data,clin.LOCRC_nonhyper@clinical.data),
                             Gene_test = 'KRAS', aachange = 'p.G12D', TMB_Status = 'Non-hypermutated')

  
  
  logstic_mafcompare_hotspot(m1 = clin.EOCRC, m2 = clin.LOCRC, Clinical.data = rt@clinical.data,
                             Gene_test = 'KRAS', aachange = 'p.G12C', TMB_Status = 'Overall')
  
  logstic_mafcompare_hotspot(m1 = clin.EOCRC_hyper, m2 = clin.LOCRC_hyper, 
                             Clinical.data = rbind(clin.EOCRC_hyper@clinical.data,clin.LOCRC_hyper@clinical.data),
                             Gene_test = 'KRAS', aachange = 'p.G12C', TMB_Status = 'Hypermutated')
  
  logstic_mafcompare_hotspot(m1 = clin.EOCRC_nonhyper, m2 = clin.LOCRC_nonhyper, 
                             Clinical.data = rbind(clin.EOCRC_nonhyper@clinical.data,clin.LOCRC_nonhyper@clinical.data),
                             Gene_test = 'KRAS', aachange = 'p.G12C', TMB_Status = 'Non-hypermutated')

  
  result_hotspot <- result_df %>% 
                    group_by(TMB_status) %>% 
                    filter(EOCRC >0 & LOCRC >0) %>% 
                    mutate(adjPval = p.adjust(p = pval, method = "fdr"),
                           mutation_rate_adjPval =  p.adjust(p = mutation_rate_pval, method = "fdr")) %>% 
                    ungroup()  %>% 
                    mutate(Hugo_Symbol = case_when(
                      TMB_status == 'Non-hypermutated' ~ paste(Hugo_Symbol, ' '),  
                      TMB_status == 'Hypermutated' ~ paste(Hugo_Symbol, '  '),
                      TRUE ~ Hugo_Symbol 
                    )) %>% 
                    mutate( Hugo_Symbol = case_when(
                      Hotspot == 'p.G12C' ~ paste(Hugo_Symbol, '    '),  
                      # Hotspot == 'p.G12C' ~ paste(Hugo_Symbol, '  '),
                      TRUE ~ Hugo_Symbol 
                              )
                    )
  
  write.csv(result_df,'./table/10.hotspot_original_data.csv')
  write.csv(result_hotspot,'./table/10.hotspot.csv')

  plot <- result_hotspot %>% 
            mutate(EOCRC_freq = format_frequency(EOCRC_freq),
                   LOCRC_freq = format_frequency(LOCRC_freq),
                   sample_size = paste0(LOCRC_freq, " v.s. ", EOCRC_freq),
                   OR_CI = paste0(format(round(or, 2), nsmall = 2), " (", 
                                  format(round(ci.low, 2), nsmall = 2), " to ", 
                                  format(round(ci.up, 2), nsmall = 2), ")"),
                   adjPval = format_adjPval(adjPval),
                   EOCRC_mutation_rate = format_frequency(EOCRC_mutation_rate),
                   LOCRC_mutation_rate = format_frequency(LOCRC_mutation_rate),
                   mutation_rate = paste0(LOCRC_mutation_rate, " v.s. ", EOCRC_mutation_rate),
                   mutation_rate_adjPval = format_adjPval(mutation_rate_adjPval)  ) %>% 
            mutate(across(everything(), as.character)) %>% 
            bind_rows(data.frame(Hugo_Symbol = "Gene", 
                                 Hotspot = "Hotspot",
                                 TMB_status = "TMB_status",
                                 sample_size = "Frequency\n(LOCRC v.s. EOCRC)", 
                                 adjPval = "adjPval", 
                                 OR_CI = "OR (95% CI)",
                                 mutation_rate = 'mutation_rate\n(LOCRC v.s. EOCRC)',
                                 mutation_rate_adjPval = 'mutation_rate_adjPval')) %>% 
            mutate(Hugo_Symbol = fct_rev(fct_relevel(Hugo_Symbol, "Gene")))


  p_mid <- result_hotspot %>% 
            ggplot(aes(y = fct_rev(Hugo_Symbol))) +
            theme_classic() +
            geom_point(aes(x = or, size = or), shape = 23, fill = "black", show.legend = F) +
            geom_errorbarh(aes(xmin = ci.low, xmax = ci.up), height = 0.3) +
            labs(x = "Odds Ratio") +
            scale_x_continuous(trans = 'log10',breaks = c(0.01,0.1,1,10),
                               labels = scales::number_format(accuracy = 0.01)) +
            coord_cartesian(ylim = c(1,10), xlim = c(0.01,10)) +
            geom_vline(xintercept = 1, linetype = 2) +  
            theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                  axis.text.y = element_blank(), axis.title.y = element_blank(),
                  axis.text.x = element_text(color = "black"))

  
  p_left <- plot %>%
            ggplot(aes(y = Hugo_Symbol)) +
            geom_text(aes(x = 0, label = Hugo_Symbol), hjust = 0, fontface = "bold", size = 4) +
            geom_text(aes(x = 0.8, label = Hotspot), hjust = 0, fontface = "bold", size = 4) +
            geom_text(aes(x = 1.6, label = TMB_status), hjust = 0, fontface = "bold", size = 4) +
            geom_text(aes(x = 3, label = sample_size), hjust = 0.5, size = 4,vjust=0.5,
                      fontface = ifelse(plot$sample_size == "Frequency\n(LOCRC v.s. EOCRC)", "bold", "plain")) +
            theme_void() + coord_cartesian(xlim = c(0, 3.5))


  p_right <- plot %>% 
            ggplot() +
            geom_text(aes(x = 0, y = Hugo_Symbol, label = OR_CI), hjust = 0,
                      fontface = ifelse(plot$OR_CI == "OR (95% CI)", "bold", "plain"), size = 4) +
            geom_text(aes(x = 1.1, y = Hugo_Symbol, label = adjPval), hjust = 0,
                      fontface = ifelse(plot$adjPval == "adjPval", "bold", "plain"), size = 4) +
            geom_text(aes(x = 2.3, y = Hugo_Symbol, label = mutation_rate), hjust = 0.5,
                      fontface = ifelse(plot$mutation_rate == 'mutation_rate\n(LOCRC v.s. EOCRC)', "bold", "plain"), size = 4) +
            geom_text(aes(x = 3, y = Hugo_Symbol, label = mutation_rate_adjPval), hjust = 0,
                      fontface = ifelse(plot$mutation_rate_adjPval == "mutation_rate_adjPval", "bold", "plain"), size = 4) +
            theme_void() + coord_cartesian(xlim = c(0, 4))


  layout <- c(
              patchwork::area(t = 0, l = 0, b = 30, r = 10), 
              patchwork::area(t = 0, l = 11, b = 30, r = 16), 
              patchwork::area(t = 0, l = 17, b = 30, r = 28))

  p_hyper_hotspot =p_left + p_mid + p_right + plot_layout(design = layout)
  
  pdf(paste(dir,'/10.hotspot_forestplot.pdf',sep = ""), 
      width = 16, height = 5)
  p_hyper_hotspot
  dev.off()


######BRAF V600E proportion calculation
  proportion_calculation_df <- data.frame(
    TMB_Status = character(),
    BRAF_EO = numeric(),
    BRAF_LO =  numeric(),
    BRAF_V600E_EO =  numeric(),
    BRAF_V600E_LO =  numeric(),
    BRAF_V600E_EO_rate =  numeric(),
    BRAF_V600E_LO_rate =  numeric(),
    Pval =  numeric()
  )

hotspot_proportion_calculation <- function(clin_data_1, clin_data_2, 
                                           gene = 'BRAF', mutation = 'p.V600E',
                                           status) {
 
  count_1_overall <- nrow(clin_data_1@data[Hugo_Symbol == gene]%>% 
                            filter(!duplicated(Tumor_Sample_Barcode)))
  count_1_mutation <- nrow(clin_data_1@data[Hugo_Symbol == gene & HGVSp_Short == mutation]%>% 
                             filter(!duplicated(Tumor_Sample_Barcode)))

  count_2_overall <- nrow(clin_data_2@data[Hugo_Symbol == gene]%>% 
                            filter(!duplicated(Tumor_Sample_Barcode)))
  count_2_mutation <- nrow(clin_data_2@data[Hugo_Symbol == gene & HGVSp_Short == mutation]%>% 
                             filter(!duplicated(Tumor_Sample_Barcode)))

  rate_1 <- count_1_mutation / count_1_overall
  rate_2 <- count_2_mutation / count_2_overall
  
  cat("Mutation rate for group 1 (", mutation, "): ", rate_1, "\n")
  cat("Mutation rate for group 2 (", mutation, "): ", rate_2, "\n")

  table <- matrix(c(
    count_1_mutation, count_1_overall - count_1_mutation,
    count_2_mutation, count_2_overall - count_2_mutation
  ), nrow = 2, byrow = TRUE)
  
  if(count_1_mutation <=5  | count_2_mutation <=5){
    chisq_result <- fisher.test(table)
    print('fisher')
  }else{
    chisq_result <- chisq.test(table)
    print('chisq')
  }
  
  print(chisq_result)
  
  df <- data.frame(
    TMB_Status = status,
    BRAF_EO = count_1_overall,
    BRAF_LO = count_2_overall,
    BRAF_V600E_EO = count_1_mutation,
    BRAF_V600E_LO = count_2_mutation,
    BRAF_V600E_EO_rate = rate_1,
    BRAF_V600E_LO_rate = rate_2,
    # Chisq_X = chisq_result$X
    Pval = chisq_result$p.value
  )
  
  proportion_calculation_df  <<- rbind(proportion_calculation_df,df)
  
}


  hotspot_proportion_calculation(clin_data_1 = clin.EOCRC, 
                                 clin_data_2 = clin.LOCRC,
                                 status = 'Overall')
  hotspot_proportion_calculation(clin_data_1 = clin.EOCRC_hyper, 
                                 clin_data_2 = clin.LOCRC_hyper,
                                 status = 'Hyper')
  hotspot_proportion_calculation(clin_data_1 = clin.EOCRC_nonhyper, 
                                 clin_data_2 = clin.LOCRC_nonhyper,
                                 status = 'Non-hyper')

  proportion_calculation_df <- proportion_calculation_df %>% 
                                mutate(adjPval = p.adjust(p = Pval, method = "fdr"))
  
  write.csv(proportion_calculation_df,'./table/10.V600E_proportion_calculation_df.csv')