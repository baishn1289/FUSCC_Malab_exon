library(ggplot2)
library(ggthemes)
library(maftools)
library(dplyr)
library(stringr)
library(data.table)
library(tidyr)
library(tidyverse)
library(cowplot)

uniquegenes <- function(m1, m2, m1Name = NULL, m2Name = NULL,
                        minMut = 5){
  
  m1.gs <- getGeneSummary(x = m1)
  m2.gs <- getGeneSummary(x = m2)
  
  if (is.null(m1Name)) {
    m1Name = "M1"
  }
  if (is.null(m2Name)) {
    m2Name = "M2"
  }
  
  m1.genes = as.character(m1.gs[AlteredSamples >= minMut, 
                                Hugo_Symbol])
  m2.genes = as.character(m2.gs[AlteredSamples >= minMut, 
                                Hugo_Symbol])
  uniquegenes =  intersect(m1.genes, m2.genes)
  
  return(uniquegenes)
  
}

logstic_country_mafcompare <- function(m1, m2, Clinical.data, Gene_test,Guojia,TMB_status) {
  
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
    print(paste(Gene_test,Guojia,TMB_status,
                "Low mutation frequencies, returning from function.",
                sep = '____'))
    return()
  }
 
  Gene_clinical_info  <-  Clinical.data[, Gene_status := ifelse(Tumor_Sample_Barcode %in% c(Alter_m1_Gene$Tumor_Sample_Barcode,
                                                                                            Alter_m2_Gene$Tumor_Sample_Barcode), 1, 0)]%>% 
                            mutate(Age_class = factor(Age_class,levels = c('LOCRC','EOCRC'))) %>% 
                            filter(Panel %in% Gene_panels)
  
  each_gene_mutation_counts = gene_mutation_counts %>% 
                              dplyr::select(all_of(Gene_test), Panel,Age_class, TMB,Tumor_Sample_Barcode) %>% 
                              filter(Tumor_Sample_Barcode %in% Gene_clinical_info$Tumor_Sample_Barcode)  %>% 
                              mutate(mutation_rate = .data[[Gene_test]] / TMB)
      
  each_gene_mutation_EO = each_gene_mutation_counts %>% 
                          filter(Age_class == 'EOCRC') %>% 
                          pull(mutation_rate)
  each_gene_mutation_LO = each_gene_mutation_counts %>% 
                          filter(Age_class == 'LOCRC') %>% 
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
                    or = OR_age_class,
                    ci.up = ci_up,
                    ci.low = ci_low,
                    EOCRC_mutation_rate = EOCRC_rate,
                    LOCRC_mutation_rate = LOCRC_rate,
                    mutation_rate_pval = mutation_rate_pvalue,
                    formula = formula_base,
                    TMB_Status = TMB_status,
                    Country = Guojia
                                     )
  
  result_df <<- rbind(result_df, df_Gene_test)  
}

generate_bubble_table <- function(guojia,Minmut){

  country_tsb<- rt@clinical.data[Country == guojia]$Tumor_Sample_Barcode
 
  if (length(intersect(clin.EOCRC@clinical.data$Tumor_Sample_Barcode, country_tsb)) > 0 &
      length(intersect(clin.LOCRC@clinical.data$Tumor_Sample_Barcode, country_tsb)) > 0){
   
      rt_country_EO <- subsetMaf(clin.EOCRC,tsb=country_tsb)
      rt_country_LO <- subsetMaf(clin.LOCRC,tsb=country_tsb)
      
      bubble_genes <- uniquegenes(m1 = rt_country_EO, m2 = rt_country_LO, minMut = Minmut)
  
      for (i in bubble_genes) {
          logstic_country_mafcompare(m1 = rt_country_EO, m2 = rt_country_LO,
                                        Clinical.data = rbind(rt_country_EO@clinical.data,
                                                                 rt_country_LO@clinical.data),  
                                           Gene_test =i,Guojia = guojia,TMB_status = 'Overall')
       
      }
  
  }
  
  if (length(intersect(clin.EOCRC_hyper@clinical.data$Tumor_Sample_Barcode, country_tsb)) > 0 &
      length(intersect(clin.LOCRC_hyper@clinical.data$Tumor_Sample_Barcode, country_tsb)) > 0) {
    
      rt_country_EO_hyper <- subsetMaf(clin.EOCRC_hyper,tsb=country_tsb)
      rt_country_LO_hyper <- subsetMaf(clin.LOCRC_hyper,tsb=country_tsb)
    
      bubble_genes_hyper <- uniquegenes(m1 = rt_country_EO_hyper, m2 = rt_country_LO_hyper,
                                      minMut = Minmut)
    
      for (i in bubble_genes_hyper) {
          logstic_country_mafcompare(m1 = rt_country_EO_hyper, m2 = rt_country_LO_hyper,
                                     Clinical.data = rbind(rt_country_EO_hyper@clinical.data,
                                                               rt_country_LO_hyper@clinical.data), 
                                      Gene_test =i,Guojia = guojia,TMB_status = 'Hyper')
            
       }
  } 
  
  if (length(intersect(clin.EOCRC_nonhyper@clinical.data$Tumor_Sample_Barcode, country_tsb)) > 0 &
      length(intersect(clin.LOCRC_nonhyper@clinical.data$Tumor_Sample_Barcode, country_tsb)) > 0) {
   
    rt_country_EO_nonhyper <- subsetMaf(clin.EOCRC_nonhyper,tsb=country_tsb)
    rt_country_LO_nonhyper <- subsetMaf(clin.LOCRC_nonhyper,tsb=country_tsb)

    bubble_genes_nonhyper <- uniquegenes(m1 = rt_country_EO_nonhyper, m2 = rt_country_LO_nonhyper,
                                      minMut = Minmut)
   
    for (i in bubble_genes_nonhyper) {
       logstic_country_mafcompare(m1 = rt_country_EO_nonhyper, m2 = rt_country_LO_nonhyper,
                                     Clinical.data = rbind(rt_country_EO_nonhyper@clinical.data,
                                                           rt_country_LO_nonhyper@clinical.data), 
                                     Gene_test =i,Guojia = guojia,TMB_status = 'Non-hyper')
    }
      
  }

}

formula_base <- "Gene_status ~ Age_class"

result_df <-  data.table::data.table(
                Hugo_Symbol = character(),
                EOCRC_total_case = numeric(),
                LOCRC_total_case = numeric(),
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
                TMB_Status = character(), 
                Country = character()
  )

  generate_bubble_table('France',Minmut = 2)
  generate_bubble_table('Nigeria',Minmut = 2)
  generate_bubble_table('Netherlands',Minmut = 2)
  generate_bubble_table('Canada',Minmut = 2)
  generate_bubble_table('Spain',Minmut = 2)
  generate_bubble_table('China',Minmut = 2)
  generate_bubble_table('US',Minmut = 2)
  
  total_bubble_table <- result_df %>% 
                        rename(Gene = Hugo_Symbol, OR =or )

  filter_gene_country <- total_bubble_table %>%
                          filter(TMB_Status == "Overall" & EOCRC_freq < 0.05)  %>%
                          dplyr::select(Gene, Country,EOCRC_freq)

  total_bubble_table_filter <- total_bubble_table %>%
                                anti_join(filter_gene_country, 
                                          by = c("Gene", "Country"))

{
  gene_country_counts <- total_bubble_table_filter %>%
                          group_by(Gene) %>%
                          summarise(Country_count = n_distinct(Country))
  
  genes_in_6_countries <- gene_country_counts %>%
                          filter(Country_count >= 6) %>%
                          dplyr::select(Gene)

}

  bubble_table = total_bubble_table_filter %>%
                filter(Gene %in%  genes_in_6_countries$Gene)%>% 
                mutate(adjPval = p.adjust(p = pval, method = "fdr")) 

  bubble_table$p_sig <- as.character(symnum(bubble_table$adjPval, 
                                          cutpoints = c(0, 0.01, 0.05, 0.25, 1), 
                                          symbols = c("***", "**", "*", "")))
  p1 <- ggplot(data = bubble_table, 
               aes(x = Gene, y = TMB_Status)) +
        geom_point(aes(color = factor(case_when(
                    OR > 1 ~ 'OR > 1',
                    OR < 1 ~ 'OR < 1',
                    OR == 1 ~ 'No significance'
                    ), levels = c('OR > 1', 'OR < 1', 'No significance'))), 
                    size = 6, alpha = 0.7, show.legend = c(size = FALSE))  +
        geom_text(aes(label = p_sig), vjust = 0.8, hjust = 0.5, size = 3, color = "black" ) +
        facet_grid(Country ~ ., scales = "free_y", space = "free") +  
        scale_y_discrete(position = "left") +  
        scale_color_manual(values = c('OR > 1' = '#d53e4f', 'OR < 1' = '#4393c3', 'No significance' = 'gray')) +
        theme_few() +
        labs(y ='TMB Category') +
        theme(
          strip.background = element_blank(), 
          strip.text.y.right = element_text(size = 13, angle = 0),
          panel.spacing = unit(0, "lines"),   
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
          axis.text.y = element_text(color = "black", size = 10, angle = 0),
          axis.text.x = element_text(color = "black", size = 10, angle = 45, vjust = 0.5, hjust = 0.5),
          legend.text = element_text(color = "black", size = 10),
          legend.title = element_blank()  ) +
        guides(color = guide_legend(override.aes = list(size = 4)))

  p1
  
  ggsave(paste(dir,'/8_intersect_gene_among_countries.pdf',sep = ''),
         width = 12,height = 6,dpi = 300)
  
  write.csv(total_bubble_table,file = '8_total_bubble_data.csv',
            quote = FALSE, row.names = FALSE)

