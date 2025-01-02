library(tidyverse)
library(patchwork)

dir= './picture'

format_pvalue <- function(p) {
  ifelse(p < 0.001, format(p, scientific = TRUE, digits = 3), 
         format(p, digits = 3, nsmall = 3))
}

custom_round <- function(x, digits) {
  posneg <- sign(x)
  z <- abs(x)*10^digits
  z <- z + 0.5
  z <- trunc(z)
  z <- z/10^digits
  z * posneg
}

format_frequency <- function(freq) {
  rounded_value <- custom_round(freq * 100, 1)
  paste0(format(rounded_value, nsmall = 1), "%")
}

format_adjPval <- function(p) {
  format(p, digits = 3, scientific = TRUE)
}



{
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
 
  mutation_data <- rt@data[, 2:30]
  dim(mutation_data)

  gene_mutation_counts <- mutation_data %>% 
                          group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% 
                          reframe(Mutation_Count = n()) %>%
                          pivot_wider(names_from = Hugo_Symbol, 
                                      values_from = Mutation_Count,
                                      values_fill = list(Mutation_Count = 0)) %>% 
                          left_join(TMB_mutload[, .(Panel, Tumor_Sample_Barcode, 
                                                    Age_class, TMB,Race)], by = 'Tumor_Sample_Barcode') 
  
  gene_panel_map <- assays_genes %>% 
                    group_by(Hugo_Symbol) %>% 
                    reframe(Genes = list(Panel)) %>%
                    deframe() 
  
logstic_mafcompare <- function(m1, m2, m1Name = NULL, 
                               m2Name = NULL, Clinical.data,Gene_test) {
    
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
    
    mut_samples_m2 <- length(unique(Alter_m2_Gene$Tumor_Sample_Barcode))
    
    panels_samples_m2 <- length(unique(m2@clinical.data[m2@clinical.data$Panel %in% Gene_panels, ]$Tumor_Sample_Barcode))
    
    m2_freq <- mut_samples_m2 / panels_samples_m2

    if ( m1_freq < 0.01 | m2_freq<0.01 ) {
          print(paste(Gene_test,
                      "Low mutation frequencies, returning from function.",
                      sep = '____'))
          return()
    }
    
    Gene_clinical_info  <-  Clinical.data[, Gene_status := ifelse(Tumor_Sample_Barcode %in% c(Alter_m1_Gene$Tumor_Sample_Barcode,
                                                                                              Alter_m2_Gene$Tumor_Sample_Barcode), 1, 0)]%>% 
                            mutate(Age_class = factor(Age_class,levels = c('LOCRC','EOCRC'))) %>% 
                            filter(Panel %in% Gene_panels)

    #mutation_rate for each gene
    each_gene_mutation_counts = gene_mutation_counts %>% 
                                dplyr::select(all_of(Gene_test), Tumor_Sample_Barcode,Age_class, TMB,Race) %>%
                                filter(Tumor_Sample_Barcode %in% Gene_clinical_info$Tumor_Sample_Barcode)  %>% 
                                mutate(mutation_rate = .data[[Gene_test]] / TMB)
    
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
                        formula = formula_base )
 
    result_df <<- rbind(result_df, df_Gene_test)  

}

  overall_genes <- uniquegenes(clin.EOCRC,clin.LOCRC,minMut = 5)
  
  problematic_genes <- c() 
  
  
  result_df_template_3.0 <- data.table::data.table(
    Hugo_Symbol = character(),
    EOCRC_total_case= numeric(),
    LOCRC_total_case= numeric(),
    EOCRC = numeric(), 
    LOCRC = numeric(), 
    EOCRC_freq = numeric(),
    LOCRC_freq = numeric(),
    pval = numeric(), 
    or = numeric(),
    ci.up =  numeric(),
    ci.low =  numeric(),
    EOCRC_mutation_rate = numeric(),
    LOCRC_mutation_rate = numeric(),
    mutation_rate_pval = numeric(),
    formula =character()
  )
  
  result_df <- result_df_template_3.0
  
  formula_base <- "Gene_status ~ Age_class"
 
for (i in overall_genes) {
     tryCatch({
      logstic_mafcompare(m1 = clin.EOCRC, m2 = clin.LOCRC,
                         Clinical.data = rt@clinical.data, Gene_test = i)
      NULL 
    }, warning = function(w) {
      problematic_genes <<- c(problematic_genes, i) 
      return(NULL) 
    })
}
   
  result_df_overall <- result_df
  
  result_logistic <- result_df_overall  %>% 
                      filter(!Hugo_Symbol %in% problematic_genes) %>%
                      mutate(mean_freq = (EOCRC_freq+LOCRC_freq)/2) %>% 
                      filter(mean_freq > 0.01,
                             #MSS：formula=='Gene_status ~ Age_class + Panel + Sex + Tumor_site + Histology_type'
                             # formula=='Gene_status ~ Age_class + Race + Panel + Country + Sex + Tumor_site + Sample_type + Histology_type'
                             formula=='Gene_status ~ Age_class + Panel + Sex + Tumor_site + Histology_type'  )%>%
                      mutate(adjPval = p.adjust(p = pval, method = "fdr"),
                             mutation_rate_adjPval =  p.adjust(p = mutation_rate_pval, method = "fdr")) 

  rownames(result_logistic) <- seq_len(nrow(result_logistic))
  
  save(result_df,problematic_genes,
       file='overall_logistic_result.Rdata')
  
  write.table(result_logistic, file='./overall_EOCRC_vs_LOCRC_freq.tsv', 
              quote=FALSE, row.names=FALSE, sep="\t")
  
} 

{   #overall_forest
    
  df_overall  <- read_tsv('overall_EOCRC_vs_LOCRC_freq.tsv')
    
  df_overall <- df_overall %>% 
                filter(!is.infinite(or) & !is.infinite(ci.low) & !is.infinite(ci.up))
   
  filtered_df_overall <- df_overall %>% 
                        filter(adjPval <= 0.05) %>%
                        mutate(Freq_average = (LOCRC_freq+EOCRC_freq )/2) %>%
                        arrange(desc(Freq_average)) 
    
  write.csv(filtered_df_overall,
            file = paste(dir,'/3.overall_differential_gene.csv',sep=''),
            row.names = F)
    
  filtered_df_overall <-filtered_df_overall %>% 
                        head(20)
    
    
  # add new labels
  plot <- filtered_df_overall %>% 
          mutate(EOCRC_freq = format_frequency(EOCRC_freq),
                 LOCRC_freq = format_frequency(LOCRC_freq),
                 sample_size = paste0(LOCRC_freq, " v.s. ", EOCRC_freq),
                 OR_CI = paste0(format(round(or, 2), nsmall = 2), " (", format(round(ci.low, 2), nsmall = 2), " to ", format(round(ci.up, 2), nsmall = 2), ")"),
                 adjPval = format_adjPval(adjPval),
                 EOCRC_mutation_rate = format_frequency(EOCRC_mutation_rate),
                 LOCRC_mutation_rate = format_frequency(LOCRC_mutation_rate),
                 mutation_rate = paste0(LOCRC_mutation_rate, " v.s. ", EOCRC_mutation_rate),
                 mutation_rate_adjPval = format_adjPval(mutation_rate_adjPval)  ) %>%
          mutate(across(everything(), as.character)) %>% 
          bind_rows(data.frame(Hugo_Symbol = "Gene", 
                               sample_size = "Frequency\n(LOCRC v.s. EOCRC)", 
                               adjPval = "adjPval", 
                               OR_CI = "OR (95% CI)",
                               mutation_rate = 'mutation_rate\n(LOCRC v.s. EOCRC)',
                               mutation_rate_adjPval = 'mutation_rate_adjPval')) %>% 
          mutate(Hugo_Symbol = fct_rev(fct_relevel(Hugo_Symbol, "Gene")))

  p_mid <- filtered_df_overall %>% 
            ggplot(aes(y = fct_rev(Hugo_Symbol))) +
            theme_classic() +
            geom_point(aes(x = or, size = or), shape = 23, fill = "black", show.legend = F) +
            geom_errorbarh(aes(xmin = ci.low, xmax = ci.up), height = 0.3) +
            labs(x = "Odds Ratio") +
            scale_x_continuous( trans = 'log10',breaks = c(0.1,1,2,4),
                                labels = scales::number_format(accuracy = 0.1)) +
            coord_cartesian(ylim = c(1, 21), xlim = c(0.1,4)) +
            geom_vline(xintercept = c(1), linetype = 1) +  
            theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                  axis.text.y = element_blank(), axis.title.y = element_blank(),
                  axis.text.x = element_text(color = "black"))
    
  p_left <- plot %>%
            ggplot(aes(y = Hugo_Symbol)) +
            geom_text(aes(x = 0, label = Hugo_Symbol), hjust = 0, fontface = "bold", size = 4) +
            geom_text(aes(x = 1, label = sample_size), hjust = 0.5, size = 4,vjust=0.5,
                      fontface = ifelse(plot$sample_size == "Frequency\n(LOCRC v.s. EOCRC)", "bold", "plain")) +
            theme_void() + coord_cartesian(xlim = c(0, 1.5))
    
 p_right <- plot %>% 
            ggplot() +
            geom_text(aes(x = 0, y = Hugo_Symbol, label = OR_CI), hjust = 0,
                      fontface = ifelse(plot$OR_CI == "OR (95% CI)", "bold", "plain"), size = 4) +
            geom_text(aes(x = 1.2, y = Hugo_Symbol, label = adjPval), hjust = 0,
                      fontface = ifelse(plot$adjPval == "adjPval", "bold", "plain"), size = 4) +
            geom_text(aes(x = 2.3, y = Hugo_Symbol, label = mutation_rate), hjust = 0.5,
                      fontface = ifelse(plot$mutation_rate == 'mutation_rate\n(LOCRC v.s. EOCRC)', "bold", "plain"), size = 4) +
            geom_text(aes(x = 3, y = Hugo_Symbol, label = mutation_rate_adjPval), hjust = 0,
                      fontface = ifelse(plot$mutation_rate_adjPval == "mutation_rate_adjPval", "bold", "plain"), size = 4) +
            theme_void() + coord_cartesian(xlim = c(0, 4))

 layout <- c(   patchwork::area(t = 0, l = 0, b = 30, r = 6), 
                patchwork::area(t = 0, l = 7, b = 30, r = 13), 
                patchwork::area(t = 0, l = 14.5, b = 30, r = 26.5))
    
  }


  pdf(paste(dir,'/3_overall_forestPlot_freq_EOLOadded20.pdf',sep = ""), 
      width = 14, height = 8)
  print(p_left + p_mid + p_right + plot_layout(design = layout))
  
  dev.off()