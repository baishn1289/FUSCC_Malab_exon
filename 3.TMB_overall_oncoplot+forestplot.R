library(tidyverse)
library(patchwork)
library(scales)
library(dplyr)
library(stringr)



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
    
    # numMutatedSamples = maf[!Variant_Type %in% 'CNV', .(MutatedSamples = length(unique(Tumor_Sample_Barcode))), by = Hugo_Symbol]
    # numAlteredSamples = maf[, .(AlteredSamples = length(unique(Tumor_Sample_Barcode))), by = Hugo_Symbol]
    # numAlteredSamples = merge(numMutatedSamples, numAlteredSamples, by = 'Hugo_Symbol', all = TRUE)
    
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
  
  
  
  cl$Age_class <-  factor(cl$Age_class,levels = c('LOCRC','EOCRC'))
  rt@clinical.data$Age_class <-  factor(rt@clinical.data$Age_class
                                        ,levels = c('LOCRC','EOCRC'))
  


 
  logstic_mafcompare <- function(m1, m2, m1Name = NULL, m2Name = NULL, Clinical.data,Gene_test) {
    
   
    Alter_m1_Gene = m1@data[Hugo_Symbol == Gene_test,] %>% 
      filter(!duplicated(Tumor_Sample_Barcode))
    
    Alter_m2_Gene = m2@data[Hugo_Symbol == Gene_test,] %>% 
      filter(!duplicated(Tumor_Sample_Barcode))
    
    
    Gene_panels <- panel_hugo_symbol[Hugo_Symbol == Gene_test] %>%
      merge(Clinical.data[, c('Tumor_Sample_Barcode', 'Age_class', 'Sex', 'Country', 'Race')], by = "Tumor_Sample_Barcode") %>%
      filter(!duplicated(Tumor_Sample_Barcode)) %>%
      group_by(Panel) %>%
      summarise(
        Age_class_count = n_distinct(Age_class)
      ) %>%
      ungroup() %>% 
      filter(Age_class_count == 2) %>% 
      pull(Panel)
  
  
    
    
    if (length(Gene_panels) == 0) {
      print(Gene_test)
     
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
    
   
    
    Gene_clinical_info = Clinical.data[, Gene_status := ifelse(Tumor_Sample_Barcode %in% c(Alter_m1_Gene$Tumor_Sample_Barcode,
                                                                                              Alter_m2_Gene$Tumor_Sample_Barcode), 1, 0)] %>% 
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
      or = OR_age_class,
      ci.up = ci_up,
      ci.low = ci_low,
      formula = formula_base
    )
    
    
    result_df <<- rbind(result_df, df_Gene_test)  
    
    
  }
  
  
  gc()
 
 
  overall_genes <- uniquegenes(clin.EOCRC,clin.LOCRC,minMut = 10)

  Hugo_clinical_info <- panel_hugo_symbol[Hugo_Symbol %in% overall_genes] %>%
    merge(rt@clinical.data[, c('Tumor_Sample_Barcode', 'Age_class','Country',"Sex","Race")], 
          by = "Tumor_Sample_Barcode") %>%
    group_by(Hugo_Symbol) %>%
    summarise(
      Panel_count = n_distinct(Panel),
      Country_count = n_distinct(Country),
      Race_count = n_distinct(Race),
      Sex_count = n_distinct(Sex)
    ) %>% ungroup()


  Hugo_clinical_info_xxxx <- Hugo_clinical_info %>% 
    pull(Hugo_Symbol)

  
  singlesex_panel_record_log <- data.frame(Gene_test = character(),
                                           stringsAsFactors = FALSE)
  
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
    or = NA,
    ci.up =  NA,
    ci.low =  NA,
    formula =NA
  )
  
  formula_base <- "Gene_status ~ Age_class"
  
  for (i in Hugo_clinical_info_xxxx) {
    result <- tryCatch({
      logstic_mafcompare(m1 = clin.EOCRC, m2 = clin.LOCRC, Clinical.data = rt@clinical.data, Gene_test = i)
      NULL  
    }, warning = function(w) {
      problematic_genes <<- c(problematic_genes, i)
      return(NULL)  
    })
  }
  
  
  result_df_overall <- result_df
  
  result_logistic <- result_df_overall[-1,] %>% 
    filter(!Hugo_Symbol %in% problematic_genes) %>% 
    filter(formula == 'Gene_status ~ Age_class + Race + Panel + Country + Sex') %>% 
    mutate(adjPval = p.adjust(p = pval, method = "fdr")) 
  
  
  rownames(result_logistic) <- seq_len(nrow(result_logistic))
  
  save(result_df,problematic_genes,singlesex_panel_record_log,
       file='overall_logistic_result.Rdata')
  
  write.table(result_logistic, file='./overall_EOCRC_vs_LOCRC_freq.tsv', 
              quote=FALSE, row.names=FALSE, sep="\t")
  
  
} 


  {#overall_forest
    
    df_overall <- read_tsv('overall_EOCRC_vs_LOCRC_freq.tsv')
    
    df_overall <- df_overall %>% 
      filter(!is.infinite(or) & !is.infinite(ci.low) & !is.infinite(ci.up))
    
   
    adjustPval_threshold <- 0.05
    
    filtered_df_overall <- df_overall %>% 
      filter(
        adjPval <= adjustPval_threshold) %>%
     
      mutate(
        Freq_average = (LOCRC_freq+EOCRC_freq )/2
      ) %>%
      arrange(desc(Freq_added)) 
    
    write.csv(filtered_df_overall,
              file = paste(dir,'/3.overall_differential_gene.csv',sep=''),
              row.names = F)
    
    filtered_df_overall <-filtered_df_overall %>% head(20)
    
    
    # add new labels
    plot <- filtered_df_overall %>% 
      mutate(EOCRC_freq = format_frequency(EOCRC_freq),
             LOCRC_freq = format_frequency(LOCRC_freq),
             sample_size = paste0(LOCRC_freq, " v.s. ", EOCRC_freq),
             OR_CI = paste0(format(round(or, 2), nsmall = 2), " (", format(round(ci.low, 2), nsmall = 2), " to ", format(round(ci.up, 2), nsmall = 2), ")"),
             adjPval = format_adjPval(adjPval)) %>% 
      mutate(across(everything(), as.character)) %>% 
      bind_rows(data.frame(Hugo_Symbol = "Gene", sample_size = "Frequency (LOCRC v.s. EOCRC)", adjPval = "adjPval", OR_CI = "OR (95% CI)")) %>% 
      mutate(Hugo_Symbol = fct_rev(fct_relevel(Hugo_Symbol, "Gene")))
    
    # middle part
    p_mid <- filtered_df_overall %>% 
      ggplot(aes(y = fct_rev(Hugo_Symbol))) +
      theme_classic() +
      geom_point(aes(x = or, size = or), shape = 23, fill = "black", show.legend = F) +
      geom_errorbarh(aes(xmin = ci.low, xmax = ci.up), height = 0.3) +
      labs(x = "Odds Ratio") +
      scale_x_continuous( trans = 'log10',breaks = c(0.5,1,2,4),
                          labels = scales::number_format(accuracy = 0.1)) +
      coord_cartesian(ylim = c(1, 21), xlim = c(0.5,4)) +
      geom_vline(xintercept = 1, linetype = "dashed") +  
      theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
            axis.text.y = element_blank(), axis.title.y = element_blank(),
            axis.text.x = element_text(color = "black"))
    
    # left part
    p_left <- plot %>%
      ggplot(aes(y = Hugo_Symbol)) +
      geom_text(aes(x = 0, label = Hugo_Symbol), hjust = 0, fontface = "bold", size = 4) +
      geom_text(aes(x = 1, label = sample_size), hjust = 0, size = 4, vjust = 0,
                fontface = ifelse(plot$sample_size == "Frequency (LOCRC v.s. EOCRC)", "bold", "plain")) +
      theme_void() + coord_cartesian(xlim = c(0, 3))
    
    # right part
    p_right <- plot %>% ggplot() +
      geom_text(aes(x = 0, y = Hugo_Symbol, label = OR_CI), hjust = 0,
                fontface = ifelse(plot$OR_CI == "OR (95% CI)", "bold", "plain"), size = 4) +
      geom_text(aes(x = 1, y = Hugo_Symbol, label = adjPval), hjust = 0,
                fontface = ifelse(plot$adjPval == "adjPval", "bold", "plain"), size = 4) +
      theme_void() + coord_cartesian(xlim = c(0, 2))
    
    
    layout <- c(
      area(t = 0, l = 0, b = 30, r = 9), 
      area(t = 0, l = 10, b = 30, r = 17), 
      area(t = 0, l = 18, b = 30, r = 25)) 
    
    
  }


pdf(paste(dir,'/3_overall_forestPlot_freq_EOLOadded20.pdf',sep = ""), 
    width = 11, height = 6)
print(p_left + p_mid + p_right + plot_layout(design = layout))

dev.off()
