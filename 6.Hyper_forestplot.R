
{
  # TMB hypermutated
  clin.EOCRC_hyper <- subsetMaf(maf=rt, tsb=TMB_EOCRC_hyper_sample, isTCGA=FALSE)
  clin.LOCRC_hyper <- subsetMaf(maf=rt, tsb=TMB_LOCRC_hyper_sample, isTCGA=FALSE)
  
  Clinical_data_hyper <- rt@clinical.data[Tumor_Sample_Barcode %in% c(TMB_EOCRC_hyper_sample,TMB_LOCRC_hyper_sample)]
 
  
   save(clin.EOCRC_hyper,clin.LOCRC_hyper,Clinical_data_hyper,
         file = 'hyper_data.Rdata')
  
  
  
  hyper_genes <- uniquegenes(clin.EOCRC_hyper,clin.LOCRC_hyper,minMut = 10)
  
  Hugo_clinical_info_hyper <- panel_hugo_symbol[Hugo_Symbol %in% hyper_genes] %>%
    merge(Clinical_data_hyper[, c('Tumor_Sample_Barcode', 'Age_class','Country',"Sex","Race")], 
          by = "Tumor_Sample_Barcode") %>%
    group_by(Hugo_Symbol) %>%
    summarise(
      Panel_count = n_distinct(Panel),
      Country_count = n_distinct(Country),
      Race_count = n_distinct(Race),
      Sex_count = n_distinct(Sex)
    ) %>% ungroup()
  
  
  Hugo_clinical_info_xxxx <- Hugo_clinical_info_hyper %>% 
   
    pull(Hugo_Symbol)
 
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
 
  problematic_genes_hyper <- c() 
  
  singlesex_panel_record_log <- data.frame(Gene_test = character(),
                                           stringsAsFactors = FALSE)
  
  for (i in Hugo_clinical_info_xxxx) {
    result <- tryCatch({
      logstic_mafcompare(m1=clin.EOCRC_hyper, m2= clin.LOCRC_hyper, 
                         Clinical.data=Clinical_data_hyper, Gene_test = i)
      NULL  
    }, warning = function(w) {
      problematic_genes_hyper <<- c(problematic_genes_hyper, i)
      return(NULL) 
    })
  }

  
  result_df_hyper <- result_df
  
  result_logistic_hyper <- result_df_hyper[-1,] %>% 
    filter(!Hugo_Symbol %in% problematic_genes_hyper) %>% 
    filter(formula == 'Gene_status ~ Age_class + Race + Panel + Country + Sex') %>% 
    mutate(adjPval = p.adjust(p = pval, method = "fdr"))
  
  
  rownames(result_logistic_hyper) <- seq_len(nrow(result_logistic_hyper))
  
  save(result_df_hyper,problematic_genes_hyper,singlesex_panel_record_log,
       file='hyper_logistic_result.Rdata')
  
  write.table(result_logistic_hyper, file='./hyper_EOCRC_vs_LOCRC_freq.tsv', 
              quote=FALSE, row.names=FALSE, sep="\t")
  
  
}



  {#EO_low_forestplot
    
    df_hyper <- read_tsv('hyper_EOCRC_vs_LOCRC_freq.tsv')
    df_hyper <- df_hyper %>% 
      filter(!is.infinite(or) & !is.infinite(ci.low) & !is.infinite(ci.up))
    
    adjustPval_threshold <- 0.05
    
    
    filtered_df_hyper <- df_hyper %>% 
      filter(
             adjPval <= adjustPval_threshold
             ) %>%
      filter(
        LOCRC_freq > EOCRC_freq
      ) %>%
      mutate(
        Freq_added = (LOCRC_freq+EOCRC_freq )
      ) %>%
      arrange(desc(Freq_added))
    
    write.csv(filtered_df_hyper,
              file = paste(dir,'/6.hyper_differential_gene_EO_low.csv',sep=''),
              row.names = F)
    
    filtered_df_hyper <-filtered_df_hyper %>% head(20)
    
    # add new labels
    plot <- filtered_df_hyper %>% 
      mutate(EOCRC_freq = format_frequency(EOCRC_freq),
             LOCRC_freq = format_frequency(LOCRC_freq),
             sample_size = paste0(LOCRC_freq, " v.s. ", EOCRC_freq),
             OR_CI = paste0(format(round(or, 2), nsmall = 2), " (", format(round(ci.low, 2), nsmall = 2), " to ", format(round(ci.up, 2), nsmall = 2), ")"),
             adjPval = format_adjPval(adjPval)) %>% 
      mutate(across(everything(), as.character)) %>% 
      bind_rows(data.frame(Hugo_Symbol = "Gene", sample_size = "Frequency (LOCRC v.s. EOCRC)", adjPval = "adjPval", OR_CI = "OR (95% CI)")) %>% 
      mutate(Hugo_Symbol = fct_rev(fct_relevel(Hugo_Symbol, "Gene")))
    
    # middle part
    p_mid <- filtered_df_hyper %>% 
      ggplot(aes(y = fct_rev(Hugo_Symbol))) +
      theme_classic() +
      geom_point(aes(x = or, size = or), shape = 23, fill = "black", show.legend = F) +
      geom_errorbarh(aes(xmin = ci.low, xmax = ci.up), height = 0.3) +
      labs(x = "Odds Ratio") +
      scale_x_continuous(breaks = c(0.1,0.5,1.0),
                         labels = scales::number_format(accuracy = 0.1)) +
      coord_cartesian(ylim = c(1,3), xlim = c(0.1,1.0)) +
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
    
    # right_part
    p_right <- plot %>% ggplot() +
      geom_text(aes(x = 0, y = Hugo_Symbol, label = OR_CI), hjust = 0,
                fontface = ifelse(plot$OR_CI == "OR (95% CI)", "bold", "plain"), size = 4) +
      geom_text(aes(x = 1, y = Hugo_Symbol, label = adjPval), hjust = 0,
                fontface = ifelse(plot$adjPval == "adjPval", "bold", "plain"), size = 4) +
      theme_void() + coord_cartesian(xlim = c(0, 2))
    
    
    layout <- c(
      area(t = 0, l = 0, b = 30, r = 10), 
      area(t = 0, l = 11, b = 30, r = 17), 
      area(t = 0, l = 18.5, b = 30, r = 25)) 
    
    
  }
 


p_hyper_low =p_left + p_mid + p_right + plot_layout(design = layout)


{#EO_high_forestplot
 
  
  filtered_df_hyper <- df_hyper %>% 
    filter(
      adjPval <= adjustPval_threshold
    ) %>%
    filter(
      LOCRC_freq < EOCRC_freq
    ) %>%
    mutate(
      Freq_added = (LOCRC_freq+EOCRC_freq )
    ) %>%
    arrange(desc(Freq_added))
  
  write.csv(filtered_df_hyper,
            file = paste(dir,'/6.hyper_differential_gene_EO_high.csv',sep=''),
            row.names = F)
  
  filtered_df_hyper <-filtered_df_hyper %>% head(20)
  
  # add new labels
  plot <- filtered_df_hyper %>% 
    mutate(EOCRC_freq = format_frequency(EOCRC_freq),
           LOCRC_freq = format_frequency(LOCRC_freq),
           sample_size = paste0(LOCRC_freq, " v.s. ", EOCRC_freq),
           OR_CI = paste0(format(round(or, 2), nsmall = 2), " (", format(round(ci.low, 2), nsmall = 2), " to ", format(round(ci.up, 2), nsmall = 2), ")"),
           adjPval = format_adjPval(adjPval)) %>% 
    mutate(across(everything(), as.character)) %>% 
    bind_rows(data.frame(Hugo_Symbol = "Gene", sample_size = "Frequency (LOCRC v.s. EOCRC)", adjPval = "adjPval", OR_CI = "OR (95% CI)")) %>% 
    mutate(Hugo_Symbol = fct_rev(fct_relevel(Hugo_Symbol, "Gene")))
  
  # middle part
  p_mid <- filtered_df_hyper %>% 
    ggplot(aes(y = fct_rev(Hugo_Symbol))) +
    theme_classic() +
    geom_point(aes(x = or, size = or), shape = 23, fill = "black", show.legend = F) +
    geom_errorbarh(aes(xmin = ci.low, xmax = ci.up), height = 0.3) +
    labs(x = "Odds Ratio") +
    scale_x_continuous(trans = 'log10',breaks = c(1,2,4,8,16),
                       labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(1,21), xlim = c(1.0,16)) +
  
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
  
  # right_part
  p_right <- plot %>% ggplot() +
    geom_text(aes(x = 0, y = Hugo_Symbol, label = OR_CI), hjust = 0,
              fontface = ifelse(plot$OR_CI == "OR (95% CI)", "bold", "plain"), size = 4) +
    geom_text(aes(x = 1, y = Hugo_Symbol, label = adjPval), hjust = 0,
              fontface = ifelse(plot$adjPval == "adjPval", "bold", "plain"), size = 4) +
    theme_void() + coord_cartesian(xlim = c(0, 2))
  
  
  layout <- c(
    area(t = 0, l = 0, b = 30, r = 10), 
    area(t = 0, l = 11, b = 30, r = 17), 
    area(t = 0, l = 18.5, b = 30, r = 25)) 
  
  
}


p_hyper_high =p_left + p_mid + p_right + plot_layout(design = layout)


