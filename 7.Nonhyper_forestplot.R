# dir = 'picture'

{
  # TMB nonhypermutated
  clin.EOCRC_nonhyper <- subsetMaf(maf=rt, tsb=TMB_EOCRC_nonhyper_sample, isTCGA=FALSE)
  clin.LOCRC_nonhyper <- subsetMaf(maf=rt, tsb=TMB_LOCRC_nonhyper_sample, isTCGA=FALSE)
  
  Clinical_data_nonhyper <- rt@clinical.data[Tumor_Sample_Barcode %in% c(TMB_EOCRC_nonhyper_sample,TMB_LOCRC_nonhyper_sample)]
  
  
 
  # fvsm <- mafCompare(m1=clin.EOCRC_nonhyper, m2=clin.LOCRC_nonhyper, 
  #                    m1Name="EOCRC", m2Name="LOCRC", minMut=5)
  
  
  save(clin.EOCRC_nonhyper,clin.LOCRC_nonhyper,Clinical_data_nonhyper,
      file = 'nonhyper_data.Rdata')
  
  
  nonhyper_genes <- uniquegenes(clin.EOCRC_nonhyper,clin.LOCRC_nonhyper,minMut = 10)
  
  #提取出nonhyper队列的基因，用对应的TSB匹配对应的临床特征信息，并统计临床特征信息
  Hugo_clinical_info_nonhyper <- panel_hugo_symbol[Hugo_Symbol %in% nonhyper_genes] %>%
    merge(Clinical_data_nonhyper[, c('Tumor_Sample_Barcode', 'Age_class','Country',"Sex","Race")], 
          by = "Tumor_Sample_Barcode") %>%
    group_by(Hugo_Symbol) %>%
    summarise(
      Panel_count = n_distinct(Panel),
      Country_count = n_distinct(Country),
      Race_count = n_distinct(Race),
      Sex_count = n_distinct(Sex)
    ) %>% ungroup()
  
  
  Hugo_clinical_info_xxxx <- Hugo_clinical_info_nonhyper %>% 
    # 函数中已经有对逻辑回归公式进行调整，不需要此处再做筛选
    # filter(Panel_count>1,Country_count>1,Race_count>1,Sex_count>1) %>% 
    pull(Hugo_Symbol)
  
  result_df <- data.table::data.table(
    Hugo_Symbol = NA,
    EOCRC_total_case= NA,
    LOCRC_total_case= NA,
    EOCRC = NA, 
    LOCRC = NA, 
    EOCRC_freq = NA,
    LOCRC_freq = NA,
    EOCRC_freq_normalized = NA,
    LOCRC_freq_normalized = NA,
    EOCRC_panel_TMB= NA,
    LOCRC_panel_TMB= NA,
    pval = NA, 
    or = NA,
    ci.up =  NA,
    ci.low =  NA,
    formula =NA
  )
  
  formula_base <- "Gene_status ~ Age_class+ TMB"
  
  singlesex_panel_record_log <- data.frame(Gene_test = character(),
                                           stringsAsFactors = FALSE)
  
 
  problematic_genes_nonhyper <- c() 
  
  for (i in Hugo_clinical_info_xxxx) {
    result <- tryCatch({
      logstic_mafcompare(m1=clin.EOCRC_nonhyper, m2= clin.LOCRC_nonhyper, 
                         Clinical.data=Clinical_data_nonhyper, Gene_test = i)
      NULL  # 如果没有错误，返回NULL
    }, warning = function(w) {
      problematic_genes_nonhyper <<- c(problematic_genes_nonhyper, i)  # 记录产生警告的基因
      return(NULL)  # 返回NULL以继续循环
    })
  }
  
  result_df_nonhyper <- result_df
  
  result_logistic_nonhyper <- result_df_nonhyper[-1,] %>% 
    filter(!Hugo_Symbol %in% problematic_genes_nonhyper) %>% 
    filter(formula == 'Gene_status ~ Age_class + Race + Panel + Country + Sex') %>% 
    mutate(adjPval = p.adjust(p = pval, method = "fdr"))
  
  
  rownames(result_logistic_nonhyper) <- seq_len(nrow(result_logistic_nonhyper))
  
  save(result_df_nonhyper,problematic_genes_nonhyper,singlesex_panel_record_log,
       file='nonhyper_logistic_result.Rdata')
  
  write.table(result_logistic_nonhyper, file='./nonhyper_EOCRC_vs_LOCRC_freq.tsv', 
              quote=FALSE, row.names=FALSE, sep="\t")
  
  
}



{
  # TMB nonhypermutated
  # clin.EOCRC_nonhyper <- subsetMaf(maf=rt, tsb=TMB_EOCRC_nonhyper_sample, isTCGA=FALSE)
  # clin.LOCRC_nonhyper <- subsetMaf(maf=rt, tsb=TMB_LOCRC_nonhyper_sample, isTCGA=FALSE)
  # 
  # fvsm <- mafCompare(m1=clin.EOCRC_nonhyper, m2=clin.LOCRC_nonhyper, m1Name="EOCRC", m2Name="LOCRC", minMut=5)
  # write.table(fvsm$results, file=paste(dir,'/nonhyperEOCRC_vs_LOCRC_counts.tsv',sep = ""), quote=FALSE, row.names=FALSE, sep="\t")
  # 
  # fvsm$SampleSummary$SampleSize[1]
  # op = fvsm$results
  # op$EOCRC_freq = round(op$EOCRC/fvsm$SampleSummary$SampleSize[1], 4)
  # op$LOCRC_freq = round(op$LOCRC/fvsm$SampleSummary$SampleSize[2], 4)
  # genes = op$Hugo_Symbol[op$pval < 0.05] # adj
  # genes
  # 
  # write.table(op, file='./nonhyper_EOCRC_vs_LOCRC_freq.tsv', 
  #             quote=FALSE, row.names=FALSE, sep="\t")
  
  
  {
    #nonhyper_EO_low
    
    df_nonhyper <- read_tsv('nonhyper_EOCRC_vs_LOCRC_freq.tsv')
    df_nonhyper <- df_nonhyper %>% 
      filter(!is.infinite(or) & !is.infinite(ci.low) & !is.infinite(ci.up))
    
    # set P_value and adjustPval if needed
    #p_value_threshold <- 0.05
    adjustPval_threshold <- 0.05
    
    filtered_df_nonhyper <- df_nonhyper %>% 
      filter(#pval <= p_value_threshold,
             adjPval <= adjustPval_threshold) %>%
      filter(
        LOCRC_freq > EOCRC_freq
      ) %>%
      mutate(
        Freq_average = (LOCRC_freq+EOCRC_freq )/2
      ) %>%
      arrange(desc(Freq_added))
    
    write.csv(filtered_df_nonhyper,
              file = paste(dir,'/7.nonhyper_differential_gene_EO_low.csv',sep=''),
              row.names = F)
    
    filtered_df_nonhyper <-filtered_df_nonhyper %>% head(20)
    
    # add new labels
    plot <- filtered_df_nonhyper %>% 
      mutate(EOCRC_freq = format_frequency(EOCRC_freq),
             LOCRC_freq = format_frequency(LOCRC_freq),
             sample_size =  paste0(LOCRC_freq, " v.s. ", EOCRC_freq),
             OR_CI = paste0(format(round(or, 2), nsmall = 2), " (", format(round(ci.low, 2), nsmall = 2), " to ", format(round(ci.up, 2), nsmall = 2), ")"),
             adjPval = format_adjPval(adjPval)) %>% 
      mutate(across(everything(), as.character)) %>% 
      bind_rows(data.frame(Hugo_Symbol = "Gene", sample_size = "Frequency (LOCRC v.s. EOCRC)", adjPval = "adjPval", OR_CI = "OR (95% CI)")) %>% 
      mutate(Hugo_Symbol = fct_rev(fct_relevel(Hugo_Symbol, "Gene")))
    
    
    p_mid <- filtered_df_nonhyper %>% 
      ggplot(aes(y = fct_rev(Hugo_Symbol))) +
      theme_classic() +
      geom_point(aes(x = or, size = or), shape = 23, fill = "black", show.legend = F) +
      geom_errorbarh(aes(xmin = ci.low, xmax = ci.up), height = 0.3) +
      labs(x = "Odds Ratio") +
       scale_x_continuous(breaks = c(0.4,0.8,1.0),
                          labels = scales::number_format(accuracy = 0.01)) +
       coord_cartesian(ylim = c(1, 8),xlim = c(0.4,1.0)) +
      # scale_x_continuous(trans = 'log10', breaks = c(1.0, 1.5, 2.0),
      #                    labels = scales::number_format(accuracy = 0.01)) +
      # coord_cartesian(ylim = c(1, 3),xlim = c(1,2)) +
      
      geom_vline(xintercept = 1, linetype = "dashed") +  
      theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
            axis.text.y = element_blank(), axis.title.y = element_blank(),
            axis.text.x = element_text(color = "black"))
    
    
    p_left <- plot %>%
      ggplot(aes(y = Hugo_Symbol)) +
      geom_text(aes(x = 0, label = Hugo_Symbol), hjust = 0, fontface = "bold", size = 4) +
      geom_text(aes(x = 1, label = sample_size), hjust = 0, size = 4, vjust = 0,
                fontface = ifelse(plot$sample_size == "Frequency (LOCRC v.s. EOCRC)", "bold", "plain")) +
      theme_void() + coord_cartesian(xlim = c(0, 3))
    
   
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
  
  
}

p_nonhyper_low =p_left + p_mid + p_right + plot_layout(design = layout)


{
  #nonhyper_EO_high
  
  # df_nonhyper <- read_tsv('nonhyper_EOCRC_vs_LOCRC_freq.tsv')
  # df_nonhyper <- df_nonhyper %>% 
  #   filter(!is.infinite(or) & !is.infinite(ci.low) & !is.infinite(ci.up))
  # 
  # set P_value and adjustPval if needed
  #p_value_threshold <- 0.05
  adjustPval_threshold <- 0.05
  
  filtered_df_nonhyper <- df_nonhyper %>% 
    filter(#pval <= p_value_threshold,
      adjPval <= adjustPval_threshold) %>%
    filter(
      LOCRC_freq < EOCRC_freq
    ) %>%
    mutate(
      Freq_average = (LOCRC_freq+EOCRC_freq )/2
    ) %>%
    arrange(desc(Freq_added))
  
  write.csv(filtered_df_nonhyper,
            file = paste(dir,'/7.nonhyper_differential_gene_EO_high.csv',sep=''),
            row.names = F)
  
  filtered_df_nonhyper <-filtered_df_nonhyper %>% head(20)
  
  # add new labels
  plot <- filtered_df_nonhyper %>% 
    mutate(EOCRC_freq = format_frequency(EOCRC_freq),
           LOCRC_freq = format_frequency(LOCRC_freq),
           sample_size =  paste0(LOCRC_freq, " v.s. ", EOCRC_freq),
           OR_CI = paste0(format(round(or, 2), nsmall = 2), " (", format(round(ci.low, 2), nsmall = 2), " to ", format(round(ci.up, 2), nsmall = 2), ")"),
           adjPval = format_adjPval(adjPval)) %>% 
    mutate(across(everything(), as.character)) %>% 
    bind_rows(data.frame(Hugo_Symbol = "Gene", sample_size = "Frequency (LOCRC v.s. EOCRC)", adjPval = "adjPval", OR_CI = "OR (95% CI)")) %>% 
    mutate(Hugo_Symbol = fct_rev(fct_relevel(Hugo_Symbol, "Gene")))
  
  
  p_mid <- filtered_df_nonhyper %>% 
    ggplot(aes(y = fct_rev(Hugo_Symbol))) +
    theme_classic() +
    geom_point(aes(x = or, size = or), shape = 23, fill = "black", show.legend = F) +
    geom_errorbarh(aes(xmin = ci.low, xmax = ci.up), height = 0.3) +
    labs(x = "Odds Ratio") +
    scale_x_continuous( breaks = c(1.0,2.0),
                       labels = scales::number_format(accuracy = 0.01)) +
    coord_cartesian(ylim = c(1, 2),xlim = c(1.0,2.0)) +
    # scale_x_continuous(trans = 'log10', breaks = c(1.0, 1.5, 2.0),
    #                    labels = scales::number_format(accuracy = 0.01)) +
    # coord_cartesian(ylim = c(1, 3),xlim = c(1,2)) +
    
    geom_vline(xintercept = 1, linetype = "dashed") +  
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(color = "black"))
  
  
  p_left <- plot %>%
    ggplot(aes(y = Hugo_Symbol)) +
    geom_text(aes(x = 0, label = Hugo_Symbol), hjust = 0, fontface = "bold", size = 4) +
    geom_text(aes(x = 1, label = sample_size), hjust = 0, size = 4, vjust = 0,
              fontface = ifelse(plot$sample_size == "Frequency (LOCRC v.s. EOCRC)", "bold", "plain")) +
    theme_void() + coord_cartesian(xlim = c(0, 3))
  
  
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


 
p_nonhyper_high =p_left + p_mid + p_right + plot_layout(design = layout)



library(cowplot)

pdf(paste(dir,'/6_ForestPlot_binded.pdf',sep = ""), 
    width = 11, height = 14)

plot_grid(p_hyper_high,p_hyper_low,
          p_nonhyper_high,p_nonhyper_low,
          
          ncol = 1,rel_heights = c(24,6,5,11))


dev.off()

# 
# pdf(paste(dir,'/4_nonhyper_forestPlot_EOhigher_freqadded20.pdf',sep = ""), 
#     width = 11, height = 1.6)
# 
# print(p_left + p_mid + p_right + plot_layout(design = layout))
# 
# dev.off()

# 
# {
#   
#   #p1 = forestPlot(mafCompareRes=fvsm, pVal=0.05, color=c("maroon", "royalblue"), geneFontSize=0.8)
#   pdf(paste(dir,'/4_nonhyper_forestPlot_mutation_p0.001_fdr0.001.pdf',sep = ""), width = 12, height = 10)
#   forestPlot(mafCompareRes = fvsm, pVal=0.001,  fdr = 0.001,
#              color=c("maroon", "royalblue"), 
#              titleSize = 2,
#              geneFontSize = 1,
#              lineWidth = 1.5
#   )
#   dev.off()
#   
#   se_gene = getGeneSummary(clin.EOCRC)[1:30]$Hugo_Symbol
#   png(paste(dir,'/5_nonhyper_coOncoplot_top30.png',sep = ""),
#       width = 1500, height =1000)
#   #pdf(paste(dir,'/5_nonhyper_coOncoplot_top30.pdf',sep = ""),width = 15, height = 10)
#   coOncoplot(m1=clin.EOCRC, m2=clin.LOCRC, 
#              m1Name="EOCRC", m2Name="LOCRC", 
#              genes = se_gene)
#   dev.off()
#   
#   png(paste(dir,'/6_nonhyper_coOncoplot_diff_top30.png',sep = ""),
#       width = 1500, height =1000)
#   #pdf(paste(dir,'/6_nonhyper_coOncoplot_diff_top30.pdf',sep = ""), width = 15, height = 10)
#   coOncoplot(m1=clin.EOCRC, m2=clin.LOCRC, 
#              m1Name="EOCRC", m2Name="LOCRC", 
#              genes = genes) # [1:30]
#   dev.off()
# }
