
{
  # TMB nonhypermutated
  clin.EOCRC_nonhyper <- subsetMaf(maf=rt, tsb=TMB_EOCRC_nonhyper_sample, isTCGA=FALSE)
  clin.LOCRC_nonhyper <- subsetMaf(maf=rt, tsb=TMB_LOCRC_nonhyper_sample, isTCGA=FALSE)
  
  fvsm <- mafCompare(m1=clin.EOCRC_nonhyper, m2=clin.LOCRC_nonhyper, m1Name="EOCRC", m2Name="LOCRC", minMut=5)
  write.table(fvsm$results, file=paste(dir,'/nonhyperEOCRC_vs_LOCRC_counts.tsv',sep = ""), quote=FALSE, row.names=FALSE, sep="\t")
  
  fvsm$SampleSummary$SampleSize[1]
  op = fvsm$results
  op$EOCRC_freq = round(op$EOCRC/fvsm$SampleSummary$SampleSize[1], 4)
  op$LOCRC_freq = round(op$LOCRC/fvsm$SampleSummary$SampleSize[2], 4)
  genes = op$Hugo_Symbol[op$pval < 0.05] # adj
  genes
  
  write.table(op, file='./nonhyper_EOCRC_vs_LOCRC_freq.tsv', 
              quote=FALSE, row.names=FALSE, sep="\t")
  
  
  {
    
    
    df_nonhyper <- read_tsv('nonhyper_EOCRC_vs_LOCRC_freq.tsv')
    df_nonhyper <- df_nonhyper %>% 
      filter(!is.infinite(or) & !is.infinite(ci.low) & !is.infinite(ci.up))
    
    # set P_value and adjustPval if needed
    p_value_threshold <- 0.05
    adjustPval_threshold <- 0.05
    
    filtered_df_nonhyper <- df_nonhyper %>% 
      filter(pval <= p_value_threshold, adjPval <= adjustPval_threshold) %>%
      arrange(or) %>% head(20)
    
    # add new labels
    plot <- filtered_df_nonhyper %>% 
      mutate(EOCRC_freq = format_frequency(EOCRC_freq),
             LOCRC_freq = format_frequency(LOCRC_freq),
             sample_size = paste0(EOCRC_freq, " vs. ", LOCRC_freq),
             OR_CI = paste0(format(round(or, 2), nsmall = 2), " (", format(round(ci.low, 2), nsmall = 2), " to ", format(round(ci.up, 2), nsmall = 2), ")"),
             adjPval = format_adjPval(adjPval)) %>% 
      mutate(across(everything(), as.character)) %>% 
      bind_rows(data.frame(Hugo_Symbol = "Gene", sample_size = "Frequency (EOCRC vs LOCRC)", adjPval = "adjPval", OR_CI = "OR (95% CI)")) %>% 
      mutate(Hugo_Symbol = fct_rev(fct_relevel(Hugo_Symbol, "Gene")))
    
    
    p_mid <- filtered_df_nonhyper %>% 
      ggplot(aes(y = fct_rev(Hugo_Symbol))) +
      theme_classic() +
      geom_point(aes(x = or, size = or), shape = 23, fill = "black", show.legend = F) +
      geom_errorbarh(aes(xmin = ci.low, xmax = ci.up), height = 0.3) +
      labs(x = "Odds Ratio") +
      scale_x_continuous(trans = 'log10', breaks = c(0.01, 0.1, 0.3, 1), 
                         labels = scales::number_format(accuracy = 0.01)) +
      coord_cartesian(ylim = c(1, 21),xlim = c(0.01, 1)) + 
      
      geom_vline(xintercept = 1, linetype = "dashed") +  
      theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
            axis.text.y = element_blank(), axis.title.y = element_blank(),
            axis.text.x = element_text(color = "black"))
    
    
    p_left <- plot %>% 
      ggplot(aes(y = Hugo_Symbol)) + 
      geom_text(aes(x = 0, label = Hugo_Symbol), hjust = 0, fontface = "bold", size = 4) +
      geom_text(aes(x = 1, label = sample_size), hjust = 0, size = 4,
                fontface = ifelse(plot$sample_size == "Frequency (EOCRC vs LOCRC)", "bold", "plain")) +
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
  
  
  pdf(paste(dir,'/4_nonhyper_forestPlot.pdf',sep = ""), 
      width = 11, height = 6)
  
  print(p_left + p_mid + p_right + plot_layout(design = layout))
  
  dev.off()
  
  
  #p1 = forestPlot(mafCompareRes=fvsm, pVal=0.05, color=c("maroon", "royalblue"), geneFontSize=0.8)
  pdf(paste(dir,'/4_nonhyper_forestPlot_mutation_p0.001_fdr0.001.pdf',sep = ""), width = 12, height = 10)
  forestPlot(mafCompareRes = fvsm, pVal=0.001,  fdr = 0.001,
             color=c("maroon", "royalblue"), 
             titleSize = 2,
             geneFontSize = 1,
             lineWidth = 1.5
  )
  dev.off()
  
  se_gene = getGeneSummary(clin.EOCRC)[1:30]$Hugo_Symbol
  png(paste(dir,'/5_nonhyper_coOncoplot_top30.png',sep = ""),
      width = 1500, height =1000)
  #pdf(paste(dir,'/5_nonhyper_coOncoplot_top30.pdf',sep = ""),width = 15, height = 10)
  coOncoplot(m1=clin.EOCRC, m2=clin.LOCRC, 
             m1Name="EOCRC", m2Name="LOCRC", 
             genes = se_gene)
  dev.off()
  
  png(paste(dir,'/6_nonhyper_coOncoplot_diff_top30.png',sep = ""),
      width = 1500, height =1000)
  #pdf(paste(dir,'/6_nonhyper_coOncoplot_diff_top30.pdf',sep = ""), width = 15, height = 10)
  coOncoplot(m1=clin.EOCRC, m2=clin.LOCRC, 
             m1Name="EOCRC", m2Name="LOCRC", 
             genes = genes) # [1:30]
  dev.off()
}
