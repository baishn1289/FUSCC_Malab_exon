
{
  # TMB hypermutated
  clin.EOCRC_hyper <- subsetMaf(maf=rt, tsb=TMB_EOCRC_hyper_sample, isTCGA=FALSE)
  clin.LOCRC_hyper <- subsetMaf(maf=rt, tsb=TMB_LOCRC_hyper_sample, isTCGA=FALSE)
  
  
  fvsm <- mafCompare(m1=clin.EOCRC_hyper, m2=clin.LOCRC_hyper, 
                     m1Name="EOCRC", m2Name="LOCRC", minMut=5)
  write.table(fvsm$results, file=paste(dir,'/hyper_EOCRC_vs_LOCRC_counts.tsv',sep = ""),
              quote=FALSE, row.names=FALSE, sep="\t")
  
  fvsm$SampleSummary$SampleSize[1]
  op = fvsm$results
  op$EOCRC_freq = round(op$EOCRC/fvsm$SampleSummary$SampleSize[1], 4)
  op$LOCRC_freq = round(op$LOCRC/fvsm$SampleSummary$SampleSize[2], 4)
  genes = op$Hugo_Symbol[op$pval < 0.05] # adj
  genes
  
  write.table(op, file='./hyper_EOCRC_vs_LOCRC_freq.tsv', 
              quote=FALSE, row.names=FALSE, sep="\t")
  
  
  {
    df_hyper <- read_tsv('hyper_EOCRC_vs_LOCRC_freq.tsv')
    df_hyper <- df_hyper %>% 
      filter(!is.infinite(or) & !is.infinite(ci.low) & !is.infinite(ci.up))
    
    #add pvalue and adjustPvalue if need
    p_value_threshold <- 0.000005
    adjustPval_threshold <- 0.000005
    
    
    filtered_df_hyper <- df_hyper %>% 
      filter(pval <= p_value_threshold, adjPval <= adjustPval_threshold) %>%
      arrange(or) %>% head(20)
    
    # add new labels
    plot <- filtered_df_hyper %>% 
      mutate(EOCRC_freq = format_frequency(EOCRC_freq),
             LOCRC_freq = format_frequency(LOCRC_freq),
             sample_size = paste0(EOCRC_freq, " vs. ", LOCRC_freq),
             OR_CI = paste0(format(round(or, 2), nsmall = 2), " (", format(round(ci.low, 2), nsmall = 2), " to ", format(round(ci.up, 2), nsmall = 2), ")"),
             adjPval = format_adjPval(adjPval)) %>% 
      mutate(across(everything(), as.character)) %>% 
      bind_rows(data.frame(Hugo_Symbol = "Gene", sample_size = "Frequency (EOCRC vs LOCRC)", adjPval = "adjPval", OR_CI = "OR (95% CI)")) %>% 
      mutate(Hugo_Symbol = fct_rev(fct_relevel(Hugo_Symbol, "Gene")))
    
    # middle part
    p_mid <- filtered_df_hyper %>% 
      ggplot(aes(y = fct_rev(Hugo_Symbol))) +
      theme_classic() +
      geom_point(aes(x = or, size = or), shape = 23, fill = "black", show.legend = F) +
      geom_errorbarh(aes(xmin = ci.low, xmax = ci.up), height = 0.3) +
      labs(x = "Odds Ratio") +
      scale_x_continuous(trans = 'log10', breaks = c(0.1, 0.3, 1, 2, 3, 5, 10),
                         labels = scales::number_format(accuracy = 0.1)) + 
      coord_cartesian(ylim = c(1, 21), xlim = c(0.3, max(filtered_df_hyper$ci.up, na.rm = TRUE) + 1)) + 
      geom_vline(xintercept = 1, linetype = "dashed") +  
      theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
            axis.text.y = element_blank(), axis.title.y = element_blank(),
            axis.text.x = element_text(color = "black"))
    
    # left part
    p_left <- plot %>% 
      ggplot(aes(y = Hugo_Symbol)) + 
      geom_text(aes(x = 0, label = Hugo_Symbol), hjust = 0, fontface = "bold", size = 4) +
      geom_text(aes(x = 1, label = sample_size), hjust = 0, size = 4,
                fontface = ifelse(plot$sample_size == "Frequency (EOCRC vs LOCRC)", "bold", "plain")) +
      theme_void() + coord_cartesian(xlim = c(0, 3))
    
    # right_part
    p_right <- plot %>% ggplot() +
      geom_text(aes(x = 0, y = Hugo_Symbol, label = OR_CI), hjust = 0,
                fontface = ifelse(plot$OR_CI == "OR (95% CI)", "bold", "plain"), size = 4) +
      geom_text(aes(x = 1, y = Hugo_Symbol, label = adjPval), hjust = 0,
                fontface = ifelse(plot$adjPval == "adjPval", "bold", "plain"), size = 4) +
      theme_void() + coord_cartesian(xlim = c(0, 2))
    
    
    layout <- c(
      area(t = 0, l = 0, b = 30, r = 8), 
      area(t = 0, l = 9, b = 30, r = 15), 
      area(t = 0, l = 16, b = 30, r = 22)) 
    
    
  }
  pdf(paste(dir,'/4_hyper_forestPlot.pdf',sep = ""), width = 11, height = 6)
  
  print(p_left + p_mid + p_right + plot_layout(design = layout))
  
  dev.off()
   
}

