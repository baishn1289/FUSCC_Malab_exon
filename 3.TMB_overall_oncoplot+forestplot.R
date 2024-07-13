library(tidyverse)
library(patchwork)
library(scales)

#load()

format_pvalue <- function(p) {
  ifelse(p < 0.001, format(p, scientific = TRUE, digits = 3), format(p, digits = 3, nsmall = 3))
}


format_frequency <- function(freq) {
  paste0(format(round(freq * 100, 1), nsmall = 1), "%")
}


format_adjPval <- function(p) {
  format(p, digits = 3, scientific = TRUE)
}

{
  fvsm <- mafCompare(m1=clin.EOCRC, m2=clin.LOCRC, 
                     m1Name="EOCRC", m2Name="LOCRC", minMut=5)
  
  write.table(fvsm$results, 
              file=paste(dir,'/overall_EOCRC_vs_LOCRC_counts.tsv',sep = ""), 
              quote=FALSE, row.names=FALSE, sep="\t")
  
  fvsm$SampleSummary$SampleSize[1]
  op = fvsm$results
  op$EOCRC_freq = round(op$EOCRC/fvsm$SampleSummary$SampleSize[1], 4)
  op$LOCRC_freq = round(op$LOCRC/fvsm$SampleSummary$SampleSize[2], 4)
  genes = op$Hugo_Symbol[op$adjPval < 0.05] # adj
  
  write.table(op, file='./overall_EOCRC_vs_LOCRC_freq.tsv', 
              quote=FALSE, row.names=FALSE, sep="\t")
  
  {#overall_forest
    
    df_overall <- read_tsv('overall_EOCRC_vs_LOCRC_freq.tsv')
    
    df_overall <- df_overall %>% 
      filter(!is.infinite(or) & !is.infinite(ci.low) & !is.infinite(ci.up))
    
    #set p_value and adjustPval if need
    p_value_threshold <- 0.005
    adjustPval_threshold <- 0.005
    
    filtered_df_overall <- df_overall %>% 
      filter(pval <= p_value_threshold, adjPval <= adjustPval_threshold) %>%
      arrange(or) %>% head(20)
    
    
    # add new labels
    plot <- filtered_df_overall %>% 
      mutate(EOCRC_freq = format_frequency(EOCRC_freq),
             LOCRC_freq = format_frequency(LOCRC_freq),
             sample_size = paste0(EOCRC_freq, " vs. ", LOCRC_freq),
             OR_CI = paste0(format(round(or, 2), nsmall = 2), " (", format(round(ci.low, 2), nsmall = 2), " to ", format(round(ci.up, 2), nsmall = 2), ")"),
             adjPval = format_adjPval(adjPval)) %>% 
      mutate(across(everything(), as.character)) %>% 
      bind_rows(data.frame(Hugo_Symbol = "Gene", sample_size = "Frequency (EOCRC vs LOCRC)", adjPval = "adjPval", OR_CI = "OR (95% CI)")) %>% 
      mutate(Hugo_Symbol = fct_rev(fct_relevel(Hugo_Symbol, "Gene")))
    
    # middle part
    p_mid <- filtered_df_overall %>% 
      ggplot(aes(y = fct_rev(Hugo_Symbol))) +
      theme_classic() +
      geom_point(aes(x = or, size = or), shape = 23, fill = "black", show.legend = F) +
      geom_errorbarh(aes(xmin = ci.low, xmax = ci.up), height = 0.3) +
      labs(x = "Odds Ratio") +
      scale_x_continuous(trans = 'log10', breaks = c(0.1, 0.3, 1, 2, 5, 10,20),
                         labels = scales::number_format(accuracy = 0.1)) + # 使用对数尺度，并设置刻度
      coord_cartesian(ylim = c(1, 21), xlim = c(0.3, max(filtered_df_overall$ci.up, na.rm = TRUE) + 1)) + # 调整 xlim 范围
      geom_vline(xintercept = 1, linetype = "dashed") +  # 修改为 xintercept = 1
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
    
    # right part
    p_right <- plot %>% ggplot() +
      geom_text(aes(x = 0, y = Hugo_Symbol, label = OR_CI), hjust = 0,
                fontface = ifelse(plot$OR_CI == "OR (95% CI)", "bold", "plain"), size = 4) +
      geom_text(aes(x = 1, y = Hugo_Symbol, label = adjPval), hjust = 0,
                fontface = ifelse(plot$adjPval == "adjPval", "bold", "plain"), size = 4) +
      theme_void() + coord_cartesian(xlim = c(0, 2))
    
    
    layout <- c(
      area(t = 0, l = 0, b = 30, r = 8), # 调整左侧标签宽度
      area(t = 0, l = 9, b = 30, r = 16), # 调整中间部分宽度
      area(t = 0, l = 17, b = 30, r = 24)) # 调整右侧标签宽度
    
    pdf(paste(dir,'/4_overall_forestPlot.pdf',sep = ""), 
        width = 11, height = 6)
    print(p_left + p_mid + p_right + plot_layout(design = layout))
    
    dev.off()
  }
}
