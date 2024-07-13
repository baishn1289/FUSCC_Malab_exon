
library(maftools)
library(data.table)
library(dplyr)
library(mclust)
library(NMF)
library(pheatmap)
library(data.table)
library(stringr)
library(Cairo)

dir= './picture'

#load()

{# MafSummary
  se_gene = getGeneSummary(clin.EOCRC)[1:50]$Hugo_Symbol
  
  CairoPDF(paste(dir,'/coOncoplot_top50_Cairo.pdf',sep = ""),
           width = 20, height = 10,)
  coOncoplot(m1=clin.EOCRC, m2=clin.LOCRC,
             m1Name="EOCRC", m2Name="LOCRC", 
             genes = se_gene)
  dev.off()
  
  pdf(paste(dir,'/1_plotmafSummary.pdf',sep = ""), width = 12, height = 5)
  plotmafSummary(maf = rt)
  dev.off()
  
  CairoPDF(paste(dir,'/2_oncoplot_Cairo.pdf',sep = ""), width = 12, height = 6)

  oncoplot(maf = rt,
           top = 30,
           fontSize = 0.6,
           showTumorSampleBarcodes = F,
           sepwd_samples = 0,
           clinicalFeatures = c('Age_class','country'),
           bgCol = "#FFFFFF",
           sortByAnnotation = TRUE,
           groupAnnotationBySize = FALSE)
  dev.off()
}


{
  captureSize = 2
  pdf(paste(dir,'/3_TMB.pdf',sep = ""), width = 4, height = 4)
  TMBdata = tmb(maf = rt, captureSize = captureSize)
  head(TMBdata)
  dev.off()
  
  pdf(paste(dir,'/3_TMB_EOCRC.pdf',sep = ""), width = 4, height = 4)
  tmb(maf = clin.EOCRC, captureSize = captureSize, logScale = T)
  dev.off()
  
  pdf(paste(dir,'/3_TMB_LOCRC.pdf',sep = ""), width = 4, height = 4)
  tmb(maf = clin.LOCRC, captureSize = captureSize, logScale = T)
  dev.off()
}


{
  maf.mutload = getSampleSummary(clin.EOCRC)[, .(Tumor_Sample_Barcode, 
                                                 total)]
  maf.mutload[, ':='(total_perMB, total/captureSize)]
  maf.mutload[, `:=`(total_perMB_log, log10(total_perMB))]
  maf.mutload = maf.mutload[order(total_perMB, decreasing = FALSE)]
  
  TMB_EOCRC = maf.mutload$total_perMB_log
  TMB_EOCRC_hyper_sample = maf.mutload$Tumor_Sample_Barcode[maf.mutload$total_perMB_log > 1]
  TMB_EOCRC_nonhyper_sample = maf.mutload$Tumor_Sample_Barcode[maf.mutload$total_perMB_log <= 1]
  
  maf.mutload = getSampleSummary(clin.LOCRC)[, .(Tumor_Sample_Barcode, 
                                                 total)]
  maf.mutload[, `:=`(total_perMB, total/captureSize)]
  maf.mutload[, `:=`(total_perMB_log, log10(total_perMB))]
  maf.mutload = maf.mutload[order(total_perMB, decreasing = FALSE)]
  TMB_LOCRC = maf.mutload$total_perMB_log
  TMB_LOCRC_sample = maf.mutload$Tumor_Sample_Barcode[maf.mutload$total_perMB_log > 1]
  TMB_LOCRC_hyper_sample = maf.mutload$Tumor_Sample_Barcode[maf.mutload$total_perMB_log > 1]
  TMB_LOCRC_nonhyper_sample = maf.mutload$Tumor_Sample_Barcode[maf.mutload$total_perMB_log <= 1]
}

{
  hist((TMB_EOCRC), breaks = 200,plot=T)
  hist((TMB_LOCRC), breaks = 200,plot=T)
  
  
  lower_bound_EOCRC = quantile(TMB_EOCRC, 0.05)
  upper_bound_EOCRC = quantile(TMB_EOCRC, 0.95)
  TMB_EOCRC_filtered = TMB_EOCRC[TMB_EOCRC >= lower_bound_EOCRC & TMB_EOCRC <= upper_bound_EOCRC]
  lower_bound_LOCRC = quantile(TMB_LOCRC, 0.05)
  upper_bound_LOCRC = quantile(TMB_LOCRC, 0.95)
  TMB_LOCRC_filtered = TMB_LOCRC[TMB_LOCRC >= lower_bound_LOCRC & TMB_LOCRC <= upper_bound_LOCRC]
  t.test(TMB_EOCRC_filtered, TMB_LOCRC_filtered)
  boxplot(TMB_EOCRC_filtered, TMB_LOCRC_filtered)
}

library(magrittr) 
library(ggsignif)
library(tidyverse)
library(readxl)
library(ggpubr)
library(rstatix)

{#***Low TMB
  table(TMB_EOCRC <= 1)
  table(TMB_LOCRC <= 1)
  x = TMB_EOCRC[TMB_EOCRC <= 1]
  y = TMB_LOCRC[TMB_LOCRC <= 1]
  t.test(x, y)
  wilcox.test(x, y)
  #x = TMB_EOCRC_filtered[TMB_EOCRC_filtered <= 1]
  #y = TMB_LOCRC_filtered[TMB_LOCRC_filtered <= 1]
  #t.test(x, y)
  #wilcox.test(x, y)
  
  
  nonhyper_TMB_EOCRC = data.frame(TMB = x,class='EOCRC')
  nonhyper_TMB_LOCRC = data.frame(TMB = y,class='LOCRC')
  nonhyper_TMB = rbind(nonhyper_TMB_EOCRC,nonhyper_TMB_LOCRC) 
  nonhyper_TMB$class <-as.factor(nonhyper_TMB$class)
  
  palette <- c('#deebf7','#3182bd')
  
  picture_nonhyper <-  
    ggplot(nonhyper_TMB, aes(x = class, y = TMB, fill = class)) +
    geom_boxplot(outlier.shape = T, width = 0.5) + 
    #geom_jitter(width = 0.03, size = 0.1) +
    #stat_boxplot(geom = "errorbar", width = 0.2) + 
    #scale_fill_npg() + 
    #ylim(-0.5,1.3)+
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1),
                       limits = c(-0.5,1.3)) +
    scale_fill_manual(values = palette) +
    theme_pubr() +
    
    labs(title = "TMB_EOCRC_vs_LOCRC", 
         subtitle = "nonhypermutated", 
         x = "",
         y = "TMB/MB(log10)"
    ) +
    geom_signif(comparisons =  list(c("EOCRC", "LOCRC")),
                y_position = 1.1,
                tip_length = 0.02,
                #map_signif_level = `p`, 
                test = t.test 
    ) +
    theme(
      plot.title    = element_text(color = "black", size   = 16, hjust = 0.5),
      plot.subtitle = element_text(color = "black", size   = 16,hjust = 0.5),
      plot.caption  = element_text(color = "black", size   = 16,face = "italic", hjust = 1),
      axis.text.x   = element_text(color = "black", size = 16, angle = 0),
      axis.text.y   = element_text(color = "black", size = 16, angle = 0),
      axis.title.x  = element_text(color = "black", size = 16, angle = 0),
      axis.title.y  = element_text(color = "black", size = 16, angle = 90),
      legend.position = "none",
      #legend.title  = element_text(color = "black", size  = 16),
      #legend.text   = element_text(color = "black", size   = 16),
      axis.line.y = element_line(color = "black", linetype = "solid"), 
      axis.line.x = element_line (color = "black",linetype = "solid"), 
      panel.border = element_rect(linetype = "solid", linewidth = 1.2,fill = NA),
      
      
    )
 
}

pdf(paste(dir,'/2_TMB_nonhyper_boxplot.pdf',sep = ""),width=4,height=5)
picture_nonhyper
dev.off()


{#*** High TMB
  table(TMB_EOCRC > 1)
  table(TMB_LOCRC > 1)
  x = TMB_EOCRC[TMB_EOCRC > 1]
  y = TMB_LOCRC[TMB_LOCRC > 1]
  #t.test(x, y)
  #wilcox.test(x, y)
  #x = TMB_EOCRC_filtered[TMB_EOCRC_filtered > 1]
  #y = TMB_LOCRC_filtered[TMB_LOCRC_filtered > 1]
  #t.test(x, y)
  #wilcox.test(x, y)
  
  hyper_TMB_EOCRC = data.frame(TMB = x,class='EOCRC')
  hyper_TMB_LOCRC = data.frame(TMB = y,class='LOCRC')
  hyper_TMB = rbind(hyper_TMB_EOCRC,hyper_TMB_LOCRC) 
  hyper_TMB$class <-as.factor(hyper_TMB$class)
  
  picture_hyper <- ggplot(hyper_TMB, aes(x = class, y = TMB, fill = class)) +
    geom_boxplot(outlier.shape = T, width =0.5) + 
    #geom_jitter(width = 0.03, size = 0.1) + 
    stat_boxplot(geom = "errorbar", width = 0.2) +
    #scale_fill_npg() + 
    scale_fill_manual(values = palette)+
    #ylim(1.0,4.5)+
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1),
                       limits = c(1.0,4.5)) +
    
    theme_pubr() +
    
    labs(title = "TMB_EOCRC_vs_LOCRC", 
         subtitle = "hypermutated", 
         x = "", 
         y = "TMB/MB(log10)"
         
    ) +
    geom_signif(comparisons =  list(c("EOCRC", "LOCRC")),
                y_position = 4.2,
                tip_length = 0.02,
                #map_signif_level = p, 
                test = t.test
    ) +
    theme(
      plot.title    = element_text(color = "black", size   = 16, hjust = 0.5),
      plot.subtitle = element_text(color = "black", size   = 16,hjust = 0.5),
      plot.caption  = element_text(color = "black", size   = 16,face = "italic", hjust = 1),
      axis.text.x   = element_text(color = "black", size = 16, angle = 0),
      axis.text.y   = element_text(color = "black", size = 16, angle = 0),
      axis.title.x  = element_text(color = "black", size = 16, angle = 0),
      axis.title.y  = element_text(color = "black", size = 16, angle = 90),
      legend.position = "none",
      #legend.title  = element_text(color = "black", size  = 16),
      #legend.text   = element_text(color = "black", size   = 16),
      axis.line.y = element_line(color = "black", linetype = "solid"), 
      axis.line.x = element_line (color = "black",linetype = "solid"), 
      panel.border = element_rect(linetype = "solid", linewidth = 1.2,fill = NA)
      
    )
}

pdf(paste(dir,'/2_TMB_hyper_boxplot.pdf',sep = ""), width=4, height=5)
picture_hyper
dev.off()
  
