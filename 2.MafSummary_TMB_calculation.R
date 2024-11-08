library(Cairo)
library(ggthemes)
library(ggpubr)
library(effsize)
library(magrittr) 
library(ggsignif)
library(readxl)
library(ggpubr)
library(rstatix)



dir= './picture'

{# MafSummary
  se_gene = getGeneSummary(clin.EOCRC)[1:50]$Hugo_Symbol
  
  # CairoPDF(paste(dir,'/coOncoplot_top50_Cairo.pdf',sep = ""),
  #          width = 20, height = 10,)
  # coOncoplot(m1=clin.EOCRC, m2=clin.LOCRC,
  #            m1Name="EOCRC", m2Name="LOCRC", 
  #            genes = se_gene)
  # dev.off()
  
  # pdf(paste(dir,'/1_plotmafSummary.pdf',sep = ""), width = 12, height = 5)
  # plotmafSummary(maf = rt)
  # dev.off()
  
  CairoPDF(paste(dir,'/2_oncoplot_Cairo.pdf',sep = ""), width = 12, height = 6)
  
  oncoplot(maf = rt_panel_common_genes,
           top = 30,
           fontSize = 0.6,
           showTumorSampleBarcodes = F,
           sepwd_samples = 0,
           clinicalFeatures = c('Age_class','Country'),
           bgCol = "#FFFFFF",
           sortByAnnotation = TRUE,
           groupAnnotationBySize = FALSE)
  dev.off()
}




{
 
  # rt@clinical.data <-  rt@clinical.data %>%
  #                       mutate(Race = ifelse(Race == 'Unknown',NA,Race),
  #                              Sex = ifelse(Sex == 'Unknown',NA,Sex),
  #                              Histology_type =  ifelse(Histology_type == 'Unknown',NA,Histology_type),
  #                              MSI_Status = ifelse(MSI_Status == 'Unknown',NA,MSI_Status),
  #                              Tumor_site = ifelse(Tumor_site == 'Unknown',NA,Tumor_site)
  #                            )

  #EOCRC_TMB_calculation
  TMB_EOCRC_mutload <-  left_join(getSampleSummary(clin.EOCRC)[, .(Tumor_Sample_Barcode,total)],
                            rt@clinical.data[,.(Number,Tumor_Sample_Barcode,Panel)],
                            by='Tumor_Sample_Barcode')
  
  TMB_EOCRC_mutload <-  left_join(TMB_EOCRC_mutload,capture_size,by = 'Panel')
  
  table(is.na(TMB_EOCRC_mutload$Panel))
  TMB_EOCRC_mutload$TMB <- (TMB_EOCRC_mutload$total) / TMB_EOCRC_mutload$total_size
  
  TMB_EOCRC_mutload$logTMB <-  log10(TMB_EOCRC_mutload$TMB)
  TMB_EOCRC_mutload = TMB_EOCRC_mutload[order(TMB, decreasing = FALSE)]
  
  TMB_EOCRC = TMB_EOCRC_mutload$logTMB
  TMB_EOCRC_hyper_sample = TMB_EOCRC_mutload$Tumor_Sample_Barcode[TMB_EOCRC_mutload$logTMB > 1]
  TMB_EOCRC_nonhyper_sample = TMB_EOCRC_mutload$Tumor_Sample_Barcode[TMB_EOCRC_mutload$logTMB <= 1]
 
  #LOCRC_TMB_calculation
  TMB_LOCRC_mutload <-  left_join(getSampleSummary(clin.LOCRC)[, .(Tumor_Sample_Barcode,total)],
                            rt@clinical.data[,.(Number,Tumor_Sample_Barcode,Panel)],
                            by='Tumor_Sample_Barcode')
  
  TMB_LOCRC_mutload <-  left_join(TMB_LOCRC_mutload,capture_size,by = 'Panel')
  
  table(is.na(TMB_LOCRC_mutload$Panel))
  TMB_LOCRC_mutload$TMB <- (TMB_LOCRC_mutload$total) / TMB_LOCRC_mutload$total_size
  
  TMB_LOCRC_mutload$logTMB <-  log10(TMB_LOCRC_mutload$TMB)
  TMB_LOCRC_mutload = TMB_LOCRC_mutload[order(TMB, decreasing = FALSE)]
  
  TMB_LOCRC = TMB_LOCRC_mutload$logTMB
  TMB_LOCRC_hyper_sample = TMB_LOCRC_mutload$Tumor_Sample_Barcode[TMB_LOCRC_mutload$logTMB > 1]
  TMB_LOCRC_nonhyper_sample = TMB_LOCRC_mutload$Tumor_Sample_Barcode[TMB_LOCRC_mutload$logTMB <= 1]
  

  # 加入临床特征构建线性回归模型
  # 此处TMB_mutload没有经过1%的首尾筛选
  TMB_mutload = rbind(TMB_LOCRC_mutload,TMB_EOCRC_mutload) %>% 
    left_join(rt@clinical.data[,.(Age,Age_class,Country,Center,
                                  Sex,Race,Tumor_site,
                                  Histology_type,Sample_type,Tumor_Sample_Barcode)]
              ,by='Tumor_Sample_Barcode')

 
  
  # model_overall = lm(TMB~total_size+Country+Race+Sex+Tumor_site+Sample_type+Histology_type,
  #                    data = TMB_mutload)
  # 
  # 
  # TMB_mutload_regression = TMB_mutload %>% 
  #   mutate(residuals = residuals(model_overall),
  #          logresiduals = log10(residuals)) 
  
    # filter(residuals > quantile(residuals, 0.01) & residuals < quantile(residuals, 0.99))
  
  #EO和LO比较时，可能会相互干扰对方组的回归，并且如果纳入Age_class变量相当于消除年龄组之间的异质性
  #因此两个年龄组做隔离单独回归，单独建模
  EO_regression = TMB_mutload[Age_class == 'EOCRC']
  LO_regression = TMB_mutload[Age_class == 'LOCRC']
  
  model_EO <- lm(TMB~Panel+Country+Race+Sex+Tumor_site+Sample_type+Histology_type,
                 data = EO_regression)
  model_LO <- lm(TMB~Panel+Country+Race+Sex+Tumor_site+Sample_type+Histology_type,
                 data = LO_regression)
  
  # library(lme4)
  # model_mixed_EO <- lmer(TMB ~ Panel + Race + Sex + Tumor_site +# Histology_type + #Sample_type + 
  #                          (1 | Center),  
  #                        data = EO_regression)
  # 
  # model_mixed_LO <- lmer(TMB ~ Panel + Race +Sex + Tumor_site + #Histology_type + #Sample_type + 
  #                          (1 | Center),  
  #                        data = LO_regression)
  # 
  # AIC(model_mixed_EO)
  # AIC(model_mixed_LO)
  
  # library(quantreg)
  # model_quantile_EO <- rq(TMB ~ Panel + Race + Sex + Tumor_site + Sample_type + Histology_type,
  #                      data = EO_regression, tau = 0.5)
  # model_quantile_LO <- rq(TMB ~ Panel,# + Race + Sex + Tumor_site + Sample_type + Histology_type,
  #                         data = LO_regression, tau = 0.5)
  
  # EO_regression = EO_regression %>% 
  #   mutate(residuals = residuals(model_EO)) %>% 
  #   filter(residuals > quantile(residuals, 0.05) & residuals < quantile(residuals, 0.95))
  # 
  # LO_regression = LO_regression %>% 
  #   mutate(residuals = residuals(model_LO)) %>%  
  #   filter(residuals > quantile(residuals, 0.05) & residuals < quantile(residuals, 0.95))
  
  EO_regression = EO_regression %>% 
    mutate(residuals = residuals(model_EO)) 
    # filter(residuals > quantile(residuals, 0.01) & residuals < quantile(residuals, 0.99))
  
  LO_regression = LO_regression %>% 
    mutate(residuals = residuals(model_LO)) 
    # filter(residuals > quantile(residuals, 0.01) & residuals < quantile(residuals, 0.99))
 
  # EO_mixed_regression = EO_regression %>% 
  #   mutate(residuals = residuals(model_mixed_EO)) %>% 
  #   filter(residuals > quantile(residuals, 0.01) & residuals < quantile(residuals, 0.99))
  # 
  # LO_mixed_regression = LO_regression %>% 
  #   mutate(residuals = residuals(model_mixed_LO)) %>%  
  #   filter(residuals > quantile(residuals, 0.01) & residuals < quantile(residuals, 0.99))
  # 

  
  summary(EO_regression[TMB>10]$residuals)
  summary(LO_regression[TMB>10]$residuals)
  summary(EO_regression[TMB<=10]$residuals)
  summary(LO_regression[TMB<=10]$residuals)

  wilcox.test(EO_regression$residuals,LO_regression$residuals)

  
  wilcox.test(EO_regression[residuals>10]$residuals,
              LO_regression[residuals>10]$residuals)
  wilcox.test(EO_regression[residuals<=10]$residuals,
              LO_regression[residuals<=10]$residuals)
  
  summary(model_EO)$r.squared
  summary(model_LO)$r.squared
  summary(model_overall)$r.squared
  
  write.csv(TMB_mutload,file = 'TMB_mutload.csv')
  write.csv(EO,file = 'TMB_EOCRC_mutload_regression.csv')
  write.csv(LO,file = 'TMB_LOCRC_mutload_regression.csv')
  
 
}



{
  hist((TMB_EOCRC), breaks = 200,plot=T)
  hist((TMB_LOCRC), breaks = 200,plot=T)
 
}


palette <- c('#deebf7','#3182bd')

{#TMB_across_panels_for_EOandLO
  TMB_panel_filtered <- TMB_mutload_regression
    
  picture_acorss_panels <-  
    ggplot(TMB_panel_filtered, aes(x = Age_class, y = residuals, fill = Age_class)) +
    geom_boxplot(outlier.shape = T, width = 0.5) + 
    facet_wrap(~ Panel, scales = "free_y",nrow = 3) +  
    theme_few() +
    scale_fill_manual(values = palette) +
    theme(
      # plot.title    = element_text(color = "black", size   = 10, hjust = 0.5),
      # plot.subtitle = element_text(color = "black", size   = 10,hjust = 0.5),
      # plot.caption  = element_text(color = "black", size   = 10,face = "italic", hjust = 1),
      axis.text.x   = element_text(color = "black", size = 10, angle = 45,hjust = 1),
      axis.text.y   = element_text(color = "black", size = 10, angle = 0),
      axis.title.x  = element_text(color = "black", size = 14, angle = 0),
      axis.title.y  = element_text(color = "black", size = 14, angle = 90),
      legend.position = 'right',
      #legend.title  = element_text(color = "black", size  = 10),
      #legend.text   = element_text(color = "black", size   = 10),
      axis.line.y = element_line(color = "black", linetype = "solid"), 
      axis.line.x = element_line (color = "black",linetype = "solid"), 
      panel.border = element_rect(linetype = "solid", linewidth = 1.2,fill = NA),
      strip.text = element_text(size = 8, color = "black", face = "bold") 
      
    )
  
  pdf(paste(dir,'/2.TMB_across_panels_for_EOandLO.pdf',sep =''),width = 22,height = 6)
  picture_acorss_panels
  dev.off()
  
  
}



{#***nonhyper original_TMB
  table(TMB_EOCRC <= 1)
  table(TMB_LOCRC <= 1)
  x = TMB_EOCRC[TMB_EOCRC <= 1]
  y = TMB_LOCRC[TMB_LOCRC <= 1]
  
  t.test(x, y)
  wilcox.test(x, y)
  # cliff.delta(x, y)
  TMB_nonhyper_FC <- mean(y)/mean(x)
  
  summary(x)
  summary(y)
  #x = TMB_EOCRC_filtered[TMB_EOCRC_filtered <= 1]
  #y = TMB_LOCRC_filtered[TMB_LOCRC_filtered <= 1]
  #t.test(x, y)
  #wilcox.test(x, y)
  
  
  nonhyper_TMB_EOCRC = data.frame(TMB = x,class='EOCRC')
  nonhyper_TMB_LOCRC = data.frame(TMB = y,class='LOCRC')
  nonhyper_TMB = rbind(nonhyper_TMB_EOCRC,nonhyper_TMB_LOCRC) 
  nonhyper_TMB$class <- factor(nonhyper_TMB$class,levels = c('LOCRC','EOCRC'))
  
  picture_nonhyper <-  
    ggplot(nonhyper_TMB, aes(x = class, y = TMB, fill = class)) +
    geom_boxplot(outlier.shape = T, width = 0.5) + 
    #geom_jitter(width = 0.03, size = 0.1) +
    #stat_boxplot(geom = "errorbar", width = 0.2) + 
    #scale_fill_npg() + 
    #ylim(-0.5,1.3)+
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1),
                       limits = c(-0.7,1.3)) +
    scale_fill_manual(values = palette) +
    theme_pubr() +
    
    labs(title = "TMB_EOCRC_vs_LOCRC", 
         subtitle = "nonhypermutated", 
         x = "",
         y = "log10(TMB)"
    ) +
    geom_signif(comparisons =  list(c("EOCRC", "LOCRC")),
                y_position = 1.1,
                tip_length = 0.02,
                #map_signif_level = `p`, 
                test = t.test, size= 0.5,textsize = 5
    ) +
    annotate("text", 
             x = 0.7,  
             y = 1.3, 
             label = paste("FC =", round(TMB_nonhyper_FC, 2)), 
             size = 5, 
             color = "black") +
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

pdf(paste(dir,'/2_original_TMB_nonhyper_boxplot.pdf',sep = ""),width=4.5,height=5)
picture_nonhyper
dev.off()


{#*** High original_TMB
  table(TMB_EOCRC > 1)
  table(TMB_LOCRC > 1)
  x = TMB_EOCRC[TMB_EOCRC > 1]
  y = TMB_LOCRC[TMB_LOCRC > 1]
  
  
  t.test(x, y)
  wilcox.test(x, y)
  # cliff.delta(x, y)
  TMB_hyper_FC <- mean(y)/mean(x)
  summary(x)
  summary(y)
  #x = TMB_EOCRC_filtered[TMB_EOCRC_filtered > 1]
  #y = TMB_LOCRC_filtered[TMB_LOCRC_filtered > 1]
  #t.test(x, y)
  #wilcox.test(x, y)
  
  hyper_TMB_EOCRC = data.frame(TMB = x,class='EOCRC')
  hyper_TMB_LOCRC = data.frame(TMB = y,class='LOCRC')
  hyper_TMB = rbind(hyper_TMB_EOCRC,hyper_TMB_LOCRC) 
  hyper_TMB$class <- factor(hyper_TMB$class,levels = c('LOCRC','EOCRC'))
  
  table(hyper_TMB$class)
  table(hyper_TMB$class)/nrow(hyper_TMB)
  
  
  picture_hyper <- ggplot(hyper_TMB, aes(x = class, y = TMB, fill = class)) +
    geom_boxplot(outlier.shape = T, width =0.5) + 
    #geom_jitter(width = 0.03, size = 0.1) + 
    #stat_boxplot(geom = "errorbar", width = 0.2) +
    #scale_fill_npg() + 
    scale_fill_manual(values = palette)+
    #ylim(1.0,4.5)+
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1),
                       limits = c(1.0,3.4)) +
    
    theme_pubr() +
    
    labs(title = "TMB_EOCRC_vs_LOCRC", 
         subtitle = "hypermutated", 
         x = "", 
         y = "log10(TMB)"
         
    ) +
    annotate("text", 
             x = 0.7,  
             y = 3.4, 
             label = paste("FC =", round(TMB_hyper_FC , 2)), 
             size = 5, 
             color = "black")+
    geom_signif(comparisons =  list(c("EOCRC", "LOCRC")),
                y_position = 3.2,
                tip_length = 0.02,
                #map_signif_level = p, 
                test = t.test,textsize = 5
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

pdf(paste(dir,'/2_original_TMB_hyper_boxplot.pdf', sep=''),width=4.5, height=5)
picture_hyper
dev.off()



{#*** overall original_TMB
  
  x = TMB_EOCRC
  y = TMB_LOCRC
 
  t.test(x, y)
  wilcox.test(x, y)
  # cliff.delta(x, y)
  
  TMB_overall_FC <- mean(y)/mean(x)
  summary(x)
  summary(y)
  # x = TMB_EOCRC_filtered[TMB_EOCRC_filtered > 1]
  # y = TMB_LOCRC_filtered[TMB_LOCRC_filtered > 1]
  # t.test(x, y)
  # wilcox.test(x, y)

  
  
  overall_TMB_EOCRC = data.frame(TMB = x,class='EOCRC')
  overall_TMB_LOCRC = data.frame(TMB = y,class='LOCRC')
  overall_TMB = rbind(overall_TMB_EOCRC,overall_TMB_LOCRC) 
  overall_TMB$class <-factor(overall_TMB$class,levels = c('LOCRC','EOCRC'))
  
  
  picture <- ggplot(overall_TMB, aes(x = class, y = TMB, fill = class)) +
    geom_boxplot(outlier.shape = T,
                 width =0.5) + 
    #geom_jitter(width = 0.03, size = 0.1) + 
    #stat_boxplot(geom = "errorbar", width = 0.2) +
    #scale_fill_npg() + 
    scale_fill_manual(values = palette)+
    #ylim(1.0,4.5)+
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1),
                       limits = c(-0.7,2.4)) +
    
    theme_pubr() +
    
    labs(title = "TMB_EOCRC_vs_LOCRC", 
         subtitle = "overall", 
         x = "", 
         y = "log10(TMB)"
         
    ) +
    geom_signif(comparisons =  list(c("EOCRC", "LOCRC")),
                y_position = 2.2,
                tip_length = 0.02,
                #map_signif_level = p, 
                test = t.test,textsize = 5
    ) +
    # 添加 Fold Change 值到图形中
    annotate("text", 
             x = 0.7,  # 放在 EOCRC 和 LOCRC 中间
             y = 2.4,  # 设置文本显示的位置
             label = paste("FC =", round(TMB_overall_FC, 2)),  # 显示 Fold Change，保留两位小数
             size = 5, 
             color = "black") +
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

pdf(paste(dir,'/2_original_TMB_overall_boxplot.pdf', sep = ''),width=4.5, height=5)
picture
dev.off() 



{
  #nonhyper regression_TMB
  nonhyper_regression_TMB = rbind(EO_regression[TMB<=10],LO_regression[TMB>10]) %>% 
    mutate(Age_class = factor(Age_class,levels = c('LOCRC','EOCRC')))
  
  # nonhyper_regression_FC <- mean(LO_regression[residuals<=5]$residuals)/mean(EO_regression[residuals<=5]$residuals)
  
  summary(LO_regression[residuals<=5]$residuals)
  summary(EO_regression[residuals<=5]$residuals)
  
  cliff.delta(LO_regression[residuals<=5]$residuals, 
              EO_regression[residuals<=5]$residuals)
  
  picture_nonhyper_regression <- ggplot(nonhyper_regression_TMB, 
                                        aes(x = Age_class, 
                                            y = residuals, fill = Age_class)) +
    geom_boxplot(outlier.shape = T,
                 width =0.5) +
    scale_fill_manual(values = palette)+
    scale_y_continuous(labels = scales::number_format(accuracy = 2),
                       limits = c(-32,8)) +
    
    theme_pubr() +
    
    labs(title = "TMB_EOCRC_vs_LOCRC", 
         subtitle = "non-hypermutated", 
         x = "", 
         y = "adjusted TMB(residuals)"
         
    ) +
    geom_signif(comparisons =  list(c("EOCRC", "LOCRC")),
                y_position = 6,
                tip_length = 0.02,
                #map_signif_level = p, 
                test =  wilcox.test,
                textsize = 5
    ) +
    # 添加 Fold Change 值到图形中
    # annotate("text", 
    #          x = 0.7,  # 放在 EOCRC 和 LOCRC 中间
    #          y = 10,  # 设置文本显示的位置
    #          label = paste("FC =", round(nonhyper_regression_FC, 2)),  # 显示 Fold Change，保留两位小数
    #          size = 5, 
    #          color = "black") +
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
  
  pdf(file = paste(dir,'/2.regression_nonhyepr_TMB.pdf',sep = ''),
      width = 4.5,height = 5.5)
  picture_nonhyper_regression
  dev.off()
  
  
}


{
  #hyper regression_TMB
  hyper_regression_TMB = rbind(EO_regression[residuals>5],LO_regression[residuals>5]) %>% 
    mutate(Age_class = factor(Age_class,levels = c('LOCRC','EOCRC')))
  
  # hyper_regression_FC <- mean(LO_regression[residuals>5]$residuals)/mean(EO_regression[residuals>5]$residuals)
  
  summary(LO_regression[residuals>5]$residuals)
  summary(EO_regression[residuals>5]$residuals)
  hist(LO_regression[residuals>5]$residuals,breaks = 200)
  hist(EO_regression[residuals>5]$residuals,breaks = 200)
 
  wilcox.test(LO_regression[residuals>5]$residuals,
              EO_regression[residuals>5]$residuals)
  
  cliff.delta( EO_regression[residuals>5]$residuals, 
               LO_regression[residuals>5]$residuals)
  
  picture_hyper_regression <- ggplot(hyper_regression_TMB, 
                                        aes(x = Age_class, 
                                            y = residuals, fill = Age_class)) +
    geom_boxplot(outlier.shape = T,
                 width =0.5) + 
    #geom_jitter(width = 0.03, size = 0.1) + 
    #stat_boxplot(geom = "errorbar", width = 0.2) +
    #scale_fill_npg() + 
    scale_fill_manual(values = palette)+
    #ylim(1.0,4.5)+
    scale_y_continuous(labels = scales::number_format(accuracy = 2),
                       limits = c(0,100)) +
    
    theme_pubr() +
    
    labs(title = "TMB_EOCRC_vs_LOCRC", 
         subtitle = "hypermutated", 
         x = "", 
         y = "adjusted TMB(residuals)"
         
    ) +
    geom_signif(comparisons =  list(c("EOCRC", "LOCRC")),
                y_position = 90,
                tip_length = 0.02,
                #map_signif_level = p, 
                test = wilcox.test,
                textsize = 5
    ) +
    # 添加 Fold Change 值到图形中
    # annotate("text", 
    #          x = 0.7,  # 放在 EOCRC 和 LOCRC 中间
    #          y = 100,  # 设置文本显示的位置
    #          label = paste("FC =", round(hyper_regression_FC, 2)),  # 显示 Fold Change，保留两位小数
    #          size = 5, 
    #          color = "black") +
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
  
  pdf(file = paste(dir,'/2.regression_hyepr_TMB.pdf',sep = ''),
      width = 4.5,height = 5.5)
  picture_hyper_regression
  dev.off()
  
  
}



{#TMB_mutload_worldwide
#此处TMB_mutload没有经过1%的首尾筛选
  
  TMB_mutload_binded = TMB_mutload %>% 
    dplyr::arrange(desc(TMB)) %>% 
    mutate(
      order = c(1:length(Tumor_Sample_Barcode))
    )
  
      # lower_TMB_mutload_binded = quantile(TMB_mutload_binded$TMB, 0.01)
      # upper_TMB_mutload_binded = quantile(TMB_mutload_binded$TMB, 0.99)
      # # TMB_EOCRC = TMB_EOCRC[TMB_EOCRC >= lower_TMB_mutload_binded & TMB_EOCRC <= upper_TMB_mutload_binded]
      # TMB_mutload_binded <- TMB_mutload_binded %>% 
      #   filter(TMB >= lower_TMB_mutload_binded) %>% 
      #   filter(TMB<= upper_TMB_mutload_binded )
  
  
  
  #quantiles <- quantile(TMB_mutload_binded$TMB, probs = c(0.05, 0.95))
  #TMB_mutload_binded <- TMB_mutload_binded %>% 
  #  filter(TMB >= quantiles[1] & TMB <= quantiles[2])
  
  palette = c('black')
  
  last_total_perMB_10 <- TMB_mutload_binded %>% 
    filter(TMB == 10) %>% 
    slice_head(n = 1)
  
  x_position <- last_total_perMB_10$order
  
  ggplot(TMB_mutload_binded,aes(x=order,y=logTMB))+
    geom_point() +
    geom_hline(yintercept = 1,
               linetype=2,
               lwd=1,
               color='red')+
    geom_vline(xintercept = x_position,
               linetype=2,
               color='red',
               lwd=1) +
    theme_few()+
    #ylim(0,1000)+
    labs(title = "Tumor Mutation Burden Distribution", 
         #subtitle = "nonhypermutated", 
         x = "",
         y = "log10(TMB)"
    ) + theme(
      plot.title    = element_text(color = "black", size   = 14, hjust = 0.5),
      #plot.subtitle = element_text(color = "black", size   = 16,hjust = 0.5),
      plot.caption  = element_text(color = "black", size   = 12,face = "italic", hjust = 1),
      axis.text.x   = element_blank(),
      axis.text.y   = element_text(color = "black", size = 12, angle = 0),
      #axis.title.x  = element_text(color = "black", size = 16, angle = 0),
      axis.title.y  = element_text(color = "black", size = 12, angle = 90),
      #legend.title  = element_text(color = "black", size  = 16),
      #legend.text   = element_text(color = "black", size   = 16),
      #strip.text = element_text(size = 14),
      axis.line.y = element_line(color = "black", linetype = "solid"), 
      axis.line.x = element_line (color = "black",linetype = "solid"), 
      panel.border = element_rect(linetype = "solid", linewidth = 1.2,fill = NA) )
  
  ggsave(filename = paste(dir,'/2_TMB_binded.pdf',sep = ''),
         width = 7,height = 4,dpi = 300)
  
}


{#regression_TMB_mutload_worldwide
  #此处TMB_mutload没有经过1%的首尾筛选
  
  regression_TMB_mutload_binded = TMB_mutload_regression %>% 
    dplyr::arrange(desc(residuals)) %>% 
    mutate(
      order = c(1:length(Tumor_Sample_Barcode))
    )
  
  
  # palette = c('black')
  
  closest_point <- regression_TMB_mutload_binded %>%
    mutate(abs_diff = abs(residuals - 5)) %>%   # 计算 residuals 与 5 的差值
    arrange(abs_diff) %>%                       # 按差值升序排列
    slice(1)                                    # 选择第一个点
  
  x_position <- closest_point$order  # 提取 order 值作为交点的 x 坐标
  
  ggplot(regression_TMB_mutload_binded, aes(x = order, y = residuals)) +
    geom_point() +
    geom_hline(yintercept = 5, linetype = 2, lwd = 1, color = 'red') +
    geom_vline(xintercept = x_position, linetype = 2, color = 'red', lwd = 1) +
    theme_few() +
    labs(
      title = "Adjusted Tumor Mutation Burden Distribution", 
      x = "",
      y = "adjusted TMB(residuals)"
    ) + theme(
      plot.title = element_text(color = "black", size = 14, hjust = 0.5),
      plot.caption = element_text(color = "black", size = 12, face = "italic", hjust = 1),
      axis.text.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 12, angle = 0),
      axis.title.y = element_text(color = "black", size = 12, angle = 90),
      axis.line.y = element_line(color = "black", linetype = "solid"), 
      axis.line.x = element_line(color = "black", linetype = "solid"), 
      panel.border = element_rect(linetype = "solid", linewidth = 1.2, fill = NA)
    )
  
  ggsave(filename = paste(dir, '/2_regression_TMB_binded.pdf', sep = ''), 
         width = 8, height = 4, dpi = 300)
  
  

  
}




{#TMB_mutload_each_country
  
  
  
  plot_tmb_distribution <- function(country_name, rt_clinical_data, TMB_mutload) {
    # Screening samples from designated countries
    country_tmb_sample <- rt_clinical_data$Tumor_Sample_Barcode[rt_clinical_data$Country == country_name]
    country_tmb <- TMB_mutload[TMB_mutload$Tumor_Sample_Barcode %in% country_tmb_sample,]
    country_tmb <- country_tmb %>% dplyr::arrange(desc(TMB))
    country_tmb$order <- 1:nrow(country_tmb)
    
    last_total_perMB_10 <- country_tmb %>% 
      filter(TMB == 10) %>% 
      slice_head(n = 1)
    
    if (nrow(last_total_perMB_10) == 0) {
      
      closest_above <- country_tmb %>% 
        filter(TMB > 10) %>% 
        slice_tail(n = 1)
      
      closest_below <- country_tmb %>% 
        filter(TMB < 10) %>% 
        slice_head(n = 1)
      
      if (nrow(closest_above) > 0 & nrow(closest_below) > 0) {
        x_position <- mean(c(closest_above$order, closest_below$order))
      } else if (nrow(closest_above) > 0) {
        x_position <- closest_above$order
      } else if (nrow(closest_below) > 0) {
        x_position <- closest_below$order
      } else {
        x_position <- NA 
      }
    } else {
      x_position <- last_total_perMB_10$order
    }
    
      # lower_country_tmb = quantile(country_tmb$TMB, 0.01)
      # upper_country_tmb = quantile(country_tmb$TMB, 0.99)
      # # TMB_EOCRC = TMB_EOCRC[TMB_EOCRC >= lower_country_tmb & TMB_EOCRC <= upper_country_tmb]
      # country_tmb <- country_tmb %>% 
      #   filter(TMB >= lower_country_tmb) %>% 
      #   filter(TMB<= upper_country_tmb )
    
    
    ggplot(country_tmb, aes(x = order, y = logTMB)) +
      geom_point() +
      geom_hline(yintercept = 1, linetype = 2, lwd = 1, color = 'red') +
      geom_vline(xintercept = x_position, linetype = 2, color = 'red', lwd = 1) +
      theme_few() +
      labs(title = paste("Tumor Mutation Burden Distribution in", country_name), 
           x = "", 
           y = "log10(TMB)"
      ) +
      theme(
        plot.title    = element_text(color = "black", size = 14, hjust = 0.5),
        plot.caption  = element_text(color = "black", size = 12, face = "italic", hjust = 1),
        axis.text.x   = element_blank(),
        axis.text.y   = element_text(color = "black", size = 12, angle = 0),
        axis.title.y  = element_text(color = "black", size = 12, angle = 90),
        axis.line.y   = element_line(color = "black", linetype = "solid"), 
        axis.line.x   = element_line(color = "black", linetype = "solid"), 
        panel.border  = element_rect(linetype = "solid", linewidth = 1.2, fill = NA)
      )
    
    ggsave(filename = paste(dir,'/2_',country_name,'_TMB_distribution.pdf',sep = ''),width = 7,height = 4,dpi = 300)
    
  }
  
  plot_tmb_distribution("Canada", rt@clinical.data, TMB_mutload)
  plot_tmb_distribution(country_name = "US", rt_clinical_data = rt@clinical.data, 
                        TMB_mutload = TMB_mutload)
  plot_tmb_distribution("China", rt@clinical.data, TMB_mutload)
  plot_tmb_distribution("Nigeria", rt@clinical.data, TMB_mutload)
  plot_tmb_distribution("Netherlands", rt@clinical.data, TMB_mutload)
  plot_tmb_distribution("Spain", rt@clinical.data, TMB_mutload)
  plot_tmb_distribution("Korea", rt@clinical.data, TMB_mutload)
  plot_tmb_distribution("France", rt@clinical.data, TMB_mutload)

  
  # plot_tmb_distribution_part <- function(country_name, rt_clinical_data, TMB_mutload) {
  #   
  #   country_tmb_sample <- rt_clinical_data$Tumor_Sample_Barcode[rt_clinical_data$Country == country_name]
  #   country_tmb <- TMB_mutload[TMB_mutload$Tumor_Sample_Barcode %in% country_tmb_sample,]
  #   country_tmb <- country_tmb %>% dplyr::arrange(desc(TMB))
  #   country_tmb$order <- 1:nrow(country_tmb)
  #   
  #   
  #   last_total_perMB_10 <- country_tmb %>% 
  #     filter(TMB == 10) %>% 
  #     slice_head(n = 1)
  #   
  #   x_position <- last_total_perMB_10$order
  #   
  #   
  #   ggplot(country_tmb, aes(x = order, y = logTMB)) +
  #     geom_point() +
  #     #geom_hline(yintercept = 1, linetype = 2, lwd = 1, color = 'red') +
  #     #geom_vline(xintercept = x_position, linetype = 2, color = 'red', lwd = 1) +
  #     theme_few() +
  #     labs(title = paste("Tumor Mutation Rate Distribution in", country_name), 
  #          x = "", 
  #          y = "log10(TMB)"
  #     ) +
  #     theme(
  #       plot.title    = element_text(color = "black", size = 14, hjust = 0.5),
  #       plot.caption  = element_text(color = "black", size = 12, face = "italic", hjust = 1),
  #       axis.text.x   = element_blank(),
  #       axis.text.y   = element_text(color = "black", size = 12, angle = 0),
  #       axis.title.y  = element_text(color = "black", size = 12, angle = 90),
  #       axis.line.y   = element_line(color = "black", linetype = "solid"), 
  #       axis.line.x   = element_line(color = "black", linetype = "solid"), 
  #       panel.border  = element_rect(linetype = "solid", linewidth = 1.2, fill = NA)
  #     )
  #   
  #   ggsave(filename = paste(dir,'/2_',country_name,'_TMB_distribution.pdf',sep = ''),width = 7,height = 4,dpi = 300)
  #   
  # }
  # 
  # plot_tmb_distribution_part("Korea", rt@clinical.data, TMB_mutload)
  # plot_tmb_distribution_part("France", rt@clinical.data, TMB_mutload)
  
  
}


{#Tumor Mutation Rate Distribution among Countries
  countries_TMB = data.frame(Tumor_Sample_Barcode='',total='',Number='',
                             Panel='',total_size='',TMB='',
                             logTMB='',order='',Country='',
                             Age='',Age_class = '',Center='',
                             Sex='',Race = '')
  for (kokka in unique(rt@clinical.data$Country)){
    #kokka = 'Korea'
    country_tmb_sample <- rt@clinical.data$Tumor_Sample_Barcode[rt@clinical.data$Country == kokka]
    country_tmb <- TMB_mutload[TMB_mutload$Tumor_Sample_Barcode %in% country_tmb_sample,]
    country_tmb <- country_tmb %>% dplyr::arrange(desc(TMB))
    country_tmb$order <- 1:nrow(country_tmb)
    
        # #去掉首尾离群值
        # lower_country_tmb = quantile(country_tmb$TMB, 0.01)
        # upper_country_tmb = quantile(country_tmb$TMB, 0.99)# TMB_EOCRC = TMB_EOCRC[TMB_EOCRC >= lower_country_tmb & TMB_EOCRC <= upper_country_tmb]
        # country_tmb <- country_tmb %>% 
        #   filter(TMB >= lower_country_tmb) %>% 
        #   filter(TMB<= upper_country_tmb )
        # 
    
    country_size = length(country_tmb_sample)
    country_tmb$order = country_tmb$order/country_size
    
    country_tmb$Country = kokka
    
    countries_TMB = rbind(countries_TMB,country_tmb)
  }
  
  countries_TMB = countries_TMB[-1,]
  countries_TMB$logTMB <- as.numeric(countries_TMB$logTMB)
  countries_TMB$order <- as.numeric(countries_TMB$order)
  
  summary(countries_TMB$logTMB)
  
  ggplot(countries_TMB, aes(x = order, y = logTMB)) +
    geom_point(size = 0.5) +  
    facet_wrap(~ Country, scales = "free_x", nrow = 1) +
    geom_hline(yintercept = 1,
               linetype = 2,
               lwd = 1,
               color = 'red') +
    ylim(-1.5, 3.5) +  
    theme_few() +
    labs(title = "Tumor Mutation Burden Distribution among Countries",
         x = "Country",
         y = "log10(TMB)") +
    theme(
      plot.title = element_text(color = "black", size = 12, hjust = 0.5),
      plot.caption = element_text(color = "black", size = 12, face = "italic", hjust = 1),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 12, angle = 0),
      axis.title.y = element_text(color = "black", size = 12, angle = 90),
      axis.line.y = element_line(color = "black", linetype = "solid"),
      axis.line.x = element_line(color = "black", linetype = "solid"),
      panel.border = element_rect(linetype = "solid", linewidth = 1.2, fill = NA),
      strip.background = element_blank(),
      panel.spacing = unit(0, "lines")
    )
  
  ggsave(filename = paste(dir,'/2_Tumor Mutation Rate Distribution among Countries.pdf',sep = ''),
         width = 12,height = 5,dpi = 300)
  
}

