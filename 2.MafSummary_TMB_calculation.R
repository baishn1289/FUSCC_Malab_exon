library(Cairo)
library(ggplot2)
library(ggthemes)
library(ggpubr)

dir= './picture'

{# MafSummary
  se_gene = getGeneSummary(clin.EOCRC)[1:50]$Hugo_Symbol
 
  
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
  
  #EOCRC_TMB_calculation
  maf.mutload <-  left_join(getSampleSummary(clin.EOCRC)[, .(Tumor_Sample_Barcode,total)],
                            rt@clinical.data[,.(Number,Tumor_Sample_Barcode,Panel)],
                            by='Tumor_Sample_Barcode')
  
  maf.mutload <-  left_join(maf.mutload,capture_size,by = 'Panel')
  
  table(is.na(maf.mutload$Panel))
  maf.mutload$total_perMB <- (maf.mutload$total) / maf.mutload$total_size
  
  maf.mutload$total_perMB_log <-  log10(maf.mutload$total_perMB)
  maf.mutload = maf.mutload[order(total_perMB, decreasing = FALSE)]
  
  TMB_EOCRC = maf.mutload$total_perMB_log
  TMB_EOCRC_hyper_sample = maf.mutload$Tumor_Sample_Barcode[maf.mutload$total_perMB_log > 1]
  TMB_EOCRC_nonhyper_sample = maf.mutload$Tumor_Sample_Barcode[maf.mutload$total_perMB_log <= 1]
  
  TMB_EOCRC_mutload = maf.mutload
  write.csv(TMB_EOCRC_mutload,file = 'TMB_EOCRC.csv')
  
  #LOCRC_TMB_calculation
  maf.mutload <-  left_join(getSampleSummary(clin.LOCRC)[, .(Tumor_Sample_Barcode,total)],
                            rt@clinical.data[,.(Number,Tumor_Sample_Barcode,Panel)],
                            by='Tumor_Sample_Barcode')
  
  maf.mutload <-  left_join(maf.mutload,capture_size,by = 'Panel')
  
  table(is.na(maf.mutload$Panel))
  maf.mutload$total_perMB <- (maf.mutload$total) / maf.mutload$total_size
  
  maf.mutload$total_perMB_log <-  log10(maf.mutload$total_perMB)
  maf.mutload = maf.mutload[order(total_perMB, decreasing = FALSE)]
  
  
  TMB_LOCRC = maf.mutload$total_perMB_log
 
  TMB_LOCRC_hyper_sample = maf.mutload$Tumor_Sample_Barcode[maf.mutload$total_perMB_log > 1]
  TMB_LOCRC_nonhyper_sample = maf.mutload$Tumor_Sample_Barcode[maf.mutload$total_perMB_log <= 1]
  
  TMB_LOCRC_mutload = maf.mutload
  write.csv(TMB_LOCRC_mutload,file = 'TMB_LOCRC.csv')
  
  
  TMB_mutload = rbind(TMB_LOCRC_mutload,TMB_EOCRC_mutload) %>% 
    left_join(rt@clinical.data[,.(Age,Age_class,Country,Center,Sex,Race,Tumor_Sample_Barcode)]
              ,by='Tumor_Sample_Barcode')
  write.csv(TMB_mutload,file = 'TMB_mutload.csv')
  
  
}


{
  hist((TMB_EOCRC), breaks = 200,plot=T)
  hist((TMB_LOCRC), breaks = 200,plot=T)
  
  
  lower_bound_EOCRC = quantile(TMB_EOCRC, 0.01)
  upper_bound_EOCRC = quantile(TMB_EOCRC, 0.99)
  TMB_EOCRC = TMB_EOCRC[TMB_EOCRC >= lower_bound_EOCRC & TMB_EOCRC <= upper_bound_EOCRC]
  
  lower_bound_LOCRC = quantile(TMB_LOCRC, 0.01)
  upper_bound_LOCRC = quantile(TMB_LOCRC, 0.99)
  TMB_LOCRC = TMB_LOCRC[TMB_LOCRC >= lower_bound_LOCRC & TMB_LOCRC <= upper_bound_LOCRC]
 
  
  summary(TMB_EOCRC)
  summary(TMB_LOCRC)
  

  
}

library(magrittr) 
library(ggsignif)
library(tidyverse)
library(readxl)
library(ggpubr)
library(rstatix)

{#nonhyper TMB
  table(TMB_EOCRC <= 1)
  table(TMB_LOCRC <= 1)
  x = TMB_EOCRC[TMB_EOCRC <= 1]
  y = TMB_LOCRC[TMB_LOCRC <= 1]
  
  t.test(x, y)
  wilcox.test(x, y)
 
  
  nonhyper_TMB_EOCRC = data.frame(TMB = x,class='EOCRC')
  nonhyper_TMB_LOCRC = data.frame(TMB = y,class='LOCRC')
  nonhyper_TMB = rbind(nonhyper_TMB_EOCRC,nonhyper_TMB_LOCRC) 
  nonhyper_TMB$class <-as.factor(nonhyper_TMB$class)
  
  table(nonhyper_TMB$class)
  table(nonhyper_TMB$class)/nrow(nonhyper_TMB)
  
  
  palette <- c('#deebf7','#3182bd')
  
  picture_nonhyper <-  
    ggplot(nonhyper_TMB, aes(x = class, y = TMB, fill = class)) +
    geom_boxplot(outlier.shape = T, width = 0.5) + 
   
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1),
                       limits = c(-1,1.3)) +
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
     
      axis.line.y = element_line(color = "black", linetype = "solid"), 
      axis.line.x = element_line (color = "black",linetype = "solid"), 
      panel.border = element_rect(linetype = "solid", linewidth = 1.2,fill = NA),
      
      
    )
  
}

pdf(paste(dir,'/2_TMB_nonhyper_boxplot.pdf',sep = ""),width=4,height=5)
picture_nonhyper
dev.off()


{#*** Hyper TMB
  table(TMB_EOCRC > 1)
  table(TMB_LOCRC > 1)
  x = TMB_EOCRC[TMB_EOCRC > 1]
  y = TMB_LOCRC[TMB_LOCRC > 1]
  
  
  t.test(x, y)
  wilcox.test(x, y)
 
  
  hyper_TMB_EOCRC = data.frame(TMB = x,class='EOCRC')
  hyper_TMB_LOCRC = data.frame(TMB = y,class='LOCRC')
  hyper_TMB = rbind(hyper_TMB_EOCRC,hyper_TMB_LOCRC) 
  hyper_TMB$class <-as.factor(hyper_TMB$class)
  
  table(hyper_TMB$class)
  table(hyper_TMB$class)/nrow(hyper_TMB)
  
  
  picture_hyper <- ggplot(hyper_TMB, aes(x = class, y = TMB, fill = class)) +
    geom_boxplot(outlier.shape = T, width =0.5) + 
    
    scale_fill_manual(values = palette)+

    scale_y_continuous(labels = scales::number_format(accuracy = 0.1),
                       limits = c(1.0,2.2)) +
    
    theme_pubr() +
    
    labs(title = "TMB_EOCRC_vs_LOCRC", 
         subtitle = "hypermutated", 
         x = "", 
         y = "log10(TMB)"
         
    ) +
    geom_signif(comparisons =  list(c("EOCRC", "LOCRC")),
                y_position = 2.1,
                tip_length = 0.02,
                
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
      axis.line.y = element_line(color = "black", linetype = "solid"), 
      axis.line.x = element_line (color = "black",linetype = "solid"), 
      panel.border = element_rect(linetype = "solid", linewidth = 1.2,fill = NA)
      
    )
}

pdf(paste(dir,'/2_TMB_hyper_boxplot.pdf', sep=''),width=4, height=5)
picture_hyper
dev.off()



{#*** overall TMB
  
  x = TMB_EOCRC
  y = TMB_LOCRC
 
  t.test(x, y)
  wilcox.test(x, y)
 
  
  overall_TMB_EOCRC = data.frame(TMB = x,class='EOCRC')
  overall_TMB_LOCRC = data.frame(TMB = y,class='LOCRC')
  overall_TMB = rbind(overall_TMB_EOCRC,overall_TMB_LOCRC) 
  overall_TMB$class <-as.factor(overall_TMB$class)
  
  table(overall_TMB$class)
  table(overall_TMB$class)/nrow(overall_TMB)
  
  picture <- ggplot(overall_TMB, aes(x = class, y = TMB, fill = class)) +
    geom_boxplot(outlier.shape = T, width =0.5) + 
    
    scale_fill_manual(values = palette)+
   
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1),
                       limits = c(-1,2.5)) +
    
    theme_pubr() +
    
    labs(title = "TMB_EOCRC_vs_LOCRC", 
         subtitle = "overall", 
         x = "", 
         y = "log10(TMB)"
         
    ) +
    geom_signif(comparisons =  list(c("EOCRC", "LOCRC")),
                y_position = 2.2,
                tip_length = 0.02,
               
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
      axis.line.y = element_line(color = "black", linetype = "solid"), 
      axis.line.x = element_line (color = "black",linetype = "solid"), 
      panel.border = element_rect(linetype = "solid", linewidth = 1.2,fill = NA)
      
    )
}

pdf(paste(dir,'/2_TMB_overall_boxplot.pdf', sep = ''),width=4, height=5)
picture
dev.off() 




{#TMB_mutload_binded
  
  TMB_mutload_binded = TMB_mutload %>% dplyr::arrange(desc(total_perMB))
  TMB_mutload_binded$order = c(1:length(TMB_mutload_binded$Tumor_Sample_Barcode))
  
  
  lower_TMB_mutload_binded = quantile(TMB_mutload_binded$total_perMB, 0.01)
  upper_TMB_mutload_binded = quantile(TMB_mutload_binded$total_perMB, 0.99)

  TMB_mutload_binded <- TMB_mutload_binded %>% 
    filter(total_perMB >= lower_TMB_mutload_binded) %>% 
    filter(total_perMB<= upper_TMB_mutload_binded )
  
  
  palette = c('black')
  
  last_total_perMB_10 <- TMB_mutload_binded %>% 
    filter(total_perMB == 10) %>% 
    slice_head(n = 1)
  
  x_position <- last_total_perMB_10$order
  
  ggplot(TMB_mutload_binded,aes(x=order,y=total_perMB_log))+
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
 
    labs(title = "Tumor Mutation Burden Distribution", 
        
         x = "",
         y = "log10(TMB)"
    ) + theme(
      plot.title    = element_text(color = "black", size   = 14, hjust = 0.5),
    
      plot.caption  = element_text(color = "black", size   = 12,face = "italic", hjust = 1),
      axis.text.x   = element_blank(),
      axis.text.y   = element_text(color = "black", size = 12, angle = 0),
    
      axis.title.y  = element_text(color = "black", size = 12, angle = 90),
      
      axis.line.y = element_line(color = "black", linetype = "solid"), 
      axis.line.x = element_line (color = "black",linetype = "solid"), 
      panel.border = element_rect(linetype = "solid", linewidth = 1.2,fill = NA) )
  
  ggsave(filename = paste(dir,'/2_TMB_binded.pdf',sep = ''),
         width = 7,height = 4,dpi = 300)
  
}


{#TMB_mutload_each_country
  
  plot_tmb_distribution <- function(country_name, rt_clinical_data, TMB_mutload) {
  
    country_tmb_sample <- rt_clinical_data$Tumor_Sample_Barcode[rt_clinical_data$Country == country_name]
    country_tmb <- TMB_mutload[TMB_mutload$Tumor_Sample_Barcode %in% country_tmb_sample,]
    country_tmb <- country_tmb %>% dplyr::arrange(desc(total_perMB))
    country_tmb$order <- 1:nrow(country_tmb)
    
    last_total_perMB_10 <- country_tmb %>% 
      filter(total_perMB == 10) %>% 
      slice_head(n = 1)
    
    if (nrow(last_total_perMB_10) == 0) {
      
      closest_above <- country_tmb %>% 
        filter(total_perMB > 10) %>% 
        slice_tail(n = 1)
      
      closest_below <- country_tmb %>% 
        filter(total_perMB < 10) %>% 
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
    
    lower_country_tmb = quantile(country_tmb$total_perMB, 0.01)
    upper_country_tmb = quantile(country_tmb$total_perMB, 0.99)
  
    country_tmb <- country_tmb %>% 
      filter(total_perMB >= lower_country_tmb) %>% 
      filter(total_perMB<= upper_country_tmb )
    
    
    ggplot(country_tmb, aes(x = order, y = total_perMB_log)) +
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
  
  
}


{#Tumor Mutation Rate Distribution among Countries
  countries_TMB = data.frame(Tumor_Sample_Barcode='',total='',Number='',
                             Panel='',total_size='',total_perMB='',
                             total_perMB_log='',order='',Country='',
                             Age='',Age_class = '',Center='',
                             Sex='',Race = '')
  for (Guojia in unique(rt@clinical.data$Country)){
   
    country_tmb_sample <- rt@clinical.data$Tumor_Sample_Barcode[rt@clinical.data$Country == Guojia]
    country_tmb <- TMB_mutload[TMB_mutload$Tumor_Sample_Barcode %in% country_tmb_sample,]
    country_tmb <- country_tmb %>% dplyr::arrange(desc(total_perMB))
    country_tmb$order <- 1:nrow(country_tmb)
    
    lower_country_tmb = quantile(country_tmb$total_perMB, 0.01)
    upper_country_tmb = quantile(country_tmb$total_perMB, 0.99)
    country_tmb <- country_tmb %>% 
      filter(total_perMB >= lower_country_tmb) %>% 
      filter(total_perMB<= upper_country_tmb )
    
    
    country_size = length(country_tmb_sample)
    country_tmb$order = country_tmb$order/country_size
    
    country_tmb$Country = Guojia
    
    countries_TMB = rbind(countries_TMB,country_tmb)
  }
  
  countries_TMB = countries_TMB[-1,]
  countries_TMB$total_perMB_log <- as.numeric(countries_TMB$total_perMB_log)
  countries_TMB$order <- as.numeric(countries_TMB$order)
  
  ggplot(countries_TMB, aes(x = order, y = total_perMB_log)) +
    geom_point(size = 0.5) +  
    facet_wrap(~ Country, scales = "free_x", nrow = 1) +
    geom_hline(yintercept = 1,
               linetype = 2,
               lwd = 1,
               color = 'red') +
    ylim(-1.5, 2.5) +  
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

