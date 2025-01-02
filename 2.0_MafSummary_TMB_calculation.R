library(Cairo)
library(ggthemes)
library(ggpubr)
library(magrittr) 
library(ggsignif)
library(readxl)
library(tidyverse)
library(rstatix)


dir= './picture'

{    # MafSummary
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
  TMB_EOCRC_mutload <-  left_join(getSampleSummary(clin.EOCRC)[, .(Tumor_Sample_Barcode,total)],
                            rt@clinical.data[,.(Number,Tumor_Sample_Barcode,Panel)],
                            by='Tumor_Sample_Barcode') 
  
  TMB_EOCRC_mutload <-  left_join(TMB_EOCRC_mutload,capture_size,by = 'Panel')
  
  table(is.na(TMB_EOCRC_mutload$Panel))
  TMB_EOCRC_mutload$TMB <- (TMB_EOCRC_mutload$total) / TMB_EOCRC_mutload$total_size
  
  TMB_EOCRC_mutload$logTMB <-  log10(TMB_EOCRC_mutload$TMB)
  TMB_EOCRC_mutload = TMB_EOCRC_mutload[order(TMB, decreasing = FALSE)]
  
  TMB_EOCRC = TMB_EOCRC_mutload$logTMB
  TMB_EOCRC_hyper_sample = TMB_EOCRC_mutload$Tumor_Sample_Barcode[TMB_EOCRC_mutload$TMB > 15]
  TMB_EOCRC_nonhyper_sample = TMB_EOCRC_mutload$Tumor_Sample_Barcode[TMB_EOCRC_mutload$TMB <= 15]
 
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
  TMB_LOCRC_hyper_sample = TMB_LOCRC_mutload$Tumor_Sample_Barcode[TMB_LOCRC_mutload$TMB > 15]
  TMB_LOCRC_nonhyper_sample = TMB_LOCRC_mutload$Tumor_Sample_Barcode[TMB_LOCRC_mutload$TMB <= 15]
  
  TMB_mutload = rbind(TMB_LOCRC_mutload,TMB_EOCRC_mutload) %>% 
    left_join(rt@clinical.data[,.(Age,Age_class,Country,Center,
                                  Sex,Race,Tumor_site,
                                  Histology_type,Sample_type,Tumor_Sample_Barcode)]
              ,by='Tumor_Sample_Barcode')
 
  #MSS
  # TMB_mutload_model =  lm(logTMB~Age_class+Panel+Sex+Tumor_site+Histology_type,
  #                         data = TMB_mutload)
  
  TMB_mutload_model =  lm(logTMB~Age_class+Panel+Race+Sex+Tumor_site+Sample_type+Histology_type,
          data = TMB_mutload)
  
  TMB_mutload_regression = TMB_mutload %>%
    mutate(residuals = residuals(TMB_mutload_model))

  write.csv(TMB_mutload,file = './table/2.0_TMB_mutload.csv')
 
}



{
  hist((TMB_EOCRC), breaks = 200,plot=T)
  hist((TMB_LOCRC), breaks = 200,plot=T)
 
}


  palette <- c('#deebf7','#3182bd')

{#TMB_across_panels_for_EOandLO
 
  picture_acorss_panels <-  
    ggplot(TMB_mutload, aes(x = Age_class, y = logTMB, fill = Age_class)) +
    geom_boxplot(outlier.shape = T, width = 0.5) + 
    facet_wrap(~ Panel, scales = "free_y",nrow = 3) +  
    theme_few() +
    scale_fill_manual(values = palette) +
    theme(
      axis.text.x   = element_text(color = "black", size = 10, angle = 45,hjust = 1),
      axis.text.y   = element_text(color = "black", size = 10, angle = 0),
      axis.title.x  = element_text(color = "black", size = 14, angle = 0),
      axis.title.y  = element_text(color = "black", size = 14, angle = 90),
      legend.position = 'right',
      axis.line.y = element_line(color = "black", linetype = "solid"), 
      axis.line.x = element_line (color = "black",linetype = "solid"), 
      panel.border = element_rect(linetype = "solid", linewidth = 1.2,fill = NA),
      strip.text = element_text(size = 8, color = "black", face = "bold") 
    )
  
}

  pdf(file = './picture/2.TMB_across_panels_for_EOandLO.pdf',width = 22,height = 6)
  picture_acorss_panels
  dev.off()

{
  #nonhyper regression_TMB
  nonhyper_regression_TMB = TMB_mutload_regression %>% 
                            filter(TMB<=15) %>% 
                            mutate(Age_class = factor(Age_class,levels = c('LOCRC','EOCRC')))
  
  picture_nonhyper_regression <- ggplot(nonhyper_regression_TMB, 
                                        aes(x = Age_class, 
                                            y = residuals, fill = Age_class)) +
                                geom_boxplot(outlier.shape = T,
                                             width =0.5) +
                                scale_fill_manual(values = palette)+
                                scale_y_continuous(labels = scales::number_format(accuracy = 0.05),
                                                   limits = c(-2,1.5)) +
                                theme_pubr() +
                                labs(title = "TMB_EOCRC_vs_LOCRC", 
                                     subtitle = "Non-hypermutated", 
                                     x = "", 
                                     y = "adjusted TMB(residuals)") +
                                geom_signif(comparisons =  list(c("EOCRC", "LOCRC")),
                                            y_position =1.3,
                                            tip_length =0.02, 
                                            test =  wilcox.test,
                                            textsize = 5) +
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

  pdf(file = './picture/2.regression_nonhyepr_TMB.pdf',
      width = 4.5,height = 5.5)
  picture_nonhyper_regression
  dev.off()


{
  #hyper regression_TMB
  hyper_regression_TMB = TMB_mutload_regression %>% 
                         filter(TMB>15) %>% 
                         mutate(Age_class = factor(Age_class,levels = c('LOCRC','EOCRC')))

  picture_hyper_regression <- ggplot(hyper_regression_TMB, 
                                      aes(x = Age_class, 
                                          y = residuals, fill = Age_class)) +
                              geom_boxplot(outlier.shape = T,
                                         width =0.5) + 
                              scale_fill_manual(values = palette)+
                              scale_y_continuous(labels = scales::number_format(accuracy = 0.05),
                                                 limits = c(-0.1,2.8)) +
                              theme_pubr() +
                              labs(title = "TMB_EOCRC_vs_LOCRC", 
                                   subtitle = "Hypermutated", 
                                   x = "", 
                                   y = "adjusted TMB(residuals)") +
                              geom_signif(comparisons =  list(c("EOCRC", "LOCRC")),
                                          y_position = 2.6,
                                          tip_length = 0.02,
                                          test = wilcox.test,
                                          textsize = 5) +
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

  pdf(file = './picture/2.regression_hyepr_TMB.pdf',
      width = 4.5,height = 5.5)
  picture_hyper_regression
  dev.off()



  regression_EO_vs_LO_data <- data.frame(
    TMB_Status = c('Hyper','Non-hyper'),
    EO_number = c(nrow(hyper_regression_TMB[Age_class == 'EOCRC']),nrow(nonhyper_regression_TMB[Age_class == 'EOCRC'])),
    LO_number = c(nrow(hyper_regression_TMB[Age_class == 'LOCRC']),nrow(nonhyper_regression_TMB[Age_class == 'LOCRC'])),
    pval = c(wilcox.test(hyper_regression_TMB[Age_class == 'LOCRC']$residuals,
                         hyper_regression_TMB[Age_class == 'EOCRC']$residuals)$p.value,
             wilcox.test(nonhyper_regression_TMB[Age_class == 'LOCRC']$residuals,
                         nonhyper_regression_TMB[Age_class == 'EOCRC']$residuals)$p.value ),
    EO_mean = c(summary(hyper_regression_TMB[Age_class == 'EOCRC']$residuals)['Mean'],
                summary(nonhyper_regression_TMB[Age_class == 'EOCRC']$residuals)['Mean']),
    LO_mean = c(summary(hyper_regression_TMB[Age_class == 'LOCRC']$residuals)['Mean'],
                summary(nonhyper_regression_TMB[Age_class == 'LOCRC']$residuals)['Mean']),
    EO_median = c(summary(hyper_regression_TMB[Age_class == 'EOCRC']$residuals)['Median'],
                summary(nonhyper_regression_TMB[Age_class == 'EOCRC']$residuals)['Median']),
    LO_median = c(summary(hyper_regression_TMB[Age_class == 'LOCRC']$residuals)['Median'],
                summary(nonhyper_regression_TMB[Age_class == 'LOCRC']$residuals)['Median']),
    Effect_size = c(cliff.delta( hyper_regression_TMB[Age_class == 'EOCRC']$residuals,
                                 hyper_regression_TMB[Age_class == 'LOCRC']$residuals)$estimate,
                    cliff.delta(nonhyper_regression_TMB[Age_class == 'EOCRC']$residuals,
                                nonhyper_regression_TMB[Age_class == 'LOCRC']$residuals)$estimate),
    Magnitude =  c( as.character(cliff.delta( hyper_regression_TMB[Age_class == 'EOCRC']$residuals,
                                              hyper_regression_TMB[Age_class == 'LOCRC']$residuals)$magnitude),
                    as.character(cliff.delta( nonhyper_regression_TMB[Age_class == 'EOCRC']$residuals,
                                              nonhyper_regression_TMB[Age_class == 'LOCRC']$residuals)$magnitude)
                  )
  )

write_csv(regression_EO_vs_LO_data,file = './table/2.0_regression_EO_vs_LO_data.csv')


{#TMB_mutload_worldwide

  TMB_mutload_binded = TMB_mutload %>%
                        dplyr::arrange(desc(TMB)) %>% 
                        mutate( order = c(1:length(Tumor_Sample_Barcode)) )

  palette = c('black')
  
  last_total_perMB_15 <- TMB_mutload_binded %>% 
                          filter(TMB == 15) %>% 
                          slice_head(n = 1)
  
  if (nrow(last_total_perMB_15) == 0) {

      closest_above <- TMB_mutload_binded %>% 
        filter(TMB > 15) %>% 
        slice_tail(n = 1)
      
      closest_below <- TMB_mutload_binded %>% 
        filter(TMB < 15) %>% 
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
      x_position <- last_total_perMB_15$order
    }
    
    ggplot(TMB_mutload_binded,aes(x=order,y=logTMB))+
      geom_point() +
      geom_hline(yintercept = log10(15),
                 linetype=5,
                 lwd=1,
                 color='red')+
      geom_vline(xintercept = x_position,
                 linetype=5,
                 color='red',
                 lwd=1) +
      theme_few()+
      labs(title = "Tumor Mutation Burden Distribution",
           x = "",
           y = "log10(TMB)") + 
      theme(
        plot.title    = element_text(color = "black", size   = 14, hjust = 0.5),
        plot.caption  = element_text(color = "black", size   = 12,face = "italic", hjust = 1),
        axis.text.y   = element_text(color = "black", size = 12, angle = 0),
        axis.title.y  = element_text(color = "black", size = 12, angle = 90),
        axis.line.y = element_line(color = "black", linetype = "solid"), 
        axis.line.x = element_line (color = "black",linetype = "solid"), 
        panel.border = element_rect(linetype = "solid", linewidth = 1.2,fill = NA) )
  
  ggsave(filename ='./picture/2_TMB_binded.pdf',
         width = 7,height = 4.7,dpi = 300)
  
}

{#TMB_mutload_each_country
  
  plot_tmb_distribution <- function(country_name, rt_clinical_data, TMB_mutload) {
    country_tmb_sample <- rt_clinical_data$Tumor_Sample_Barcode[rt_clinical_data$Country == country_name]
    country_tmb <- TMB_mutload[TMB_mutload$Tumor_Sample_Barcode %in% country_tmb_sample,]
    country_tmb <- country_tmb %>% dplyr::arrange(desc(TMB))
    country_tmb$order <- 1:nrow(country_tmb)
    
    last_total_perMB_15 <- country_tmb %>% 
      filter(TMB == 15) %>% 
      slice_head(n = 1)
    
    if (nrow(last_total_perMB_15) == 0) {
      
      closest_above <- country_tmb %>% 
        filter(TMB >15) %>% 
        slice_tail(n = 1)
      
      closest_below <- country_tmb %>% 
        filter(TMB < 15) %>% 
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
      x_position <- last_total_perMB_15$order
    }
    
   
    ggplot(country_tmb, aes(x = order, y = logTMB)) +
      geom_point() +
      geom_hline(yintercept = log10(15), linetype = 5, lwd = 1, color = 'red') +
      geom_vline(xintercept = x_position, linetype = 5, color = 'red', lwd = 1) +
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
    
    ggsave(filename = paste(dir,'/2_',country_name,'_TMB_distribution.pdf',sep = ''),width = 6,height = 4,dpi = 300)
    
  }
  
  plot_tmb_distribution("Canada", rt@clinical.data, TMB_mutload)
  plot_tmb_distribution("US", rt@clinical.data, TMB_mutload)
  plot_tmb_distribution("China", rt@clinical.data, TMB_mutload)
  plot_tmb_distribution("Nigeria", rt@clinical.data, TMB_mutload)
  plot_tmb_distribution("Netherlands", rt@clinical.data, TMB_mutload)
  plot_tmb_distribution("Spain", rt@clinical.data, TMB_mutload)
  plot_tmb_distribution("Korea", rt@clinical.data, TMB_mutload)
  plot_tmb_distribution("France", rt@clinical.data, TMB_mutload)
 
}





