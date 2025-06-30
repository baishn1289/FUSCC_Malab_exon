library(maftools)
library(tidyverse)
library(data.table)
library(ggthemes)
library(ggpubr)
library(ggsignif)
library(readxl)
library(car)
library(e1071)
library(lmtest)
library(mclust)
library(NMF)
library(mice)
library(Cairo)
library(patchwork)
library(parallel)
library(foreach)
library(doParallel)
library(ggsci)
library(treeio)
library(ape)
library(ggnewscale)

# set your path
# dir= ""
# dir_data = ''

# Load the MAF files
# read.maf()

# Multiple imputation
clin_tab <- rt@clinical.data
str(clin_tab)
colnames(clin_tab)

var <- names(clin_tab)
tab <- list()
total <- nrow(clin_tab)

for (i in 1:length(var)){
  tab[i]=data.frame(Variables=var[i],
                    Total=total,
                    Freq=total-sum(is.na(clin_tab[,get(var[i])])),
                    Missing=sum(is.na(clin_tab[,get(var[i])])),
                    miss_p=sprintf("%0.2f",sum(is.na(clin_tab[,get(var[i])]))/total*100)) %>% list()}
na<-do.call(rbind, tab)
view(na)

write.csv(na,file=paste0(dir,'table/1_na_data_percentage.csv'))

md.pattern(clin_tab)
md.pattern(clin_tab %>% select(Tumor_site,Race,Sample_type,Histology_type))


na_patient_number <- clin_tab %>% 
  filter(is.na(Tumor_site)|is.na(Histology_type)|is.na(Sample_type))
na_patient_percentage <- nrow(na_patient_number)/(nrow(clin_tab))

write.csv(data.frame(na_patient_number = nrow(na_patient_number),
                     na_patient_percentage = na_patient_percentage),
          file =paste0(dir,'table/1_na_patient_percentage.csv'))

mice_plot <- aggr(clin_tab, col=c('navyblue','yellow'),
                  numbers=TRUE, sortVars=TRUE,
                  labels=names(clin_tab), cex.axis=.7,
                  gap=5, ylab=c("Missing data","Pattern"),)

v <-c( "Sex","Tumor_site","Histology_type","Age_class","Country",
       "Panel","Sample_type" )
clin_tab =  clin_tab[, (v) := lapply(.SD, factor), .SDcols = v]

str(clin_tab)

imputed <- mice(clin_tab,maxit = 0)

names(clin_tab)
method_manual <- c(
  'Tumor_Sample_Barcode' = '',
  'Patient_ID' = '',
  'Sex' = '',
  'Age' = '',
  'MSI_Status'  = '',
  "Tumor_site" = "logreg",
  "Histology_type" = "logreg",
  'Stage' = '',
  "Age_class" = "",
  "Country" = "",
  'Number' = '',
  'Center' = '',
  "Panel" = "",
  "Sample_type" = "logreg",
  "Race" = ""
)

imputed<-mice(clin_tab, m=10, maxit = 20, seed =123,method = method_manual,
              parallel = "multicore",cores = 30)

summary(imputed)

imp_data <- complete(imputed, "all")
imp_list <- lapply(1:10, function(i) complete(imputed, i)) %>% 
  lapply(function(df){
    df %>% mutate(Age_class = factor(Age_class,levels = c('LOCRC','EOCRC')))
  })


format_Pval <- function(p) {
  ifelse(p < 0.0001,
         "<0.0001",
         case_when(
           p<0.001 ~ sprintf("%.5f", round(p*100000)/100000),
           p<0.01 ~sprintf("%.4f", round(p*10000)/10000),
           p<0.1 ~ sprintf("%.3f", round(p*1000)/1000),
           TRUE ~ sprintf("%.2f", round(p*100)/100)
         )
  )
}


#Waterfall plot
{
  CairoPDF(paste0(dir,'picture/2_oncoplot_Cairo.pdf'), width = 12, height = 6)
  
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


#TMB_calculation
{   
  
  #EOCRC_TMB_calculation
  TMB_EOCRC_mutload <-  left_join(getSampleSummary(clin.EOCRC)[, .(Tumor_Sample_Barcode,total)],
                                  rt@clinical.data[,.(Number,Tumor_Sample_Barcode,Panel)],
                                  by='Tumor_Sample_Barcode') 
  
  TMB_EOCRC_mutload <-  left_join(TMB_EOCRC_mutload,capture_size,by = 'Panel')
  
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
  
  TMB_LOCRC_mutload$TMB <- (TMB_LOCRC_mutload$total) / TMB_LOCRC_mutload$total_size
  
  TMB_LOCRC_mutload$logTMB <-  log10(TMB_LOCRC_mutload$TMB)
  TMB_LOCRC_mutload = TMB_LOCRC_mutload[order(TMB, decreasing = FALSE)]
  
  TMB_LOCRC = TMB_LOCRC_mutload$logTMB
  TMB_LOCRC_hyper_sample = TMB_LOCRC_mutload$Tumor_Sample_Barcode[TMB_LOCRC_mutload$TMB > 15]
  TMB_LOCRC_nonhyper_sample = TMB_LOCRC_mutload$Tumor_Sample_Barcode[TMB_LOCRC_mutload$TMB <= 15]
  
  TMB_mutload <-  rbind(TMB_LOCRC_mutload,TMB_EOCRC_mutload) %>% 
    dplyr::select(-Number,-Panel)  
  
  
  TMB_mutload_for_regression  <-  TMB_mutload %>% 
    filter(TMB < quantile(TMB,0.99) & TMB > quantile(TMB,0.01))
  
  write.csv(TMB_mutload,file = paste0(dir,'table/2.0_TMB_mutload.csv'))
}


#TMB distribution worldwide
palette <- c('#deebf7','#3182bd')

{
  TMB_mutload_binded = TMB_mutload %>%
    dplyr::arrange(desc(TMB)) %>% 
    mutate(
      order = c(1:length(Tumor_Sample_Barcode))
    )
  
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
  
  ggsave(filename =paste0(dir,'picture/2_TMB_binded.pdf'),
         width = 7,height = 4.7,dpi = 300)
  
}


#TMB_mutload_each_country
{
  plot_tmb_distribution <- function(country_name, rt_clinical_data, TMB_mutload) {
    # Screening samples from designated countries
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
    
    ggsave(filename = paste0(dir,'picture/2_',country_name,
                             '_TMB_distribution.pdf'),
           width = 6,height = 4,dpi = 300)
    
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


#TMB comparison between EOCRC and LOCRC
{
  
  imp_list_with_tmb <- lapply(imp_list, function(df) {
    left_join(df, TMB_mutload_for_regression, by = "Tumor_Sample_Barcode") %>% 
      filter(!is.na(TMB)) 
  })
  
  imp_list_with_tmb_hyper <- lapply(imp_list, function(df) {
    left_join(df, TMB_mutload_for_regression,by = "Tumor_Sample_Barcode") %>% 
      filter(TMB>15) %>% 
      filter(!is.na(TMB))
  })
  
  imp_list_with_tmb_nonhyper <- lapply(imp_list, function(df) {
    left_join(df, TMB_mutload_for_regression, by = "Tumor_Sample_Barcode") %>% 
      filter(TMB<=15) %>% 
      filter(!is.na(TMB)) %>% 
      mutate(adjTMB =  (1.31-logTMB) )   
  })
  
  
  fit_models_hyper <- lapply(imp_list_with_tmb_hyper, function(df) {
    model <- glm(formula = TMB ~ Age_class + Sex + Panel + Tumor_site + Histology_type,
                 family = Gamma(link = 'log'), data = df) 
    return(model)
  })
  
  class(fit_models_hyper) <- c("mira", "list")
  pooled_models_hyper <- pool(fit_models_hyper)
  summary_pooled_hyper <- summary(pooled_models_hyper, conf.int = TRUE)
  age_effect_hyper <- summary_pooled_hyper %>%
    filter(term == "Age_classEOCRC") %>%
    mutate(MR =  round(exp(estimate),digits = 2),
           `95% upper limit` = round(exp(conf.low),digits = 2),
           `95% lower limit` = round(exp(conf.high),digits = 2),
           p.value = format_Pval(as.numeric(p.value))
    ) %>%
    dplyr::select(term, MR, `95% upper limit`, 
                  `95% lower limit`,estimate, p.value)
  print(age_effect_hyper)

  
  fit_models_nohyper <- lapply(imp_list_with_tmb_nonhyper, function(df) {
    model <- glm(formula = adjTMB ~ Age_class + Sex + Panel + Tumor_site +Histology_type,
                 family = Gamma(link = 'log'), data = df)
    
    return(model)
  })
  
  class(fit_models_nohyper) <- c("mira", "list")
  pooled_models_nohyper <- pool(fit_models_nohyper)
  summary_pooled_nohyper <- summary(pooled_models_nohyper, conf.int = TRUE, exponentiate = TRUE)
  age_effect_nohyper <- summary_pooled_nohyper %>%
    filter(term == "Age_classEOCRC") %>%
    mutate(MR =  round(exp(estimate),digits = 2),
           `95% upper limit` = round(exp(conf.low),digits = 2),
           `95% lower limit` = round(exp(conf.high),digits = 2),
           p.value = format_Pval(as.numeric(p.value))
    ) %>%
    dplyr::select(term, MR, `95% upper limit`, 
                  `95% lower limit`,estimate, p.value)
  print(age_effect_nohyper)
  
  
  regression_data <- rbind(age_effect_hyper,age_effect_nohyper) %>% 
    mutate(Status = c('hyper','nonhyper'),
           Dependant = c(as.character(fit_models_hyper[[1]]$formula)[[2]],
                         as.character(fit_models_nohyper[[1]]$formula)[[2]]),
           Formula = c(as.character(fit_models_hyper[[1]]$formula)[[3]],
                       as.character(fit_models_nohyper[[1]]$formula)[[3]])) 
  write.csv(regression_data,file = paste0(dir,'table/2.0_regression_data.csv'))
  
}

{ 
  hyper_TMB = TMB_mutload_for_regression %>% 
    left_join(rt@clinical.data[,.(Tumor_Sample_Barcode,Age_class)],
              by='Tumor_Sample_Barcode') %>% 
    filter(TMB>15) %>% 
    mutate(Age_class = factor(Age_class,levels = c('LOCRC','EOCRC')))
  
  summary(hyper_TMB$logTMB)
  
  picture_hyper_regression <- ggplot(hyper_TMB, 
                                     aes(x = Age_class, 
                                         y = logTMB, fill = Age_class)) +
    geom_boxplot(outlier.shape = T,
                 width =0.5) + 
    scale_fill_manual(values = palette)+
    scale_y_continuous(labels = scales::number_format(accuracy = 0.05),
                       limits = c(1.16,2.1)) +
    theme_pubr() +
    labs(title = "Hypermutated", 
         x = "", 
         y = "logTMB"
         
    ) +
    annotate("text",
             x = 1.5,  
             y = 2.08,  
             label = paste("Mean ratio = ", round(regression_data[regression_data$Status =='hyper',]$MR,
                                                  2)),  
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
      axis.line.y = element_line(color = "black", linetype = "solid"), 
      axis.line.x = element_line (color = "black",linetype = "solid"), 
      panel.border = element_rect(linetype = "solid", linewidth = 1.2,fill = NA)
      
    )
}
pdf(file = paste0(dir,'picture/2.boxplot_hyepr_TMB.pdf'),
    width = 4.5,height = 5.5)
picture_hyper_regression
dev.off()


{
  nonhyper_TMB = TMB_mutload_for_regression %>% 
    left_join(rt@clinical.data[,.(Tumor_Sample_Barcode,Age_class)],
              by='Tumor_Sample_Barcode') %>% 
    filter(TMB<=15) %>% 
    mutate(Age_class = factor(Age_class,levels = c('LOCRC','EOCRC')),
           adjTMB = 1.31 - logTMB)
  
  summary(nonhyper_TMB$adjTMB)
  
  picture_nonhyper_regression <- ggplot(nonhyper_TMB, 
                                        aes(x = Age_class, 
                                            y = adjTMB, fill = Age_class)) +
    geom_boxplot(outlier.shape = T,
                 width =0.5) + 
    scale_fill_manual(values = palette)+
    scale_y_continuous(labels = scales::number_format(accuracy = 0.05),
                       limits = c(0.1,2.1)) +
    
    theme_pubr() +
    labs(title = "Non-hypermutated", 
         x = "", 
         y = "adjTMB"
         
    ) +
    annotate("text",
             x = 1.5,  
             y = 2.05, 
             label = paste("Mean ratio = ", round(regression_data[regression_data$Status =='nonhyper',]$MR,
                                                  2)), 
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
      axis.line.y = element_line(color = "black", linetype = "solid"), 
      axis.line.x = element_line (color = "black",linetype = "solid"), 
      panel.border = element_rect(linetype = "solid", linewidth = 1.2,fill = NA)
      
    )
  
}
pdf(file = paste0(dir,'picture/2.boxplot_nonhyepr_TMB.pdf'),
    width = 4.5,height = 5.5)
picture_nonhyper_regression
dev.off()


