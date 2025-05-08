
races_names <- c("Asian Pacific lslander", "White", "Black")

Race_regression <- function(race_input){
  
  race_samples <- rt@clinical.data[Race == race_input]$Tumor_Sample_Barcode
  
  imp_list_with_tmb_hyper_race <- lapply(imp_list, function(df) {
    left_join(df, TMB_mutload_for_regression,by = "Tumor_Sample_Barcode") %>% 
      filter(TMB>15) %>% 
      filter(Tumor_Sample_Barcode %in% race_samples,!is.na(TMB))
  })
  
  imp_list_with_tmb_nonhyper_race <- lapply(imp_list, function(df) {
    left_join(df, TMB_mutload_for_regression,by = "Tumor_Sample_Barcode") %>% 
      filter(TMB<=15) %>% 
      filter(Tumor_Sample_Barcode %in% race_samples,!is.na(TMB)) %>% 
      mutate(adjTMB =  (1.31-logTMB) )    
  })
  
  fit_models_hyper_race <- lapply(imp_list_with_tmb_hyper_race, function(df) {
    model <- glm(formula = TMB ~ Age_class + Sex + Panel + Tumor_site + Histology_type,
                 family = Gamma(link = 'log'), data = df)   
    return(model)
  })
  
  class(fit_models_hyper_race) <- c("mira", "list")
  pooled_models_hyper <- pool(fit_models_hyper_race)
  summary_pooled_hyper <- summary(pooled_models_hyper, conf.int = TRUE)
  age_effect_hyper <- summary_pooled_hyper %>%
    filter(term == "Age_classEOCRC") %>%
    mutate(MR =  round(exp(estimate),digits = 2),
           `95% upper limit` = round(exp(conf.low),digits = 2),
           `95% lower limit` = round(exp(conf.high),digits = 2),
           p.value = format_Pval(as.numeric(p.value))
    ) %>%
    dplyr::select(term, MR, `95% upper limit`, `95% lower limit`,estimate, p.value)
  
  fit_models_nonhyper_race <- lapply(imp_list_with_tmb_nonhyper_race, function(df) {
    model <- glm(formula = adjTMB ~ Age_class + Sex + Panel + Tumor_site + Histology_type,
                 family = Gamma(link = 'log'), data = df)   
    return(model)
  })
  
  
  skewness(imp_list_with_tmb_nonhyper_race[[1]]$adjTMB)
  
  class(fit_models_nonhyper_race) <- c("mira", "list")
  pooled_models_nonhyper <- pool(fit_models_nonhyper_race)
  summary_pooled_nonhyper <- summary(pooled_models_nonhyper, conf.int = TRUE)
  age_effect_nonhyper <- summary_pooled_nonhyper %>%
    filter(term == "Age_classEOCRC") %>%
    mutate(MR =  round(exp(estimate),digits = 2),
           `95% upper limit` = round(exp(conf.low),digits = 2),
           `95% lower limit` = round(exp(conf.high),digits = 2),
           p.value = format_Pval(as.numeric(p.value))
    ) %>%
    dplyr::select(term, MR, `95% upper limit`, `95% lower limit`,estimate, p.value)
  
  race_TMB_info <- rbind(age_effect_hyper,age_effect_nonhyper) %>% 
    mutate(Status = c('hyper','nonhyper'),
           Dependant = c(as.character(fit_models_hyper_race[[1]]$formula)[[2]],
                         as.character(fit_models_nonhyper_race[[1]]$formula)[[2]]),
           Race = race_input
    )
  
  return(race_TMB_info)
  
}

race_TMB <- lapply(c('White','Black',"Asian Pacific lslander"),Race_regression)

race_TMB_result <- rbindlist(race_TMB) 

write.csv(race_TMB_result,
          file = paste0(dir,'table/2.1.1_race_TMB_comprare.csv'))

palette <- c('#deebf7','#3182bd')

{
  
  hyper_race_data = TMB_mutload_for_regression[TMB>15] %>% 
    left_join(rt@clinical.data[,.(Tumor_Sample_Barcode,Race,Age_class)],
              by = 'Tumor_Sample_Barcode') %>% 
    filter(Race %in% races_names) %>% 
    mutate(Age_class = factor(Age_class,levels = c('LOCRC','EOCRC')))
  
  summary(hyper_race_data$logTMB)
  
  race_annotation_hyper <- race_TMB_result %>% 
    filter(Status == "hyper") %>%  
    mutate(label = paste0("MR = ", round(MR, 2)),
           x_pos = 1.5,   
           y_pos = 2.15)  
  
  picture_race_hyper <- ggplot(hyper_race_data, 
                               aes(x = Age_class, y = logTMB, fill = Age_class)) +
    geom_boxplot(outlier.shape = T,
                 width =0.5) +
    facet_wrap('Race', scales = "free",nrow =1)+
    scale_fill_manual(values = palette)+
    theme_pubr() +
    labs(title = "Hypermutated CRC",
         x = "",
         y = "logTMB"
    ) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.05),
                       limits = c(1.16,2.2)) +
    geom_text(data = race_annotation_hyper,
              aes(x = x_pos, y = y_pos, label = label),
              inherit.aes = FALSE,
              size = 5, color = "black") +
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
      strip.text = element_text(size = 14)
    )
  
  picture_race_hyper
  
}

pdf(paste0(dir,'picture/2.1_race_hyper_TMB_compare.pdf'),
    width=9,height =5)
picture_race_hyper
dev.off()


{
  
  nonhyper_race_data = TMB_mutload_for_regression[TMB<=15] %>% 
    left_join(rt@clinical.data[,.(Tumor_Sample_Barcode,Race,Age_class)],
              by = 'Tumor_Sample_Barcode') %>% 
    filter(Race %in% races_names) %>% 
    mutate(Age_class = factor(Age_class,levels = c('LOCRC','EOCRC')),
           adjTMB = 1.31 - logTMB)
  
  summary(nonhyper_race_data$adjTMB)
  
  race_annotation_nonhyper <- race_TMB_result %>% 
    filter(Status == "nonhyper") %>%  
    mutate(label = paste0("MR = ", round(MR, 2)),
           x_pos = 1.5,   
           y_pos = 2) 
  
  picture_race_nonhyper <- ggplot(nonhyper_race_data, 
                                  aes(x = Age_class, y = adjTMB, fill = Age_class)) +
    geom_boxplot(outlier.shape = T,
                 width =0.5) +
    facet_wrap('Race', scales = "free",nrow =1)+
    scale_fill_manual(values = palette)+
    theme_pubr() +
    labs(title = "Non-hypermutated CRC",
         x = "",
         y = "adjTMB"
    ) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.05),
                       limits = c(0.1,2.1)) +
    geom_text(data = race_annotation_nonhyper,
              aes(x = x_pos, y = y_pos, label = label),
              inherit.aes = FALSE,
              size = 5, color = "black") +
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
      strip.text = element_text(size = 14)
    )
  
  picture_race_nonhyper

  
}

pdf(paste0(dir,'picture/2.1_race_nonhyper_TMB_compare.pdf'),
    width=9,height =5)
picture_race_nonhyper
dev.off()


