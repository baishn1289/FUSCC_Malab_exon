#race stratified

#TMB_calculation

#单独建模，原理参照第二段代码

library(tidyverse)
library(ggpubr)
library(ggthemes)

dir = './picture'


regression_model <-  function(regression_data, race,age_class){
  
  race_data = regression_data %>% filter(Race  == race)
  
  model = lm(TMB~Panel+Country+Sex+Tumor_site+Sample_type+Histology_type,
             data = race_data)
  
  race_data = race_data %>% 
    mutate(residuals = residuals(model)) %>% 
    filter(residuals > quantile(residuals, 0.01) & residuals < quantile(residuals, 0.99))
  
  assign(paste0(race, "_", age_class), race_data, envir = .GlobalEnv)
}


for (i in c("Asian Pacific lslander","Black","White")){
  
  regression_model(regression_data = EO_regression,race = i,age_class = 'EO')
  regression_model(regression_data = LO_regression,race = i,age_class = 'LO')
  
}


race_data_collective  = rbind(
  `Asian Pacific lslander_EO`,`Asian Pacific lslander_LO`,
  White_EO,White_LO,Black_EO,Black_LO
)

palette <- c('#deebf7','#3182bd')


{
  #hyper_race_TMB_compare
  
  hyper_race_data = race_data_collective[TMB>10]
  
  summary(hyper_race_data$residuals)
  
  picture_race_hyper <- ggplot(hyper_race_data, aes(x = Age_class, y = residuals, fill = Age_class)) +
    geom_boxplot(outlier.shape = T,
                 width =0.5) + 
    facet_wrap('Race', scales = "free",nrow =1)+
    scale_fill_manual(values = palette)+
    theme_pubr() +
    
    labs(title = "adjusted_TMB_EOCRC_vs_LOCRC", 
         subtitle = "hypermutated_CRC", 
         x = "", 
         y = "adjusted_TMB(residuals)"
         
    ) +
    geom_signif(comparisons =  list(c("EOCRC", "LOCRC")),
                # y_position = 270,
                tip_length = 0.02,
                #map_signif_level = p, 
                test = wilcox.test,textsize = 5
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
      strip.text = element_text(size = 14) 
    )
  
  picture_race_hyper
  
  pdf(paste(dir,'/2.1_race_hyper_TMB_compare.pdf',sep = ''),
      width=8,height =6)
  picture_race_hyper
  dev.off()
  
}



{
  #nonhyper_race_TMB_compare
  
  nonhyper_race_data = race_data_collective[TMB<=10]
  
  summary(nonhyper_race_data$residuals)
  
  nonhyper_race_data =left_join(nonhyper_race_data,y_position_data,
                                by = "Race")
  
  picture_race_nonhyper <- ggplot(nonhyper_race_data, aes(x = Age_class, y = residuals, fill = Age_class)) +
    geom_boxplot(outlier.shape = T, width = 0.5) + 
    facet_wrap('Race', scales = "free",nrow =1) +  
    scale_fill_manual(values = palette) +
    theme_pubr() +
    
    labs(title = "adjusted_TMB_EOCRC_vs_LOCRC", 
         subtitle = "nonhypermutated_CRC", 
         x = "", 
         y = "adjusted_TMB(residuals)") +
    
    
    geom_signif(
      comparisons = list(c("EOCRC", "LOCRC")),
      tip_length = 0.02,
      test = wilcox.test,
      textsize = 5
    ) +
    
    theme(
      plot.title    = element_text(color = "black", size   = 16, hjust = 0.5),
      plot.subtitle = element_text(color = "black", size   = 16, hjust = 0.5),
      plot.caption  = element_text(color = "black", size   = 16, face = "italic", hjust = 1),
      axis.text.x   = element_text(color = "black", size = 16, angle = 0),
      axis.text.y   = element_text(color = "black", size = 16, angle = 0),
      axis.title.x  = element_text(color = "black", size = 16, angle = 0),
      axis.title.y  = element_text(color = "black", size = 16, angle = 90),
      legend.position = "none",
      axis.line.y = element_line(color = "black", linetype = "solid"), 
      axis.line.x = element_line(color = "black", linetype = "solid"), 
      panel.border = element_rect(linetype = "solid", linewidth = 1.2, fill = NA),
      strip.text = element_text(size = 14) 
    )
  
  picture_race_nonhyper
  
  
  pdf(paste(dir,'/2.1_race_nonhyper_TMB_compare.pdf',sep = ''),
      width=8,height =6)
  picture_race_nonhyper
  dev.off()
  
}


