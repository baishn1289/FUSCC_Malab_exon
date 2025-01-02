#race stratified
#TMB_calculation

library(tidyverse)
library(ggpubr)
library(ggthemes)
library(car)
dir = './picture'

regression_race_model <-  function(regression_data, zhongzu){
  
  race_data = regression_data %>%  filter(Race  == zhongzu)
  
  model = lm(logTMB ~ Age_class+Panel+Sex+Tumor_site+Sample_type+Histology_type+Country,
             data = race_data)
  
  race_data = race_data %>%
              mutate(residuals = residuals(model))
              
  assign(paste0(zhongzu, "_data"), race_data, envir = .GlobalEnv)
}

for (i in c("Asian Pacific lslander","Black","White")){
  regression_race_model(regression_data = TMB_mutload,zhongzu = i)
}

  race_data_collective  = rbind(`Asian Pacific lslander_data`,
                                 Black_data,White_data)
  
  write.csv(race_data_collective,file = './table/2.1.1_across_race_data_all.csv')
  
  palette <- c('#deebf7','#3182bd')


{
  #hyper_race_TMB_compare
  
  hyper_race_data = race_data_collective[TMB>15]

  picture_race_hyper <- ggplot(hyper_race_data, 
                               aes(x = Age_class, y = residuals, fill = Age_class)) +
                        geom_boxplot(outlier.shape = T,width =0.5) + 
                        facet_wrap('Race', scales = "free",nrow =1)+
                        scale_fill_manual(values = palette)+
                        theme_pubr() +
                        labs(title = "adjusted_TMB_EOCRC_vs_LOCRC", 
                             subtitle = "Hypermutated_CRC", 
                             x = "", 
                             y = "adjusted_TMB(residuals)" ) +
                        geom_signif(comparisons =  list(c("EOCRC", "LOCRC")),
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
                          panel.border = element_rect(linetype = "solid", linewidth = 1.2,fill = NA),
                          strip.text = element_text(size = 14) 
                        )

}

  pdf(paste(dir,'/2.1_race_hyper_TMB_compare.pdf',sep = ''),
      width=8,height =6)
  picture_race_hyper
  dev.off()


{
  #nonhyper_race_TMB_compare
  
  nonhyper_race_data = race_data_collective[TMB<= 15]
  
  picture_race_nonhyper <- ggplot(nonhyper_race_data, aes(x = Age_class, y = residuals, fill = Age_class)) +
                            geom_boxplot(outlier.shape = T, width = 0.5) + 
                            facet_wrap('Race', scales = "free",nrow =1) +  
                            scale_fill_manual(values = palette) +
                            theme_pubr() +
                            labs(title = "adjusted_TMB_EOCRC_vs_LOCRC", 
                                 subtitle = "Non-hypermutated_CRC", 
                                 x = "", 
                                 y = "adjusted_TMB(residuals)") +
                            geom_signif(comparisons = list(c("EOCRC", "LOCRC")),
                                        tip_length = 0.02,
                                        test = wilcox.test,
                                        textsize = 5 ) +
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

}

  pdf(paste(dir,'/2.1_race_nonhyper_TMB_compare.pdf',sep = ''),
      width=8,height =6)
  picture_race_nonhyper
  dev.off()

regression_data_cal <- function(data,tmb,race){
  
  regression_race_data = data.table(
                            TMB_Status = tmb,
                            Race = race,
                            EO_number = nrow(data %>%  filter(Race == race,Age_class == 'EOCRC') ),
                            LO_number = nrow(data %>%  filter(Race == race,Age_class == 'LOCRC') ),
                            pval =   wilcox.test(data %>% 
                                                   filter(Race == race,Age_class == 'EOCRC') %>% pull(residuals),
                                                 data %>% 
                                                   filter(Race == race,Age_class == 'LOCRC') %>% pull(residuals))$p.value,
                            effect_size = cliff.delta(data %>% 
                                                        filter(Race == race,Age_class == 'EOCRC') %>% pull(residuals),
                                                      data %>% 
                                                        filter(Race == race,Age_class == 'LOCRC') %>% pull(residuals))$estimate,
                            Magnitude =  as.character(cliff.delta(data %>% 
                                                                       filter(Race == race,Age_class == 'EOCRC') %>% pull(residuals),
                                                                  data %>% 
                                                                       filter(Race == race,Age_class == 'LOCRC') %>% pull(residuals))$magnitude),
                            EO_mean = as.numeric(summary(data %>%  filter(Race == race,Age_class == 'EOCRC')%>% pull(residuals))['Mean']),
                            LO_mean = as.numeric(summary(data %>%  filter(Race == race,Age_class == 'LOCRC')%>% pull(residuals))['Mean']),
                            EO_median = as.numeric(summary(data %>%  filter(Race == race,Age_class == 'EOCRC')%>% pull(residuals))['Median']),
                            LO_median = as.numeric(summary(data %>%  filter(Race == race,Age_class == 'LOCRC')%>% pull(residuals))['Median']) )
                          
  regression_EO_vs_LO_race_data <<-  rbind(regression_race_data,regression_EO_vs_LO_race_data)
}

  regression_EO_vs_LO_race_data <- data.frame(
                                TMB_Status = character(),
                                Race = character(),
                                EO_number = numeric(),
                                LO_number = numeric(),
                                pval =   numeric(),
                                effect_size =numeric(),
                                Magnitude = character(),
                                EO_mean = numeric(),
                                LO_mean = numeric(),
                                EO_median = numeric(),
                                LO_median = numeric() )
                                                        

for ( zhongzu in c("Asian Pacific lslander","Black","White") ){
  
  regression_data_cal(data = hyper_race_data ,tmb = 'Hyper',race = zhongzu)
  regression_data_cal(data = nonhyper_race_data ,tmb = 'Non-hyper',race = zhongzu)
  
}

write.csv(
  regression_EO_vs_LO_race_data,file = './table/2.1_regression_EO_vs_LO_race_data.csv'
)
