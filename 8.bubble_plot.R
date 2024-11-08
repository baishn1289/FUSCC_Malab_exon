library(ggplot2)
library(ggthemes)
library(maftools)
library(dplyr)
library(stringr)
library(data.table)
library(tidyr)
library(tidyverse)
library(cowplot)

# dir = './picture/table_9.30'

barplot_tmb_status_hyper_EO <- clin.EOCRC_hyper@clinical.data %>% 
  group_by(Country) %>% 
  summarise(
    hyper = n_distinct(Tumor_Sample_Barcode),
  ) %>% 
  ungroup() %>%
 mutate(Age_class = 'EOCRC') 

barplot_tmb_status_hyper_LO <- clin.LOCRC_hyper@clinical.data %>% 
  group_by(Country) %>% 
  summarise(
    hyper = n_distinct(Tumor_Sample_Barcode),
  ) %>% 
  ungroup() %>% 
  mutate(Age_class = 'LOCRC') 


barplot_tmb_status_nonhyper_EO <- clin.EOCRC_nonhyper@clinical.data %>% 
  group_by(Country) %>% 
  summarise(
    nonhyper = n_distinct(Tumor_Sample_Barcode),
  ) %>% 
  ungroup() %>% 
  mutate(Age_class = 'EOCRC') 


barplot_tmb_status_nonhyper_LO <- clin.LOCRC_nonhyper@clinical.data %>% 
  group_by(Country) %>% 
  summarise(
    nonhyper = n_distinct(Tumor_Sample_Barcode),
  ) %>% 
  ungroup() %>% 
  mutate(Age_class = 'LOCRC') 

barplot_tmb_status_EO <- full_join(barplot_tmb_status_nonhyper_EO,
                                   barplot_tmb_status_hyper_EO,
                                   by =c("Country", "Age_class") )%>%
  gather(key = "Type", value = "Value", -Country, -Age_class) %>%
  group_by(Country, Age_class) %>%
  mutate(Percentage = Value / sum(Value) * 100)

barplot_tmb_status_LO <- full_join(barplot_tmb_status_nonhyper_LO,
                                   barplot_tmb_status_hyper_LO,
                                   by =c("Country", "Age_class") )%>%
  gather(key = "Type", value = "Value", -Country, -Age_class) %>%
  group_by(Country, Age_class) %>%
  mutate(Percentage = Value / sum(Value) * 100)


t.test(barplot_tmb_status_LO[barplot_tmb_status_LO$Type == 'hyper',]$Percentage,
       barplot_tmb_status_EO[barplot_tmb_status_EO$Type == 'hyper'& barplot_tmb_status_EO$Country != "Korea",]$Percentage)

t.test(barplot_tmb_status_LO[barplot_tmb_status_LO$Type == 'nonhyper',]$Percentage,
       barplot_tmb_status_EO[barplot_tmb_status_EO$Type == 'nonhyper'& barplot_tmb_status_EO$Country != "Korea",]$Percentage)




p_EO <- ggplot(barplot_tmb_status_EO, aes(x=Country, y=Percentage, fill=Type)) +
  geom_bar(stat="identity") +
  facet_wrap(~Age_class, nrow=1) +
  scale_fill_manual(values=c("hyper"="#d8bfd8", "nonhyper"="#6495ed")) +
  labs(y="TMB_Status Proportion(%)", x="Country", fill="Age_class") +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill="white", colour="black", linewidth=1),
    strip.text = element_text(size=12, face="bold"),
    axis.text.x = element_text(angle=45, hjust=1),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p_LO <- ggplot(barplot_tmb_status_LO, aes(x=Country, y=Percentage, fill=Type)) +
  geom_bar(stat="identity") +
  facet_wrap(~Age_class, nrow=1) +
  scale_fill_manual(values=c("hyper"="#d8bfd8", "nonhyper"="#6495ed")) +
  labs(y=NULL, x="Country", fill="Age_class") +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill="white", colour="black", linewidth=1),
    strip.text = element_text(size=12, face="bold"),
    axis.text.x = element_text(angle=45, hjust=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )



print(p_EO)
print(p_LO)

library(gridExtra)
pdf(paste(dir,"/8_TMB_Status_distribution.pdf", sep = ""), width = 9, height = 6)
grid.arrange(p_EO, p_LO, nrow = 1,widths = c(1, 1.25))
dev.off()


{
  # US_sample <- rt@clinical.data[Country =='US']$Tumor_Sample_Barcode

  barplot_tmb_status_hyper_EO_US <- clin.EOCRC_hyper@clinical.data %>%
    # filter(Tumor_Sample_Barcode%in% US_sample) %>%
    group_by(Panel) %>%
    summarise(
      hyper = n_distinct(Tumor_Sample_Barcode),
    ) %>%
    ungroup() %>%
    mutate(Age_class = 'EOCRC')

  barplot_tmb_status_hyper_LO_US <- clin.LOCRC_hyper@clinical.data %>%
    # filter(Tumor_Sample_Barcode%in% US_sample) %>%
    group_by(Panel) %>%
    summarise(
      hyper = n_distinct(Tumor_Sample_Barcode),
    ) %>%
    ungroup() %>%
    mutate(Age_class = 'LOCRC')



  barplot_tmb_status_nonhyper_EO_US <- clin.EOCRC_nonhyper@clinical.data %>%
    # filter(Tumor_Sample_Barcode%in% US_sample) %>%
    group_by(Panel) %>%
    summarise(
      nonhyper = n_distinct(Tumor_Sample_Barcode),
    ) %>%
    ungroup() %>%
    mutate(Age_class = 'EOCRC')


  barplot_tmb_status_nonhyper_LO_US <- clin.LOCRC_nonhyper@clinical.data %>%
    # filter(Tumor_Sample_Barcode%in% US_sample) %>%
    group_by(Panel) %>%
    summarise(
      nonhyper = n_distinct(Tumor_Sample_Barcode),
    ) %>%
    ungroup() %>%
    mutate(Age_class = 'LOCRC')

  barplot_tmb_status_EO_US <- full_join(barplot_tmb_status_nonhyper_EO_US,
                                        barplot_tmb_status_hyper_EO_US,
                                        by =c("Panel", "Age_class") ) %>%
    gather(key = "Type", value = "Value", -Panel, -Age_class) %>%

    # 将NA值替换为0
    mutate(Value = replace_na(Value, 0)) %>%

    group_by(Panel, Age_class) %>%
    mutate(Percentage = Value / sum(Value) * 100)

  barplot_tmb_status_LO_US <- full_join(barplot_tmb_status_nonhyper_LO_US,
                                        barplot_tmb_status_hyper_LO_US,
                                        by =c("Panel", "Age_class") ) %>%
    gather(key = "Type", value = "Value", -Panel, -Age_class) %>%

    # 将NA值替换为0
    mutate(Value = replace_na(Value, 0)) %>%

    group_by(Panel, Age_class) %>%
    mutate(Percentage = Value / sum(Value) * 100)

  sum(is.na(barplot_tmb_status_EO_US$Percentage))
  sum(is.na(barplot_tmb_status_LO_US$Percentage))
  summary(barplot_tmb_status_EO_US$Percentage)
  summary(barplot_tmb_status_LO_US$Percentage)


  p_EO_US <- ggplot(barplot_tmb_status_EO_US, aes(x=Panel, y=Percentage, fill=Type)) +
    geom_bar(stat="identity") +
    facet_wrap(~Age_class, nrow=1) +
    scale_fill_manual(values=c("hyper"="#d8bfd8", "nonhyper"="#6495ed")) +
    labs(y="TMB_Status Proportion(%)", x="Panel", fill="Age_class") +
    theme_minimal() +
    theme(
      strip.background = element_rect(fill="white", colour="black", linewidth=1),
      strip.text = element_text(size=12, face="bold"),
      axis.text.x = element_text(angle=45, hjust=1),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  p_LO_US <- ggplot(barplot_tmb_status_LO_US, aes(x=Panel, y=Percentage, fill=Type)) +
    geom_bar(stat="identity") +
    facet_wrap(~Age_class, nrow=1) +
    scale_fill_manual(values=c("hyper"="#d8bfd8", "nonhyper"="#6495ed")) +
    labs(y=NULL, x="Panel", fill="Age_class") +
    theme_minimal() +
    theme(
      strip.background = element_rect(fill="white", colour="black", linewidth=1),
      strip.text = element_text(size=12, face="bold"),
      axis.text.x = element_text(angle=45, hjust=1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
    )

  print(p_EO_US)
  print(p_LO_US)

  library(gridExtra)
  pdf(paste(dir,"/8_TMB_Status_distribution_allpanels.pdf", sep = ""), width = 30, height = 6)
  grid.arrange(p_EO_US, p_LO_US, nrow = 1,widths = c(1, 1.25))
  dev.off()



  }




#bubble plot
# setwd('C:/Users/baish/Desktop/TMB')
# dir = 'picture'

# load('2024-09-28_overall_data.Rdata')
# load('hyper_data.Rdata')
# load('nonhyper_data.Rdata')
# load('2024-09-28_filtered_panel_data.Rdata')

#筛选三个TMB范围共同基因
# filtered_df_overall = read_tsv(paste(dir,'/overall_EOCRC_vs_LOCRC_freq.tsv',sep='')) %>% 
#   filter(adjPval <0.05)
# filtered_df_nonhyper = read_tsv(paste(dir,'/nonhyper_EOCRC_vs_LOCRC_freq.tsv',sep=''))%>% 
#   filter(adjPval <0.05)
# filtered_df_hyper = read_tsv(paste(dir,'/hyper_EOCRC_vs_LOCRC_freq.tsv',sep=''))%>% 
#   filter(adjPval <0.05)
# 
# insect_gene = Reduce(intersect, list(filtered_df_hyper$Hugo_Symbol,
#                       filtered_df_nonhyper$Hugo_Symbol,
#                     filtered_df_overall$Hugo_Symbol))



uniquegenes <- function(m1, m2, m1Name = NULL, m2Name = NULL,
                        minMut = 5){
  
  # numMutatedSamples = maf[!Variant_Type %in% 'CNV', .(MutatedSamples = length(unique(Tumor_Sample_Barcode))), by = Hugo_Symbol]
  # numAlteredSamples = maf[, .(AlteredSamples = length(unique(Tumor_Sample_Barcode))), by = Hugo_Symbol]
  # numAlteredSamples = merge(numMutatedSamples, numAlteredSamples, by = 'Hugo_Symbol', all = TRUE)
  
  m1.gs <- getGeneSummary(x = m1)
  m2.gs <- getGeneSummary(x = m2)
  
  if (is.null(m1Name)) {
    m1Name = "M1"
  }
  if (is.null(m2Name)) {
    m2Name = "M2"
  }
  
  m1.genes = as.character(m1.gs[AlteredSamples >= minMut, 
                                Hugo_Symbol])
  m2.genes = as.character(m2.gs[AlteredSamples >= minMut, 
                                Hugo_Symbol])
  uniquegenes =  intersect(m1.genes, m2.genes)
  
  return(uniquegenes)
  
}

# m1 = rt_country_EO
# m2 = rt_country_LO
# country='Spain'
# Clinical.data = rt_clinicaldata_country
# Gene_test = 'TP53'

logstic_country_mafcompare <- function(m1, m2, Clinical.data, Gene_test,Guojia,TMB_status) {
  
  
  
  #从maf@data中提取出来的Tumor_Sample_Barcode要去重，因为一个样本可能在一个基因上有多个突变外显子，这导致了一个TSB可能会被多次记录，不能代表样本数
  Alter_m1_Gene = m1@data[Hugo_Symbol == Gene_test,] %>% 
    filter(!duplicated(Tumor_Sample_Barcode))
  
  Alter_m2_Gene = m2@data[Hugo_Symbol == Gene_test,] %>% 
    filter(!duplicated(Tumor_Sample_Barcode))
  
  # ！！！！！检查函数时使用
  # Gene_panels_information <- panel_hugo_symbol[Hugo_Symbol == Gene_test] %>%
  #   merge(Clinical.data[, c('Tumor_Sample_Barcode', 'Age_class', 'Sex', 'Country', 'Race')], by = "Tumor_Sample_Barcode") %>%
  #   filter(!duplicated(Tumor_Sample_Barcode)) %>%
  #   group_by(Panel) %>%
  #   summarise(
  #     Age_class_count = n_distinct(Age_class),
  #     Sex_count = n_distinct(Sex),
  #     Race_count = n_distinct(Race),
  #     Country_count = n_distinct(Country)
  #   ) %>%
  #   ungroup()
  # 
  
  
  Gene_panels <- panel_hugo_symbol[Hugo_Symbol == Gene_test] %>%
    merge(Clinical.data[, c('Tumor_Sample_Barcode', 'Age_class', 'Sex', 'Country', 'Race')], by = "Tumor_Sample_Barcode") %>%
    filter(!duplicated(Tumor_Sample_Barcode)) %>%
    group_by(Panel) %>%
    summarise(
      Age_class_count = n_distinct(Age_class)
    ) %>%
    ungroup() %>% 
    filter(Age_class_count == 2) %>% 
    pull(Panel)
  
  # 如果 Gene_panels 为空，记录信息并返回
  if (length(Gene_panels) == 0) {
    print(paste(Gene_test, '没有一个中心有两种性别', sep = ' '))
    
    # 将当前的记录追加到全局记录日志中
    singlesex_panel_record_log <<- rbind(singlesex_panel_record_log, 
                                         data.frame(Gene_test = Gene_test,
                                                    Country = Guojia,
                                                    TMB_status = TMB_status,
                                                    stringsAsFactors = FALSE))
    
    return()
  }
  
 
  #基因突变频率计算
  #突变频率按照患者突变频率计算：突变患者人数/突变基因所在panel的患者总人数
  #只要患者有一个样本突变了，就认定该患者突变
  # table(is.na(cl$Patient_ID))
  #rt为所有测序的样本对应的患者临床信息
  # table(is.na(rt@clinical.data$Patient_ID))
  # table(duplicated(rt@clinical.data$Patient_ID))
  
  #提出该基因变异的样本数量
  mut_samples_m1 <- length(unique(Alter_m1_Gene$Tumor_Sample_Barcode))
  
  
  #提出该基因所在的panels的所有样本数量
  panels_samples_m1 <- length(unique(m1@clinical.data[m1@clinical.data$Panel %in% Gene_panels, ]$Tumor_Sample_Barcode))
  
  
  m1_freq <- mut_samples_m1 / panels_samples_m1
  
  # ------------
  mut_samples_m2 <- length(unique(Alter_m2_Gene$Tumor_Sample_Barcode))
  
  
  panels_samples_m2 <- length(unique(m2@clinical.data[m2@clinical.data$Panel %in% Gene_panels, ]$Tumor_Sample_Barcode))
  
  
  m2_freq <- mut_samples_m2 / panels_samples_m2
  
  
  #TMB_normalize
  # 以panels中所有样本的TMB平均值做标准化
  m1_TMB <-  left_join(m1@clinical.data[m1@clinical.data$Panel %in% Gene_panels, ],
                       TMB_mutload,by = 'Tumor_Sample_Barcode') %>% 
    summarise(
      TMB_Panel_mean = mean(TMB)
    ) %>% pull(TMB_Panel_mean)
  
  m1_freq_normalized <- m1_freq/m1_TMB
  
  
  m2_TMB <-  left_join(m2@clinical.data[m2@clinical.data$Panel %in% Gene_panels, ],
                       TMB_mutload,by = 'Tumor_Sample_Barcode') %>% 
    summarise(
      TMB_Panel_mean = mean(TMB)
    ) %>% pull(TMB_Panel_mean)
  
  m2_freq_normalized <- m2_freq/m2_TMB
  
  
  Gene_clinical_info = Clinical.data[, Gene_status := ifelse(Tumor_Sample_Barcode %in% c(Alter_m1_Gene$Tumor_Sample_Barcode,
                                                                                         Alter_m2_Gene$Tumor_Sample_Barcode), 1, 0)] %>% 
    left_join(TMB_mutload[,.(Tumor_Sample_Barcode,TMB)],by = 'Tumor_Sample_Barcode') %>% 
    filter(Panel %in% Gene_panels)
 
  
  if (length(unique(Gene_clinical_info$Race)) > 1) {
    formula_base <- paste(formula_base, "+ Race")
  }
  if (length(unique(Gene_clinical_info$Panel)) > 1) {
    formula_base <- paste(formula_base, "+ Panel")
  }
  if (length(unique(Gene_clinical_info$Country)) > 1) {
    formula_base <- paste(formula_base, "+ Country")
  }
  if (length(unique(Gene_clinical_info$Sex)) > 1) {
    formula_base <- paste(formula_base, "+ Sex")
  }
  
  
  fit <- glm(as.formula(formula_base), family = binomial(link = "logit"), data = Gene_clinical_info)

  OR_age_class <- exp(coef(fit)["Age_classEOCRC"])

  pval_age_class <- summary(fit)$coefficients["Age_classEOCRC", "Pr(>|z|)"]
  
 
  coef_age_class <- coef(fit)["Age_classEOCRC"]
  se_age_class <- summary(fit)$coefficients["Age_classEOCRC", "Std. Error"]
  ci_low <- exp(coef_age_class - 1.96 * se_age_class)
  ci_up  <- exp(coef_age_class + 1.96 * se_age_class)
  
  df_Gene_test <- data.table::data.table(
    Hugo_Symbol = Gene_test,
    EOCRC_total_case= panels_samples_m1,
    LOCRC_total_case= panels_samples_m2,
    EOCRC_mutated_case = mut_samples_m1, 
    LOCRC_mutated_case = mut_samples_m2, 
    EOCRC_freq = m1_freq,
    LOCRC_freq = m2_freq,
    EOCRC_freq_normalized = m1_freq_normalized,
    LOCRC_freq_normalized = m2_freq_normalized,
    EOCRC_panel_TMB= m1_TMB,
    LOCRC_panel_TMB= m2_TMB,
    pval = pval_age_class, 
    or = OR_age_class,
    ci.up = ci_up,
    ci.low = ci_low,
    formula = formula_base,
    TMB_Status = TMB_status,
    Country = Guojia
  )
  
  result_df <<- rbind(result_df, df_Gene_test)  
}




rt@clinical.data$Age_class <-  factor(rt@clinical.data$Age_class
                                      ,levels = c('LOCRC','EOCRC'))


# guojia = 'US'

generate_bubble_table <- function(guojia,Minmut){
  
  #挑选出当前国家的样本
  country_tsb<- rt@clinical.data[Country == guojia]$Tumor_Sample_Barcode
  
  
 
  # 确保当前国家在总和的老年和青年队列中都有对应样本
  if (length(intersect(clin.EOCRC@clinical.data$Tumor_Sample_Barcode, country_tsb)) > 0 &
      length(intersect(clin.LOCRC@clinical.data$Tumor_Sample_Barcode, country_tsb)) > 0){
      
      #挑出该国家在青年和老年中的对应样本
      rt_country_EO <- subsetMaf(clin.EOCRC,tsb=country_tsb)
      rt_country_LO <- subsetMaf(clin.LOCRC,tsb=country_tsb)
      rt_clinicaldata_country <- rt@clinical.data[Country == guojia]
    
      bubble_genes <- uniquegenes(m1 = rt_country_EO, m2 = rt_country_LO, minMut = Minmut)
      
      if(length(bubble_genes) == 0 ){
        return(NULL)
      }
     
      for (i in bubble_genes) {
        result <- tryCatch({
          logstic_country_mafcompare(m1 = rt_country_EO, m2 = rt_country_LO,
                                     Clinical.data = rt_clinicaldata_country, 
                                     Gene_test =i,Guojia = guojia,TMB_status = 'Overall')
          NULL  # 如果没有错误，返回NULL
        }, warning = function(w) {
          problematic_genes <<- 
            rbind(problematic_genes,
                  data.frame(problematic_gene = i,
                            country =guojia,
                            TMB_Staus = 'Overall'))# 记录产生警告的基因
          return(NULL)  # 返回NULL以继续循环
        })
      }
      
      country_comp <- result_df[-1,]
     
    
    
  
  }else{
    
    country_comp <- data.frame(Hugo_Symbol=NA, 
                               EOCRC_total_case= NA,
                               LOCRC_total_case= NA,
                               EOCRC_mutated_case = NA, 
                               LOCRC_mutated_case = NA,
                               EOCRC_freq =NA, 
                               LOCRC_freq =NA, 
                               EOCRC_freq_normalized = NA,
                               LOCRC_freq_normalized = NA,
                               EOCRC_panel_TMB= NA,
                               LOCRC_panel_TMB= NA,
                               pval =NA,
                               or =NA, 
                               ci.up =NA, 
                               ci.low =NA,
                               formula = NA, TMB_Status =NA, Country = NA)
  }
  
  
  if (length(intersect(clin.EOCRC_hyper@clinical.data$Tumor_Sample_Barcode, country_tsb)) > 0 &
      length(intersect(clin.LOCRC_hyper@clinical.data$Tumor_Sample_Barcode, country_tsb)) > 0) {
    
    rt_country_EO_hyper <- subsetMaf(clin.EOCRC_hyper,tsb=country_tsb)
    rt_country_LO_hyper <- subsetMaf(clin.LOCRC_hyper,tsb=country_tsb)
    Clinical_data_hyper <- Clinical_data_hyper[Country == guojia]
    
    
    
    bubble_genes_hyper <- uniquegenes(m1 = rt_country_EO_hyper, m2 = rt_country_LO_hyper,
                                      minMut = Minmut)
    
    if(length(bubble_genes_hyper) == 0 ){
      return(NULL)
    }
   
       
    for (i in bubble_genes_hyper) {
      result <- tryCatch({
        logstic_country_mafcompare(m1 = rt_country_EO_hyper, m2 = rt_country_LO_hyper,
                                   Clinical.data = Clinical_data_hyper, 
                                   Gene_test =i,Guojia = guojia,TMB_status = 'hyper')
        NULL  # 如果没有错误，返回NULL
      }, warning = function(w) {
        problematic_genes <<- 
          rbind(problematic_genes,
                data.frame(problematic_gene = i,
                           country =guojia,
                           TMB_Staus = 'hyper'))
        return(NULL)  # 返回NULL以继续循环
      })
    }
    
   country_comp_hyper <- result_df[-1,]
   
    #   mutate(
    #     Fold = EOCRC_freq/ LOCRC_freq,
    #     TMB = 'Hyper',
    #     Country = guojia
    #   )
  
  }else{
    
    country_comp_hyper <- data.frame(Hugo_Symbol=NA, 
                                     EOCRC_total_case= NA,
                                     LOCRC_total_case= NA,
                                     EOCRC_mutated_case = NA, 
                                     LOCRC_mutated_case = NA,
                                     EOCRC_freq =NA, 
                                     LOCRC_freq =NA, 
                                     EOCRC_freq_normalized = NA,
                                     LOCRC_freq_normalized = NA,
                                     EOCRC_panel_TMB= NA,
                                     LOCRC_panel_TMB= NA,
                                     pval =NA,
                                     or =NA, 
                                     ci.up =NA, 
                                     ci.low =NA,
                                     formula = NA, TMB_Status =NA, Country = NA)
    
  } 
  
  if (length(intersect(clin.EOCRC_nonhyper@clinical.data$Tumor_Sample_Barcode, country_tsb)) > 0 &
      length(intersect(clin.LOCRC_nonhyper@clinical.data$Tumor_Sample_Barcode, country_tsb)) > 0) {
   
    rt_country_EO_nonhyper <- subsetMaf(clin.EOCRC_nonhyper,tsb=country_tsb)
    rt_country_LO_nonhyper <- subsetMaf(clin.LOCRC_nonhyper,tsb=country_tsb)
    Clinical_data_nonhyper <- Clinical_data_nonhyper[Country == guojia]
    
    
    
    bubble_genes_nonhyper <- uniquegenes(m1 = rt_country_EO_nonhyper, m2 = rt_country_LO_nonhyper,
                                      minMut = Minmut)
    
    if(length(bubble_genes_nonhyper) == 0 ){
      return(NULL)
    }
    
  
    
    for (i in bubble_genes_nonhyper) {
      result <- tryCatch({
        logstic_country_mafcompare(m1 = rt_country_EO_nonhyper, m2 = rt_country_LO_nonhyper,
                                   Clinical.data = Clinical_data_nonhyper, 
                                   Gene_test =i,Guojia = guojia,TMB_status = 'nonhyper')
        NULL  # 如果没有错误，返回NULL
      }, warning = function(w) {
        problematic_genes <<- 
          rbind(problematic_genes,
                data.frame(problematic_gene = i,
                           country =guojia,
                           TMB_Staus = 'nonhyper'))
        return(NULL)  # 返回NULL以继续循环
      })
    }
    
    country_comp_nonhyper <- result_df[-1] 
    
    #   mutate(
    #     Fold = EOCRC_freq/ LOCRC_freq,
    #     TMB = 'nonhyper',
    #     Country = guojia
    #   )
    
    
  } else {
    country_comp_nonhyper <- data.frame(Hugo_Symbol=NA, 
                                        EOCRC_total_case= NA,
                                        LOCRC_total_case= NA,
                                        EOCRC_mutated_case = NA, 
                                        LOCRC_mutated_case = NA,
                                        EOCRC_freq =NA, 
                                        LOCRC_freq =NA, 
                                        EOCRC_freq_normalized = NA,
                                        LOCRC_freq_normalized = NA,
                                        EOCRC_panel_TMB= NA,
                                        LOCRC_panel_TMB= NA,
                                        pval =NA,
                                        or =NA, 
                                        ci.up =NA, 
                                        ci.low =NA,
                                        formula = NA, TMB_Status =NA, Country = NA)
  }
  
  
  bubble_table = rbind(bubble_table,country_comp_nonhyper,
                           country_comp_hyper,country_comp) 
  
  
  return(bubble_table)
  
}



formula_base <- "Gene_status ~ Age_class"

singlesex_panel_record_log <- data.frame(Gene_test = character(),
                          Country = character(),
                          TMB_status = character(),
                          stringsAsFactors = FALSE)


problematic_genes <- data.frame(
  problematic_gene = NA,
  country =NA,
  TMB_Staus = NA
)


result_df <- data.table::data.table(
  Hugo_Symbol=NA, 
  EOCRC_total_case= NA,
  LOCRC_total_case= NA,
  EOCRC_mutated_case = NA, 
  LOCRC_mutated_case = NA,
  EOCRC_freq =NA, 
  LOCRC_freq =NA, 
  EOCRC_freq_normalized = NA,
  LOCRC_freq_normalized = NA,
  EOCRC_panel_TMB= NA,
  LOCRC_panel_TMB= NA,
  pval =NA,
  or =NA, 
  ci.up =NA, 
  ci.low =NA,
  formula =NA,
  TMB_Status = NA,
  Country = NA
)


#根据国家和TMB范畴分队列，生成对应EO和LO比较数据
bubble_table = data.frame(Hugo_Symbol=NA, 
                          EOCRC_total_case= NA,
                          LOCRC_total_case= NA,
                          EOCRC_mutated_case = NA, 
                          LOCRC_mutated_case = NA,
                          EOCRC_freq =NA, 
                          LOCRC_freq =NA, 
                          EOCRC_freq_normalized = NA,
                          LOCRC_freq_normalized = NA,
                          EOCRC_panel_TMB= NA,
                          LOCRC_panel_TMB= NA,
                          pval =NA,
                          or =NA, 
                          ci.up =NA, 
                          ci.low =NA,
                          formula = NA, TMB_Status =NA, Country = NA)

total_bubble_table <-  data.frame(Hugo_Symbol=NA, 
                                  EOCRC_total_case= NA,
                                  LOCRC_total_case= NA,
                                  EOCRC_mutated_case = NA, 
                                  LOCRC_mutated_case = NA,
                                  EOCRC_freq =NA, 
                                  LOCRC_freq =NA, 
                                  EOCRC_freq_normalized = NA,
                                  LOCRC_freq_normalized = NA,
                                  EOCRC_panel_TMB= NA,
                                  LOCRC_panel_TMB= NA,
                                  pval =NA,
                                  or =NA, 
                                  ci.up =NA, 
                                  ci.low =NA,
                                  formula = NA, TMB_Status =NA, Country = NA)


  unique(rt@clinical.data$Country)

  France_bubble_table <- generate_bubble_table('France',Minmut = 2)
  Netherlands_bubble_table  <- generate_bubble_table('Netherlands',Minmut = 2)
  Canada_bubble_table  <- generate_bubble_table('Canada',Minmut = 2)
  Spain_bubble_table  <- generate_bubble_table('Spain',Minmut = 2)
  Nigeria_bubble_table  <- generate_bubble_table('Nigeria',Minmut = 2)
  China_bubble_table  <- generate_bubble_table('China',Minmut = 2)
  US_bubble_table  <- generate_bubble_table('US',Minmut = 2)
  
  
  
  total_bubble_table <- US_bubble_table[-1,]
  
  names(total_bubble_table)
  names(total_bubble_table)[names(total_bubble_table) =='Hugo_Symbol'] = 'Gene'
  names(total_bubble_table)[names(total_bubble_table) =='or'] = 'OR'


#### 筛除标准：overall范畴 EOCRC_Freq 小于 0.05 的基因和国家组合
#先拉出overall范畴的所有国家所有基因列表
filter_gene_country <- total_bubble_table %>%
  filter(TMB_Status == "Overall" & EOCRC_freq < 0.05)  %>%
  select(Gene, Country,EOCRC_freq)



###检查一下total_bubble_table_filter最后的结果中是否包含problem_genes

# 从 total_bubble_table 中删除这些基因和国家组合的行
total_bubble_table_filter <- total_bubble_table %>%
  anti_join(filter_gene_country, by = c("Gene", "Country"))%>% 
  mutate(adjPval = p.adjust(p = pval, method = "fdr")) 
  
total_bubble_table_filter$p_sig <- as.character(symnum(total_bubble_table_filter$adjPval, 
                                              cutpoints = c(0, 0.01, 0.05, 0.25, 1), 
                                              symbols = c("***", "**", "*", "")))



{
  # 筛选出出现在至少 6 个国家的基因
  # 统计基因对应的国家数
  gene_country_counts <- total_bubble_table_filter %>%
    group_by(Gene) %>%
    summarise(Country_count = n_distinct(Country))
  
  genes_in_6_countries <- gene_country_counts %>%
    filter(Country_count >= 6) %>%
    select(Gene)
  }

bubble_table = total_bubble_table_filter %>%
  filter(Gene %in%  genes_in_6_countries$Gene)


# p1 <- ggplot(data = bubble_table, 
#              aes(x = Gene, y = TMB_Status)) +
#   geom_point(aes(color = factor(case_when(
#     OR > 1 ~ 'OR > 1',
#     OR < 1 ~ 'OR < 1',
#     OR == 1 ~ 'No significance'
#   ), levels = c('OR > 1', 'OR < 1', 'No significance'))), 
#   size = 4, alpha = 0.7, show.legend = c(size = FALSE))  +
#   facet_grid(Country ~ ., scales = "free_y", space = "free") +  
#   scale_y_discrete(position = "left") +  
#   scale_color_manual(values = c('OR > 1' = '#d53e4f', 'OR < 1' = '#4393c3', 'No significance' = 'gray')) +
#   theme_few() +
#   labs(y ='TMB Category')+
#   theme(
#     strip.background = element_blank(), 
#     strip.text.y.right = element_text(size = 10, angle = 0),
#     panel.spacing = unit(0, "lines"),   
#     panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
#     axis.text.y = element_text(color = "black", size = 10, angle = 0),
#     axis.text.x = element_text(color = "black", size = 8, angle = 45, vjust = 0.5, hjust = 0.5),
#     legend.text = element_text(color = "black", size = 10),
#     legend.title = element_blank()
#   ) +
#   guides(color = guide_legend(override.aes = list(size = 4)))



p1 <- ggplot(data = bubble_table, 
             aes(x = Gene, y = TMB_Status)) +
  geom_point(aes(color = factor(case_when(
    OR > 1 ~ 'OR > 1',
    OR < 1 ~ 'OR < 1',
    OR == 1 ~ 'No significance'
  ), levels = c('OR > 1', 'OR < 1', 'No significance'))), 
  size = 6, alpha = 0.7, show.legend = c(size = FALSE))  +
  
  # 添加显著性标记
  geom_text(aes(label = p_sig), vjust = 0.8, hjust = 0.5, size = 3, color = "black"
            #vjust调整上下
            ) +
  facet_grid(Country ~ ., scales = "free_y", space = "free") +  
  scale_y_discrete(position = "left") +  
  scale_color_manual(values = c('OR > 1' = '#d53e4f', 'OR < 1' = '#4393c3', 'No significance' = 'gray')) +
  theme_few() +
  labs(y ='TMB Category') +
  theme(
    strip.background = element_blank(), 
    strip.text.y.right = element_text(size = 13, angle = 0),
    panel.spacing = unit(0, "lines"),   
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    axis.text.y = element_text(color = "black", size = 10, angle = 0),
    axis.text.x = element_text(color = "black", size = 10, angle = 45, vjust = 0.5, hjust = 0.5),
    legend.text = element_text(color = "black", size = 10),
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

p1

ggsave(paste(dir,'/8_intersect_gene_among_countries.pdf',sep = ''),
       width = 15,height = 6.5,dpi = 300)

write.csv(total_bubble_table,file = '8_total_bubble_data.csv',
          quote = FALSE, row.names = FALSE)

