library(maftools)
library(mclust)
library(NMF)
library(data.table)
library(tidyverse)

# fangwuchu()

setwd('./Rcode')
# rm(list=ls())
gc()

dir= './picture'

{
  for (i in list.files(pattern = "rawdata_[A-Z][0-9]+\\.Rdata$")) {
    load(i)
  }
  
  load('FUSCC_2020to2024.Rdata')
  
  #cl_A7$`Family history of cancer`=as.character(cl_A7$`Family history of cancer`)
  
  cl <- dplyr::bind_rows(cl_FUSCC_2020to2024, 
                         #cl_A7,cl_C9,
                         cl_A8,cl_A10,
                         cl_C8, cl_C10,
                         cl_D1,cl_D7
  )
  
  
  table(cl$Country)
  
  cl = cl[!duplicated(cl$Tumor_Sample_Barcode),]
  table(cl$Country)
  
  rt = merge_mafs(list(rt_FUSCC_2020to2024,
                       #rt_A7,rt_C9,
                       rt_A8,rt_A10,
                       rt_C8,rt_C10,
                       rt_D1,rt_D7
  ))
  
  rt_back = rt
# rt=rt_back
  

# patient_country_count = bind_rows(rt_A7@clinical.data,
  #                                   rt_A8@clinical.data,
  #                                   rt_A10@clinical.data,
  #                                   rt_C8@clinical.data,
  #                                   rt_C9@clinical.data,
  #                                   rt_C10@clinical.data,
  #                                   rt_D1@clinical.data,
  #                                   rt_D7@clinical.data,
  #                                   rt_FUSCC_2020to2024@clinical.data
  #                                   )
  # patient_MSKCC =  bind_rows(rt_A7@clinical.data,
  #                            rt_A8@clinical.data,
  #                            rt_C9@clinical.data)
  # 
  # patient_country_count = patient_country_count[!duplicated(patient_country_count$Patient_ID),]
  # 
  
  # 
  # rm(cl_A7,
  #    cl_A8,cl_A10,
  #    cl_C8,cl_C9,cl_C10,cl_D1,cl_D7
  #    ,cl_FUSCC_2020to2024
  #    
  # )
  # 
  # rm(rt_A7,
  #    rt_A8,rt_A10,
  #    rt_C8,rt_C9,rt_C10,rt_D1,rt_D7
  #    ,rt_FUSCC_2020to2024
  #    
  # )
  
  rt@clinical.data$Country <- as.factor(rt@clinical.data$Country)
  rt@clinical.data$Age_class <- as.factor(rt@clinical.data$Age_class)
  
  # rt@data <- rt@data[,c(1:31)]
  # rt@maf.silent <-rt@maf.silent[,c(1:31)]
  
  
  
  { #Race数据整理
    
    table(cl$Number,cl$Race)
    
    cl<- cl%>% 
      mutate(Race= ifelse(is.na(Race), 'Unknown', Race)) %>% 
      mutate(Race= ifelse(Race %in% c('Not Applicable','Not Collected','Other',
                                      'UNKNOWN_OTHER','Native American',
                                      'NATIVE AMERICAN-AM IND/ALASKA'), 
                          'Unknown', Race)) %>% 
      mutate(Race= ifelse(Race %in% c('Asian','ASIAN-FAR EAST/INDIAN SUBCONT',
                                      'East Asian','Pacific Islander'), 
                          'Asian Pacific lslander', Race)) %>% 
      mutate(Race= ifelse(Race %in% c('BLACK OR AFRICAN AMERICAN'),
                          'Black',Race) )%>% 
      mutate(Race= ifelse(Race %in% c('WHITE'),
                          'White',Race) )
    
    table(rt@clinical.data$Tumor_Sample_Barcode %in% cl$Tumor_Sample_Barcode)
    # 以下代码错误，并不是一一匹配
    # merged_data <- merge(rt@clinical.data[, c('Tumor_Sample_Barcode')],
    #                      cl[, c('Tumor_Sample_Barcode', 'Race')],
    #                      by = 'Tumor_Sample_Barcode')
    
    rt@clinical.data <- left_join(rt@clinical.data[,-c('Race')],cl[, c('Tumor_Sample_Barcode', 'Race')],
                                  by ='Tumor_Sample_Barcode')
  }
  
  
  {#数据整理
    rt@clinical.data <- rt@clinical.data %>%
      mutate(
        Sex = case_when(
          is.na(Sex) ~ "Unknown", 
          Sex == "female" ~ "Female",  
          Sex == "male" ~ "Male",  
          TRUE ~ Sex  # 其他情况保持不变
                        ),
        Tumor_site = case_when(
          Tumor_site %in% c('Left hemicolon','Right hemicolon')~'Colon',  
          TRUE ~ Tumor_site
                              ),
        Histology_type =  case_when(
          Histology_type %in% c('Colon Adenocarcinoma In Situ')  ~'Colon Adenocarcinoma',
          Histology_type %in% c('Colorectal Adenocarcinoma')& Tumor_site == 'Colon'~'Colon Adenocarcinoma',
          Histology_type %in% c('Colorectal Adenocarcinoma')& Tumor_site == 'Rectum'~'Rectal Adenocarcinoma',
          Histology_type %in% c('Mucinous adenocarcinoma')& Tumor_site == 'Colon'~'Colon Mucinous adenocarcinoma',
          Histology_type %in% c('Mucinous adenocarcinoma')& Tumor_site == 'Rectum'~'Rectal Mucinous adenocarcinoma',
          Histology_type %in% c('Mucinous adenocarcinoma')& Tumor_site == 'Unknown'~'Colorectal Mucinous adenocarcinoma',
          Histology_type %in% c('Mucinous Adenocarcinoma of the Colon and Rectum')~'Colorectal Mucinous adenocarcinoma',
          TRUE ~ Histology_type)
      
        
      )
   
  }
  
  
  
  {#capture_size
    unique(rt@clinical.data$Panel)
    panel_size_manual <- data.table(
      SEQ_ASSAY_ID = c('FUSCC_ClinSeq','WESPlus gene panel','OncoPanel_AMC_version3','Sidra-LUMC'),
      total_size = c(2.3,65.755218,60,40)
    )
    
    capture_size <- rbind(panel_size_manual,panel_size) %>% 
      rename( Panel=SEQ_ASSAY_ID )
    
    # capture_size_test <- panel_size_test
    
    #GENIE_16.0版本数据计算的VHIO-300 靶向测序范围不具有生物学合理性，予以替换
    #UCSF-IDTV5-TN过大，暂时不替换
    #UCSF-IDTV5-TO过大,没有找到合适的capture size，暂时不替换
    capture_size[Panel == 'VHIO-300']$total_size =  1.3
    
    # capture_size[Panel == 'UCSF-IDTV5-TN']$total_size =  6.14
    # capture_size[Panel == 'UCSF-IDTV5-TO']$total_size = 
   
  }
  
  

  {#筛除有多个样本的患者
    
    rt_clinical_data_filter <- rt@clinical.data  %>% 
      mutate(
        
        Patient_ID = ifelse(Number == 'C10' & Center == 'MSK',
                            gsub("GENIE-MSK-", "", Patient_ID),
                            Patient_ID)
      )  
    
    #筛选至少有10个患者的panel
    rt_clinical_panel_filter <- rt_clinical_data_filter %>% 
      group_by(Panel) %>% 
      reframe(
        Panel_count = n_distinct(Patient_ID)
      ) %>% filter(Panel_count >=10)
    
    #先找出有多个样本的患者
    #在多个样本年龄不匹配的患者中，除去年龄大的样本
    #在每个患者最小年龄的样本中，保留原发的样本(除去转移或者不确定的样本)
    #如果仍然有多个原发样本，随机挑选1个原发样本
    
    
    
    multi_samples <- rt_clinical_data_filter[duplicated(rt_clinical_data_filter$Patient_ID) | 
                                        duplicated(rt_clinical_data_filter$Patient_ID, fromLast = TRUE), ] 
    
    multi_samples_info  <-  multi_samples %>%  
    group_by(Patient_ID) %>% 
      reframe(
          # MSI_same = n_distinct(MSI_Status) == 1,
          # Tumor_site_same = n_distinct(Tumor_site) == 1,
          # Stage_same = n_distinct(Stage) == 1,
          Age_same = n_distinct(Age) == 1,
          Number_same = n_distinct(Number) == 1,
          Min_Age = min(Age),
          # 筛出大龄样本
          Tumor_Sample_Barcode_diffage = list(Tumor_Sample_Barcode[Age != min(Age)]),
          # 统计 Min_Age 的样本数
          Min_Age_Sample_Count = sum(Age == min(Age)),
          
          Tumor_Sample_Barcode_sameage = list(Tumor_Sample_Barcode[Age == min(Age)]) # 确保返回空列表，而不是 NULL
        
    ) 
    
    sameage_tsb <- unlist(multi_samples_info$Tumor_Sample_Barcode_sameage)
    diffage_tsb <- unlist(multi_samples_info$Tumor_Sample_Barcode_diffage)
    
    
    #挑出有不同sample_type的患者下的样本
    multi_samples_same_age <-multi_samples %>% 
      filter(Tumor_Sample_Barcode %in% sameage_tsb) %>% 
      group_by(Patient_ID) %>% 
      reframe(
      sample_type_count = n_distinct(Sample_type)) %>% 
      filter(sample_type_count ==2) %>% 
      pull(Patient_ID)
      
    
    #对sample_type进一步筛选
   #有一西班牙患者sample type为转移和不确定两个样本，单独保留该患者的转移样本
    #转移：GENIE-VHIO-1266-001；不确定：GENIE-VHIO-1266-002
    multi_samples_diff_sampletype <- multi_samples %>% 
      filter(Patient_ID %in% multi_samples_same_age,
             Tumor_Sample_Barcode %in% sameage_tsb,
             Sample_type == 'Metastasis') %>% mutate(Tumor_Sample_Barcode = ifelse(Tumor_Sample_Barcode == 'GENIE-VHIO-1266-001',
                                                                                                                  'GENIE-VHIO-1266-002',
                                                                                                                  Tumor_Sample_Barcode)) %>% pull(Tumor_Sample_Barcode)
      
    barcode_discard <- c( multi_samples_diff_sampletype, diffage_tsb)
    
    multi_samples_tobe_filtered <- multi_samples %>% 
      filter(!Tumor_Sample_Barcode %in% barcode_discard) %>% 
      filter(duplicated(Patient_ID)) %>% pull(Tumor_Sample_Barcode)
   
    #删去所有不符合要求的样本
    barcode_new <- setdiff(rt@clinical.data$Tumor_Sample_Barcode, 
                           c(multi_samples_tobe_filtered,barcode_discard))
    
    
    rt@clinical.data <- rt@clinical.data[Tumor_Sample_Barcode %in% barcode_new] %>% 
      filter(Panel %in% rt_clinical_panel_filter$Panel) %>% 
      ##C10整理的时候已经筛除了小于0.2的panel，所以这边没有变化是正常的
      #这个代码只是做个兜底和检查
      filter(Panel %in% capture_size$Panel)
    
    
  } 
  
  rt <- subsetMaf(rt,tsb =rt@clinical.data$Tumor_Sample_Barcode )
  
  table(rt@clinical.data$Country)
  table(is.na(rt@data$Variant_Classification))
  table(rt@data$Variant_Classification)
  table(is.na(rt@clinical.data$Panel))
}

# 
# {#test
#   #Center_panel_information
#   table(rt@clinical.data$Panel,rt@clinical.data$Center)
#   table(is.na(rt@clinical.data$Panel))
#   table(rt@clinical.data$Number,rt@clinical.data$Panel)
#   table(is.na(rt@clinical.data$Center))
#   table(is.na(rt@clinical.data$Continent))
#   
#   
#   test.info = rt@clinical.data %>% 
#     filter(is.na(rt@clinical.data$Panel))
#   
#   unique(test.info$Number)
#   
#   table(is.na(rt@clinical.data$Number))
#   
#   test.info = rt@clinical.data %>% 
#     filter(Number =='NA')
#   
#   
#   table(duplicated(names(rt@clinical.data)))
#   
#   #   
#   # rt@clinical.data$Panel <- ifelse(rt@clinical.data$Number %in% c('A8','C8','C9'),
#   #                                  paste0('MSK-',rt@clinical.data$Panel),
#   #                                  rt@clinical.data$Panel)
#   #   
#   
# }


{#upset plot
  table(rt@clinical.data$Center)
  unique(rt@clinical.data$Panel)
  
  
  # rt@clinical.data <-  rt@clinical.data %>% 
  #   mutate(
  #     `Panel` = ifelse(Panel %in% c("IMPACT410", "IMPACT468", "IMPACT341",'MSK_US'), "MSKCC-IMPACT", `Panel`))
  
  
  
  
  #给每个Hugo_symbol加上center和country信息
  
  #有一个问题，rt@data中的数据不一定代表panel覆盖基因的全部，有可能只是没有测到被遗漏了
  #panel_hugo_symbol主要是创造一个数据框，每个样本（tsb）对应的panel以及这个样本测到的突变基因
  panel_hugo_symbol = left_join(rt@data[,c('Hugo_Symbol','Tumor_Sample_Barcode')],
                                rt@clinical.data[,c('Tumor_Sample_Barcode','Panel')],
                                by = 'Tumor_Sample_Barcode')
  
  #计算各个panel对应的基因数量
  panel_genenumber = panel_hugo_symbol %>%
    group_by(Panel) %>%        # 按Center和Country分组
    distinct(Hugo_Symbol) %>%            # 去除重复的Hugo_Symbol
    summarise(non_duplicate_count = n()) %>% as.data.frame()
  
  df_unique <- unique(panel_hugo_symbol[, c("Hugo_Symbol", "Panel")])
  
  df_wide <- df_unique %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = Panel, values_from = value, values_fill = 0) %>% 
    as.data.frame()
  
  rownames(df_wide) <- df_wide$Hugo_Symbol
  df_wide <- df_wide[, -1]
  
  
  #靶向测序取交集剩下的基因太少，计划放弃法国
  #尽量减少panel数量，筛选掉基因数量少于300的panel
  #美国人数足够，计划放弃美国患者人数少于300的panel
  filter_panel <-  panel_genenumber %>% 
    filter(non_duplicate_count <300)
  
  df_wide_filtered <- df_wide[, !colnames(df_wide) %in% filter_panel[,1]]
  
  table(rt@clinical.data$Center)
  
  #查看panel，center以及患者数量
  filtered_panel_country <- rt@clinical.data[,.(patients=.N),by=c('Panel','Center','Country')] %>% 
    filter(Panel %in% colnames(df_wide_filtered))  %>% 
    left_join(panel_genenumber,by = 'Panel')
  
  
  
  names(df_wide_filtered)
  df_wide_filtered <- df_wide_filtered[, -c(11:13)]
  
  country_patient_sum <- filtered_panel_country[, .(total_patients = sum(patients)), by = Country]
  
  library(UpSetR)
  pdf(file =  'targeted_panel.pdf',width = 12,height = 8)
  upset(df_wide_filtered, sets = colnames(df_wide_filtered), keep.order = TRUE, order.by = "freq", 
        main.bar.color = "dodgerblue", sets.bar.color = "orange")
  dev.off()
  
  
  pdf(file = 'all_panel.pdf',width =12,height = 15)
  upset(df_wide, sets = colnames(df_wide), keep.order = TRUE, order.by = "freq", 
        main.bar.color = "dodgerblue", sets.bar.color = "orange")
  dev.off()
  

  common_genes <- panel_hugo_symbol[Panel %in% colnames(df_wide_filtered)] %>%  
    group_by(Hugo_Symbol) %>%
    filter(n_distinct(Panel) == ncol(df_wide_filtered)) %>%
    distinct(Hugo_Symbol)
  
  
  
  
  clinical_filtered <- rt@clinical.data[Panel %in% colnames(df_wide_filtered)]
  #rt_panel_filtered:取了筛选后的panel的rt数据（注意，这里的rt@data仍然包含了筛选后panel的所有基因（相当于并集），并不是common_genes)
  rt_panel_filtered <- subsetMaf(rt,tsb = clinical_filtered$Tumor_Sample_Barcode)
  #rt_panel_common_genes取了panel的基因交集（commongenes）
  rt_panel_common_genes <- subsetMaf(rt_panel_filtered, genes = common_genes$Hugo_Symbol)
  
  
}




{
  table(cl$Age)
  
  # cl$Age_class <-  factor(cl$Age_class,levels = c('LOCRC','EOCRC'))
  # rt@clinical.data$Age_class <-  factor(rt@clinical.data$Age_class
  #                                       ,levels = c('LOCRC','EOCRC'))
  
  clin.EOCRC.sample <- subset(cl, Age <= 50)$Tumor_Sample_Barcode # %>% na.omit()
  clin.LOCRC.sample <- subset(cl, Age > 50)$Tumor_Sample_Barcode # %>% na.omit()
  
  length(clin.EOCRC.sample)
  length(clin.LOCRC.sample)
  
  
  clin.EOCRC <- subsetMaf(maf=rt, tsb=clin.EOCRC.sample, isTCGA=FALSE)
  clin.LOCRC <- subsetMaf(maf=rt, tsb=clin.LOCRC.sample, isTCGA=FALSE)
  
  #rt_panel_filtered
  clin.EOCRC_panel_filtered <- subsetMaf(maf=rt_panel_filtered, tsb=clin.EOCRC.sample, isTCGA=FALSE)
  clin.LOCRC_panel_filtered <- subsetMaf(maf=rt_panel_filtered, tsb=clin.LOCRC.sample, isTCGA=FALSE)
  
  
  #rt_panel_common_genes
  clin.EOCRC_panel_common_genes <- subsetMaf(maf=rt_panel_common_genes, tsb=clin.EOCRC.sample, isTCGA=FALSE)
  clin.LOCRC_panel_common_genes <- subsetMaf(maf=rt_panel_common_genes, tsb=clin.LOCRC.sample, isTCGA=FALSE)
  
  
  current_data = format(Sys.time(),"%Y-%m-%d")
  save(cl,rt,clin.EOCRC,clin.LOCRC,
       file=paste(current_data,"overall_data.Rdata",sep = '_'))
  
  save(filter_panel,panel_hugo_symbol,capture_size,
       file=paste(current_data,"filtered_panel_data.Rdata",sep = '_'))
  
  save(clinical_filtered,rt_panel_filtered,clin.EOCRC_panel_filtered,clin.LOCRC_panel_filtered,
       file=paste(current_data,"panel_filtered_data.Rdata",sep = '_'))
  
  save(rt_panel_common_genes,clin.EOCRC_panel_common_genes,clin.LOCRC_panel_common_genes,
       file=paste(current_data,"common_genes_data.Rdata",sep = '_'))
  
}

