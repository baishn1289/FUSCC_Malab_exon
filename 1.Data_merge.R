library(maftools)
library(mclust)
library(NMF)
library(pheatmap)
library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(UpSetR)

dir= './picture'

{
  for (i in list.files(pattern = "rawdata_[A-Z][0-9]+\\.Rdata$")) {
    load(i)
  }
  
  load('FUSCC_2020to2024.Rdata')
  
  
  cl <- dplyr::bind_rows(cl_FUSCC_2020to2024, 
                         cl_A8,cl_A10,
                         cl_C8, cl_C9,cl_C10,
                         cl_D1,cl_D7
  )
  
 
  
  rt = merge_mafs(list(rt_FUSCC_2020to2024,
                       rt_A8,rt_A10,
                       rt_C8,rt_C9,rt_C10,
                       rt_D1,rt_D7
  ))
  
  
  rt@clinical.data$Country <- as.factor(rt@clinical.data$Country)
  rt@clinical.data$Age_class <- as.factor(rt@clinical.data$Age_class)
 
  
  { 
    
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
   
    
    rt@clinical.data <- left_join(rt@clinical.data[,-c('Race')],cl[, c('Tumor_Sample_Barcode', 'Race')],
                                  by ='Tumor_Sample_Barcode')
  }
  
  {
    
    table(rt@clinical.data$Sex)
    rt@clinical.data <- rt@clinical.data %>%
      mutate(Sex = case_when(
        is.na(Sex) ~ "Unknown", 
        Sex == "female" ~ "Female",  
        Sex == "male" ~ "Male",  
        TRUE ~ Sex  
      ))
    
  }
  
  
  rt_backup <- rt@clinical.data 
  
  rt_clinical_data_test <- rt@clinical.data %>%
    filter(Center %in% c('MSK','MSKCC')) %>% 
    mutate(
      Tumor_Sample_Barcode_backup = Tumor_Sample_Barcode
    ) %>% 
    mutate(Patient_ID = ifelse(Number == 'C10' & Center == 'MSK',
                          gsub("GENIE-MSK-", "", Patient_ID),
                          Patient_ID),
      Tumor_Sample_Barcode = ifelse(Number == 'C10' & Center == 'MSK',
                                    gsub("GENIE-MSK-|\\-IM[0-9]+", "", Tumor_Sample_Barcode),
                                    Tumor_Sample_Barcode)
          ) 
  

  rt_patient_dup_num = rt_clinical_data_test$Patient_ID[duplicated(rt_clinical_data_test$Patient_ID)]
  
  rt_info_dup_patients <- rt_clinical_data_test %>%
    filter(Patient_ID %in% rt_patient_dup_num)
  
   
   
  {
     
      Number_dup_patient_A8_C9_C10 <- rt_info_dup_patients %>% 
        filter(
               Number %in% c('A8','C9','C10')) %>% 
        group_by(Number) %>% 
        reframe(
          Number_Tumor_Sample_Barcode = list(Tumor_Sample_Barcode))%>%
        ungroup()
      
      intersect(gsub("-IM[0-9]+", "",unlist(Number_dup_patient_A8_C9_C10[1,2])),
                gsub("-IM[0-9]+", "",unlist(Number_dup_patient_A8_C9_C10[3,2]))
                
                )
      
      dup_tumor_sample_barcodes_to_be_deleted <- c('P-0003046-T01-IM3')
      
      }
    
    
      {
      
        dup_patient_ID_A8_C10 <-  intersect(rt_clinical_data_test[Number == 'A8']$Patient_ID,
                                          rt_clinical_data_test[Number == 'C10']$Patient_ID)
      
      
        for (i in dup_patient_ID_A8_C10){
          
          C10_tsb = rt_clinical_data_test[Number == 'C10' & Patient_ID == i ]$Tumor_Sample_Barcode_backup
          
          if(length(C10_tsb) != 1){
            print(i)
          }
         
          rt_backup[Number == 'A8' & Patient_ID == i ]$Race  =  unique(rt_backup[Tumor_Sample_Barcode %in% C10_tsb]$Race)
        
          dup_tumor_sample_barcodes_to_be_deleted <-  c(C10_tsb,
                                                      dup_tumor_sample_barcodes_to_be_deleted)
      }
      
      {
        dup_patient_ID_C9_C10 <-  intersect(rt_clinical_data_test[Number == 'C10']$Patient_ID,
                                            rt_clinical_data_test[Number == 'C9']$Patient_ID)
        
        for (i in dup_patient_ID_C9_C10){
          
          C10_tsb = rt_clinical_data_test[Number == 'C10' & Patient_ID == i ]$Tumor_Sample_Barcode_backup
          
          if(length(C10_tsb) != 1){
            print(i)
          }
          
          if(rt_backup[Number == 'C9' & Patient_ID == i ]$Race == 'Unknown'){
            rt_backup[Number == 'C9' & Patient_ID == i ]$Race  =  unique(rt_backup[Tumor_Sample_Barcode %in% C10_tsb]$Race)
          }
          
          dup_tumor_sample_barcodes_to_be_deleted <-  c(C10_tsb,
                                                        dup_tumor_sample_barcodes_to_be_deleted)
      }
      
     
  }
      
   
  }
  
  
  rt@clinical.data <- rt_backup
  
  barcode_new <- setdiff(rt@clinical.data$Tumor_Sample_Barcode, 
                         dup_tumor_sample_barcodes_to_be_deleted)
  
  rt@clinical.data <- rt@clinical.data[Tumor_Sample_Barcode %in% barcode_new]
  
  table(rt@clinical.data$Race)
 
  unique(rt@clinical.data$Panel)
  panel_size_manual <- data.table(
    SEQ_ASSAY_ID = c('FUSCC_ClinSeq','WESPlus gene panel','OncoPanel_AMC_version3','Sidra-LUMC'),
    total_size = c(2,65.755218,60,60)
  )
  
  capture_size <- rbind(panel_size_manual,panel_size) %>% 
    rename( Panel=SEQ_ASSAY_ID )
  
  capture_size[Panel == 'VHIO-300']$total_size =  1.3

  rt@clinical.data <- rt@clinical.data[Panel %in% capture_size$Panel]

  rt <- subsetMaf(rt,tsb = rt@clinical.data$Tumor_Sample_Barcode)
  
 
}



{#upset plot
  table(rt@clinical.data$Center)
  unique(rt@clinical.data$Panel)
  
  
  panel_hugo_symbol = left_join(rt@data[,c('Hugo_Symbol','Tumor_Sample_Barcode')],
                                rt@clinical.data[,c('Tumor_Sample_Barcode','Panel')],
                                by = 'Tumor_Sample_Barcode')
 
  panel_genenumber = panel_hugo_symbol %>%
    group_by(Panel) %>%       
    distinct(Hugo_Symbol) %>%           
    summarise(non_duplicate_count = n()) %>% as.data.frame()
  
  df_unique <- unique(panel_hugo_symbol[, c("Hugo_Symbol", "Panel")])
  
  df_wide <- df_unique %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = Panel, values_from = value, values_fill = 0) %>% 
    as.data.frame()
  
  rownames(df_wide) <- df_wide$Hugo_Symbol
  df_wide <- df_wide[, -1]
  
  
 
  filter_panel <-  panel_genenumber %>% 
    filter(non_duplicate_count <300)
  
  df_wide_filtered <- df_wide[, !colnames(df_wide) %in% filter_panel[,1]]
  
  table(rt@clinical.data$Center)
  
  filtered_panel_country <- rt@clinical.data[,.(patients=.N),by=c('Panel','Center','Country')] %>% 
    filter(Panel %in% colnames(df_wide_filtered))  %>% 
    left_join(panel_genenumber,by = 'Panel')
  
  
  names(df_wide_filtered)
  df_wide_filtered <- df_wide_filtered[, -c(11:13)]
  
  country_patient_sum <- filtered_panel_country[, .(total_patients = sum(patients)), by = Country]
  

  pdf(file =  'targeted_panel.pdf',width = 10,height = 8)
  upset(df_wide_filtered, sets = colnames(df_wide_filtered), keep.order = TRUE, order.by = "freq", 
        main.bar.color = "dodgerblue", sets.bar.color = "orange")
  dev.off()
  
  
  pdf(file = 'all_panel.pdf',width = 10,height = 8)
  upset(df_wide, sets = colnames(df_wide), keep.order = TRUE, order.by = "freq", 
        main.bar.color = "dodgerblue", sets.bar.color = "orange")
  dev.off()
  

  common_genes <- panel_hugo_symbol[Panel %in% colnames(df_wide_filtered)] %>%  
    group_by(Hugo_Symbol) %>%
    filter(n_distinct(Panel) == ncol(df_wide_filtered)) %>%
    distinct(Hugo_Symbol)
  
  
  clinical_filtered <- rt@clinical.data[Panel %in% colnames(df_wide_filtered)]
  rt_panel_filtered <- subsetMaf(rt,tsb = clinical_filtered$Tumor_Sample_Barcode)
  rt_panel_common_genes <- subsetMaf(rt_panel_filtered, genes = common_genes$Hugo_Symbol)
  
  
}


{
  table(cl$Age)
  
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
