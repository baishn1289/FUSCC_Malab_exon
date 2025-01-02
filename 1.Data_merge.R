library(maftools)
library(mclust)
library(NMF)
library(data.table)
library(tidyverse)
library(readxl)
library(UpSetR)

dir= './picture'

{    #MAF_data_importing
  
  for (i in list.files(pattern = "rawdata_[A-Z][0-9]+\\.Rdata$")) {
      load(i) }
    
    
  load('FUSCC_2020to2024.Rdata')
    
  cl <- dplyr::bind_rows(cl_FUSCC_2020to2024, 
                           cl_A8,cl_A10,
                           cl_C8, cl_C10,
                           cl_D1,cl_D7)
    
    
    
  cl = cl[!duplicated(cl$Tumor_Sample_Barcode),]
   
  rt = merge_mafs(list(rt_FUSCC_2020to2024,
                         rt_A8,rt_A10,
                         rt_C8,rt_C10,
                         rt_D1,rt_D7))
  
  rt@clinical.data$Country <- as.factor(rt@clinical.data$Country)
  rt@clinical.data$Age_class <- as.factor(rt@clinical.data$Age_class)
    
  { #Race_data_sorting
   
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
       
    rt@clinical.data <- left_join(rt@clinical.data[,-c('Race')],cl[, c('Tumor_Sample_Barcode', 'Race')],
                                      by ='Tumor_Sample_Barcode')
  
  }
    
    
  {    #Sex_Histology_type_Tumor_site_data_sorting
  
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
      
      
    {    #capture_size
        
    panel_size_manual <- data.table(
          SEQ_ASSAY_ID = c('FUSCC_ClinSeq','WESPlus','SureSelect_Exon_V6','Sidra-LUMC'),
          total_size = c(2.3,65.755218,60,60) )
        
        
    capture_size <- rbind(panel_size_manual,panel_size) %>% 
          rename( Panel=SEQ_ASSAY_ID )
      
    capture_size[Panel == 'VHIO-300']$total_size =  1.3
     
  }
    
    
  {    #modify patient ID of patients from MKSCC in AACR GEINE Project
  
    rt_clinical_data_filter <- rt@clinical.data  %>% 
                                mutate(Patient_ID = ifelse(Number == 'C10' & Center == 'MSK',
                                                            gsub("GENIE-MSK-", "", Patient_ID),
                                                                  Patient_ID)  )
                                                   
    #retain panels more than 10 patients
    rt_clinical_panel_filter <- rt_clinical_data_filter %>% 
                                  group_by(Panel) %>% 
                                  reframe(
                                    Panel_count = n_distinct(Patient_ID)
                                  ) %>% filter(Panel_count >=10)
        
    #find patients with multiple samples and retain only one sample
    
    #For the same patients, the race information of the samples in the C10 dataset was assigned to the A8 sample
    #and subsequent duplicate sample screening was performed
        
    race_info_C10 <- rt_clinical_data_filter %>%
                      filter(Number == 'C10') %>%
                      dplyr::select(Patient_ID, Race) %>% 
                      filter(!duplicated(Patient_ID))
        
        
    rt_clinical_data_filter <- rt_clinical_data_filter %>%
                                left_join(race_info_C10, by = "Patient_ID", suffix = c("", "_C10")) %>%
                                mutate(
                                  Race = if_else(Number == 'A8' & !is.na(Race_C10), Race_C10, Race)
                                ) %>%
                                dplyr::select(-Race_C10)
     
        
    # First identify patients with multiple samples
    # In patients with multiple age-mismatched samples, remove the older samples
    # In each patient's youngest age sample, retain the original sample (excluding metastatic or uncertain samples)
    # If there are still multiple primary samples, pick one primary sample at random
        
    multi_samples <- rt_clinical_data_filter[duplicated(rt_clinical_data_filter$Patient_ID) | 
                                            duplicated(rt_clinical_data_filter$Patient_ID, fromLast = TRUE), ] 
        
    multi_samples_info  <-  multi_samples %>%  
                            group_by(Patient_ID) %>% 
                            reframe(
                                Age_same = n_distinct(Age) == 1,
                                Number_same = n_distinct(Number) == 1,
                                Min_Age = min(Age),
                                Tumor_Sample_Barcode_diffage = list(Tumor_Sample_Barcode[Age != min(Age)]),
                                Min_Age_Sample_Count = sum(Age == min(Age)),
                                Tumor_Sample_Barcode_sameage = list(Tumor_Sample_Barcode[Age == min(Age)]) )
                               
        
    sameage_tsb <- unlist(multi_samples_info$Tumor_Sample_Barcode_sameage)
    diffage_tsb <- unlist(multi_samples_info$Tumor_Sample_Barcode_diffage)
        
        
    multi_samples_same_age <-multi_samples %>% 
                            filter(Tumor_Sample_Barcode %in% sameage_tsb) %>% 
                            group_by(Patient_ID) %>% 
                            reframe(
                            sample_type_count = n_distinct(Sample_type)) %>% 
                            filter(sample_type_count ==2) %>% 
                            pull(Patient_ID)
          
        
    #One patient whose sample types ware metastatic and uncertain
    #the metastatic sample of this patient was retained separately
    multi_samples_diff_sampletype <- multi_samples %>% 
                                    filter(Patient_ID %in% multi_samples_same_age,
                                           Tumor_Sample_Barcode %in% sameage_tsb,
                                           Sample_type == 'Metastasis') %>% 
                                    mutate(Tumor_Sample_Barcode = ifelse(Tumor_Sample_Barcode == 'GENIE-VHIO-1266-001',
                                                                          'GENIE-VHIO-1266-002',
                                                                          Tumor_Sample_Barcode)) %>% 
                                    pull(Tumor_Sample_Barcode)
          
    barcode_discard <- c( multi_samples_diff_sampletype, diffage_tsb)
        
    multi_samples_tobe_filtered <- multi_samples %>% 
                                    filter(!Tumor_Sample_Barcode %in% barcode_discard) %>% 
                                    filter(duplicated(Patient_ID)) %>% 
                                    pull(Tumor_Sample_Barcode)
    
    barcode_new <- setdiff(rt@clinical.data$Tumor_Sample_Barcode, 
                          c(multi_samples_tobe_filtered,barcode_discard))
        
    rt_clinical_data_filter <- rt_clinical_data_filter[Tumor_Sample_Barcode %in% barcode_new] %>% 
                              filter(Panel %in% rt_clinical_panel_filter$Panel) %>% 
                              filter(Panel %in% capture_size$Panel)
      
  }
    
  rt@clinical.data <- rt_clinical_data_filter
    
  rt@data <- rt@data[,c(2:31,33,62:66,75,76)]
    
  rt <- subsetMaf(rt,tsb =rt@clinical.data$Tumor_Sample_Barcode )
}


{    #Remove any genes in the MAF file that are not in corresponding panel or hg19
  
  assays_genes <- read_xlsx("panel_genelist_modified.xlsx")%>% 
                    rename(Panel = SEQ_ASSAY_ID) %>% 
                    na.omit() %>% 
                    filter(Panel %in% rt@clinical.data$Panel) %>% 
                    group_by(Panel) %>% 
                    filter(!duplicated(Hugo_Symbol)) %>% 
                    ungroup() %>% 
                    as.data.frame()
     
  diffgenes_bind <-  data.frame(Hugo_Symbol = character(),
                                 Panel = character(),
                                 Number=character())
    
  for (i in unique(rt@clinical.data$Panel)){
      
    rt_inspection = left_join(rt@data,rt@clinical.data,
                                by = 'Tumor_Sample_Barcode') %>% 
                    filter(Panel == i)
      
    assays_genes_inspection = assays_genes %>% filter(Panel == i)
      
    diff_genes = setdiff(rt_inspection$Hugo_Symbol,
                         assays_genes_inspection$Hugo_Symbol)
      
    if(length(diff_genes) != 0){
        
        df = rt_inspection %>% 
             filter(Hugo_Symbol %in% diff_genes) %>% 
             dplyr::select(Hugo_Symbol,Panel,Number)
            
       diffgenes_bind = rbind(diffgenes_bind,df)
        
    }
}
  
  diffgenes_bind_list <- diffgenes_bind %>% 
                          group_by(Panel) %>% 
                          filter(!duplicated(Hugo_Symbol)) %>% 
                          reframe(diff_gene_list = list(Hugo_Symbol)) 
    
  
  rt_data <-  left_join(rt@data,
                        rt@clinical.data[,.(Tumor_Sample_Barcode,Panel)],
                        by='Tumor_Sample_Barcode')%>%
              left_join(diffgenes_bind_list, by = "Panel") %>%
              group_by(Panel) %>%
              filter(!Hugo_Symbol %in% unlist(diff_gene_list)) %>%
              ungroup() %>%
              dplyr::select(-c(diff_gene_list, Panel)) %>%
              setDT() 
      
   
  rt <- read.maf(maf = rt_data,clinicalData = rt@clinical.data)
}

{    #upset plot

  panel_genenumber <- assays_genes %>% 
                      group_by(Panel) %>% 
                      reframe(non_duplicate_count = n()) %>% 
                      as.data.frame()
    
  
  #In order to obtain more shared genes in subsequent oncoplot, panels with fewer than 300 genes were removed here
  filter_panel <-  panel_genenumber %>% 
                   filter(non_duplicate_count <300)
   
  
  df_wide <- assays_genes %>%
              mutate(value = 1) %>%
              pivot_wider(names_from = Panel, values_from = value, values_fill = 0) %>% 
              as.data.frame()
    
  rownames(df_wide) <- df_wide$Hugo_Symbol
  df_wide <- df_wide[, -1]
    
  df_wide_filtered <- df_wide[, !colnames(df_wide) %in% filter_panel[,1]]
    
   
  filtered_panel_country <- rt@clinical.data[,.(patients=.N),by=c('Panel','Center','Country')] %>% 
                            filter(Panel %in% colnames(df_wide_filtered))  %>% 
                            left_join(panel_genenumber,by = 'Panel')
  
  country_patient_sum <- filtered_panel_country[, .(total_patients = sum(patients)), by = Country]
  
  }
  
    
  pdf(file =  './picture/1_targeted_panel.pdf',width = 14,height = 15)
  upset(df_wide_filtered, sets = colnames(df_wide_filtered), keep.order = TRUE, order.by = "freq", 
          main.bar.color = "dodgerblue", sets.bar.color = "orange")
  dev.off()
  
  pdf(file = './picture/1_all_panel.pdf',width =14,height = 12)
  upset(df_wide, sets = colnames(df_wide), keep.order = TRUE, order.by = "freq", 
          main.bar.color = "dodgerblue", sets.bar.color = "orange")
  dev.off()
    
  common_genes <- assays_genes[assays_genes$Panel %in% colnames(df_wide_filtered),] %>%  
                  group_by(Hugo_Symbol) %>%
                  filter(n_distinct(Panel) == ncol(df_wide_filtered)) %>%
                  distinct(Hugo_Symbol)
    
  clinical_filtered <- rt@clinical.data[Panel %in% colnames(df_wide_filtered)]
  rt_panel_filtered <- subsetMaf(rt,tsb = clinical_filtered$Tumor_Sample_Barcode)
  rt_panel_common_genes <- subsetMaf(rt_panel_filtered, genes = common_genes$Hugo_Symbol)


{
  
  # rt <- subsetMaf(rt,tsb = rt@clinical.data[Country == 'China']$Tumor_Sample_Barcode)
  # rt <- subsetMaf(rt,tsb = rt@clinical.data[Country == 'China' & MSI_Status %in% c('MSS','pMMR')]$Tumor_Sample_Barcode)
  clin.EOCRC <- subsetMaf(maf=rt, tsb=rt@clinical.data[Age <= 50]$Tumor_Sample_Barcode)
  clin.LOCRC <- subsetMaf(maf=rt, tsb=rt@clinical.data[Age > 50]$Tumor_Sample_Barcode)

  #rt_panel_filtered
  clin.EOCRC_panel_filtered <- subsetMaf(maf=rt_panel_filtered, tsb=clin.EOCRC.sample, isTCGA=FALSE)
  clin.LOCRC_panel_filtered <- subsetMaf(maf=rt_panel_filtered, tsb=clin.LOCRC.sample, isTCGA=FALSE)
  
  
  #rt_panel_common_genes
  clin.EOCRC_panel_common_genes <- subsetMaf(maf=rt_panel_common_genes, tsb=clin.EOCRC.sample, isTCGA=FALSE)
  clin.LOCRC_panel_common_genes <- subsetMaf(maf=rt_panel_common_genes, tsb=clin.LOCRC.sample, isTCGA=FALSE)
  
  
  current_data = format(Sys.time(),"%Y-%m-%d")
  save(rt,clin.EOCRC,clin.LOCRC,
       file=paste(current_data,"overall_data.Rdata",sep = '_'))
  
  save(filter_panel,assays_genes,capture_size,
       file=paste(current_data,"filtered_panel_data.Rdata",sep = '_'))
  
  save(clinical_filtered,rt_panel_filtered,clin.EOCRC_panel_filtered,clin.LOCRC_panel_filtered,
       file=paste(current_data,"panel_filtered_data.Rdata",sep = '_'))
  
  save(rt_panel_common_genes,clin.EOCRC_panel_common_genes,clin.LOCRC_panel_common_genes,
       file=paste(current_data,"common_genes_data.Rdata",sep = '_'))
  
}
