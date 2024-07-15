
######Calculate EOCRC_freq and LOCRC_freq under the TMB classification of each country

generate_Country_specific_df <- function(cl_data, TMB_status, rt.EOCRC, rt.LOCRC, kokka,XOCRC) {
  
  # Extract sample ids for specific countries
  Country_samples <- cl_data %>% filter(Country == kokka) 
  Country_samples_ids <- Country_samples$Tumor_Sample_Barcode
  
  # Check whether Country_samples_ids exists in rt.EOCRC and rt.LOCRC
  common_ids_eocrc <- intersect(Country_samples_ids, rt.EOCRC@clinical.data$Tumor_Sample_Barcode)
  common_ids_locrc <- intersect(Country_samples_ids, rt.LOCRC@clinical.data$Tumor_Sample_Barcode)
  
  if (length(common_ids_eocrc) == 0 | length(common_ids_locrc) == 0) {
    return(data.frame())
  }
  
 
  clin.EOCRC_Country <- subsetMaf(maf = rt.EOCRC, tsb = common_ids_eocrc)
  clin.LOCRC_Country <- subsetMaf(maf = rt.LOCRC, tsb = common_ids_locrc)
  

  fvsm_Country <- mafCompare(m1 = clin.EOCRC_Country, m2 = clin.LOCRC_Country, 
                             m1Name = "EOCRC", m2Name = "LOCRC", minMut = 5)
  
  Country_df <- fvsm_Country$results
  Country_df$EOCRC_freq <- round(Country_df$EOCRC / fvsm_Country$SampleSummary$SampleSize[1], 4)
  Country_df$LOCRC_freq <- round(Country_df$LOCRC / fvsm_Country$SampleSummary$SampleSize[2], 4)
  Country_df$Country <- kokka
  Country_df$TMB_Status <- TMB_status
  
  Country_df <- Country_df %>% filter(EOCRC != 0 & LOCRC != 0)
  Country_df <- Country_df %>% filter(is.finite(or))
  
  Country_df <-  Country_df %>%
    arrange(or) %>% 
    group_by(or) %>%
    mutate(pval = pval[order(pval)]) %>%
    ungroup() %>% slice(1:10)
  
  if(XOCRC == 'EOCRC'){
  Country_EOCRC_df <- fvsm_Country$results
  Country_EOCRC_df$EOCRC_freq <- round(Country_EOCRC_df$EOCRC / fvsm_Country$SampleSummary$SampleSize[1], 4)
  Country_EOCRC_df$TMB_Status <- TMB_status
  names(Country_EOCRC_df)[names(Country_EOCRC_df) =="EOCRC_freq"] <- paste(kokka,"_EOCRC_freq",sep = '')
  
  Country_EOCRC_df <-Country_EOCRC_df %>% filter(EOCRC != 0 & LOCRC != 0)
  Country_EOCRC_df <-Country_EOCRC_df %>% filter(is.finite(or))
  
  Country_EOCRC_df <-   Country_EOCRC_df%>%
    arrange(or) %>% 
    group_by(or) %>%
    mutate(pval = pval[order(pval)]) %>%
    ungroup() %>% slice(1:10)
  
  Country_EOCRC_df <- Country_EOCRC_df[,c(1,9,10)]
  return(Country_EOCRC_df)
  
  }
  
  
  if(XOCRC == 'LOCRC'){
    Country_LOCRC_df <- fvsm_Country$results
    Country_LOCRC_df$LOCRC_freq <- round(Country_LOCRC_df$LOCRC / fvsm_Country$SampleSummary$SampleSize[2], 4)
    Country_LOCRC_df$TMB_Status <- TMB_status
    names(Country_LOCRC_df)[names(Country_LOCRC_df) =="LOCRC_freq"] <- paste(kokka,"_LOCRC_freq",sep = '')
    
    Country_LOCRC_df <-Country_LOCRC_df %>% filter(EOCRC != 0 & LOCRC != 0)
    Country_LOCRC_df <-Country_LOCRC_df %>% filter(is.finite(or))
    
    Country_LOCRC_df <-   Country_LOCRC_df%>%
      arrange(or) %>% 
      group_by(or) %>%
      mutate(pval = pval[order(pval)]) %>%
      ungroup() %>% slice(1:10)
    
    Country_LOCRC_df <- Country_LOCRC_df[,c(1,9,10)]
    return(Country_LOCRC_df)
    
  }
  
  
  # save
  #write.table(Country_df, file = paste('./', TMB_status, '_', kokka, '_EOCRC_vs_LOCRC_freq.tsv', sep = ''), 
         #     quote = FALSE, row.names = FALSE, sep = "\t")
 
}



countries <- unique(rt@clinical.data$Country)
clinical_data <- rt@clinical.data
df_countries <-data.frame()
clinical_data <- clinical_data[,-c(18,38)]
Country_df_backup<- data.frame(Hugo_Symbol ='',TMB_Status = '')

#Country_alone_EOCRC
df_countries_EOCRC <-data.frame(Hugo_Symbol ='',TMB_Status = '')
for (Country in countries){
  overall = generate_Country_specific_df(cl_data = clinical_data, TMB_status = "overall", rt.EOCRC = clin.EOCRC, rt.LOCRC = clin.LOCRC, kokka = Country,XOCRC = "EOCRC")
  hyper = generate_Country_specific_df(cl_data = clinical_data,TMB_status = "hyper",rt.EOCRC = clin.EOCRC_hyper,rt.LOCRC = clin.LOCRC_hyper,kokka = Country,XOCRC = "EOCRC")
  nonhyper = generate_Country_specific_df(cl_data = clinical_data,TMB_status = "nonhyper",rt.EOCRC = clin.EOCRC_nonhyper,rt.LOCRC = clin.LOCRC_nonhyper,kokka = Country,XOCRC = "EOCRC")
  
  
  if(length(overall)=='0'|length(hyper)=='0'|length(nonhyper)=='0' )
  {return(Country_df_backup)
    }
    
  df_countries_EOCRC <- bind_rows(df_countries_EOCRC, overall , hyper , nonhyper)
 
}

names(df_countries_EOCRC)[names(df_countries_EOCRC) =="Hugo_Symbol"] <- "Gene"
df_countries_EOCRC <- df_countries_EOCRC[-1,]
write.csv(x = df_countries_EOCRC,file = 'df_countries_alone_tree_EO.csv',row.names = F)

#
#Country_alone_LOCRC
df_countries_LOCRC <-data.frame(Hugo_Symbol ='',TMB_Status = '')
for (Country in countries){
  overall = generate_Country_specific_df(cl_data = clinical_data, TMB_status = "overall", rt.EOCRC = clin.EOCRC, rt.LOCRC = clin.LOCRC, kokka = Country,XOCRC = "LOCRC")
  hyper = generate_Country_specific_df(cl_data = clinical_data,TMB_status = "hyper",rt.EOCRC = clin.EOCRC_hyper,rt.LOCRC = clin.LOCRC_hyper,kokka = Country,XOCRC = "LOCRC")
  nonhyper = generate_Country_specific_df(cl_data = clinical_data,TMB_status = "nonhyper",rt.EOCRC = clin.EOCRC_nonhyper,rt.LOCRC = clin.LOCRC_nonhyper,kokka = Country,XOCRC = "LOCRC")
  
  if(length(overall)=='0'|length(hyper)=='0'|length(nonhyper)=='0' )
  {return(Country_df_backup)
  }
  
  df_countries_LOCRC <- bind_rows(df_countries_LOCRC, overall , hyper , nonhyper)
 
}

names(df_countries_LOCRC)[names(df_countries_LOCRC) =="Hugo_Symbol"] <- "Gene"
df_countries_LOCRC <- df_countries_LOCRC[-1,]
write.csv(x = df_countries_LOCRC,file = 'df_countries_alone_tree_LO.csv',row.names = F)

