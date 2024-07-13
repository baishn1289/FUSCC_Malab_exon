library(maftools)
library(data.table)
library(dplyr)
library(mclust)
library(NMF)
library(pheatmap)
library(data.table)
library(stringr)
library(Cairo)

setwd('~/R/exome 7.12')
rm(list=ls())
gc()

dir= './picture'

{
  for (i in list.files(pattern = "rawdata_[A-Z][0-9]+\\.Rdata$")) {
    load(i)
  }
  
  load('FUSCC_2020to2024.Rdata')

  cl_C9$Hypertension = as.numeric(cl_C9$Hypertension)
  cl_C9$`Diabetes mellitus` = as.numeric(cl_C9$`Diabetes mellitus`)
  cl_A7$`Family history of cancer`=as.character(cl_A7$`Family history of cancer`)
  
  cl <- dplyr::bind_rows(cl_FUSCC_2020to2024, 
                         cl_A7                         ,cl_A8,cl_A10,
                         cl_C8, cl_C9,cl_C10,
                         cl_D1,cl_D7
                          )
  
  table(cl$Country)
  cl = cl[!duplicated(cl$Tumor_Sample_Barcode),]
  table(cl$Country)
  
  rt = merge_mafs(list(rt_FUSCC_2020to2024,
                       rt_A7,rt_A8,rt_A10,
                       rt_C8,rt_C9,rt_C10,
                       rt_D1,rt_D7
  ))

  patient_country_count = bind_rows(rt_A7@clinical.data,
                                    rt_A8@clinical.data,
                                    rt_A10@clinical.data,
                                    rt_C8@clinical.data,
                                    rt_C9@clinical.data,
                                    rt_C10@clinical.data,
                                    rt_D1@clinical.data,
                                    rt_D7@clinical.data,
                                    rt_FUSCC_2020to2024@clinical.data
                                    )
  patient_MSKCC =  bind_rows(rt_A7@clinical.data,
                             rt_A8@clinical.data,
                             rt_C9@clinical.data)
  
  patient_country_count = patient_country_count[!duplicated(patient_country_count$Patient_ID),]
  
  
  
  rm(cl_A7,cl_A8,cl_A10,
     cl_C8,cl_C9,cl_C10,cl_D1,cl_D7
     ,cl_FUSCC_2020to2024
     
  )
  
  rm(rt_A7,rt_A8,rt_A10,
     rt_C8,rt_C9,rt_C10,rt_D1,rt_D7
     ,rt_FUSCC_2020to2024
   
  )
 
  rt@clinical.data$Country <- as.factor(rt@clinical.data$Country)
  rt@clinical.data$Age_class <- as.factor(rt@clinical.data$Age_class)
  
  rt@data <- rt@data[,c(1:31)]
  rt@maf.silent <-rt@maf.silent[,c(1:31)]
  
  rt <- subsetMaf(rt,tsb = rt@clinical.data$Tumor_Sample_Barcode)

}


{
  table(cl$Age)
 
  clin.EOCRC.sample <- subset(cl, Age <= 50)$Tumor_Sample_Barcode # %>% na.omit()
  clin.LOCRC.sample <- subset(cl, Age > 50)$Tumor_Sample_Barcode # %>% na.omit()
  
  length(clin.EOCRC.sample)
  length(clin.LOCRC.sample)
  

  clin.EOCRC <- subsetMaf(maf=rt, tsb=clin.EOCRC.sample, isTCGA=FALSE)
  clin.LOCRC <- subsetMaf(maf=rt, tsb=clin.LOCRC.sample, isTCGA=FALSE)
  
  current_data = format(Sys.time(),"%Y-%m-%d")
  save(cl,clin.EOCRC,clin.LOCRC,rt,file=paste("rt_Ver2.Rdata",current_data))
  
}

